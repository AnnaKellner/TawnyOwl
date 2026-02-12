# Packages ----------------------------------------------------------------
library(dplyr)             # For data processing functions
library(tidyr)             # For data processing functions
library(sp)                # For spatial data
library(sf)                # For spatial data
library(stargazer)         # for model tables
library(lubridate)         # For temporal data
library(ggplot2)           # For plotting
library(lme4)              # For GLMMs
library(broom.mixed)       # For extracting model summaries
library(gt)                # For creating tables
library(dplyr)             # For data manipulation
library(purrr)             # For efficient model handling
library(patchwork)         # For plotting
library(ggeffects)         # For plotting model fit
library(lme4)              # For GLMM fitting
library(ggeffects)         # For prediction figures
library(performance)       # For model diagnostics
library(see)               # For model diagnostics (in combo with performance package)
library(MuMIn)             # For model averaging
library(AICcmodavg)        # For model averaging
library(lmerTest)          # For P-values
library(ggspatial)         # For OSM basemap

# ggplot theme ------------------------------------------------------------
theme_set(
  theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 15) # Facet labels
    )
)

# Master data -------------------------------------------------------------
df <- read.csv("C:\\Users\\r01ak21\\OneDrive - University of Aberdeen\\OneDrive_1_10-01-2022\\Career\\PhD\\Reading\\My PhD\\Chapters\\Chapter 2\\Data\\processed_fulldataset_1993_2023.csv",
               row.names = NULL, 
               header = TRUE,
               stringsAsFactors = FALSE)

# REMOVE two out of range coordinates 
df <- df %>%
  filter(Latitude >= 54.5 & Latitude <= 56)

cat("Observations after initial coordinate filtering:", nrow(df), "\n")

# remove path and polygon_id
df <- df %>% 
  dplyr::select(-path, -polygon_id)

summary(df)

# remove where eggs is =12 since this is a mandarin and has mistakenly been kept
df <- df[df$initialeggs != 12, ]

# Create a corrected version of the dataset
# biologically impossible for NFledged to be more than initialeggs
df <- df %>%
  mutate(initialeggs = ifelse(NFledged > initialeggs, NFledged, initialeggs))

# Load polygons -------
polygons <- st_read("C:\\Users\\r01ak21\\OneDrive - University of Aberdeen\\OneDrive_1_10-01-2022\\Career\\PhD\\Reading\\My PhD\\Chapters\\Chapter 2\\Data\\KR.shp")

# Assign a unique polygon ID if not already present
polygons <- polygons %>% mutate(polygonID = row_number())

# Ensure polygons have the correct CRS
polygons <- st_transform(polygons, crs = 4326)

# Load supergrid -------
supergrid <- st_read("C:\\Users\\r01ak21\\OneDrive - University of Aberdeen\\OneDrive_1_10-01-2022\\Career\\PhD\\Reading\\My PhD\\Chapters\\Chapter 2\\Data\\GB_GridSquares\\gb-grids_5503583\\1km_grid_region.shp")

# Assign a unique GridID if not already present
supergrid <- supergrid %>% mutate(GridID = row_number())

# Ensure supergrid has the correct CRS
supergrid <- st_transform(supergrid, crs = 4326)

# SPATIAL ASSIGNMENT
# Remove rows with missing coordinates
df <- df %>%
  filter(!is.na(Longitude), !is.na(Latitude))

cat("Observations after removing NA coordinates:", nrow(df), "\n")

# Convert df to sf
df_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

# Ensure CRS consistency (defensive)
polygons  <- st_transform(polygons, 4326)
supergrid <- st_transform(supergrid, 4326)

# Ensure IDs exist (defensive)
polygons  <- polygons  %>% mutate(polygonID = row_number())
supergrid <- supergrid %>% mutate(GridID = row_number())

# ---- POLYGON ASSIGNMENT (DISTANCE – GUARANTEED) ----
distances <- st_distance(df_sf, polygons)
closest_polygon <- apply(distances, 1, which.min)

df_sf$polygonID <- polygons$polygonID[closest_polygon]

# Hard guarantee
stopifnot(!any(is.na(df_sf$polygonID)))

# ---- GRID ASSIGNMENT (SPATIAL JOIN) ----
df_sf <- st_join(
  df_sf,
  supergrid %>% dplyr::select(GridID),
  left = TRUE
)

stopifnot(!any(is.na(df_sf$GridID)))

# Drop geometry for downstream analysis
df <- df_sf %>% st_drop_geometry()

cat("Observations after spatial assignment:", nrow(df), "\n")

cat("Observations after renaming:", nrow(df), "\n")
str(df)

# CLEAN UP NAs -------
df <- df %>%
  mutate(
    NHatched = ifelse(is.na(NHatched), 0, NHatched),  # NAs in NHatched become 0
    psi = ifelse(is.na(psi), 0, psi)  # NAs in psi become 0 (no martens)
  )

cat("Observations after cleaning NAs:", nrow(df), "\n")

# Data processing -------
df <- df |> 
  mutate(year = round(year),
         year_f = factor(year)) |>
  mutate(
    initialeggs = ifelse(initialeggs == 0 & (NFledged > 0 | NHatched > 0), 
                         pmax(NFledged, NHatched, na.rm = TRUE), initialeggs),
    fledge_fail = ifelse(NFledged == 0 | is.na(NFledged), 1, 0),
    NHatched = ifelse(initialeggs == 0 & NFledged == 0, 0, NHatched)
  ) |> 
  filter(!is.na(fledge_fail)) |>
  filter(year >= 1993 & year <= 2023) |> 
  mutate(psi = ifelse(is.na(psi), 0, psi))

df$forest <- trimws(df$Forest)

df <- df %>%
  mutate(
    GridID = factor(GridID.x),
    polygonID = factor(polygonID)
  )

cat("Observations after data processing:", nrow(df), "\n")

# Calculating mean clutch size per polygon and VAR -------
mean_eggs <- df %>%
  filter(initialeggs > 0) %>%
  group_by(polygonID, year_f) %>%
  summarise(mean_eggs = mean(initialeggs, na.rm = TRUE), .groups = "drop")

mean_fledged <- df %>%
  filter(NFledged > 0) %>%
  group_by(polygonID, year_f) %>%
  summarise(mean_fledged = mean(NFledged, na.rm = TRUE), .groups = "drop")

df <- df %>%
  left_join(mean_eggs, by = c("polygonID", "year_f")) %>%
  left_join(mean_fledged, by = c("polygonID", "year_f"))

df <- df %>%
  mutate(VAR = mean_fledged / mean_eggs)

cat("Observations before scaling:", nrow(df), "\n")

# Scale numeric variables used in model
df <- df %>%
  mutate(
    scaled_eggs = scale(mean_eggs)[,1],
    scaled_meanfledge = scale(mean_fledged)[,1],
    scaled_VAR = scale(VAR)[,1],
    scaled_psi = scale(psi)[,1]
  )

# Removing duplicates -------
df <- df %>%
  group_by(year, nestbox_ID2, forest) %>%
  filter(!(forest == "Kielder" & year < 2012 & duplicated(across(everything()), fromLast = TRUE))) %>%
  ungroup()

cat("Observations after removing duplicates:", nrow(df), "\n")

# Filter for complete cases for modeling -------
model_vars1 <- c("scaled_psi", "scaled_VAR", "scaled_eggs", "fledge_fail", "year_f", "GridID", "polygonID")

df <- df %>%
  filter(complete.cases(dplyr::select(., all_of(model_vars1))))

cat("Observations with complete cases for modeling:", nrow(df), "\n")
cat("Failures:", sum(df$fledge_fail == 1, na.rm = TRUE), "\n")
cat("Failure rate %:", round(sum(df$fledge_fail == 1, na.rm = TRUE) / nrow(df) * 100, 1), "\n")


cat("NA polygonIDs before modelling:", sum(is.na(df$polygonID)), "\n")
stopifnot(!any(is.na(df$polygonID)))

cat("NA GridIDs before modelling:", sum(is.na(df$GridID)), "\n")
stopifnot(!any(is.na(df$GridID)))

# Ensure year is numeric
df <- df %>%
  mutate(year = as.numeric(year))

# Polygons of interest
polys_check <- c(4, 5, 6, 7)

# 7. FIT GLMM MODELS ####
m1 <- glmer(
  fledge_fail ~ scaled_psi + scaled_VAR * scaled_eggs +
    (1|year_f) + (1|GridID) + (1|polygonID),
  data = df,
  family = binomial
)
summary(m1)

# 2. failure ~ psi + mean clutch size
m2 <- glmer(
  fledge_fail ~ scaled_psi + scaled_eggs +
    (1 | year_f) + (1 | GridID),
  data = df,
  family = binomial
)
summary(m2)
# 3. failure ~ psi * mean clutch size
m3 <- glmer(
  fledge_fail ~ scaled_psi * scaled_eggs +
    (1 | year_f) + (1 | GridID),
  data = df,
  family = binomial
)
summary(m3)
# 4. failure ~ psi + VAR + mean clutch size
m4 <- glmer(
  fledge_fail ~ scaled_psi + scaled_VAR + scaled_eggs +
    (1 | year_f) + (1 | GridID),
  data = df,
  family = binomial
)
summary(m4)
# 5. failure ~ psi * mean clutch size + VAR
m5 <- glmer(
  fledge_fail ~ scaled_psi * scaled_eggs + scaled_VAR +
    (1 | year_f) + (1 | GridID),
  data = df,
  family = binomial
)
summary(m5)
# 6. failure ~ VAR + mean clutch size
m6 <- glmer(
  fledge_fail ~ scaled_VAR + scaled_eggs +
    (1 | year_f) + (1 | GridID),
  data = df,
  family = binomial
)
summary(m6)
# 7. failure ~ psi + VAR
m7 <- glmer(
fledge_fail ~ scaled_psi + scaled_VAR +
    (1 | year_f) + (1 | GridID),
  data = df,
  family = binomial
)
summary(m7)
# AICc MODEL SELECTION ####
# Sample size: n = 1979 breeding attempts

models <- list(m1, m2, m3, m4, m5, m6, m7)

model.names <- c(
  "failure ~ psi + VAR * mean clutch size",
  "failure ~ psi + mean clutch size",
  "failure ~ psi * mean clutch size",
  "failure ~ psi + VAR + mean clutch size",
  "failure ~ psi * mean clutch size + VAR",
  "failure ~ VAR + mean clutch size",
  "failure ~ psi + VAR"
)

aictab(
  cand.set = models,
  modnames = model.names,
  sort = TRUE
)

# AIC COMPARISON (VALID) ####
models <- list(m2, m3, m4, m5, m6, m7)
model.names <- paste0("m", 2:7)

aictab(
  cand.set = models,
  modnames = model.names,
  sort = TRUE
)

# Based on this table, the best model is Model 7
# Model 7: 
# Shows there is a significant relationship with psi. Wow I was right all along. 
# Load package
install.packages("stargazer")

# Extract AIC values and rank models
model_info <- tibble(
  Model_No = seq_along(models),
  Model = models,
  AIC = map_dbl(models, AIC)
) %>%
  arrange(AIC) %>%  # Rank models by AIC (ascending)
  mutate(Rank = row_number()) %>%
  select(Rank, Model_No, AIC)
# TAWNY OWL NEST FAILURE PLOTS ####
cat("Total observations:", nrow(df), "\n")
cat("Failures (fledge_fail == 1):", sum(df$fledge_fail == 1, na.rm = TRUE), "\n")
cat("Successes (fledge_fail == 0):", sum(df$fledge_fail == 0, na.rm = TRUE), "\n")
cat("Failure rate %:", round(sum(df$fledge_fail == 1, na.rm = TRUE) / nrow(df) * 100, 1), "\n")

# EXTRACT MODEL PREDICTIONS FOR RESULTS TEXT =================================

cat("\n=== EXTRACTING PREDICTIONS FOR RESULTS TEXT ===\n")

# Pine marten occupancy effect
psi_preds <- ggpredict(m1, terms = "scaled_psi [all]")
mean_psi <- mean(df$psi, na.rm = TRUE)
sd_psi <- sd(df$psi, na.rm = TRUE)
psi_min_scaled <- (min(df$psi, na.rm = TRUE) - mean_psi) / sd_psi
psi_max_scaled <- (max(df$psi, na.rm = TRUE) - mean_psi) / sd_psi
psi_pred_min <- psi_preds$predicted[which.min(abs(psi_preds$x - psi_min_scaled))]
psi_pred_max <- psi_preds$predicted[which.max(abs(psi_preds$x - psi_max_scaled))]

cat("Pine marten occupancy predictions:\n")
cat("Min PSI predicted failure probability:", round(psi_pred_min * 100, 1), "%\n")
cat("Max PSI predicted failure probability:", round(psi_pred_max * 100, 1), "%\n")

# Clutch size effect
eggs_preds_low <- ggpredict(m1, terms = c("scaled_eggs [all]", "scaled_VAR [-1]"))
eggs_preds_high <- ggpredict(m1, terms = c("scaled_eggs [all]", "scaled_VAR [1]"))

mean_eggs_val <- mean(df$mean_eggs, na.rm = TRUE)
sd_eggs_val <- sd(df$mean_eggs, na.rm = TRUE)

eggs_pred_low_min <- min(eggs_preds_low$predicted)
eggs_pred_low_max <- max(eggs_preds_low$predicted)
eggs_pred_high_min <- min(eggs_preds_high$predicted)
eggs_pred_high_max <- max(eggs_preds_high$predicted)

cat("\nClutch size predictions (low food availability):\n")
cat("Min clutch size failure probability:", round(eggs_pred_low_max * 100, 1), "%\n")
cat("Max clutch size failure probability:", round(eggs_pred_low_min * 100, 1), "%\n")

cat("\nClutch size predictions (high food availability):\n")
cat("Min clutch size failure probability:", round(eggs_pred_high_max * 100, 1), "%\n")
cat("Max clutch size failure probability:", round(eggs_pred_high_min * 100, 1), "%\n")

cat("\n=== END OF CRITICAL OUTPUT ===\n")
cat("Use these percentages to update your Results text\n")


# figures ######
# FIGURE 3: TAWNY OWL NEST FAILURE
# Using Model 1 (best model)
# PANEL 4A: Interaction plot (clutch size × VAR)
Low <- ggpredict(m1, terms = c("scaled_eggs [all]", "scaled_VAR [-1]"))
High <- ggpredict(m1, terms = c("scaled_eggs [all]", "scaled_VAR [1]"))

Low$VAR_Group <- "Low food availability"
High$VAR_Group <- "High food availability"

interaction_df <- rbind(Low, High)

eggs_mean <- mean(df$mean_eggs, na.rm = TRUE)
eggs_sd   <- sd(df$mean_eggs, na.rm = TRUE)
interaction_df$clutch_size <- interaction_df$x * eggs_sd + eggs_mean

# Calculate density for bottom area
dens_df <- density(df$mean_eggs, na.rm = TRUE)
dens_df <- data.frame(
  x = dens_df$x,
  y = dens_df$y / max(dens_df$y) * 0.1  # scale down to sit at bottom (0–0.1)
)

plot3a <- ggplot(interaction_df, aes(x = clutch_size, y = predicted, col = VAR_Group, fill = VAR_Group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_area(data = dens_df, aes(x = x, y = y), inherit.aes = FALSE, fill = "grey80", alpha = 0.5) +
  scale_fill_manual(values = c("Low food availability" = "#E69F00",
                               "High food availability" = "#0072B2")) +
  scale_colour_manual(values = c("Low food availability" = "#E69F00",
                                 "High food availability" = "#0072B2")) +
  scale_y_continuous(
    name = "Probability of nest failure",
    sec.axis = sec_axis(~ . / 0.1, name = "Density of mean clutch size")  # secondary axis
  ) +
  labs(x = "Mean clutch size", fill = "Food availability", col = "Food availability") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))

# PANEL 3B: Psi effect with ####
psi_preds <- ggpredict(m1, terms = "scaled_psi [all]")
psi_mean <- mean(df$psi, na.rm = TRUE)
psi_sd   <- sd(df$psi, na.rm = TRUE)
psi_preds$psi_unscaled <- psi_preds$x * psi_sd + psi_mean

dens_psi <- density(df$psi, na.rm = TRUE)
dens_psi <- data.frame(
  x = dens_psi$x,
  y = dens_psi$y / max(dens_psi$y) * 0.1
)

plot3b <- ggplot(psi_preds, aes(x = psi_unscaled, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#0072B2", alpha = 0.3) +
  geom_line(color = "#0072B2", linewidth = 1.2) +
  geom_area(data = dens_psi, aes(x = x, y = y), inherit.aes = FALSE, fill = "grey80", alpha = 0.5) +
  scale_y_continuous(
    name = "Probability of nest failure",
    sec.axis = sec_axis(~ . / 0.1, name = "Density of pine marten occupancy")
  ) +
  labs(x = "Pine marten occupancy probability") +
  theme_minimal(base_size = 16)

# Combine panels
combined_figure3 <- plot3a / plot3b + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')

# Save
ggsave("Figure3_combined_density_bottom.jpg", plot = combined_figure3, width = 12, height = 6, dpi = 300)

# changing font size####
# PANEL 3A: Interaction plot (clutch size × VAR)
plot3a <- ggplot(interaction_df, aes(x = clutch_size, y = predicted, col = VAR_Group, fill = VAR_Group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_area(data = dens_df, aes(x = x, y = y), inherit.aes = FALSE, fill = "grey80", alpha = 0.5) +
  scale_fill_manual(values = c("Low food availability" = "#E69F00",
                               "High food availability" = "#0072B2")) +
  scale_colour_manual(values = c("Low food availability" = "#E69F00",
                                 "High food availability" = "#0072B2")) +
  scale_y_continuous(
    name = "Probability of nest failure",
    sec.axis = sec_axis(~ . / 0.1, name = "Density of mean clutch size")
  ) +
  labs(x = "Mean clutch size", fill = "Food availability", col = "Food availability") +
  theme_minimal(base_size = 14) +  # slightly smaller base font
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 12),     # smaller tick labels
    axis.title = element_text(size = 13)     # slightly smaller axis titles
  )

# PANEL 3B: Psi effect with density
plot3b <- ggplot(psi_preds, aes(x = psi_unscaled, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#0072B2", alpha = 0.3) +
  geom_line(color = "#0072B2", linewidth = 1.2) +
  geom_area(data = dens_psi, aes(x = x, y = y), inherit.aes = FALSE, fill = "grey80", alpha = 0.5) +
  scale_y_continuous(
    name = "Probability of nest failure",
    sec.axis = sec_axis(~ . / 0.1, name = "Density of pine marten occupancy")
  ) +
  labs(x = "Pine marten occupancy probability") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

# Combine panels
combined_figure3 <- plot3a / plot3b + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')

# Save
ggsave("Figure3_combined_density_bottom.jpg", plot = combined_figure3, width = 12, height = 6, dpi = 300)

# TAWNY OWL NEST PREDATION & PINE MARTEN OCCUPANCY ANALYSIS ####
# --- FILE PATHS  ---
paths <- list(
  nest = "C:/Users/r01ak21/OneDrive - University of Aberdeen/OneDrive_1_10-01-2022/Career/PhD/Reading/My PhD/Chapters/Chapter 2/Data/tawnyboxlocationsallyears_1979_2024.csv",
  experiment = "C:/Users/r01ak21/OneDrive - University of Aberdeen/OneDrive_1_10-01-2022/Career/PhD/Reading/My PhD/Chapters/Chapter 2/Data/UpdatedPredationNewRules07_07_2025.csv",
  pm = "C:/Users/r01ak21/OneDrive - University of Aberdeen/OneDrive_1_10-01-2022/Career/PhD/Reading/My PhD/Chapters/Chapter 2/Data/pinemartenoccupancy.csv",
  df = "C:/Users/r01ak21/OneDrive - University of Aberdeen/OneDrive_1_10-01-2022/Career/PhD/Reading/My PhD/Chapters/Chapter 2/Data/processed_fulldataset_1993_2023.csv",
  supergrid = "C:/Users/r01ak21/OneDrive - University of Aberdeen/OneDrive_1_10-01-2022/Career/PhD/Reading/My PhD/Chapters/Chapter 2/Data/GB_GridSquares/gb-grids_5503583/1km_grid_region.shp",
  polygons = "C:/Users/r01ak21/OneDrive - University of Aberdeen/OneDrive_1_10-01-2022/Career/PhD/Reading/My PhD/Chapters/Chapter 2/Data/KR.shp"
)

# --- LOAD DATA --- ##
cat("Loading data...\n")
tawny_boxes <- read.csv(paths$nest, stringsAsFactors = FALSE)
experiment <- read.csv(paths$experiment, stringsAsFactors = FALSE)
pm <- read.csv(paths$pm, na.strings = "")
df <- read.csv(paths$df, stringsAsFactors = FALSE)  # Load from paths list
supergrid <- st_read(paths$supergrid, quiet = TRUE)
polygons <- st_read(paths$polygons, quiet = TRUE)

cat("Data loaded successfully.\n")

# PART 2: PREPARE SPATIAL DATA ####
cat("Preparing spatial data...\n")

# Clean tawny box locations
tawny_clean <- tawny_boxes %>%
  rename(BOX = nestbox_ID2, Year = year, Easting = East, Northing = North) %>%
  mutate(
    BOX = toupper(trimws(as.character(BOX))),
    Year = as.integer(Year),
    Easting = as.numeric(Easting),
    Northing = as.numeric(Northing)
  ) %>%
  filter(!is.na(BOX), !is.na(Year), !is.na(Easting), !is.na(Northing))

cat("Tawny boxes cleaned:", nrow(tawny_clean), "rows\n")

# Clean experiment data
experiment_clean <- experiment %>%
  mutate(
    BOX = toupper(trimws(as.character(BOX))),
    Year = as.integer(Year)
  ) %>%
  filter(!is.na(BOX), !is.na(Year))

cat("Experiment cleaned:", nrow(experiment_clean), "rows\n")

# Join experiment + tawny box locations
combined_df <- left_join(experiment_clean, tawny_clean, by = c("BOX", "Year")) %>%
  filter(!is.na(Easting), !is.na(Northing))

cat("After spatial join:", nrow(combined_df), "rows\n")

# Convert to sf
combined_sf <- st_as_sf(combined_df, coords = c("Easting", "Northing"), crs = 27700)

# Ensure CRS match before spatial joins
if (st_crs(combined_sf) != st_crs(supergrid)) {
  supergrid <- st_transform(supergrid, st_crs(combined_sf))
}
if (st_crs(combined_sf) != st_crs(polygons)) {
  polygons <- st_transform(polygons, st_crs(combined_sf))
}

# Add polygonID to polygons if missing
if (!"polygonID" %in% names(polygons)) {
  polygons <- polygons %>% mutate(polygonID = row_number())
}

# Spatial joins
combined_sf <- st_join(combined_sf, supergrid %>% dplyr::select(PLAN_NO))
# ---- POLYGON ASSIGNMENT (DISTANCE SAFE) ----
distances <- st_distance(combined_sf, polygons)
closest_polygon <- apply(distances, 1, which.min)

combined_sf$polygonID <- polygons$polygonID[closest_polygon]

# GUARANTEE
stopifnot(!any(is.na(combined_sf$polygonID)))

cat("Spatial joins complete\n")

# PART 3: PINE MARTEN OCCUPANCY ####

cat("Processing pine marten occupancy...\n")

# Prepare pine marten data
pm_sf <- st_as_sf(pm, coords = c("East", "North"), crs = 27700) %>%
  st_transform(st_crs(supergrid)) %>%
  st_join(supergrid %>% dplyr::select(PLAN_NO)) %>%
  mutate(
    year = occ.1.20 + 2011,
    join_id = paste(PLAN_NO, year, sep = "_")
  )

# Join psi to combined data
combined_sf <- combined_sf %>%
  mutate(join_id = paste(PLAN_NO, Year, sep = "_"))

combined_final <- left_join(
  combined_sf,
  st_drop_geometry(pm_sf),
  by = "join_id"
)

cat("Pine marten data joined\n")

# PART 4: POLYGON-LEVEL CLUTCH SIZE (from main df) ####

cat("Calculating polygon-level clutch size...\n")

# Calculate mean eggs by polygon from main df
# Note: df has column 'initialeggs' or similar - adjust if needed
polygon_means <- df %>%
  mutate(polygonID = as.integer(as.character(polygonID))) %>%
  filter(!is.na(initialeggs), initialeggs > 0) %>%
  group_by(polygonID) %>%
  summarise(polygon_mean_eggs = mean(initialeggs, na.rm = TRUE), .groups = "drop")

cat("Polygon means calculated for", nrow(polygon_means), "polygons\n")

# Add to combined data
combined_final <- combined_final %>%
  mutate(polygonID = as.integer(as.character(polygonID))) %>%
  left_join(polygon_means, by = "polygonID")

# PART 5: FINAL DATA CLEANING ####

cat("Final data cleaning...\n")

cat("Checking available columns in combined_final:\n")
print(names(combined_final))

final_dfSpatial <- combined_final %>%
  filter(Experiment %in% c("0", "1")) %>%
  mutate(
    Puncture = ifelse(
      Experiment == 1 & (is.na(Puncture) | Puncture == "" | Puncture != 1),
      0,
      Puncture
    ),
    Predated = as.integer(Predated)
  ) %>%
  filter(Experiment != 0) %>%
  # Select only columns that exist
  dplyr::select(any_of(c("BOX", "PLAN_NO", "Experiment", "Year", "Date_added", 
                         "Predated", "Puncture", "MissingEgg", "polygonID", 
                         "psi", "polygon_mean_eggs", "geometry")))


data_year <- final_dfSpatial %>% 
  filter(!is.na(Predated)) %>%  # Add this line
  st_drop_geometry() %>%
  mutate(
    Predated = as.integer(Predated),
    pred_label = ifelse(Predated == 1, "Predated", "Non-Predated")
  )

# Convert to regular dataframe with coordinates
final_df <- final_dfSpatial %>%
  mutate(
    Longitude = st_coordinates(.)[, 1],
    Latitude = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

# Add scaled variables for modeling
final_df <- final_df %>%
  mutate(
    scaled_psi = scale(psi)[, 1],
    scaled_clutch = scale(polygon_mean_eggs)[, 1]
  )

# PART 6: DESCRIPTIVE STATISTICS ####

cat("\n=== DESCRIPTIVE STATISTICS ===\n\n")
# Predation summary
cat("Overall Predation (Puncture):\n")
print(table(final_df$Puncture))

cat("\nOverall Predation (Predated):\n")
print(table(final_df$Predated))

# By year
cat("\nPredation by Year:\n")
predation_by_year <- final_df %>%
  group_by(Year) %>%
  summarise(
    n = n(),
    pred_n = sum(Puncture == 1, na.rm = TRUE),
    pred_rate = round(pred_n / n, 3),
    .groups = "drop"
  )
print(predation_by_year)

# By psi bins
cat("\nPredation by PSI bins:\n")
predation_by_psi <- final_df %>%
  mutate(psi_bin = cut(psi, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf))) %>%
  group_by(psi_bin) %>%
  summarise(
    n = n(),
    pred_n = sum(Puncture == 1, na.rm = TRUE),
    pred_rate = round(pred_n / n, 3),
    .groups = "drop"
  )
print(predation_by_psi)

# PSI summary
cat("\nPSI Summary by Predation Status:\n")
psi_summary <- final_df %>%
  group_by(Predated) %>%
  summarise(
    N = n(),
    psi_min = min(psi, na.rm = TRUE),
    psi_max = max(psi, na.rm = TRUE),
    psi_mean = round(mean(psi, na.rm = TRUE), 3),
    .groups = "drop"
  )
print(psi_summary)

# PART 7: MIXED MODELS ####

cat("\n=== FITTING MIXED MODELS ===\n\n")

# Fit GLMM with random effects
# Note: Using only polygonID since PLAN_NO not available in final data
pred_model <- glmer(
  Predated ~ scaled_psi + scaled_clutch + (1 | polygonID),
  family = binomial,
  data = final_df
)

cat("Main model (polygonID random effect):\n")
print(summary(pred_model))
cat("\nCreating prediction plots...\n")

# PSI predictions
psi_preds <- ggpredict(m1, terms = "scaled_psi [seq(-0.33, 5.45, by = 0.5)]")
mean_psi <- mean(df$psi, na.rm = TRUE)
sd_psi <- sd(df$psi, na.rm = TRUE)
psi_preds$x_unscaled <- psi_preds$x * sd_psi + mean_psi

# Plot with density distribution
p_psi <- ggplot() +
  # Main predicted line with ribbon
  geom_ribbon(data = psi_preds, aes(x = x_unscaled, ymin = conf.low, ymax = conf.high),
              fill = "skyblue", alpha = 0.3) +
  geom_line(data = psi_preds, aes(x = x_unscaled, y = predicted),
            size = 1.2, color = "blue4") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  
  # Density of psi scaled to fit
  geom_density(data = df, aes(x = psi, y = ..scaled.. * 0.2), fill = "grey40", alpha = 0.3) +
  
  # Axes and theme
  scale_y_continuous(
    name = "Predicted probability of nest failure",
    sec.axis = sec_axis(~. / 0.2, name = "Distribution of pine marten occupancy")
  ) +
  labs(x = "Pine Marten Occupancy (PSI)") +
  theme_minimal(base_size = 14)

ggsave("nest_failure_psi.jpg", p_psi, width = 7, height = 5, dpi = 300)
print(p_psi)

# PART 9: Experiment MODEL TABLE ###
# PSI predictions - EXTENDED RANGE ####
cat("\nCreating prediction plots with extended psi range...\n")

# Create a sequence from 0 to 1 psi (unscaled)
psi_sequence <- seq(0, 1, by = 0.05)

# Scale this sequence using the model's mean and SD
mean_psi <- mean(df$psi, na.rm = TRUE)
sd_psi <- sd(df$psi, na.rm = TRUE)
psi_sequence_scaled <- (psi_sequence - mean_psi) / sd_psi

# Get predictions across full range
psi_preds <- ggpredict(m1, terms = paste0("scaled_psi [", 
                                          paste(psi_sequence_scaled, collapse = ","), 
                                          "]"))
psi_preds$x_unscaled <- psi_preds$x * sd_psi + mean_psi

# Plot with extended range
p_psi <- ggplot() +
  # Main predicted line with ribbon
  geom_ribbon(data = psi_preds, aes(x = x_unscaled, ymin = conf.low, ymax = conf.high),
              fill = "skyblue", alpha = 0.3) +
  geom_line(data = psi_preds, aes(x = x_unscaled, y = predicted),
            size = 1.2, color = "blue4") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  
  # Density of psi scaled to fit
  geom_density(data = df, aes(x = psi, y = ..scaled.. * 0.2), fill = "grey40", alpha = 0.3) +
  
  # Set x-axis to go from 0 to 1
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  # Axes and theme
  scale_y_continuous(
    name = "Predicted probability of nest failure",
    sec.axis = sec_axis(~. / 0.2, name = "Distribution of pine marten occupancy")
  ) +
  labs(x = "Pine Marten Occupancy (PSI)",
       title = "Tawny Owl Nest Failure vs Pine Marten Occupancy") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave("nest_failure_psi_extended.jpg", p_psi, width = 8, height = 6, dpi = 300)
print(p_psi)

cat("Saved: nest_failure_psi_extended.jpg\n")

cat("\n=== MODEL TABLE ===\n\n")

tab_model(pred_model,
          show.re.var = TRUE,
          show.icc = FALSE,
          show.aic = TRUE,
          pred.labels = c("Intercept", "Scaled pine marten occupancy", "Scaled mean clutch size"),
          dv.labels = "Nest predation (binomial GLMM)")
# Predation Experiment Yearly Maps#######
# PART 10: MAPS - PREDATION BY YEAR (2022, 2023, 2024) 

cat("\n=== CREATING PREDATION MAPS ===\n\n")

# Transform to WGS84
final_dfSpatial_wgs <- st_transform(final_dfSpatial, crs = 4326)

# Extract coordinates
final_dfSpatial_wgs <- final_dfSpatial_wgs %>%
  mutate(
    Longitude = st_coordinates(.)[, 1],
    Latitude = st_coordinates(.)[, 2]
  )

# Get coordinate limits for consistent scaling
coord_range <- st_coordinates(final_dfSpatial_wgs)
lon_limits <- range(coord_range[, 1], na.rm = TRUE)
lat_limits <- range(coord_range[, 2], na.rm = TRUE)

# Add 5% buffer
lon_buffer <- (lon_limits[2] - lon_limits[1]) * 0.05
lat_buffer <- (lat_limits[2] - lat_limits[1]) * 0.05
lon_limits <- c(lon_limits[1] - lon_buffer, lon_limits[2] + lon_buffer)
lat_limits <- c(lat_limits[1] - lat_buffer, lat_limits[2] + lat_buffer)

# Create maps for each year
map_list <- list()

for (yr in sort(unique(final_df$Year))) {
  
  data_year <- final_dfSpatial_wgs %>% 
    filter(Year == yr) %>%
    filter(!is.na(Predated)) %>%
    filter(!is.na(Longitude), !is.na(Latitude)) %>%
    st_drop_geometry() %>%
    mutate(
      Predated = as.integer(Predated),
      pred_label = ifelse(Predated == 1, "Predated", "Non-Predated")
    )
  # Create plot with OSM background
  p <- ggplot(data_year, aes(x = Longitude, y = Latitude)) +
    annotation_map_tile(type = "osm", zoom = 12) +
    geom_point(aes(color = pred_label, fill = pred_label), size = 2, alpha = 0.8, shape = 21, stroke = 1.5) +
    scale_color_manual(values = c("Predated" = "black", "Non-Predated" = "black")) +
    scale_fill_manual(values = c("Predated" = "black", "Non-Predated" = "white")) +
    coord_sf(xlim = lon_limits, ylim = lat_limits, crs = 4326) +
    labs(title = paste(yr),
         x = "Longitude",
         y = "Latitude",
         color = "Outcome",
         fill = "Outcome") +
    theme_minimal(base_size = 14) +
    theme(aspect.ratio = 1,
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = "right")
  
  map_list[[as.character(yr)]] <- p
  
  fname <- paste0("predation_map_", yr, ".png")
  ggsave(fname, plot = p, width = 8, height = 8, dpi = 300)
  cat("Saved:", fname, "\n")
  print(p)
  
}

# Combine all maps side by side
combined_maps <- map_list[[1]] + map_list[[2]] + map_list[[3]] + 
  plot_layout(ncol = 3, guides = "collect")

ggsave("predation_maps_combined.jpg", plot = combined_maps, width = 16, height = 5, dpi = 300)
cat("Saved: predation_maps_combined.png\n")
str(df)

# SUPPLEMENTARY TABLE: Breeding attempts by polygon & year ####
  supplementary_breeding <- df %>%
  filter(!is.na(polygonID)) %>%  # Add this
  group_by(year, polygonID) %>%
  summarise(
    n_breeding_attempts = n(),
    n_successful = sum(fledge_fail == 0),
    n_failed = sum(fledge_fail == 1),
    proportion_failure = round(mean(fledge_fail), 3),
    .groups = 'drop'
  ) %>%
  arrange(year, polygonID)

# View
print(supplementary_breeding)

# Save as CSV for supplementary material
write.csv(supplementary_breeding, "Supplementary_Table_S1_BreedingAttempts.csv", row.names = FALSE)

# Optional: Create a nicely formatted version for the paper
# Pivot to wide format if you want polygon columns
supplementary_breeding_wide <- supplementary_breeding %>%
  pivot_wider(
    names_from = polygonID,
    values_from = c(n_breeding_attempts, n_successful, n_failed),
    names_sep = "_Poly"
  )

# SUPPLEMENTARY TABLE: Breeding attempts by polygon & year ###########
supplementary_breeding <- df %>%
  group_by(year, polygonID) %>%
  filter(!is.na(polygonID)) %>%  # Add this
  summarise(
    n_breeding_attempts = n(),
    n_successful = sum(fledge_fail == 0),
    n_failed = sum(fledge_fail == 1),
    proportion_failure = round(mean(fledge_fail), 3),
    .groups = 'drop'
  ) %>%
  arrange(year, polygonID)

# View
print(supplementary_breeding, n = Inf)

write.csv(supplementary_breeding, "Supplementary_Table_S1_BreedingAttempts_ByPolygonYear.csv", row.names = FALSE)

# Convert polygon object to sf
polygons_sf <- st_as_sf(polygons)

# Calculate area in km²
polygons_area <- polygons_sf %>%
  mutate(
    area_m2 = st_area(.),
    area_km2 = as.numeric(area_m2) / 1e6
  ) %>%
  st_drop_geometry() %>%
  summarise(
    mean_area = mean(area_km2, na.rm = TRUE),
    median_area = median(area_km2, na.rm = TRUE),
    min_area = min(area_km2, na.rm = TRUE),
    max_area = max(area_km2, na.rm = TRUE)
  )

print(polygons_area)

df_2023 <- df %>%
  filter(year == 2023)

median_psi_2023 <- median(df_2023$psi, na.rm = TRUE)
print(median_psi_2023)
# Scale it using your model's scaling parameters
mean_psi <- mean(df$psi, na.rm = TRUE)
sd_psi <- sd(df$psi, na.rm = TRUE)
psi_median_scaled <- (median_psi_2023 - mean_psi) / sd_psi
# Get prediction at median psi (observed)
pred_median_psi <- ggpredict(m1, terms = paste0("scaled_psi [", psi_median_scaled, "]"))
print(pred_median_psi)

# Get prediction at psi = 0 (pre-marten)
psi_zero_scaled <- (0 - mean_psi) / sd_psi
pred_zero_psi <- ggpredict(m1, terms = paste0("scaled_psi [", psi_zero_scaled, "]"))
print(pred_zero_psi)

# Get prediction at psi = 1 (extrapolated, for discussion only)
psi_one_scaled <- (1 - mean_psi) / sd_psi
pred_one_psi <- ggpredict(m1, terms = paste0("scaled_psi [", psi_one_scaled, "]"))
print(pred_one_psi)

### calculate coef 
# Using model-fitted probabilities
cat("\n=== SPATIAL HETEROGENEITY (CV) ANALYSIS ===\n\n")

# Extract exact data used in the model

model_df <- model.frame(m1)

model_df <- model_df %>%
  mutate(
    year_numeric = as.numeric(as.character(year_f))
  )

cat("Rows used in model:", nrow(model_df), "\n")

# Population-level predictions (no random effects)
model_df$pred_prob <- predict(
  m1,
  newdata = model_df,
  type = "response",
  re.form = NA
)

stopifnot(nrow(model_df) == length(model_df$pred_prob))


# CV BY YEAR (across polygons)

cv_by_year <- model_df %>%
  group_by(year_numeric, polygonID) %>%
  summarise(
    mean_pred_failure = mean(pred_prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(year_numeric) %>%
  summarise(
    cv = 100 * sd(mean_pred_failure, na.rm = TRUE) /
      mean(mean_pred_failure, na.rm = TRUE),
    n_polygons = n(),
    .groups = "drop"
  )

print(cv_by_year)


# Define periods for shading

periods_df <- tibble(
  period = c("Period 1 (1993–2003)",
             "Period 2 (2004–2013)",
             "Period 3 (2014–2023)"),
  start  = c(1993, 2004, 2014),
  end    = c(2003.99, 2013.99, 2023.99)
)


# Plot CV over time with LOESS
#     (with and without 2019)

p_cv <- ggplot(cv_by_year, aes(x = year_numeric, y = cv)) +
  
  # Period shading
  geom_rect(
    data = periods_df,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = period),
    inherit.aes = FALSE,
    alpha = 0.12
  ) +
  
  # Annual CV points
  geom_point(size = 2.2, colour = "black") +
  
  # LOESS including all years
  geom_smooth(
    method = "loess",
    span = 0.6,
    se = FALSE,
    colour = "steelblue",
    linewidth = 1.3
  ) +
  
  # LOESS excluding 2019
  geom_smooth(
    data = cv_by_year %>% filter(year_numeric != 2019),
    method = "loess",
    span = 0.6,
    se = FALSE,
    colour = "firebrick",
    linetype = "dashed",
    linewidth = 1.3
  ) +
  
  labs(
    x = "Year",
    y = "Coefficient of Variation (%)",
    title = "Temporal change in spatial heterogeneity of predicted nest failure",
    subtitle = "Solid line: LOESS including all years; dashed line: LOESS excluding 2019"
  ) +
  
  scale_fill_manual(
    values = c("Period 1 (1993–2003)" = "grey80",
               "Period 2 (2004–2013)" = "grey65",
               "Period 3 (2014–2023)" = "grey50")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_cv)

ggsave(
  "CV_spatial_heterogeneity_over_time_LOESS.jpg",
  plot = p_cv,
  width = 11,
  height = 6,
  dpi = 300
)

cat("Saved: CV_spatial_heterogeneity_over_time_LOESS.jpg\n")

# Save as RDS
# WITH 2019

# Ensure all years appear on x-axis
all_years <- seq(min(cv_by_year$year_numeric),
                 max(cv_by_year$year_numeric), 1)

# Period shading dataframe (as before)
periods_df <- data.frame(
  start = c(1993, 2004, 2014),
  end   = c(2003.99, 2013.99, 2023.99),
  period = c("Period 1", "Period 2", "Period 3")
)

p_cv <- ggplot(cv_by_year, aes(x = year_numeric, y = cv)) +
  
  # Period shading
  geom_rect(
    data = periods_df,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.15
  ) +
  
  # LOESS smooth with shaded SE (ALL years)
  geom_smooth(
    method = "loess",
    span = 0.6,
    se = TRUE,
    colour = "steelblue",
    fill = "steelblue",
    alpha = 0.25,
    linewidth = 1.6
  ) +
  
  # Annual CV points
  geom_point(
    size = 2.5,
    colour = "black"
  ) +
  
  # Vertical dashed lines at period boundaries
  geom_vline(
    xintercept = c(2003.99, 2013.99),
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  # Axes
  scale_x_continuous(
    breaks = all_years,
    limits = c(1993, 2023.5)
  ) +
  
  labs(
    x = "Year",
    y = "Coefficient of Variation (%)",
    title = "Spatial heterogeneity in predicted nest failure"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = 8
    ),
    plot.title = element_text(face = "bold")
  )

print(p_cv)

ggsave(
  "CV_spatial_heterogeneity_LOESS_shaded.jpg",
  plot = p_cv,
  width = 12,
  height = 6.5,
  dpi = 300
)

cat("Saved: CV_spatial_heterogeneity_LOESS_shaded.jpg\n")

# without 2019 ############
# CV over time with LOESS smooth + shaded uncertainty

# Ensure all years appear on x-axis
all_years <- seq(min(cv_by_year$year_numeric),
                 max(cv_by_year$year_numeric), 1)

# Period shading dataframe (as before)
periods_df <- data.frame(
  start = c(1993, 2004, 2014),
  end   = c(2003.99, 2013.99, 2023.99),
  period = c("Period 1", "Period 2", "Period 3")
)
# Ensure ALL years exist (even if CV is NA)
cv_by_year_full <- cv_by_year %>%
  complete(year_numeric = 1993:2023)

# Period shading dataframe
periods_df <- data.frame(
  start = c(1993, 2004, 2014),
  end   = c(2003.99, 2013.99, 2023.99),
  period = c("Period 1", "Period 2", "Period 3")
)

# Plot
p_cv <- ggplot(cv_by_year_full, aes(x = year_numeric, y = cv)) +
  
  # Period shading
  geom_rect(
    data = periods_df,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.15
  ) +
  
  # LOESS smooth WITH uncertainty (all years)
  geom_smooth(
    method = "loess",
    span = 0.6,
    se = TRUE,
    colour = "steelblue",
    fill = "steelblue",
    alpha = 0.25,
    linewidth = 1.6,
    na.rm = TRUE
  ) +
  
  # LOESS smooth WITHOUT 2019 (sensitivity)
  geom_smooth(
    data = cv_by_year_full %>% filter(year_numeric != 2019),
    method = "loess",
    span = 0.6,
    se = FALSE,
    colour = "firebrick",
    linetype = "dashed",
    linewidth = 1.2,
    na.rm = TRUE
  ) +
  
  # Annual CV points
  geom_point(
    size = 2.5,
    colour = "black",
    na.rm = TRUE
  ) +
  
  # Period boundaries
  geom_vline(
    xintercept = c(2003.99, 2013.99),
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  # Axes
  scale_x_continuous(
    breaks = 1993:2023,
    limits = c(1993, 2023.5)
  ) +
  
  labs(
    x = "Year",
    y = "Coefficient of Variation (%)",
    title = "Spatial heterogeneity in predicted nest failure",
    subtitle = "LOESS smooth ± SE (solid); dashed line excludes 2019"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    plot.title = element_text(face = "bold")
  )

print(p_cv)

ggsave(
  "CV_spatial_heterogeneity_LOESS_shaded.jpg",
  plot = p_cv,
  width = 12,
  height = 6.5,
  dpi = 300
)


# calculating cv without 2019
cv_by_year <- cv_by_year %>%
  mutate(
    period = case_when(
      year_numeric >= 1993 & year_numeric <= 2003 ~ "Early (1993–2003)",
      year_numeric >= 2004 & year_numeric <= 2013 ~ "Middle (2004–2013)",
      year_numeric >= 2014 & year_numeric <= 2023 ~ "Late (2014–2023)"
    )
  )

cv_period_with2019 <- cv_by_year %>%
  group_by(period) %>%
  summarise(
    mean_cv = mean(cv, na.rm = TRUE),
    .groups = "drop"
  )

print(cv_period_with2019)

cv_period_without2019 <- cv_by_year %>%
  filter(year_numeric != 2019) %>%   # ⬅️ THIS is the key change
  group_by(period) %>%
  summarise(
    mean_cv = mean(cv, na.rm = TRUE),
    .groups = "drop"
  )

print(cv_period_without2019)

# new fig 3 ###########
# Predictions at low and high VAR

Low  <- ggpredict(m1, terms = c("scaled_eggs [all]", "scaled_VAR [-1]"))
High <- ggpredict(m1, terms = c("scaled_eggs [all]", "scaled_VAR [1]"))

Low$VAR_Group  <- "Low food availability"
High$VAR_Group <- "High food availability"

interaction_df <- bind_rows(Low, High)


# Unscale eggs to original units

eggs_mean <- mean(df$mean_eggs, na.rm = TRUE)
eggs_sd   <- sd(df$mean_eggs, na.rm = TRUE)

interaction_df$clutch_size <-
  interaction_df$x * eggs_sd + eggs_mean


# Density of observed clutch size

dens_eggs <- density(df$mean_eggs, na.rm = TRUE)
dens_eggs <- data.frame(
  x = dens_eggs$x,
  y = dens_eggs$y / max(dens_eggs$y) * 0.25
)


# Plot 3A

plot3a <- ggplot(
  interaction_df,
  aes(x = clutch_size,
      y = predicted,
      colour = VAR_Group,
      fill   = VAR_Group)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.3, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_area(data = dens_eggs,
            aes(x = x, y = y),
            inherit.aes = FALSE,
            fill = "grey30",
            alpha = 0.7) +
  scale_colour_manual(values = c(
    "Low food availability"  = "#E69F00",
    "High food availability" = "#0072B2"
  )) +
  scale_fill_manual(values = c(
    "Low food availability"  = "#E69F00",
    "High food availability" = "#0072B2"
  )) +
  scale_y_continuous(
    name = "Probability of nest failure",
    sec.axis = sec_axis(~ . / 0.25,
                        name = "Density of mean clutch size")
  ) +
  labs(
    x = "Mean clutch size",
    colour = "Food availability",
    fill   = "Food availability"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

# Predictions across psi gradient

psi_preds <- ggpredict(
  m1,
  terms = "scaled_psi [all]"
)


# Unscale psi to 0–1 scale

psi_mean <- mean(df$psi, na.rm = TRUE)
psi_sd   <- sd(df$psi, na.rm = TRUE)

psi_preds$psi_unscaled <-
  psi_preds$x * psi_sd + psi_mean


# Density of observed psi

dens_psi <- density(df$psi, na.rm = TRUE)
dens_psi <- data.frame(
  x = dens_psi$x,
  y = dens_psi$y / max(dens_psi$y) * 0.25
)


# Plot 3B

plot3b <- ggplot(
  psi_preds,
  aes(x = psi_unscaled, y = predicted)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "#0072B2",
              alpha = 0.3) +
  geom_line(color = "#0072B2", linewidth = 1.2) +
  geom_area(data = dens_psi,
            aes(x = x, y = y),
            inherit.aes = FALSE,
            fill = "grey30",
            alpha = 0.7) +
  scale_x_continuous(
    name = "Pine marten occupancy probability (ψ)",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_y_continuous(
    name = "Probability of nest failure",
    sec.axis = sec_axis(~ . / 0.25,
                        name = "Density of pine marten occupancy")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 13)
  )
library(patchwork)


# Combine panels

combined_figure3 <- plot3a / plot3b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")


# Save figure

ggsave(
  filename = "Figure3_combined.jpg",
  plot     = combined_figure3,
  width    = 12,
  height   = 6,
  dpi      = 300
)

# =]


# figure 1 #########
# TAWNY OWL MAP (2023)

# Load data
df <- read.csv("C:\\Users\\r01ak21\\OneDrive - University of Aberdeen\\OneDrive_1_10-01-2022\\Career\\PhD\\Reading\\My PhD\\Chapters\\Chapter 2\\Data\\processed_fulldataset_1993_2023.csv",
               row.names = NULL, 
               header = TRUE,
               stringsAsFactors = FALSE)

# Load polygons
polygons <- st_read("C:\\Users\\r01ak21\\OneDrive - University of Aberdeen\\OneDrive_1_10-01-2022\\Career\\PhD\\Reading\\My PhD\\Chapters\\Chapter 2\\Data\\KR.shp")

# PROCESS DATA
# Remove rows with missing coordinates
df <- df %>%
  filter(!is.na(Longitude), !is.na(Latitude))

# Convert to spatial
df_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

# Ensure CRS consistency
polygons <- st_transform(polygons, 4326)
# Assign polygons by nearest distance
distances <- st_distance(df_sf, polygons)
closest_polygon <- apply(distances, 1, which.min)
df_sf$polygonID <- row_number(polygons)[closest_polygon]

# Remap polygon IDs: old -> new
id_remap <- c(1, 4, 3, 5, 2, 7, 6)
df_sf$polygonID <- id_remap[df_sf$polygonID]

# Now swap just 4 and 6
df_sf$polygonID <- ifelse(df_sf$polygonID == 4, 6, 
                          ifelse(df_sf$polygonID == 6, 4, df_sf$polygonID))

# MAP (2023)

# Filter to 2023
owl_2023 <- df_sf %>%
  filter(year == 2023)

# Rebuild polygon boundaries as concave hulls around assigned points
polygons_rebuilt <- owl_2023 %>%
  group_by(polygonID) %>%
  summarise(geometry = st_union(geometry)) %>%
  mutate(geometry = st_concave_hull(geometry, ratio = 0.9)) %>%
  mutate(geometry = st_buffer(geometry, dist = 300))

# Plot
ggplot() +
  annotation_map_tile(type = "osm") +
  
  # Rebuilt polygons
  geom_sf(
    data = polygons_rebuilt,
    fill = NA,
    colour = "black",
    linewidth = 0.7
  ) +
  
  # Owl boxes
  geom_sf(
    data = owl_2023,
    colour = "black",
    size = 4,
    alpha = 0.8
  ) +
  
  # Polygon labels
  geom_sf_text(
    data = polygons_rebuilt %>%
      mutate(label_geom = st_centroid(geometry)) %>%
      st_set_geometry("label_geom"),
    aes(label = polygonID),
    size = 6,
    fontface = "bold", 
    colour = "black"
  ) +
  
  coord_sf(expand = FALSE) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

ggsave("tawny_owl_map_2023.jpg", 
       width = 16,
       height = 10, 
       dpi = 300)

# Save as RDS to share
saveRDS(df_sf, "owl_spatial_data.rds")
saveRDS(polygons_rebuilt, "polygons_rebuilt_2023.rds")

cat("Map saved: tawny_owl_map_2023.jpg\n")
cat("Data saved as RDS files for sharing\n")
# save as RDS to share 

oldK <- read.csv("C:\\Users\\r01ak21\\OneDrive - University of Aberdeen\\OneDrive_1_10-01-2022\\Career\\PhD\\Reading\\My PhD\\Chapters\\Chapter 2\\Data\\SteveLocations.csv",
               row.names = NULL, 
               header = TRUE,
               stringsAsFactors = FALSE)
str(oldK)

# Extract numeric box number (remove letter suffix)
oldK$BoxNum <- gsub("[A-Z]$", "", oldK$BOX)

# Parse dates - handle empty strings
parse_date <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA
  as.Date(paste0(x, "-01"), format = "%b-%y-%d")
}

oldK$StartDate <- parse_date(oldK$DER1)
oldK$EndDate <- parse_date(oldK$DRM2)

# Get unique boxes with their overall ranges
box_ranges <- oldK %>%
  group_by(BoxNum) %>%
  summarise(
    Earliest_Start = min(StartDate, na.rm = TRUE),
    Latest_End = max(EndDate, na.rm = TRUE),
    Num_Deployments = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    Earliest_Start = if_else(is.infinite(Earliest_Start), as.Date(NA), Earliest_Start),
    Latest_End = if_else(is.infinite(Latest_End), as.Date(NA), Latest_End)
  ) %>%
  arrange(BoxNum)

print(box_ranges)
nrow(box_ranges)
# Create a year sequence for each box's active period
box_years <- box_ranges %>%
  rowwise() %>%
  mutate(
    # Handle NA end dates by using a cutoff year (e.g., 1998)
    Latest_End = if_else(is.na(Latest_End), as.Date("1998-12-31"), Latest_End),
    years = list(seq(year(Earliest_Start), year(Latest_End), by = 1))
  ) %>%
  unnest(years) %>%
  select(BoxNum, years)

# Count unique boxes per year
boxes_per_year <- box_years %>%
  group_by(years) %>%
  summarise(
    Num_Boxes = n_distinct(BoxNum),
    .groups = 'drop'
  ) %>%
  arrange(years)

print(boxes_per_year)

