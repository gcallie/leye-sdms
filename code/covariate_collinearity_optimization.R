# Check collinearity of raster layers, select candidate model covariates, determine appropriate scale of effect for leye SDMs

# Load packages
library(tidyverse)
library(terra)
library(usdm)
library(car)
#library(corrplot)

# Read in reprojected, cropped data
az <- vect("data/processed/az_boundary.shp")
elev <- rast("data/processed/elev_az.tif")
slope <- rast("data/processed/slope_az.tif")
aspect <- rast("data/processed/aspect_az.tif")
tri <- rast("data/processed/tri_az.tif")
tpi <- rast("data/processed/tpi_az.tif")
karst_density <- rast("data/processed/karst_density_az.tif")
karst_presence <- rast("data/processed/karst_presence_az.tif")
karst_dist <- rast("data/processed/karst_dist_az.tif")
bioclim_covs <- rast("data/processed/bioclim_az.tif")
leye <- read_csv("data/processed/leye_occ.csv")
leye_used <- project(vect(leye,
                 geom = c("longitude", "latitude"),
                 crs = "epsg:4326"), crs(crs(bioclim_covs))) # Read in leye occurrence points and reproject to Lambert Equal Area

# Generate available points using 50x multiplier
system.time(leye_avail <- spatSample(az, 
                         size = nrow(leye) * 50,
                         method = "random"))

# Preview objects
# AZ Polygon
az

# Covariate rasters
elev
slope
aspect
tri
tpi
karst_density
karst_dist
karst_presence
bioclim_covs

# Create multi-layer raster for all covariates
all_paths <- c("data/processed/elev_az.tif",
               "data/processed/slope_az.tif",
               "data/processed/aspect_az.tif",
               "data/processed/tri_az.tif",
               "data/processed/tpi_az.tif",
               "data/processed/karst_density_az.tif",
               "data/processed/karst_dist_az.tif", 
               "data/processed/bioclim_az.tif")
all_covs <- rast(all_paths)
all_covs
plot(all_covs)

# LEYE occurrence data
leye_used
leye_avail

# Test by plotting one raster with AZ boundary and occurrence points overlaid 
plot(elev) # Plot test raster
plot(az, add = TRUE) # Plot AZ boundary
points(leye_avail, # Plot available points
       col ="white", 
       pch = 1,
       cex = 0.75)
points(leye_used, # Plot used points
       col = "red")

# Create data frame of available point coordinates
leye_avail <- crds(leye_avail)
leye_avail <- as.data.frame(leye_avail)
leye_avail <- leye_avail |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 0) |> 
  mutate(weight = 5000) # weight available points as 5000
head(leye_avail) # Preview object

# Create data frame of used point coordinates
leye_used <- crds(leye_used)
leye_used <- as.data.frame(leye_used)
leye_used <-  leye_used |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 1) |> 
  mutate(weight = 1) # weight used points as 1 
head(leye_used) # Preview object

# Join used and available points into a single data frame
leye_pts <- rbind(leye_used, leye_avail)
head(leye_pts) # Preview object
table(leye_pts$pa) # Check that number of used and available points looks correct

# Pull out values for points from rasters

## All covariates
covs_leye <- extract(x = all_covs, 
                      y = leye_pts[, c("longitude", "latitude")],
                      ID = FALSE)
head(covs_leye) # Preview object


# CHECK RASTER COLLINEARITY #

# Check collinearity across AZ using LEYE used/available points

# VIF: keep covariates with VIF ~ 5 or less
vif_test <- vifstep(covs_leye, th = 5, keep = c("PPT_wt", "PPT_sm", "DD_0"))
vif_test
covs_list <- vif_test@results$Variables # Preview object
covs_list

# Pairwise (Pearson's correlation)
all_cor <- cor(covs_leye) # Create correlation matrix

## Summarize highly correlated (cor > 0.7) covariates in a table
all_cor_df <- as.data.frame(all_cor) # Create data frame object of correlation values
all_cor_df <- all_cor_df |> 
  rownames_to_column("cov1") |> 
  pivot_longer(DD_0:Tave_wt, # Update this based on covariates included
               names_to = "cov2",
               values_to = "cor") |> 
  filter(cov1 != cov2) |> 
  filter(abs(cor) >= 0.7) |> 
  mutate(cov_order = paste(cov1, cov2) |> # Start workflow to remove duplicates
           strsplit(split = ' ') |> 
           map_chr( ~ sort(.x) |>  
                      paste(collapse = ' '))) |> 
  mutate(count = 1) |> 
  group_by(cov_order) |> 
  mutate(num = row_number()) |> 
  filter(num != 2) |> 
  ungroup() |> 
  select(cov1, cov2, cor)
all_cor_df # Preview object

# Compare binomial GLMs with MWMT vs DD_0 vs Eref

## Add covariate point values to data frame with LEYE lat/long, presence/absence
leye_data <- bind_cols(covs_leye, leye_pts)
head(leye_data) # Preview object

## Model with MWMT
mwmt_mod <- glm(pa ~ slope + aspect + TPI + PPT_sm + PPT_wt + karst_dist_m + karst_density + MWMT, data = leye_data,
                weight = weight,
                family = "binomial")
summary(mwmt_mod)

## Model with DD_0
dd0_mod <- glm(pa ~ slope + aspect + TRI + TPI + PPT_sm + PPT_wt + karst_dist_m + karst_density + DD_0, data = leye_data,
                weight = weight,
                family = "binomial")
summary(dd0_mod)

## Model with Eref
eref_mod <- glm(pa ~ slope + aspect + TRI + TPI + PPT_sm + PPT_wt + karst_dist_m + karst_density + Eref, data = leye_data,
               weight = weight,
               family = "binomial")
summary(eref_mod)

## Compare model AIC
aic <- data.frame(
  mod = c("mwmt_mod", "dd0_mod", "eref_mod"),
  AIC = c(AIC(mwmt_mod), AIC(dd0_mod), AIC(eref_mod)))
aic <-  aic |> 
  arrange(AIC) |>  
  mutate(dAIC = AIC - min(AIC)) # Calculate delta AIC
aic # Preview object

### DD_0 is strongest predictor

# Covariate list for forage species SDMs: slope, aspect, TPI, TRI, DD_0, PPT_sm, PPT_wt
covs_keep_forage <- c("slope", "aspect", "TPI", "TRI", "DD_0", "PPT_sm", "PPT_wt")

# Covariate list for LEYE SDMs (pre-pseudooptimization, not including forage SDM rasters): slope, aspect, TPI, DD_0, PPT_sm, PPT_wt, karst_dist_m, karst_density
covs_keep_leye <- c("slope", "aspect", "TPI", "TRI", "DD_0", "PPT_sm", "PPT_wt", "karst_dist_m", "karst_density")

covs <- subset(all_covs, covs_keep_leye) # Keep subset of covs for leye
names(covs) # Check that correct layers included
covs_forage <- subset(all_covs, covs_keep_forage)  # Keep subset of covs for forage species
names(covs_forage) # Check that correct layers included

## PREP DATA FOR FORAGE SPECIES SDMs ##

# Saguaros
cagi <- read_csv("data/processed/cagi_occ.csv")
cagi_used <- project(vect(cagi,
                          geom = c("longitude", "latitude"),
                          crs = "epsg:4326"), crs(crs(bioclim_covs))) # Read in cagi occurrence points and reproject to Lambert Equal Area

# Generate available points using 1x multiplier
system.time(cagi_avail <- spatSample(az, 
                         size = nrow(cagi),
                         method = "random"))

# Create data frame of available point coordinates
cagi_avail <- crds(cagi_avail)
cagi_avail <- as.data.frame(cagi_avail)
cagi_avail <- cagi_avail |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 0) |> 
  mutate(weight = 5000) # weight available points as 5000
head(cagi_avail) # Preview object

# Create data frame of used point coordinates
cagi_used <- crds(cagi_used)
cagi_used <- as.data.frame(cagi_used)
cagi_used <-  cagi_used |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 1) |> 
  mutate(weight = 1) # weight used points as 1 
head(cagi_used) # Preview object

# Join used and available points into a single data frame
cagi_pts <- rbind(cagi_used, cagi_avail)
head(cagi_pts) # Preview object
table(cagi_pts$pa) # Check that number of used and available points looks correct

# Pull out values for points from rasters
covs_cagi <- extract(x = covs_forage,
                     y = cagi_pts[, c("longitude", "latitude")],
                     ID = FALSE)
head(covs_cagi) # Preview object

# Test by plotting one raster with AZ boundary and occurrence points overlaid
plot(elev) # Plot test raster
plot(az, add = TRUE) # Plot AZ boundary
points(cagi_avail, # Plot available points
       col ="white",
       pch = 1,
       cex = 0.75)
points(cagi_used, # Plot used points
       col = "red")

# Organ Pipes
stth <- read_csv("data/processed/stth_occ.csv")
stth_used <- project(vect(stth,
                          geom = c("longitude", "latitude"),
                          crs = "epsg:4326"), crs(crs(bioclim_covs))) # Read in stth occurrence points and reproject to Lambert Equal Area

# Generate available points using 25x multiplier
system.time(stth_avail <- spatSample(az, 
                                     size = nrow(stth) * 25,
                                     method = "random"))

# Create data frame of available point coordinates
stth_avail <- crds(stth_avail)
stth_avail <- as.data.frame(stth_avail)
stth_avail <- stth_avail |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 0) |> 
  mutate(weight = 5000) # weight available points as 5000
head(stth_avail) # Preview object

# Create data frame of used point coordinates
stth_used <- crds(stth_used)
stth_used <- as.data.frame(stth_used)
stth_used <-  stth_used |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 1) |> 
  mutate(weight = 1) # weight used points as 1 
head(stth_used) # Preview object

# Join used and available points into a single data frame
stth_pts <- rbind(stth_used, stth_avail)
head(stth_pts) # Preview object
table(stth_pts$pa) # Check that number of used and available points looks correct

# Pull out values for points from rasters
covs_stth <- extract(x = covs_forage,
                     y = stth_pts[, c("longitude", "latitude")],
                     ID = FALSE)
head(covs_stth) # Preview object

# Test by plotting one raster with AZ boundary and occurrence points overlaid
plot(elev) # Plot test raster
plot(az, add = TRUE) # Plot AZ boundary
points(stth_avail, # Plot available points
       col ="white",
       pch = 1,
       cex = 0.75)
points(stth_used, # Plot used points
       col = "red")

# Agaves
agave <- read_csv("data/processed/agave_occ.csv")
agave_used <- project(vect(agave,
                          geom = c("longitude", "latitude"),
                          crs = "epsg:4326"), crs(crs(bioclim_covs))) # Read in agave occurrence points and reproject to Lambert Equal Area

# Generate available points using 10x multiplier
system.time(agave_avail <- spatSample(az, 
                                     size = nrow(agave) * 10,
                                     method = "random"))

# Create data frame of available point coordinates
agave_avail <- crds(agave_avail)
agave_avail <- as.data.frame(agave_avail)
agave_avail <- agave_avail |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 0) |> 
  mutate(weight = 5000) # weight available points as 5000
head(agave_avail) # Preview object

# Create data frame of used point coordinates
agave_used <- crds(agave_used)
agave_used <- as.data.frame(agave_used)
agave_used <-  agave_used |> 
  rename("longitude" = "x", "latitude" = "y") |> 
  mutate(pa = 1) |> 
  mutate(weight = 1) # weight used points as 1 
head(agave_used) # Preview object

# Join used and available points into a single data frame
agave_pts <- rbind(agave_used, agave_avail)
head(agave_pts) # Preview object
table(agave_pts$pa) # Check that number of used and available points looks correct

# Pull out values for points from rasters
covs_agave <- extract(x = covs_forage,
                     y = agave_pts[, c("longitude", "latitude")],
                     ID = FALSE)
head(covs_agave) # Preview object

# Test by plotting one raster with AZ boundary and occurrence points overlaid
plot(elev) # Plot test raster
plot(az, add = TRUE) # Plot AZ boundary
points(agave_avail, # Plot available points
       col ="white",
       pch = 1,
       cex = 0.75)
points(agave_used, # Plot used points
       col = "red")

# Bind lat/long and presence/absence to the extracted covariate values for forage species
cagi_data <- bind_cols(covs_cagi, cagi_pts)
head(cagi_data) # Preview object
stth_data <- bind_cols(covs_stth, stth_pts)
head(stth_data) # Preview object
agave_data <- bind_cols(covs_agave, agave_pts)
head(agave_data) # Preview object

# Save as csv files
write_csv(cagi_data, "data/processed/cagi_data.csv")
write_csv(stth_data, "data/processed/stth_data.csv")
write_csv(agave_data, "data/processed/agave_data.csv")
