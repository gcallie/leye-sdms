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

# Standardize naming conventions
all_covs <- all_covs |> 
  rename(tri = TRI, tpi = TPI, dd0 = DD_0, eref = Eref, mcmt = MCMT, mwmt = MWMT, pptat = PPT_at, pptsp = PPT_sp, pptsm = PPT_sm, pptwt = PPT_wt, tavesm = Tave_sm, tavewt = Tave_wt, tavesp = Tave_sp, taveat = Tave_at, karstdist = karst_dist_m, karstdens = karst_density)
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
vif_test <- vifstep(covs_leye, th = 5, keep = c("pptwt", "pptsm", "dd0"))
vif_test
covs_list <- vif_test@results$Variables # Preview object
covs_list

# Pairwise (Pearson's correlation)
all_cor <- cor(covs_leye) # Create correlation matrix

## Summarize highly correlated (cor > 0.7) covariates in a table
all_cor_df <- as.data.frame(all_cor) # Create data frame object of correlation values
all_cor_df <- all_cor_df |> 
  rownames_to_column("cov1") |> 
  pivot_longer(dd0:tavewt, # Update this based on covariates included
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

# Compare binomial GLMs with mwmt vs dd0 vs eref // compare slope vs tri

## Add covariate point values to data frame with LEYE lat/long, presence/absence
leye_data <- bind_cols(covs_leye, leye_pts)
head(leye_data) # Preview object

## Model with MWMT
mwmt_mod <- glm(pa ~ slope + aspect + tpi + pptsm + pptwt + karstdist + karstdens + mwmt, data = leye_data,
                weight = weight,
                family = "binomial")
summary(mwmt_mod)

## Model with DD_0
dd0_mod <- glm(pa ~ slope + aspect + tpi + pptsm + pptwt + karstdist + karstdens + dd0, data = leye_data,
                weight = weight,
                family = "binomial")
summary(dd0_mod)

## Model with Eref
eref_mod <- glm(pa ~ slope + aspect + tpi + pptsm + pptwt + karstdist + karstdens + eref, data = leye_data,
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

## Model with slope
slope_mod <- glm(pa ~ slope + aspect + tpi + pptsm + pptwt + karstdist + karstdens + dd0, data = leye_data,
               weight = weight,
               family = "binomial")
summary(slope_mod)

## Model with tri
tri_mod <- glm(pa ~ tri + aspect + tpi + pptsm + pptwt + karstdist + karstdens + dd0, data = leye_data,
                weight = weight,
                family = "binomial")
summary(tri_mod)

## Compare model AIC
aic <- data.frame(
  mod = c("slope_mod", "tri_mod"),
  AIC = c(AIC(slope_mod), AIC(tri_mod)))
aic <-  aic |> 
  arrange(AIC) |>  
  mutate(dAIC = AIC - min(AIC)) # Calculate delta AIC
aic # Preview object

# Slope is barely better

# Covariate list for forage species SDMs: slope, aspect, tpi, tri, dd0, pptsm, pptwt
covs_keep_forage <- c("slope", "aspect", "tpi", "dd0", "pptsm", "pptwt")

# Covariate list for LEYE SDMs (pre-pseudooptimization, not including forage SDM rasters): slope, aspect, tpi, tri, dd0, pptsm, pptwt, karstdist, karstdens
covs_keep_leye <- c("slope", "aspect", "tpi", "dd0", "pptsm", "pptwt", "karstdist", "karstdens")

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


## PSEUDOOPTIMIZATION ##

# Pull out values for points from rasters for only non-collinear covariates
covs_leye <- extract(x = covs,
                     y = leye_pts[, c("longitude", "latitude")],
                     ID = FALSE)
head(covs_leye) # Preview object

# Bind extracted covariate values to leye lat/long, presence/absence
leye_data <- bind_cols(covs_leye, leye_pts)
head(leye_data)

# Create SpatVector object of leye used/available points
leye_pts_sv <- vect(leye_data,
      geom = c("longitude", "latitude"),
      crs = crs(bioclim_covs))

# Create copy of leye_data data frame for storing new values extracted at buffers
leye_pseudo <- leye_data

# Create vector object of buffer radii (2.5km, 5km, 10km, 25km)
buffers <- c(2500, 5000, 10000, 25000) 

# Extract mean values for each covariate at different scales of effect and store in leye_pseudo object 

for(i in 1:length(buffers)){
  
  buffer_i <- buffers[i]
  
  cat("Running buffer size =", buffer_i,"\n")
  
  pt_buffers <- buffer(leye_pts_sv,
                           width = buffer_i)
  
# Extract values for all covariates
cov_buffer_i <- extract(covs,
                          pt_buffers,
                          fun = "mean")
  
# Add buffer size to column names
names(cov_buffer_i) <- paste0(names(cov_buffer_i), "_", buffer_i)

# Add the columns from cov_buffer_i to leye_pseudo
leye_pseudo <- cbind(leye_pseudo, cov_buffer_i) 

}

# Remove ID columns
leye_pseudo <- leye_pseudo |> 
  select(-c(ID_2500, ID_5000, ID_10000, ID_25000)) |> 
  rename(slope_pt = slope, aspect_pt = aspect, tpi_pt = tpi, dd0_pt = dd0, pptsm_pt = pptsm, pptwt_pt = pptwt, karstdist_pt = karstdist, karstdens_pt = karstdens)

summary(leye_pseudo) # View summary

# Compare and select scales of effect

## Prepare input variables
buffers <- c("pt", 2500, 5000, 10000, 25000) # Add point scale to list of buffers
buffers # Preview object

## Prepare storage data frame
col_names <- c("Covariate","Scale","AIC","Delta_AIC","AICw") # Create list of columns
scales <- as.data.frame(matrix(NA, nrow = length(covs_keep_leye)*length(buffers),
                               ncol = length(col_names)))
colnames(scales) <- col_names # Rename columns
i.row <- 1

##  Generate models and AIC values 
for(i in 1:length(covs_keep_leye)){
  cov_i <- covs_keep_leye[i]
  cat("Starting", cov_i, "\n")
  data_i <- leye_pseudo[, c("pa", "weight", paste0(cov_i,"_", buffers))]
  cov_cols <- paste0(cov_i, "_", buffers)
  
  for(j in cov_cols){
    tmp_col <- j
    tmp_data <- data_i[,c("pa", "weight", tmp_col)]
    tmp_model <- glm(tmp_data[, 1] ~ tmp_data[, 3], weights = tmp_data[, 2], family = binomial)
    scales$Covariate[i.row] <- cov_i
    scales$AIC[i.row] <- AIC(tmp_model)
    scales$Scale[i.row] <- strsplit(j, "_")[[1]][2]
    cat("Finished scale", strsplit(j, "_")[[1]][2], "\n")
    i.row <- i.row + 1
  }
  min_i <- min(scales$AIC[which(scales$Covariate == cov_i)])
  delta_i <- scales$AIC[which(scales$Covariate == cov_i)] - min_i
  w_i <- Weights(scales$AIC[which(scales$Covariate == cov_i)])
  scales$AICw[which(scales$Covariate == cov_i)] <- w_i
  scales$Delta_AIC[which(scales$Covariate == cov_i)] <- delta_i
  
}

# Remove distance to karst covariate since it does not need to be pseudooptimized
scales <- scales |> 
  filter(!str_detect(Covariate, "karstdist"))


scale_plot <- ggplot(data = scales,
                     aes(x = Scale, y = Delta_AIC)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = ''), color = "black") +
  geom_point(color = "black")+
  facet_wrap(~ Covariate)+
  theme_bw()
scale_plot

# Pull out scale of effect for each covariate
scale_of_effect <- scales |> 
  group_by(Covariate) |> 
  slice_min(Delta_AIC)
scale_of_effect # Preview object

# aspect: pt
# dd0: 25000
# karstdens: pt
# pptsm: 5000
# pptwt: 5000
# slope: 5000
# tpi: pt 




