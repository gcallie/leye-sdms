# Visualize data partitions, build MaxEnt models, and save raw model output

# Load packages
library(tidyverse)
library(terra)
library(ENMeval)
library(ecospat)

# Define species
sp <- c("cagi", "agave", "stth")

for(species in sp){
  
# Create list of forage species
if (species == "cagi" | species == "stth" | species == "agave") {
  forage <- TRUE 
  } else {
  forage <- FALSE
}
  
# Read in data

## Species points (used/available)
occs_pa <- read_csv(paste0("data/processed/", species, "_used_thin.csv"))
occs <- occs_pa |> 
  select(-pa, -weight) # Create smaller df with just lat/long
bg_pa <- read_csv(paste0("data/processed/", species, "_avail.csv"))
bg <- bg_pa |> # Create smaller df with just lat/long
  select(-pa, -weight)

## Study area polygon/environmental covariate rasters
az <- vect("data/processed/az_boundary.shp")
az_buffered <- vect("data/processed/az_buffered.shp")

if(forage == TRUE) {
  envs <- rast("data/processed/envs_forage.tif")
} else if (species == "leye") {
  envs <- rast("data/processed/envs_leye.tif")
} else {
  warning("Species not defined.")
}
proj <- crs(envs) # Store projection
envs <- crop(envs, az_buffered, mask = TRUE) # Crop environmental raster to study area

# Partition used/avaialable data into validation and training bins for k-fold cross-validation (spatial block)

# ## Generate SpatVector of used/available points
# used_sv <- vect(occs,
#               geom = c("longitude", "latitude"),
#               crs = proj)
# used_sv # Preview object
# avail_sv <- vect(bg,
#             geom = c("longitude", "latitude"),
#             crs = proj)
# avail_sv # Preview object
# 
# plot(envs[[1]])
# plot(avail_sv, add = TRUE, col = "white", cex = 0.5)
# plot(used_sv, add = TRUE, col = "red", cex = 0.5)
# 
# ## Create blocks
# block <- get.block(occs, bg, orientation = "lat_lon")
# table(block$occs.grp) # Preview block sizes
# 
# # ## Visualize blocks
# evalplot.grps(pts = occs, pts.grp = block$occs.grp, envs = envs) +
#    ggplot2::ggtitle("Spatial block partitions: occurrences")
# 
#  evalplot.grps(pts = bg, pts.grp = block$bg.grp, envs = envs) +
#    ggplot2::ggtitle("Spatial block partitions: background")
# 
# ## Create histogram to compare how the different the environment in each partition is from the others by generating MESS predictions with each partition as the reference
#  occs.z <- cbind(occs, terra::extract(envs, occs, ID = FALSE))
#  bg.z <- cbind(bg, terra::extract(envs, bg, ID = FALSE))
#  evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = block$occs.grp,
#                       bg.grp = block$bg.grp, ref.data = "occs")

# Build candidate models

## Define tuning parameters
tune.args <- list(fc = c("L", "LQ"), rm = 1:3)

## Fit models
system.time(mod <- ENMevaluate(occs = occs, 
                   envs = envs, 
                   bg = bg, 
                   algorithm = 'maxnet', 
                   partitions = 'block',
                   tune.args = tune.args,
                   parallel = TRUE))
  
# Save raw model output
saveENMevaluation(mod, paste0("outputs/", species, "_maxent_mod.rds"))

}
