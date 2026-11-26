# Visualize data partitions, build MaxEnt models, and save raw model output

# Load packages
library(tidyverse)
library(terra)
library(ENMeval)
library(ecospat)

# Define species
#sp <- "cagi"
#sp <- "stth"
#sp <- "agave"
sp <- "leye"
  
# Create list of forage species
if (sp == "cagi" | sp == "stth" | sp == "agave") {
  forage <- TRUE 
  } else {
  forage <- FALSE
}
  
# Read in data

## Species points (used/available)
occs_pa <- read_csv(paste0("data/processed/", sp, "_used.csv"))
occs <- occs_pa |> 
  select(-pa, -weight) # Create smaller df with just lat/long
bg_pa <- read_csv(paste0("data/processed/", sp, "_avail.csv"))
bg <- bg_pa |> # Create smaller df with just lat/long
  select(-pa, -weight)

## Study area polygon/environmental covariate rasters
az <- vect("data/processed/az_boundary.shp")
az_buffered <- vect("data/processed/az_buffered.shp")

if(forage == TRUE) {
  envs <- rast("data/processed/covs4foragesdms.tif")
} else if (sp == "leye") {
  envs <- rast("data/processed/covs4leyesdm_z.tif")
} else {
  warning("Species not defined.")
}

proj <- crs(envs) # Store projection
envs <- crop(envs, az_buffered, mask = TRUE) # Crop environmental rasters to study area

# Thin data to retain only one occurrence per grid cell
occs.cells <- terra::extract(envs[[1]], occs, cellnumbers = TRUE, ID = FALSE) # Number raster cells
occs.cell.dups <- duplicated(occs.cells[, 1]) # Find cases where cells contain more than one used point
occs <- occs[!occs.cell.dups, ] # Remove extra points
occs

# Partition used/available data into validation and training bins for k-fold cross-validation (spatial block)

## Generate SpatVector of used points
used_sv <- vect(occs,
              geom = c("longitude", "latitude"),
              crs = proj)
used_sv # Preview object

plot(envs[["slope_5000"]]) # Plot test raster
plot(used_sv, add = TRUE, col = "white", pch = 1) # Overlay used points

## Create blocks
block <- get.block(occs, bg, orientation = "lat_lon")
table(block$occs.grp) # Preview block sizes

## Visualize blocks
evalplot.grps(pts = occs, pts.grp = block$occs.grp, envs = envs) +
   ggplot2::ggtitle("Spatial block partitions: occurrences")
ggsave(paste0("outputs/", sp, "_occ_blocks.png"), width = 3, height = 4, dpi = 300)

evalplot.grps(pts = bg, pts.grp = block$bg.grp, envs = envs) +
   ggplot2::ggtitle("Spatial block partitions: background")
ggsave(paste0("outputs/", sp, "_bg_blocks.png"), width = 3, height = 4, dpi = 300)

## Create histogram to compare how the different the environment in each partition is from the others by generating MESS predictions with each partition as the reference
occs.z <- cbind(occs, terra::extract(envs, occs, ID = FALSE))
bg.z <- cbind(bg, terra::extract(envs, bg, ID = FALSE))
evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = block$occs.grp,
                      bg.grp = block$bg.grp, ref.data = "occs")
ggsave(paste0("outputs/", sp, "_blocks_mess.png"), width = 3, height = 4, dpi = 300)

# Build candidate models

## Define tuning parameters
tune.args <- list(fc = c("L", "LQ"), rm = 1:4) # Linear, quadratic, and hinge functional forms with regularization multipliers of 1, 2, 3, 4

## Fit models
system.time(mod <- ENMevaluate(occs = occs, 
                   envs = envs, 
                   bg = bg, 
                   doClamp = TRUE,
                   algorithm = 'maxnet', 
                   partitions = 'block',
                   tune.args = tune.args,
                   parallel = FALSE))
  
# Save raw model output
saveENMevaluation(mod, paste0("outputs/", sp, "_maxent_mod_v2.rds"))


