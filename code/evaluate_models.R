# Evaluate MaxEnt model outputs, select optimal models

# Load packages
library(terra)
library(tidyterra)
library(ENMeval)
library(ecospat)
library(tidyverse)

# Define species
sp <- "agave"

# Read in environmental rasters
envs <- rast("data/processed/envs_forage.tif")
proj <- crs(envs) # Store projection
proj

# Read in study area polygon
az <- vect("data/processed/az_boundary.shp")

# Read in model file
mod <- loadENMevaluation(paste0("outputs/", sp, "_maxent_mod.rds"))

# Preview model
mod

# Pull results into a data frame
results <- eval.results(mod) 

# Pull the model with the lowest AICc
opt.aicc <- results |> filter(delta.AICc == 0) 
opt.aicc

# Look at beta values 
mod.seq <- eval.models(mod)[[opt.aicc$tune.args]]
mod.seq$betas

# Plot response curves
data <- eval.occs(mod)[, names(envs)]
evalplot.respCurves(mod.seq, data)
ggsave(paste0("outputs/", sp, "_response_curves.png"))

# Map prediction based on model with lowest AICc
pred.seq <- eval.predictions(mod)[[as.character(opt.aicc$tune.args)]]
plot(pred.seq, col = terrain.colors(100))

# Pull SpatRaster for best model
best_rast <- mod@predictions[[as.character(opt.aicc$tune.args)]]
crs(best_rast) <- proj
names(best_rast) <- sp
best_rast
plot(best_rast)
plot(az, add = TRUE)
ggsave(paste0("outputs/", sp, "_selection_surface.png"))

# Save output
writeRaster(best_rast, paste0("outputs/", sp, "_selection_surface.tif"), overwrite = TRUE)

