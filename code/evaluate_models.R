# Evaluate MaxEnt model outputs, select optimal models

# Load packages
library(terra)
library(tidyterra)
library(ENMeval)
library(ecospat)
library(tidyverse)

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

# Read in environmental rasters
if(forage == TRUE) {
  envs <- rast("data/processed/covs4foragesdms.tif")
} else if (sp == "leye") {
  envs <- rast("data/processed/covs4leyesdm_z.tif")
} else {
  warning("Species not defined.")
}

proj <- crs(envs) # Store projection
proj
plot(envs) # Plot to make sure correct covariates appear

# Read in study area polygon
az <- vect("data/processed/az_boundary.shp")

# Read in model file
mod <- loadENMevaluation(paste0("outputs/", sp, "_maxent_mod.rds"))

# Preview model
mod

# Create table of model evaluation for different tuning arguments
results <- mod@results
results

# Save model table output
write_csv(results, paste0("outputs/", sp, "_mod_results.csv"))

# Rank models based on CBI, deltaAIC, or.10p.avg
## IDEAL: 
### cbi.val.avg: closest to 1 (0.5 = moderate predictive power, above 0.7 = strong predictive power)
### or.10p.avg: closest to or below 0.1 (above 0.1 indicates overfitting)
### delta.AICc: close to 0
results_ranked <- results |> 
  relocate(delta.AICc, .before = auc.train) |> 
  relocate(or.10p.avg, .before = delta.AICc) |> 
  relocate(cbi.val.avg, .before = or.10p.avg) |> 
  filter(or.10p.avg <= 0.2) |>  # Keep models only with low omission rates (to penalize overfitting)
  filter(cbi.val.avg >= 0.7) |> # Keep models with moderate to strong predictive power |> 
  arrange(delta.AICc) # |> # Show lowest delta.AICc values at the top
  #arrange(desc(cbi.val.avg)) #|> # Show highest cbi.val.avg values at the top
  #arrange(or.10p.avg) # Show lowest remaining or.10p.avg values at the top

results_ranked # Preview object

# Save ranked model table output
write_csv(results_ranked, paste0("outputs/", sp, "_mod_results_ranked.csv"))

opt.tune.args <- results_ranked[1, "tune.args"] |> 
  as.character()# Pull out the tune.args parameter for the best model from the above workflow
opt.tune.args # Preview

# Select the optimal model from the ENMevaluation object
opt.mod <- eval.models(mod)[[opt.tune.args]]

# Store beta values for optimal model in a data frame with the model tune.args
betas <- opt.mod$betas |> 
  as.data.frame() |> 
  rename(beta = `opt.mod$betas`) |> 
  rownames_to_column(var = "covariate")

betas # Preview object

# Save beta table output
write_csv(betas, paste0("outputs/", sp, "_betas_", opt.tune.args, ".csv"))

# Plot response curves
data <- eval.occs(mod)[, names(envs)]
evalplot.respCurves(opt.mod, data)
ggsave(paste0("outputs/", sp, "_response_curves_", opt.tune.args, ".png"), width = 8, height = 7, dpi = 300)

# Map prediction based on optimal model selected above
pred.seq <- eval.predictions(mod)[[as.character(opt.tune.args)]]
plot(pred.seq, col = terrain.colors(100))

# Pull SpatRaster for best model
best_rast <- mod@predictions[[as.character(opt.tune.args)]]
crs(best_rast) <- proj
names(best_rast) <- sp

if(forage == FALSE) {
  best_rast <- crop(best_rast, az, mask = TRUE)
  best_rast
  plot(best_rast)
} else {
  best_rast
  plot(best_rast)
}

# Save output
writeRaster(best_rast, paste0("outputs/", sp, "_selection_surface.tif"), overwrite = TRUE)
