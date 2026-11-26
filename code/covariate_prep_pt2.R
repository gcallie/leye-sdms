# Prepare used/available points, test environmental raster collinearity, export raster layers for use in forage species SDMs, and pseudooptimize covariates for lesser long-nosed bat SDMs

# Load packages
library(tidyterra)
library(tidyverse)
library(terra)
library(usdm)
library(MuMIn)
library(spThin)

# Define species
sp <- "leye"
#sp <- "cagi"
#sp <- "stth"
#sp <- "agave"

### PREPARE USED/AVAILABLE POINTS ###

# Read in data
  
## Study Area Polygons
az <- vect("data/processed/az_boundary.shp")
az_buffered <- vect("data/processed/az_buffered.shp")
  
## Species occurrence points
sp_raw <- read_csv(paste0("data/processed/", sp, "_occ.csv")) # species occurrence points

## Create multi-layer SpatRaster of environmental covariate data
envs_paths_all <- c("data/processed/elev_az_buffered.tif",
                      "data/processed/slope_az_buffered.tif",
                      "data/processed/aspect_az_buffered.tif",
                      "data/processed/tri_az_buffered.tif",
                      "data/processed/tpi_az_buffered.tif",
                      "data/processed/karst_density_az_buffered.tif",
                      "data/processed/karst_dist_az_buffered.tif", 
                      "data/processed/bioclim_az_buffered.tif")
  
envs_all <- rast(envs_paths_all) # Create SpatRaster object
names(envs_all) # Check covariate list
proj <- crs(envs_all) # Store CRS/projection
proj # Preview object
  
## Plot and check to make sure raster extents look correct
plot(envs_all)
plot(envs_all[[1]])
plot(az, add = TRUE)
plot(az_buffered, add = TRUE)
  
# Prepare data
  
## Standardize naming conventions to all lowercase, no underscores
envs_all <- envs_all |> 
    rename(tri = TRI, tpi = TPI, dd0 = DD_0, eref = Eref, mcmt = MCMT, mwmt = MWMT, pptat = PPT_at, pptsp = PPT_sp, pptsm = PPT_sm, pptwt = PPT_wt, tavesm = Tave_sm, tavewt = Tave_wt, tavesp = Tave_sp, taveat = Tave_at, karstdist = karst_dist_m, karstdens = karst_density)
  plot(envs_all) # Preview object
  
## Create a SpatRaster of covariates masked to study area
envs_az <- crop(envs_all, az, mask = TRUE)
plot(envs_az[[1]])
envs_az_buffered <- crop(envs_all, az_buffered, mask = TRUE)
plot(envs_az_buffered[[1]])
plot(az, add = TRUE)
  
## Create SpatVector object and reproject occurrence points
sp_used_sv <-  project(vect(sp_raw,
                              geom = c("longitude", "latitude"),
                              crs = "epsg:4326"), proj) 
sp_used_sv <- crop(sp_used_sv, az)
plot(envs_all[[1]]) # Plot one environmental raster
plot(az, add = TRUE) # Overlay study area polygon
plot(sp_used_sv, add = TRUE, pch = 1, cex = 0.75, col = "black") # Plot used points
  
## Convert points to data frame objects
sp_used <- crds(sp_used_sv)
sp_used <- as.data.frame(sp_used) |> 
  rename("longitude" = "x", "latitude" = "y") |> 
           mutate(pa = 1) |>  # Populate presence/absence column 
           mutate(weight = 1)# Weight used points as 1
head(sp_used) # Preview object
  
## Set seed for random number generator
set.seed(1404)

## Read in background points (if already generated)
#sp_avail <- read.csv(paste0("data/processed/", sp, "_avail.csv"))

## Generate available (i.e., background) points
system.time(sp_avail_sv <- spatSample(az,
                                        size = 100000,
                                        method = "random"))

plot(envs_all[[1]]) # Plot one environmental raster
plot(az, add = TRUE) # Overlay study area polygon
plot(sp_avail_sv, add = TRUE, pch = 1, cex = 0.75, col = "white") # Add available points
plot(sp_used_sv, add = TRUE, cex = 0.75, col = "red") # Add thinned used points

sp_avail <- crds(sp_avail_sv)
sp_avail <- as.data.frame(sp_avail)
sp_avail <- sp_avail |>
    rename("longitude" = "x", "latitude" = "y") |>
    mutate(pa = 0) |>
    mutate(weight = 5000) # weight available points as 5000
head(sp_avail) # Preview object
  
## Join used and available points into a single data frame
sp_pts <- rbind(sp_used, sp_avail)
head(sp_pts) # Preview object
table(sp_pts$pa) # Compare used/available point counts
  
## Extract covariate values at used available point locations
envs_pts <- terra::extract(x = envs_all,
                             y = sp_pts[, c("longitude", "latitude")],
                             ID = FALSE)
head(envs_pts) # Preview object
  
## Save output
write_csv(sp_used, paste0("data/processed/", sp, "_used.csv")) 
write_csv(sp_avail, paste0("data/processed/", sp, "_avail.csv")) 
write_csv(sp_pts, paste0("data/processed/", sp, "_pts.csv")) 

### TEST ENVIRONMENTAL RASTER COLLINEARITY ###

# VIF (across AZ, threshold = 5, force retention of mean summer and winter precipitation)
vif_test <- vifstep(envs_az, th = 5, keep = c("pptsm", "pptwt"))
vif_test
envs <- usdm::exclude(envs_az_buffered, vif_test) # Get subset of non-collinear covariates for buffered study area rasters
plot(envs) # Preview object
envs_subset <- names(envs) # Check covariate list
envs_subset

# Pairwise (Pearson's correlation)
cor <- layerCor(envs_all, fun = "cor") |> 
  as.data.frame() |>  # Create correlation matrix
  select(starts_with("correlation.")) |>  # Clean up columns
  rename_with(~ str_remove(., "^correlation\\."))
head(cor)

## Summarize highly correlated (cor > 0.7) covariates 
cor_df <- cor |> 
  rownames_to_column("cov1") |> 
  pivot_longer(elevation:tavewt, # Update this based on covariates included
               names_to = "cov2",
               values_to = "cor") |> 
  filter(cov1 != cov2) |> 
  mutate(cor = abs(cor)) |> 
  filter(cor >= 0.7) |> 
  mutate(cov_order = paste(cov1, cov2) |> # Start workflow to remove duplicates
           strsplit(split = ' ') |> 
           map_chr( ~ sort(.x) |>  
                      paste(collapse = ' '))) |> 
  mutate(count = 1) |> 
  group_by(cov_order) |> 
  mutate(num = row_number()) |> 
  filter(num != 2) |> 
  ungroup() |> 
  select(cov1, cov2, cor) |> 
  arrange(desc(cor))
cor_df # Preview object

# Save raster stack of non-collinear environmental covariates
names(envs) # Check list
writeRaster(envs, "data/processed/envs_noncollinear.tif", overwrite = TRUE)
# Save raster stack of environmental covariates for use in running forage species SDMS (non-collinear covariates minus karstdist and karstdens)
envs_forage_sdms <- envs |> 
  select(-karstdist, -karstdens)
names(envs_forage_sdms)
writeRaster(envs_forage_sdms, "data/processed/covs4foragesdms.tif", overwrite = TRUE)
