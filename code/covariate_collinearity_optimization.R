# Prepare used/available points, test environmental raster collinearity, export raster layers for use in forage species SDMs, and pseudooptimize covariates for lesser long-nosed bat SDMs

# Install spThin from CRAN read-only GitHub mirror
# remotes::install_github("cran/package_name")

# Load packages
library(tidyterra)
library(tidyverse)
library(terra)
library(usdm)
library(MuMIn)
library(spThin)

# Define species
sp <- "leye"
#sp <- c("cagi", "stth", "agave")

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
sp_used <- as.data.frame(sp_used)
sp_used <- sp_used |> 
    mutate(species = sp) # Create species column for use in thin()  
head(sp_used) # Preview object
  
## Thin used points to 1 occurrence per 2.5km for leye and 1 occurrence per 1km for forage species
if(sp == "leye") {
  sp_used_thin <- spThin::thin(sp_used,
                               lat.col = "y",
                               long.col = "x",
                               spec.col = "species",
                               thin.par = 2.5, 
                               reps = 1,
                               locs.thinned.list.return = TRUE, 
                               write.files = FALSE, 
                               write.log.file = FALSE) |> 
    as.data.frame() |> 
    rename("longitude" = "Longitude", "latitude" = "Latitude") |> 
    mutate(pa = 1) |>  # Populate presence/absence column 
    mutate(weight = 1)# Weight used points as 1
  } else {
      sp_used_thin <- spThin::thin(sp_used,
                                   lat.col = "y",
                                   long.col = "x",
                                   spec.col = "species",
                                   thin.par = 1, 
                                   reps = 1,
                                   locs.thinned.list.return = TRUE, 
                                   write.files = FALSE, 
                                   write.log.file = FALSE) |> 
        as.data.frame() |> 
        rename("longitude" = "Longitude", "latitude" = "Latitude") |> 
        mutate(pa = 1) |>  # Populate presence/absence column 
        mutate(weight = 1)# Weight used points as 1
  }
  
## Create SpatVector of thinned used points
sp_used_thin_sv <- vect(sp_used_thin,
                          geom = c("longitude", "latitude"),
                          crs = proj)
  
## Set seed for random number generator
set.seed(1404)
  
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
sp_pts <- rbind(sp_used_thin, sp_avail)
head(sp_pts) # Preview object
table(sp_pts$pa) # Compare used/available point counts
  
## Extract covariate values at used available point locations
envs_pts <- terra::extract(x = envs_all,
                             y = sp_pts[, c("longitude", "latitude")],
                             ID = FALSE)
head(envs_pts) # Preview object
  
## Save output
write_csv(sp_used_thin, paste0("data/processed/", sp, "_used_thin.csv")) 
write_csv(sp_avail, paste0("data/processed/", sp, "_avail.csv")) 
write_csv(sp_pts, paste0("data/processed/", sp, "_pts.csv")) 
 

### TEST ENVIRONMENTAL RASTER COLLINEARITY ###

# VIF (across AZ, threshold = 5, force retention of mean summer and winter precipitation)
vif_test <- vifstep(envs_az, th = 5, keep = c("pptsm", "pptwt"))
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

# Compare binomial GLMs with select covariates

# Save raster stack of environmental covariates for use in running forage species SDMS (non-collinear covariates minus karstdist and karstdens)
envs_forage_sdms <- envs |> 
  select(-karstdist, -karstdens)
names(envs_forage_sdms)
writeRaster(envs_forage_sdms, "data/processed/covs4foragesdms.tif", overwrite = TRUE)


### PSEUDOOPTIMIZE COVARIATES ###
# Only for lesser long-nosed bats (i.e., leye)

# Read in selection surfaces from forage SDMs generated in build_models.R and evaluate_models.R
forage_results_paths <- c("outputs/cagi_selection_surface.tif",
                       "outputs/stth_selection_surface.tif",
                       "outputs/agave_selection_surface.tif")
forage_results <- rast(forage_results_paths) # Create SpatRaster object
plot(forage_results)

# Create composite raster
forage_sum <- sum(forage_results, na.rm = TRUE) # Sum selection surfaces for cagi, stth, agave
plot(forage_sum)
names(forage_sum) <- "forage" # Rename layer

forage_results <- c(forage_results, forage_sum) # Combine individual species and summed selection surface into single SpatRaster 
plot(forage_results)
writeRaster(forage_results, "data/processed/forage_sdm_results.tif", overwrite = TRUE) # Save output

# Extract forage SDM values at leye used/available points
forage_results_pts <- terra::extract(x = forage_results,
                                  y = sp_pts[, c("longitude", "latitude")],
                                  ID = FALSE)
head(forage_results_pts) # Preview object
summary(forage_results_pts) # View summary

# Combine environmental covariates into multi-layer raster with forage SDM result layers

## Check that extents and resolutions match
ext(envs) == ext(forage_results) 
res(envs) == res(forage_results)

## Deal with minor extent/alignment issue
envs_align <- crop(envs, forage_results)
envs_align <- extend(envs_align, forage_results)

## Recheck extents
ext(envs_align) == ext(forage_results)

## Combine into multi-layer SpatRaster
leye_covs <- c(envs_align, forage_results)
leye_covs # Preview object
plot(leye_covs)
leye_covs_list <- names(leye_covs)

# Extract covariate values at different scales of effect

## Bind extracted covariate values to lat/long, presence/absence, weight, and keep only noncollinear covariates
envs_pts <- envs_pts |> 
  select(all_of(envs_subset))
head(envs_pts)

covs_pts <- bind_cols(envs_pts, forage_results_pts)
sp_data <- bind_cols(covs_pts, sp_pts) |> 
  select(-latitude, -longitude)
head(sp_data)

## Create SpatVector object of used/available points
sp_pts_sv <- vect(sp_pts,
                  geom = c("longitude", "latitude"),
                  crs = proj)
sp_pts_sv # Preview object

## Create copy of sp_data data frame for storing new values extracted at buffers
sp_pseudo <- sp_data
head(sp_pseudo)

## Create vector object of buffer radii (2.5km, 5km, 10km, 25km, 50km)
buffers <- c(2500, 5000, 10000, 25000, 50000) 

## Extract mean values for each covariate at different scales of effect and store in sp_pseudo object 

for(i in 1:length(buffers)){
  
  buffer_i <- buffers[i]
  
  cat("Running buffer size =", buffer_i,"\n")
  
  pt_buffers <- buffer(sp_pts_sv,
                       width = buffer_i)
  
  # Extract values for all covariates
  cov_buffer_i <- terra::extract(leye_covs,
                                 pt_buffers,
                                 fun = "mean")
  
  # Add buffer size to column names
  names(cov_buffer_i) <- paste0(names(cov_buffer_i), "_", buffer_i)
  
  # Add the columns from cov_buffer_i to sp_pseudo
  sp_pseudo <- cbind(sp_pseudo, cov_buffer_i) 
  
}

head(sp_pseudo) # Preview output

## Remove ID columns and rename point scale columns to _1000
sp_pseudo <- sp_pseudo |> 
  select(-c(ID_2500, ID_5000, ID_10000, ID_25000, ID_50000)) |> 
  rename(slope_1000 = slope, aspect_1000 = aspect, tpi_1000 = tpi, tri_1000 = tri, dd0_1000 = dd0, pptsm_1000 = pptsm, pptwt_1000 = pptwt, karstdist_1000 = karstdist, karstdens_1000 = karstdens, forage_1000 = forage, cagi_1000 = cagi, stth_1000 = stth, agave_1000 = agave)

names(sp_pseudo) # Check that correct covariates/scales were evaluated
summary(sp_pseudo) # View summary

# Compare and select scales of effect

## Prepare input variables
buffers <- c(1000, 2500, 5000, 10000, 25000, 50000) # Add point scale (1km) to list of buffers
buffers # Preview object

## Prepare storage data frame
col_names <- c("covariate","scale","AIC","deltaAIC","AICw") # Create list of columns
scales <- as.data.frame(matrix(NA, nrow = ncol(sp_pseudo) * length(buffers),
                               ncol = length(col_names)))
colnames(scales) <- col_names # Rename columns
head(scales)
i.row <- 1

##  Generate models and AIC values 
for(i in 1:length(leye_covs_list)){
  cov <- leye_covs_list[i]
  cat("Starting", cov, "\n")
  data <- sp_pseudo[, c("pa", "weight", paste0(cov,"_", buffers))]
  cov_cols <- paste0(cov, "_", buffers)
  
  for(j in cov_cols){
    temp_col <- j
    temp_data <- data[,c("pa", "weight", temp_col)]
    temp_model <- glm(temp_data[, 1] ~ temp_data[, 3], weights = temp_data[, 2], family = binomial)
    scales$covariate[i.row] <- cov
    scales$AIC[i.row] <- AIC(temp_model)
    scales$scale[i.row] <- as.numeric(strsplit(j, "_")[[1]][2])
    cat("Finished scale", as.numeric(strsplit(j, "_")[[1]][2], "\n"))
    i.row <- i.row + 1
  }
  min <- min(scales$AIC[which(scales$covariate == cov)])
  delta <- scales$AIC[which(scales$covariate == cov)] - min
  w <- Weights(scales$AIC[which(scales$covariate == cov)])
  scales$AICw[which(scales$covariate == cov)] <- w
  scales$deltaAIC[which(scales$covariate == cov)] <- delta
  
}

# Remove distance to karst covariate since it does not need to be pseudooptimized
scales <- scales |> 
  filter(!str_detect(covariate, "karstdist"))
scales

scale_plot <- ggplot(data = scales,
                     aes(x = scale, y = deltaAIC)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = ''), color = "black") +
  geom_point(color = "black")+
  facet_wrap(~ covariate)+
  theme_bw()
scale_plot

ggsave("outputs/leye_scales_of_effect.png") # Save output

# Pull out scale of effect for each covariate
optimized_scales <- scales |> 
  group_by(covariate) |> 
  slice_min(deltaAIC)
optimized_scales # Preview object

# Scales of effect
## agave: 5000
## aspect: 1000
## cagi: 25000
## dd0: 50000
## forage: 5000
## karstdens: 1000
## pptsm: 5000
## pptwt: 5000
## slope: 5000
## stth: 1000
## tpi: 1000
## tri: 5000

# Smooth and export scale-optimized rasters 

## Create focal matrices at different scales (5km, 25km, 50km)
buffer_5000 <- focalMat(leye_covs,
                        d = 5000,
                        type = c('circle'))

buffer_25000 <- focalMat(leye_covs,
                        d = 25000,
                        type = c('circle'))


buffer_50000 <- focalMat(leye_covs,
                         d = 50000,
                         type = c('circle'))

# NEED TO UPDATE THIS SECTION #

## Smooth pptsm, pptwt, slope, agave, forage, tri rasters to 5km scale

### pptsm
envs$pptsm_5000 <- focal(envs[["pptsm"]],
                        w = buffer_5000,
                        fun = "mean",
                        na.rm = TRUE)
### pptwt
envs$pptwt_5000 <- focal(envs[["pptwt"]],
                         w = buffer_5000,
                         fun = "mean",
                         na.rm = TRUE)

### slope
envs$slope_5000 <- focal(envs[["slope"]],
                         w = buffer_5000,
                         fun = "mean",
                         na.rm = TRUE)

## Plot to compare original rasters vs smoothed rasters
plot(c(envs[["pptsm_5000"]], envs[["pptsm"]],
       envs[["pptwt_5000"]], envs[["pptwt"]],
       envs[["slope_5000"]], envs[["slope"]]))

# Smooth cagi raster to 25km scale

# Plot to compare original vs smoothed raster

# Smooth dd0 raster to 50km scale
envs$dd0_50000 <- focal(envs[["dd0"]],
                        w = buffer_50000,
                        fun = "mean",
                        na.rm = TRUE)

# Plot to compare original vs smoothed raster
plot(c(envs[["dd0_50000"]], envs[["dd0"]]))

# Save output
envs_leye <- envs |> 
  select(-dd0, -pptsm, -pptwt, -slope)
names(envs_leye) # Check that correct covariates are included at their scales of effect
writeRaster(envs_leye, "data/processed/envs_leye.tif", overwrite = TRUE) # Save output




