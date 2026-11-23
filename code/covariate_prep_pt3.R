# Pseudooptimize leye covariates, retest collinearity, and prepare multi-scale raster layers for use in leye SDM
## This script relies on the forage species SDMs being generated as well as outputs from covariate_prep_pt1.R and covariate_prep_pt2.R

# Load packages
library(terra)
library(tidyterra)
library(usdm)
library(tidyverse)

# Read in processed data 
sp_pts <- read_csv("data/processed/leye_pts.csv")
envs <- rast("data/processed/envs_noncollinear.tif")
az <- vect("data/processed/az_boundary.shp")

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

## Extract environmental covariate values at used available point locations
envs_pts <- terra::extract(x = envs,
                           y = sp_pts[, c("longitude", "latitude")],
                           ID = FALSE)
head(envs_pts) # Preview object

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

## Point scale (1km): aspect, karstdens, stth, tpi, karstdist

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

## Smooth pptsm, pptwt, slope, agave, forage, tri rasters to 5km scale

### pptsm
leye_covs$pptsm_5000 <- focal(leye_covs[["pptsm"]],
                         w = buffer_5000,
                         fun = "mean",
                         na.rm = TRUE)
### pptwt
leye_covs$pptwt_5000 <- focal(leye_covs[["pptwt"]],
                         w = buffer_5000,
                         fun = "mean",
                         na.rm = TRUE)

### slope
leye_covs$slope_5000 <- focal(leye_covs[["slope"]],
                         w = buffer_5000,
                         fun = "mean",
                         na.rm = TRUE)

### agave
leye_covs$agave_5000 <- focal(leye_covs[["agave"]],
                         w = buffer_5000,
                         fun = "mean",
                         na.rm = TRUE)

### forage
leye_covs$forage_5000 <- focal(leye_covs[["forage"]],
                         w = buffer_5000,
                         fun = "mean",
                         na.rm = TRUE)

### tri
leye_covs$tri_5000 <- focal(leye_covs[["tri"]],
                               w = buffer_5000,
                               fun = "mean",
                               na.rm = TRUE)

## Plot to compare original rasters vs smoothed rasters
plot(c(leye_covs[["pptsm_5000"]], leye_covs[["pptsm"]],
       leye_covs[["pptwt_5000"]], leye_covs[["pptwt"]],
       leye_covs[["slope_5000"]], leye_covs[["slope"]],
       leye_covs[["agave_5000"]], leye_covs[["agave"]],
       leye_covs[["forage_5000"]], leye_covs[["forage"]],
       leye_covs[["tri_5000"]], leye_covs[["tri"]]))

## Smooth cagi raster to 25km scale

### cagi
leye_covs$cagi_25000 <- focal(leye_covs[["cagi"]],
                            w = buffer_25000,
                            fun = "mean",
                            na.rm = TRUE)

# Plot to compare original vs smoothed raster
plot(c(leye_covs[["cagi_25000"]], leye_covs[["cagi"]]))

## Smooth dd0 raster to 50km scale

### dd0
leye_covs$dd0_50000 <- focal(leye_covs[["dd0"]],
                        w = buffer_50000,
                        fun = "mean",
                        na.rm = TRUE)

# Plot to compare original vs smoothed raster
plot(c(leye_covs[["dd0_50000"]], leye_covs[["dd0"]]))

# Remove extra covariates at point scales
covs4leyesdm <- leye_covs |> 
  select(-pptsm, -pptwt, -slope, -agave, -forage, -tri, -cagi, -dd0)
names(covs4leyesdm) # Check that correct covariates are included at their scales of effect
covs4leyesdm # Preview object
covs4leyesdm <- crop(covs4leyesdm, az, mask = TRUE) # Crop to AZ boundary (formerly were cropped to buffered AZ boundary)
plot(covs4leyesdm)

# Recheck raster collinearity
vif(covs4leyesdm)
vifstep(covs4leyesdm, th = 5, keep = c("stth", "cagi_25000", "agave_5000", "pptsm_5000"))

# Compare binomial GLM with individual forage species vs composite forage raster 
## Extract point values for covariates at respective scales of effect
scaled_pts <- terra::extract(x = covs4leyesdm,
                                     y = sp_pts[, c("longitude", "latitude")],
                                     ID = FALSE)
head(scaled_pts)

## Add covariate point values to data frame with LEYE lat/long, presence/absence
leye_data <- bind_cols(scaled_pts, sp_pts)
head(leye_data) # Preview object

## Model with individual rasters for cagi, stth, agave
indiv_mod <- glm(pa ~ aspect + tpi + karstdens + karstdist + stth + pptsm_5000 + pptwt_5000 + slope_5000 + agave_5000 + tri_5000 + cagi_25000 + dd0_50000, data = leye_data,
                weight = weight,
                family = "binomial")
summary(indiv_mod)

## Model with composite forage raster
comp_mod <- glm(pa ~ aspect + tpi + karstdens + karstdist + pptsm_5000 + pptwt_5000 + slope_5000 + tri_5000 + forage_5000 + dd0_50000, data = leye_data,
                 weight = weight,
                 family = "binomial")
summary(comp_mod)

## Look at AIC values
aic <- data.frame(
  mod = c("indiv_mod", "comp_mod"),
  AIC = c(AIC(indiv_mod), AIC(comp_mod)))
aic <-  aic |> 
  arrange(AIC) |>  
  mutate(dAIC = AIC - min(AIC)) # Calculate delta AIC
aic # Preview object

### Composite forage raster has a lower delta AIC (~1.87 units lower)
### However, I think the individual species retains more granularity about bat behavior/habitat selection

### Remove composite forage raster
covs4leyesdm <- covs4leyesdm |> 
  select(-forage_5000)
names(covs4leyesdm)

# Retest collinearity
vif(covs4leyesdm) # Slope and tri also appear very collinear

# Compare binomial GLM with tri vs slope

## Model with tri and individual forage rasters (no slope or composite forage)
tri_mod <- glm(pa ~ aspect + tpi + karstdens + karstdist + stth + pptsm_5000 + pptwt_5000 + agave_5000 + tri_5000 + cagi_25000 + dd0_50000, data = leye_data,
                 weight = weight,
                 family = "binomial")
summary(tri_mod)

## Model with slope and individual forage rasters (no tri or composite forage)
slope_mod <- glm(pa ~ aspect + tpi + karstdens + karstdist + stth + pptsm_5000 + pptwt_5000 + agave_5000 + slope_5000 + cagi_25000 + dd0_50000, data = leye_data,
               weight = weight,
               family = "binomial")
summary(slope_mod)

## Look at AIC values
aic <- data.frame(
  mod = c("tri_mod", "slope_mod"),
  AIC = c(AIC(tri_mod), AIC(slope_mod)))
aic <-  aic |> 
  arrange(AIC) |>  
  mutate(dAIC = AIC - min(AIC)) # Calculate delta AIC
aic # Preview object

### Slope mod has a lower delta AIC (~6.48 units lower)

# Remove composite forage raster and tri from leye covariate SpatRaster
covs4leyesdm <- covs4leyesdm |> 
  select(-tri_5000)
  names(covs4leyesdm)
  
# Reextract covariate values at points 
scaled_pts <- terra::extract(x = covs4leyesdm,
                               y = sp_pts[, c("longitude", "latitude")],
                               ID = FALSE)
head(scaled_pts)
  
# Add covariate point values to data frame with LEYE lat/long, presence/absence
leye_data <- bind_cols(scaled_pts, sp_pts)
head(leye_data) # Preview object

# Recheck collinearity

## VIF
vif(covs4leyesdm) # Agave and pptsm_5000 appear to be semi-collinear (VIF ~ 6.5)

## Pairwise (Pearson's correlation)
cor <- layerCor(covs4leyesdm, fun = "cor") |> 
  as.data.frame() |>  # Create correlation matrix
  select(starts_with("correlation.")) |>  # Clean up columns
  rename_with(~ str_remove(., "^correlation\\."))
head(cor)

## Summarize highly correlated (cor > 0.7) covariates 
cor_df <- cor |> 
  rownames_to_column("cov1") |> 
  pivot_longer(aspect:dd0_50000, # Update this based on covariates included
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

# agave and pptsm (0.833) and karstdens and dd0_50000 (0.770) are correlated at r > 0.7; however, I am going to retain both covariates in the model as all of my VIF values are close to my 5 threshold(max = 5.87; karstdens and dd0 < 4)

# Plot final covariate raster stack
plot(covs4leyesdm)
names(covs4leyesdm) # Check list of names
summary(covs4leyesdm) # Note: NA values are the background surrounding AZ

# Z-score standardize covariates
covs4leyesdm_z <- scale(covs4leyesdm)
global(covs4leyesdm_z, c("mean", "sd"), na.rm = TRUE) # Verify that means are 0 and SD are 1
plot(covs4leyesdm_z)
summary(covs4leyesdm_z)

# Save output
writeRaster(covs4leyesdm_z, "data/processed/covs4leyesdm_z.tif", overwrite = TRUE)


