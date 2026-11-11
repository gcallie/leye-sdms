# Prepare used/available points, test environmental raster collinearity, and pseudooptimize covariates for lesser long-nosed bat SDMs

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
# sp <- "leye"
sp <- c("cagi", "stth", "agave")

# Read in data

## Polygons
az <- vect("data/processed/az_boundary.shp")
mcp <- vect("data/processed/mcp.shp")
mcp_buffered <- vect("data/processed/mcp_buffered.shp")

for(species in sp){
  
  ## Species occurrence points
  sp_raw <- read_csv(paste0("data/processed/", sp, "_occ.csv")) # species occurrence points
  
  ## Create multi-layer SpatRaster of environmental covariate data
  envs_paths_all <- c("data/processed/elev_mcp_buffered.tif",
                      "data/processed/slope_mcp_buffered.tif",
                      "data/processed/aspect_mcp_buffered.tif",
                      "data/processed/tri_mcp_buffered.tif",
                      "data/processed/tpi_mcp_buffered.tif",
                      "data/processed/karst_density_mcp_buffered.tif",
                      "data/processed/karst_dist_mcp_buffered.tif", 
                      "data/processed/bioclim_mcp_buffered.tif")
  
  envs_all <- rast(envs_paths_all) # Create SpatRaster object
  names(envs_all) # Check covariate list
  proj <- crs(envs_all) # Store CRS/projection
  proj # Preview object
  
  ## Plot and check to make sure raster extents look correct
  plot(envs_all)
  plot(envs_all[[1]])
  plot(mcp, add = TRUE)
  plot(mcp_buffered, add = TRUE)
  
  # Prepare data
  
  ## Standardize naming conventions to all lowercase, no underscores
  envs_all <- envs_all |> 
    rename(tri = TRI, tpi = TPI, dd0 = DD_0, eref = Eref, mcmt = MCMT, mwmt = MWMT, pptat = PPT_at, pptsp = PPT_sp, pptsm = PPT_sm, pptwt = PPT_wt, tavesm = Tave_sm, tavewt = Tave_wt, tavesp = Tave_sp, taveat = Tave_at, karstdist = karst_dist_m, karstdens = karst_density)
  plot(envs_all) # Preview object
  
  ## Create a SpatRaster of covariates masked to study area MCP
  envs_mcp <- crop(envs_all, mcp, mask = TRUE)
  plot(envs_mcp[[1]])
  envs_mcp_buffered <- crop(envs_all, mcp_buffered, mask = TRUE)
  plot(envs_mcp_buffered[[1]])
  plot(mcp, add = TRUE)
  
  ## Create SpatVector object and reproject occurrence points
  sp_used_sv <-  project(vect(sp_raw,
                              geom = c("longitude", "latitude"),
                              crs = "epsg:4326"), proj) 
  sp_used_sv <- crop(sp_used_sv, mcp)
  plot(envs_all[[1]]) # Plot one environmental raster
  plot(mcp, add = TRUE) # Overlay study area polygon
  plot(sp_used_sv, add = TRUE, pch = 1, cex = 0.75, col = "black") # Plot used points
  
  ## Convert points to data frame objects
  sp_used <- crds(sp_used_sv)
  sp_used <- as.data.frame(sp_used)
  sp_used <- sp_used |> 
    mutate(species = species) # Create species column for use in thin()  
  head(sp_used) # Preview object
  
  ## Thin used points to 1 occurrence per 1km grid  cell
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
    mutate(weight = 1) 
  
  ## Create SpatVector of thinned used points
  sp_used_thin_sv <- vect(sp_used_thin,
                          geom = c("longitude", "latitude"),
                          crs = proj)
  
  ## Create conditional multiplier object
  multiplier <- case_when(nrow(sp_used_thin) < 100 ~ 100,
                          nrow(sp_used_thin) >= 100 & nrow(sp_used_thin) < 20000 ~ 50,
                          .default = 25)
  
  ## Set seed for random number generator
  set.seed(1404)
  
  ## Generate available (i.e., background) points
  system.time(sp_avail_sv <- spatSample(mcp, 
                                        size = nrow(sp_used_thin) * multiplier,
                                        method = "random"))
  
  plot(envs_all[[1]]) # Plot one environmental raster
  plot(mcp, add = TRUE) # Overlay study area polygon
  plot(sp_avail_sv, add = TRUE, pch = 1, cex = 0.75, col = "white") # Add available points
  plot(sp_used_sv, add = TRUE, cex = 0.75, col = "red") # Add thinned used points
  
  sp_avail <- crds(sp_avail_sv)
  sp_avail <- as.data.frame(sp_avail)
  sp_avail <- sp_avail |> 
    rename("longitude" = "x", "latitude" = "y") |> 
    mutate(pa = 0) |> 
    mutate(weight = 1) # weight used points as 1 
  head(sp_avail) # Preview object
  
  ## Join used and available points into a single data frame
  sp_pts <- rbind(sp_used_thin, sp_avail)
  head(sp_pts) # Preview object
  table(sp_pts$pa) # Compare used/available point counts
  
  ## Extract covariate values at used available point locations
  covs_pts <- terra::extract(x = envs_all,
                             y = sp_pts[, c("longitude", "latitude")],
                             ID = FALSE)
  head(covs_pts) # Preview object
  
  ## Save output
  write_csv(sp_used_thin, paste0("data/processed/", species, "_used_thin.csv")) 
  write_csv(sp_avail, paste0("data/processed/", species, "_avail.csv")) 
  write_csv(sp_pts, paste0("data/processed/", species, "_pts.csv")) 
  write_csv(covs_pts, paste0("data/processed/", species, "_covs_pts.csv")) 
}

# TEST COLLINEARITY #

## VIF (across MCP, threshold = 5, force retention of mean summer and winter precipitation)
vif_test <- vifstep(envs_mcp, th = 5, keep = c("pptsm", "pptwt"))
envs <- usdm::exclude(envs_mcp_buffered, vif_test) # Get subset of non-collinear covariates for buffered study area rasters
plot(envs) # Preview object
envs_subset <- names(envs) # Check covariate list
envs_subset

## Pairwise (Pearson's correlation)
cor <- layerCor(envs_all, fun = "cor") |> 
  as.data.frame() |>  # Create correlation matrix
  select(starts_with("correlation.")) |>  # Clean up columns
  rename_with(~ str_remove(., "^correlation\\."))
head(cor)

# Summarize highly correlated (cor > 0.7) covariates 
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

# Save raster stack of forage (cagi, stth, agave) covariates
# Noncollinear minus karstdist and karstdens
envs_forage <- envs |> 
  select(-karstdist, -karstdens)
names(envs_forage)
writeRaster(envs_forage, "envs_forage.tif", overwrite = TRUE)

# PSEUDOOPTIMIZE COVARIATES # 
# Only for lesser long-nosed bats (i.e., leye)

# Extract covariate values at different scales of effect

## Bind extracted covariate values to lat/long, presence/absence, weight and keep only noncollinear covariates
sp_data <- bind_cols(covs_pts, sp_pts) |> 
  select(any_of(c(envs_subset, "pa", "weight")))
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
  cov_buffer_i <- terra::extract(envs,
                                 pt_buffers,
                                 fun = "mean")
  
  # Add buffer size to column names
  names(cov_buffer_i) <- paste0(names(cov_buffer_i), "_", buffer_i)
  
  # Add the columns from cov_buffer_i to sp_pseudo
  sp_pseudo <- cbind(sp_pseudo, cov_buffer_i) 
  
}

## Remove ID columns and rename point scale columns to _1000
sp_pseudo <- sp_pseudo |> 
  select(-c(ID_2500, ID_5000, ID_10000, ID_25000, ID_50000)) |> 
  rename(slope_1000 = slope, aspect_1000 = aspect, tpi_1000 = tpi, dd0_1000 = dd0, pptsm_1000 = pptsm, pptwt_1000 = pptwt, karstdist_1000 = karstdist, karstdens_1000 = karstdens)

summary(sp_pseudo) # View summary

# Compare and select scales of effect

## Prepare input variables
buffers <- c(1000, 2500, 5000, 10000, 25000, 50000) # Add point scale (1km) to list of buffers
buffers # Preview object

## Prepare storage data frame
col_names <- c("covariate","scale","AIC","deltaAIC","AICw") # Create list of columns
scales <- as.data.frame(matrix(NA, nrow = length(envs_subset) * length(buffers),
                               ncol = length(col_names)))
colnames(scales) <- col_names # Rename columns
i.row <- 1

##  Generate models and AIC values 
for(i in 1:length(envs_subset)){
  cov <- envs_subset[i]
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

# Pull out scale of effect for each covariate
optimized_scales <- scales |> 
  group_by(covariate) |> 
  slice_min(deltaAIC)
optimized_scales # Preview object

# Scales of effect
## aspect: 1000
## dd0: 25000
## karstdens: 50000
## pptsm: 5000
## pptwt: 5000
## slope: 5000
## tpi: 1000
# Forage: ? 




