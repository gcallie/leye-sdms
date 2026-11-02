# Check collinearity of raster layers, select candidate model covariates, determine appropriate scale of effect for leye SDMs

# Load packages
library(tidyverse)
library(terra)
library(usdm)
library(corrplot)

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
bioclim <- rast("data/processed/bioclim_az.tif")
leye_data <- read_csv("data/processed/leye_occ.csv")
leye_pts <- project(vect(leye_data,
                 geom = c("decimal_longitude", "decimal_latitude"),
                 crs = "epsg:4326"), crs(crs(bioclim))) # Read in leye occurrence points and reproject to Lambert Equal Area

# Preview objects
az
elev
slope
aspect
tri
tpi
karst_density
karst_dist
karst_presence
bioclim
leye_data

# Test by plotting one raster with AZ boundary and occurrence points overlaid 
plot(elev)
plot(az, add = TRUE)
plot(leye_pts, add = TRUE)

# Check that 50km buffers of points fall within the buffered AZ rasters
pt_50km <- buffer(leye_pts, width = 50000)
plot(pt_50km, add = T)


# CHECK RASTER COLLINEARITY #

# Bioclim covariates
## Calculate correlation values between bioclim layers (Pierson's correlation)
cor <- layerCor(bioclim, 
                fun = "cor",
                use = "pairwise.complete")
cor_matrix <- as.matrix(cor$correlation) # Convert correlation values to matrix
corrplot(cor_matrix, method = "number") # Plot correlation matrix

## Summarize highly correlated (cor > 0.7) bioclim variables in a table
cor_df <- as.data.frame(cor$correlation) # Create data frame object of correlation values
cor_df <- cor_df |> 
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
  select(-count, -num, -cov_order)
cor_df # Preview object

## Check VIF collinearity
vif(bioclim)

# DEM-derived covariates
## Check VIF collinearity 
vif(c(slope, aspect, tri, elev, tpi))

# Karst-derived covariates
cor2 <- layerCor(c(karst_density, karst_dist), 
                 fun = "cor",
                 use = "pairwise.complete")
cor2

# DEM- and karst-derived covariates
vif(c(slope, aspect, tri, elev, tpi, karst_density, karst_dist)) # Need to look at why karst_dist_m is coming back with VIF = NaN

