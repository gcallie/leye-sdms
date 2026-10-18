# Prepare covariate data for use in SDMs

# Load packages
library(tidyverse)
library(terra)
library(usdm)


# READ IN, REPROJECT, AND CROP GIS DATA #

# Read in and rename bioclim normals rasters
bioclim_files <- list.files("data/raw/bioclim_normals") # Create list of bioclimatic normal rasters 
subset <- c("DD_0", "PPT_wt", "PPT_sm", "Eref") # Define subset of bioclim variables to read in
subset_regex <- str_c(subset, collapse = "|") # Collapse into regex
bioclim_subset <- str_subset(bioclim_files, subset_regex)
bioclim_subset # Preview result
bioclim_labels <- str_extract(bioclim_subset, "(?<=2020_).*?(?=\\.)") # Pull covariate names from list of file name subset
bioclim_labels # Preview result

# Create multi-layer SpatRaster object for bioclim variables
bioclim_raw <- rast(paste0("data/raw/bioclim_normals/", bioclim_subset), 
                    bioclim_labels)
bioclim_raw # Preview object
names(bioclim_raw) <- bioclim_labels # Rename raster layers with covariate names
bioclim_raw # Check that names are correct
writeRaster(bioclim_raw, "data/raw/bioclim_raw.tif") # Save output

# Read in DEM raster and create SpatRaster object
dem_raw <- rast("data/raw/elevation.tif")
dem_raw # Preview object
plot(dem_raw) # Plot DEM
dem_project <- project(dem_raw, crs(bioclim_raw)) # Reproject to same CRS as bioclim variables
dem_project # Preview object
plot(dem_project) # Plot projected DEM

# Read in karst and pseudokarst polygons and create SpatVector object
karst_raw <- vect("data/raw/karst.geojson")
karst_raw # Preview object
plot(karst_raw) # Plot karst polygons
karst_project <- project(karst_raw, crs(bioclim_raw)) # Reproject to same CRS as bioclim variables
karst_project # Preview object
plot(karst_project) # Plot projected karst polygons
 
# Read in AZ boundary and create SpatVector object
az_raw <- vect("data/raw/az_boundary_raw.geojson")
az_raw # Preview object
plot(az_raw) # Plot boundary
az_project <- project(az_raw, crs(bioclim_raw)) # Reproject to same CRS as bioclim variables
az_project # Preview object
plot(az_project) # Plot reprojected boundary
#writeVector(karst_crop, "data/final/az_boundary.geojson") # Save output

# Crop bioclim rasters to study area
bioclim_crop <- crop(bioclim_raw, az_project, mask = TRUE)
bioclim_crop # Preview object
plot(bioclim_crop) # Plot cropped bioclim rasters
#writeRaster(bioclim_crop, "data/final/bioclim_az.tif") # Save output

# Crop DEM raster
dem_crop <- crop(dem_project, az_project, mask = TRUE) 
dem_crop # Preview object
plot(dem_crop) # Plot cropped DEM raster
plot(az_project, add = TRUE) # Overlay AZ boundary
#writeRaster(dem_crop, "data/final/dem_az.tif") # Save output

# Crop karst polygons
karst_crop <- crop(karst_project, az_project)
karst_crop # Preview object
plot(karst_crop) # Plot cropped karst polygons
plot(az_project, add = TRUE) # Overlay AZ boundary
#writeVector(karst_crop, "data/final/karst_az.geojson") # Save output


# CALCULATE ADDITIONAL COVARIATES #

# Generate covariates from DEM
## Slope
slope <- terrain(dem_crop, 
                  v = "slope", 
                  unit = "degrees",
                  neighbors = 8)
slope # Preview object
plot(slope) # Plot slope raster
writeRaster(slope, "data/final/slope_az.tif") # Save output

## Aspect
aspect <- terrain(dem_crop, 
                  v = "aspect", 
                  unit = "degrees",
                  neighbors = 8)
aspect # Preview object
plot(aspect) # Plot aspect raster
writeRaster(slope, "data/final/aspect_az.tif") # Save output

## Terrain Ruggedness Index
tri <- terrain(dem_crop, 
               v = "TRI")
tri # Preview object
plot(tri) # Plot TRI raster
writeRaster(slope, "data/final/tri_az.tif") # Save output

# Wrangle karst data

# CHECK FOR COLLINEARITY #

# DEM-derived covariates
## Check VIF collinearity 
vif(c(slope, aspect, tri))

# Bioclim covariates
## Calculate correlation values between bioclim layers (Pierson's correlation)
cor <- layerCor(bioclim_crop, 
                fun = "cor",
                use = "pairwise.complete")
cor_matrix <- as.matrix(cor$correlation) # Convert correlation values to matrix
corrplot(cor_matrix, method = "number") # Plot correlation matrix

## Column names to input into first argument of pivot_longer
length <- length(subset)
subset[1]
subset[length]

## Summarize bioclim correlation values in a table
cor_df <- as.data.frame(cor$correlation) # Create data frame object of correlation values
cor_df <- cor_df |> 
  rownames_to_column("cov1") |> 
  pivot_longer(DD_0:Tave_sm, # Update this based on code above
               names_to = "cov2",
               values_to = "cor") |> 
  filter(cov1 != cov2) |> 
  #filter(abs(cor) <= 0.6) |> 
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
usdm::vif(bioclim_crop)

