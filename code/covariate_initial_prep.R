# Prepare covariate data for use in SDMs

# Load packages
library(tidyverse)
library(terra)

# READ IN, REPROJECT, AND CROP RASTER DATA #

# Read in and rename bioclim normals rasters
bioclim_files <- list.files("data/raw/bioclim_normals") # Create list of bioclimatic normal rasters 
subset <- c("MWMT", "MCMT", "DD_0", "Eref",
            "Tave_wt", "Tave_sp", "Tave_sm", "Tave_at",
            "PPT_wt", "PPT_sp", "PPT_sm", "PPT_at") # Define subset of bioclim variables to read in
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
writeRaster(bioclim_raw, "data/raw/bioclim_raw.tif", overwrite = TRUE) # Save output

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
#plot(karst_raw) # Plot karst polygons
karst_project <- project(karst_raw, crs(bioclim_raw)) # Reproject to same CRS as bioclim variables
karst_project # Preview object
#plot(karst_project) # Plot projected karst polygons
 
# Read in AZ boundary and create SpatVector object
az_project <- vect("data/raw/az_boundary_lambert.shp")
az_project # Preview object
plot(az_project) # Plot boundary
writeVector(az_project, "data/processed/az_boundary.shp", overwrite = TRUE) # Save output

# Create study area extent that is 12.5% larger than study area to ensure pseudooptimization buffers for occurrence points do not fall outside of rasters
study_extent <- ext(az_project) * 1.125
study_extent # Preview extent

# Crop bioclim rasters to study area
bioclim_crop <- crop(bioclim_raw, study_extent)
bioclim_crop # Preview object
plot(bioclim_crop) # Plot cropped bioclim rasters
writeRaster(bioclim_crop, "data/processed/bioclim_az.tif", overwrite = TRUE) # Save output

# Crop DEM raster
dem_crop <- crop(dem_project, study_extent) 
dem_crop # Preview object
plot(dem_crop) # Plot cropped DEM raster
plot(az_project, add = TRUE) # Overlay AZ boundary
writeRaster(dem_crop, "data/processed/elev_az.tif", overwrite = TRUE) # Save output

# Crop karst polygons
karst_crop <- crop(karst_project, az_project)
karst_crop # Preview object
plot(karst_crop) # Plot cropped karst polygons
plot(az_project, add = TRUE) # Overlay AZ boundary
writeVector(karst_crop, "data/processed/karst_az.geojson", overwrite = TRUE) # Save output


# CALCULATE ADDITIONAL COVARIATES #

# Generate covariates from DEM
## Slope
slope <- terrain(dem_crop, 
                  v = "slope", 
                  unit = "degrees",
                  neighbors = 8)
slope # Preview object
plot(slope) # Plot slope raster
writeRaster(slope, "data/processed/slope_az.tif", overwrite = TRUE) # Save output

## Aspect
aspect <- terrain(dem_crop, 
                  v = "aspect", 
                  unit = "degrees",
                  neighbors = 8)
aspect # Preview object
plot(aspect) # Plot aspect raster
writeRaster(aspect, "data/processed/aspect_az.tif", overwrite = TRUE) # Save output

## Terrain Ruggedness Index
tri <- terrain(dem_crop, 
               v = "TRI")
tri # Preview object
plot(tri) # Plot TRI raster
writeRaster(tri, "data/processed/tri_az.tif", overwrite = TRUE) # Save output

## Topographic Positioning Index
tpi <- terrain(dem_crop, 
               v = "TPI")
tpi # Preview object
plot(tpi) # Plot TRI raster
writeRaster(tpi, "data/processed/tpi_az.tif", overwrite = TRUE) # Save output


# Convert karst polygons to SpatRaster object
## Generate overlay grid
karst_grid <- rast(nrows = 685, ncols = 587, 
                   extent = ext(bioclim_crop),
                   crs = crs(bioclim_crop),
                   resolution = 1000)
karst_grid # Preview object

## Convert to density raster
karst_rast <- rasterize(karst_crop, karst_grid, 
                        field = "OBJECTID",
                        cover = TRUE)
karst_rast # Preview object
plot(karst_rast) # Plot raster
karst_rast[karst_rast == 0] <- NA # Reclassify karst raster so 0s are NAs
summary(karst_rast)
plot(karst_rast, colNA = "gray") # Replot raster 
karst_rast <- crop(karst_rast, study_extent) # Crop to AZ boundary
names(karst_rast) <- "karst_density" # Rename raster layer with covariate names
writeRaster(karst_rast, "data/processed/karst_density_az.tif", overwrite = TRUE) # Save output

## Reclassify as binary raster
karst_presence <- karst_rast # Create new object
karst_presence[karst_presence > 0] <- 1 # Convert values > 0 to 1s
karst_presence # Preview object
summary(karst_presence) # Check that # of NA values matches that of karst_rast
plot(karst_presence, colNA = "gray") # Plot raster
names(karst_presence) <- "karst_presence"
writeRaster(karst_presence, "data/processed/karst_presence_az.tif", overwrite = TRUE) # Save output

# Create distance to karst raster
karst_dist <- terra::distance(karst_presence)
karst_dist <- crop(karst_dist, study_extent) # Recrop to AZ boundary
karst_dist # Preview object
plot(karst_dist) # Replot raster
plot(az_project, add = TRUE) # Overlay AZ boundary
names(karst_dist) <- "karst_dist_m"
writeRaster(karst_dist, "data/processed/karst_dist_az.tif", overwrite = TRUE) # Save output
