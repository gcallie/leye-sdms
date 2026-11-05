# Clean and visualize raw occurrence data downloaded via gbif_extraction.R

# Load packages
library(tidyverse)
library(terra)
library(CoordinateCleaner)
library(sf)

# List of species codes 
# leye (lesser long-nosed bats)
# cagi (saguaros)
# stth (organ pipes)
# agave (agave sp.)

# Define species code
code <- "stth"

# Read in raw GBIF download data from gbif_extraction.R
gbif_raw <- readRDS(paste0("data/raw/", code, "_occ_raw.rds"))

# PRELIMINARY DATA CLEANING

# Create sf object of AZ boundary (study area) for use in ggplots
az_sf <- st_read("data/raw/az_boundary_4326.shp")

# Create SpatVector object of study area
az_sv <- vect("data/raw/az_boundary_4326.shp")

# Plot study area polygon
plot(az_sv)

# Check for and remove data recorded based on fossil or living specimens
clean_step1 <- gbif_raw |>
  as_tibble() |>
  filter(!basisOfRecord %in% c("FOSSIL_SPECIMEN",
                               "LIVING_SPECIMEN"))
print(paste0(nrow(gbif_raw) - nrow(clean_step1), " records deleted; ",
             nrow(clean_step1), " records remaining."))

# Plot initial raw vs clean_step1 data 
# Red points have been removed from the data set
ggplot() +
  geom_sf(data = az_sf) +
  geom_point(data = gbif_raw,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "red", size = 1.5) +
  geom_point(data = clean_step1,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "black", size = 1.5) +
  theme_bw()

# Check for and remove records with problematic occurrence information using CoordinateCleaner 
clean_step2 <- clean_step1 |>
  filter(!is.na(decimalLatitude),
         !is.na(decimalLongitude)) |> # remove NA lat/long values
  cc_dupl() |> # remove distinct records with identical coordinates 
  cc_zero() |> # remove records with zeros in coordinates
  cc_equ() |> # remove records where the lat/long are equal
  cc_val() |> # remove records outside lat/long coordinate system (i.e., data entry errors)
  cc_sea() |> # remove records within oceans
  cc_cen(buffer = 2000) # remove records from around country and province centroids (i.e., imprecise georeferencing)

print(paste0(nrow(gbif_raw) - nrow(clean_step2), " records deleted; ",
             nrow(clean_step2), " records remaining."))

# Plot initial clean_step1 vs clean_step2 data
# Red points have been removed from the data set
ggplot() +
  geom_sf(data = az_sf) +
  geom_point(data = clean_step1,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "red", size = 1.5) +
  geom_point(data = clean_step2,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "black", size = 1.5, alpha = 0.25) + # Lower alpha to see when duplicates were removed
  theme_bw()

# Remove records with coordinate uncertainty >= 500m
clean_step3 <- clean_step2 |> 
  filter(is.na(coordinateUncertaintyInMeters) |
           coordinateUncertaintyInMeters < 500)

print(paste0(nrow(gbif_raw) - nrow(clean_step3), " records deleted; ",
             nrow(clean_step3), " records remaining." ))

# Plot clean_step2 vs clean_step3 data
# Red points have been removed from the data set
ggplot() +
  geom_sf(data = az_sf) +
  geom_point(data = clean_step2,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "red", size = 1.5) +
  geom_point(data = clean_step3,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "black", size = 1.5) +
  theme_bw()

# Plot initial raw data vs clean_step3 data
# Red points have been removed from the data set
ggplot() +
  geom_sf(data = az_sf) +
  geom_point(data = gbif_raw,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "red", size = 2) +
  geom_point(data = clean_step3,
             aes(x = decimalLongitude,
                 y = decimalLatitude),
             color = "black", size = 2) +
  theme_bw()

# Histogram of species occurrences plotted by year
clean_step3 |> 
  mutate(year = as.numeric(year)) |> 
  ggplot(aes(x = year)) + 
  geom_histogram() +
  stat_count(binwidth = 1, 
             geom="text", aes(label=..count..), vjust=-0.5, hjust = 0.3) +
  theme_bw()

# Convert year column to numeric (for histogram) and standardize column naming conventions
gbif_clean <- clean_step3 |>
  rename(longitude = decimalLongitude, latitude = decimalLatitude) |> 
  select(longitude, latitude) 
head(gbif_clean) # Preview object

# Save cleaned occurrence data as .rds
saveRDS(gbif_clean, paste0("data/processed/", code, "_occ.rds"))
write.csv(gbif_clean, paste0("data/processed/", code, "_occ.csv"))

# Plot the AZ boundary SpatVector object
terra::plot(az_sv)

# Build a SpatVector object of cleaned occurrence points
occ_map <- vect(gbif_clean, 
                geom = c("longitude", "latitude"),
                crs = "epsg:4326")
points(occ_map, cex = 0.7) # Overlay points over AZ map
