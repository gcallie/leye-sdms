# Extract species occurrence data from GBIF using the rgbif package

# Load packages
library(rgbif)
library(tidyverse)

# Input GBIF user information
user <- ""
pwd <- ""
email <- ""

# Set name and pull taxa key
taxa <- name_backbone("Leptonycteris yerbabuenae")
taxa
taxa$usageKey

# Build four-letter species code
code <- toupper(paste0(substr(taxa$genericName, 1, 2), substr(taxa$specificEpithet, 1, 2)))

# Set parameters and initiate download
# Set to LEYE occurrences in Arizona, 1991-2020
gbif_download <- occ_download(
  pred("taxonKey", taxa$usageKey),
  pred_in("country", "US"),
  pred("stateProvince", "Arizona"), 
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus", "PRESENT"),
  pred_gte("year", 1991), # greater than or equal to
  pred_lte("year", 2020), # less than or equal to
  user = user, pwd = pwd, email = email,
  format = "SIMPLE_CSV")

# Check status
occ_download_wait(gbif_download)

# Import download
gbif_download <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Save raw occurrence data as .rds and .csv
saveRDS(gbif_download, paste0("data/raw/", code, "_occ_raw.rds"))
write_csv(gbif_download, paste0("data/raw/", code, "_occ_raw.csv"))

