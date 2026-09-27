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
gbif_download <- occ_download(
  pred("taxonKey", taxa$usageKey),
  pred_in("country", c("US", "MX")),
  #pred("stateProvince", "Arizona"),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus", "PRESENT"),
  pred_gte("year", 1991), # greater than or equal to
  pred_lte("year", 2020), # less than or equal to
  pred_within( # 2000-2025 GBIF saguaro occurrence minimum convex polygon
    "POLYGON((
     -114.89871399999998 33.15554700000007,
     -114.62999999999994 32.700000000000045,
     -112.14844499999998 28.978038000000026,
     -111.68316699999997 28.451944000000026,
     -111.23698199999995 28.054960000000051,
     -111.08930799999996 27.938528000000076,
     -109.28783299999998 26.740488000000028,
     -109.00083299999994 30.775000000000034,
     -109.76818799999995 33.037445000000048,
     -110.61082499999998 33.823405000000037,
     -113.74968699999994 35.717705000000024,
     -114.17593299999999 34.906513000000075,
     -114.37024699999995 34.452586000000053,
     -114.89871399999998 33.15554700000007
   ))"
  ),
  user = user, pwd = pwd, email = email,
  format = "SIMPLE_CSV")


# Check status
occ_download_wait(gbif_download)

# Import download
gbif_download <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Save raw occurrence data as .rds and .csv
saveRDS(gbif_download, paste0("data/", code, "_occ_raw.rds"))
write_csv(gbif_download, paste0("data/", code, "_occ_raw.csv"))
