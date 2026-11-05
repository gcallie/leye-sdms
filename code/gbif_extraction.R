# Extract species occurrence data from GBIF using the rgbif package

# Load packages
library(rgbif)
library(tidyverse)

# Input GBIF user information
user <- ""
pwd <- ""
email <- ""

# Species list:
## Leptonycteris yerbabuenae
## Carnegiea gigantea
## Stenocereus thurberi
## Agave palmeri
## Agave parryi
## Agave schottii
## Agave deserti

# Set name and pull taxa key
taxa <- name_backbone("Agave parryi")
taxa
taxa$usageKey

# Build four-letter species code
# Used three letters for species to different agaves
code <- tolower(paste0(substr(taxa$genericName, 1, 2), substr(taxa$specificEpithet, 1, 3)))

# Set parameters and initiate download
# Set to species occurrences in Arizona, 1991-2020
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


# Combine data for agave species into a single csv and rds file
agpar <- readRDS("data/raw/agpar_occ_raw.rds")
agpar <- agpar |> 
  mutate(eventDate = as.character(eventDate)) # resolve type issue preventing row bind
agpal <- readRDS("data/raw/agpal_occ_raw.rds")
agdes <- readRDS("data/raw/agdes_occ_raw.rds")
agdes <- agdes |> 
  mutate(eventDate = as.character(eventDate)) # resolve type issue preventing row bind
agsch <- readRDS("data/raw/agsch_occ_raw.rds")

agave <- bind_rows(agpar, agpal, agdes, agsch)
write_rds(agave, "data/raw/agave_occ_raw.rds")
write_csv(agave, "data/raw/agave_occ_raw.csv")
