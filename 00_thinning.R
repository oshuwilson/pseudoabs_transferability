#-------------------------------------------------------------------------------
# Apply spatiotemporal thinning
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(tidyverse)
library(GeoThinneR)
library(terra)

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

#for each case study
for(z in 1:21){
  
  # define species site and stage
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  
  # read in data
  data <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
  
  # format data
  data <- data %>%
    select(individual_id, date, x, y) %>%
    mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"))
  
  # read in depth raster to spatially thin data to the same grid
  depth <- rast("~/OneDrive - University of Southampton/Documents/Predictor Data/processing/dShelf/depth.nc")
  
  # thin spatially to one point per grid cell
  quick_thin <- thin_points(
    data = data,
    lon_col = "x",
    lat_col = "y",
    method = "grid",
    raster_obj = depth
  )
  
  # get thinned data
  thinned <- largest(quick_thin)
  
  # temporal thinning - 1 point per day for most species, 1 per 2 days for SOES, SUFS, and ANFS
  if(this.species %in% c("ANFS", "SOES", "SUFS")){
    thinned <- thinned %>%
      mutate(period = floor_date(date, "2 days")) %>%
      group_by(individual_id, period) %>%
      slice_sample(n = 1) %>%
      ungroup() %>%
      select(-period)
  } else {
    thinned <- thinned %>%
      mutate(period = floor_date(date, "days")) %>%
      group_by(individual_id, period) %>%
      slice_sample(n = 1) %>%
      ungroup() %>% 
      select(-period)
  }
  
  # export thinned tracks
  saveRDS(thinned, 
          file = paste0("output/thinned_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_thinned.rds"))

  # print progress
  print(paste("Thinned data for", this.species, this.site, this.stage))
  
  # remove objects
  rm(data, thinned, quick_thin, this.species, this.stage, this.site)
  
}


#-------------------------------------------------------------------------------
# Repeat for Spatial Testing Data
#-------------------------------------------------------------------------------

# clear up
rm(list=ls())

#read in table with info for each species, site and stage
meta <- read.csv("data/spatial_site_metadata.csv")

#for each case study
for(z in 1:25){
  
  # define species site and stage
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  
  # read in data
  data <- read.csv(paste0("output_old/spatial/", this.species, "/", this.stage, "/", this.site, ".csv"))
  
  # format data
  data <- data %>%
    select(individual_id, date, x, y) %>%
    mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"))
  
  # read in depth raster to spatially thin data to the same grid
  depth <- rast("~/OneDrive - University of Southampton/Documents/Predictor Data/processing/dShelf/depth.nc")
  
  # thin spatially to one point per grid cell
  quick_thin <- thin_points(
    data = data,
    lon_col = "x",
    lat_col = "y",
    method = "grid",
    raster_obj = depth
  )
  
  # get thinned data
  thinned <- largest(quick_thin)
  
  # temporal thinning - 1 point per day for most species, 1 per 2 days for SOES, SUFS, and ANFS
  if(this.species %in% c("ANFS", "SOES", "SUFS")){
    thinned <- thinned %>%
      mutate(period = floor_date(date, "2 days")) %>%
      group_by(individual_id, period) %>%
      slice_sample(n = 1) %>%
      ungroup() %>%
      select(-period)
  } else {
    thinned <- thinned %>%
      mutate(period = floor_date(date, "days")) %>%
      group_by(individual_id, period) %>%
      slice_sample(n = 1) %>%
      ungroup() %>% 
      select(-period)
  }
  
  # export thinned tracks
  saveRDS(thinned, 
          file = paste0("output/thinned_tracks/spatial/", this.species, "_", this.site, "_", this.stage, "_tracks_thinned.rds"))
  
  # print progress
  print(paste("Thinned data for", this.species, this.site, this.stage))
  
  # remove objects
  rm(data, thinned, quick_thin, this.species, this.stage, this.site)
  
}
