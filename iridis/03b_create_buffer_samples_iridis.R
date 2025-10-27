#-------------------------------------------------------------------------------
# Create buffer pseudo-absences from tracking data
#-------------------------------------------------------------------------------

#clear workspace and set working directory
rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/dsdm")

#load required passages
{
  library(sf)
  library(terra)
  library(tidyterra)
  library(dplyr)
  library(lubridate)
}

#read in coast file for masking
coast <- readRDS("data/coastline.RDS")

#read in species, site, stage, and season info
meta <- read.csv("data/species_site_stage_metadata.csv")

#loop over each case study
for(z in 18){
  
  #clear out variables
  rm(list=setdiff(ls(), c("coast", "meta", "z")))
  
  #change species, site, breeding stage, and season variable
  this.species <- meta$Species[z]
  this.site <- meta$Site[z]
  this.stage <- meta$Stage[z]
  season <- meta$Season[z] #TRUE if breeding stage overlaps the new year, FALSE if not
  
  #get buffer value from 03a script
  buff.value <- read.csv("output/step_lengths.csv") %>%
    filter(Species == this.species) %>%
    pull(X75th.Percentile)
  
#-------------------------------------------------------------------------------
# 1. Format data
#-------------------------------------------------------------------------------
  
  #read in tracks
  tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
  tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
  tracks$individual_id <- as.factor(tracks$individual_id)
  tracks <- tracks %>% select(date, x, y, individual_id)
  
  # convert to terra
  tracks_terra <- tracks %>%
    vect(geom = c("x", "y"), crs = "epsg:4326") %>%
    project("epsg:6932")
  
  #crop coast to track extent
  e <- ext(tracks_terra) + c(buff.value, buff.value, buff.value, buff.value)
  crop_coast <- crop(coast, e) %>%
    aggregate()
  
  #mask out tracks to exclude land-based points
  tracks_terra <- mask(tracks_terra, crop_coast, inverse=T)
  
  
#-------------------------------------------------------------------------------
# 2. Run buffer sample generation
#-------------------------------------------------------------------------------
  
  # set max iterations to number of 15000 or number of tracks
  if(nrow(tracks_terra) < 15000){
    max_n <- nrow(tracks_terra)
  } else {
    max_n <- 15000
  }
  
  # loop over each track location
  for(i in 1:max_n){ 
    
    # get track location
    if(max_n < 15000) {
      loc <- tracks_terra[i]
    } else {
      loc <- sample(tracks_terra, 1)
    }
    
    # buffer location
    loc_buff <- buffer(loc, buff.value)
    
    # erase inner circle to ensure buffer sample falls in different cell
    inner_lim <- buffer(loc, 11320)
    buff <- erase(loc_buff, inner_lim)
    
    # erase crop_coast
    buff <- erase(buff, crop_coast)
    
    # sample buffer point
    pt <- spatSample(buff, 1)
    
    # assign date and individual_id from location
    pt$date <- loc$date
    pt$individual_id <- loc$individual_id
    
    # combine with all points
    if(i == 1){
      buffers <- pt
    } else {
      buffers <- bind_spat_rows(buffers, pt)
    }
    
    # print progress
    if(i %% 100 == 0){
      print(paste0(i, "/", max_n, " complete"))
    }
  }
  
  
#-------------------------------------------------------------------------------
# 3. Export 
#-------------------------------------------------------------------------------
  
  # convert to dataframe
  buff_df <- buffers %>%
    project("epsg:4326") %>%
    as.data.frame(geom = "XY")
  
  #export
  write.csv(buff_df, paste0("output/buffers/", this.species, "_", this.site, "_", this.stage, ".csv"))
  
  # print completion of whole case study
  print(paste0(this.species, " ", this.site, " ", this.stage, " completed"))
  
}
