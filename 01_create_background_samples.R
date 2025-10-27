#-------------------------------------------------------------------------------
# Create background samples from tracking data
#-------------------------------------------------------------------------------

#clear workspace and set working directory
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

#load required packages
{
  library(sf)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(CCAMLRGIS)
  library(ggplot2)
}

#-------------------------------------------------------------------------------
# 1. Data Preparation
#-------------------------------------------------------------------------------

#import coast file from CCAMLRGIS - used to mask out land in background sampling
coast <- readRDS("data/coastline.RDS")

#read in species, site, stage, and season info
meta <- read.csv("data/species_site_stage_metadata.csv")

#loop over each case study
for(z in 1:nrow(meta)){
  
  #clear out variables
  rm(list=setdiff(ls(), c("coast", "meta", "z")))
  
  #change species, site, breeding stage, and season variable
  this.species <- meta$Species[z]
  this.site <- meta$Site[z]
  this.stage <- meta$Stage[z]
  season <- meta$Season[z] #TRUE if breeding stage overlaps the new year, FALSE if not
  
  #import and format tracks
  tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
  tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
  tracks <- tracks %>% select(individual_id, date, x, y)
  
  
#-------------------------------------------------------------------------------
# 2. Run background sampling for each test run
#-------------------------------------------------------------------------------
  
  #if season = FALSE, separate by year
  if(season == FALSE){
    tracks$season <- year(tracks$date)
  }
  
  #if season = TRUE, separate by season
  if(season == TRUE){
    tracks$season <- year(round_date(tracks$date, unit="year"))
  }
  
  #set up loop to run over each season
  tracks$season <- as.factor(tracks$season)
  seasons <- levels(tracks$season)
  
  for(i in seasons){
    
    #extract training data
    train_tracks <- tracks %>% filter(season != i)
    
    #convert tracks to terra
    tracks_terra <- vect(train_tracks,
                         geom = c("x", "y"),
                         crs = "epsg:4326")
    tracks_terra <- project(tracks_terra, "EPSG:6932")
    
    #create minimum convex hull
    mch <- convHull(tracks_terra)
    
    #crop coast to land within minimum convex hull
    e <- ext(tracks_terra)
    crop_coast <- crop(coast,e) 
    
    # create background points
    bg <- spatSample(mch, 40000)
    
    # erase those that fall within coastline
    bg <- erase(bg, crop_coast)
    
    # subsample background points to 20000
    if(nrow(bg) > 20000) {
      bg <- sample(bg, 20000)
    } else {
      stop("Fewer than 20000 background points remaining")
    }
    
    # sample dates and individual_ids from tracks
    bg$date <- sample(train_tracks$date, replace = T)
    bg$individual_id <- sample(train_tracks$individual_id, replace = T)
    
    #format dataframe
    background <- bg %>%
      project("epsg:4326") %>%
      as.data.frame(geom = "XY")
    
    #export
    write.csv(background, paste0("output/background/", this.species, "_", this.site, "_", this.stage, "_", i, ".csv")) 
    
    #update progress
    print(i)
  }
  
  
#-------------------------------------------------------------------------------
# 3. Repeat the process for all tracks for spatial transfer
#-------------------------------------------------------------------------------
  
  #convert tracks to terra
  tracks_terra <- vect(tracks,
                       geom = c("x", "y"),
                       crs = "epsg:4326")
  tracks_terra <- project(tracks_terra, "EPSG:6932")
  
  #create minimum convex hull
  mch <- convHull(tracks_terra)
  
  #crop coast to land within minimum convex hull
  e <- ext(tracks_terra)
  crop_coast <- crop(coast,e)
  
  # create background points
  bg <- spatSample(mch, 40000)
  
  # erase those that fall within coastline
  bg <- erase(bg, crop_coast)
  
  # subsample background points to 20000
  if(nrow(bg) > 20000) {
    bg <- sample(bg, 20000)
  } else {
    stop("Fewer than 20000 background points remaining")
  }
  
  # sample dates and individual_ids from tracks
  bg$date <- sample(tracks$date, replace = T)
  bg$individual_id <- sample(tracks$individual_id, replace = T)
  
  #format dataframe
  background <- bg %>%
    project("epsg:4326") %>%
    as.data.frame(geom = "XY")
  
  #export
  write.csv(background, paste0("output/background/", this.species, "_", this.site, "_", this.stage, ".csv")) 
  
  # print completion of whole case study
  print(paste0(this.species, " ", this.site, " ", this.stage, " completed"))
  
}


#-------------------------------------------------------------------------------
# 4. Create Background Samples for the Spatial Testing Populations
#-------------------------------------------------------------------------------

# clean up
rm(list = ls())

#import coast file from CCAMLRGIS - used to mask out land in background sampling
coast <- readRDS("data/coastline.RDS")

#read in species, site, stage, and season info
meta <- read.csv("data/spatial_site_metadata.csv")

#loop over each case study
for(z in 1:25){
  
  #clear out variables
  rm(list=setdiff(ls(), c("coast", "meta", "z")))
  
  #change species, site, breeding stage, and season variable
  this.species <- meta$Species[z]
  this.site <- meta$Site[z]
  this.stage <- meta$Stage[z]
  
  #import and format tracks
  tracks <- read.csv(paste0("output_old/spatial/", this.species, "/", this.stage, "/", this.site, ".csv"))
  tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
  tracks <- tracks %>% select(individual_id, date, x, y)
  
  # remove tracks from before 1993
  tracks <- tracks %>%
    filter(date >= as_date("1993-01-01"))
  
  #convert tracks to terra
  tracks_terra <- vect(tracks,
                       geom = c("x", "y"),
                       crs = "epsg:4326")
  tracks_terra <- project(tracks_terra, "EPSG:6932")
  
  #create minimum convex hull
  mch <- convHull(tracks_terra)
  
  #crop coast to land within minimum convex hull
  e <- ext(tracks_terra)
  crop_coast <- crop(coast,e)
  
  # create background points
  bg <- spatSample(mch, 40000)
  
  # erase those that fall within coastline
  bg <- erase(bg, crop_coast)
  
  # subsample background points to 20000
  if(nrow(bg) > 20000) {
    bg <- sample(bg, 20000)
  } else {
    stop("Fewer than 20000 background points remaining")
  }
  
  # sample dates and individual_ids from tracks
  bg$date <- sample(tracks$date, replace = T)
  bg$individual_id <- sample(tracks$individual_id, replace = T)
  
  #format dataframe
  background <- bg %>%
    project("epsg:4326") %>%
    as.data.frame(geom = "XY")
  
  #export
  write.csv(background, paste0("output/background/spatial/", this.species, "_", this.site, "_", this.stage, ".csv")) 
  
  # print completion of whole case study
  print(paste0(this.species, " ", this.site, " ", this.stage, " completed"))
  
}
