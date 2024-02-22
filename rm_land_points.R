#script to remove land points from background samples

#working directory on home laptop - may need to change
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(dplyr)
  library(lubridate)
  library(terra)
  library(sf)
  library(CCAMLRGIS)
}

#load coastfile
coast <- load_Coastline()
coast_v <- vect(coast)

#refresh from here
rm(list=setdiff(ls(), c("coast", "coast_v")))

#set parameters
this.species <- "SOES"
this.site <- "WAP"
this.stage <- "post-moult"

#read in tracks and background
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
background <- read.csv(paste0("output/background/", this.species, "/", this.site, "/", this.stage, ".csv"))
buffers <- read.csv(paste0("output/buffers/", this.species, "/", this.site, "/", this.stage, ".csv"))

#convert tracks to terra and EPSG:6932
tracks_terra <- vect(tracks,
                     geom = c("x", "y"),
                     crs = "epsg:4326")
tracks_terra <- project(tracks_terra, "EPSG:6932")
plot(tracks_terra, pch=".")

#create extent of tracks for smaller coast file
e <- ext(tracks_terra)
crop_coast <- crop(coast_v, e)
plot(crop_coast, add=T)


#mask out land points
tracks_terra <- mask(tracks_terra, crop_coast, inverse=T)

#check
plot(tracks_terra, pch=".")
plot(crop_coast, add=T)

#convert to lat/lon
tracks_terra <- project(tracks_terra, "EPSG:4326")

#convert to dataframe
tracks_crop <- as.data.frame(tracks_terra, geom="XY")

#create column for land/sea
tracks_crop$land <- F

#new dataframe to isolate land-based tracks
tracks2 <- tracks %>% left_join(tracks_crop, by=c("individual_id", "date")) %>%
  select(individual_id, date, land)

tracks2[is.na(tracks2)] <- T

#join land column to buffer dataset
#sort both by date first
tracks2$date <- as.POSIXct(tracks2$date, format="%Y-%m-%d %H:%M:%S")
background$date <- as.POSIXct(background$date, format="%Y-%m-%d")

tracks2 <- tracks2 %>% arrange(date)
background <- background %>% arrange(date)

background <- cbind(background, tracks2$land)

#filter to sea-points
background <- background %>% filter(`tracks2$land` == FALSE)
background <- background %>% select(date, x, y)

#remove columns
tracks_crop <- tracks_crop %>% select(-land)


#export at-sea tracks and background
write.csv(tracks_crop, paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, "_at_sea.csv"))
write.csv(background, paste0("output/background/", this.species, "/", this.site, "/", this.stage, "_at_sea.csv"))
