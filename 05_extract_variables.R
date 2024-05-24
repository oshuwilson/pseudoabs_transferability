#script for extracting environmental variables to tracks
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(terra)
  library(dplyr)
  library(lubridate)
  library(tidyterra)
}

#setup
rm(list=ls())
this.species <- "SOES"
this.site <- "WAP"
this.stage <- "post-moult"

#read in tracks
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))

#only keep relevant columns
tracks <- tracks %>% select(individual_id, date, stage, x, y)

#format date (for temporal extraction later)
tracks$date <- as_datetime(tracks$date)
min_date <- min(tracks$date) #important for wind and chlorophyll

#---------------
#Static Variables

###Depth###
depth <- rast("D:/Satellite_Data/static/depth/depth.nc")

#create SpatVector for tracks
tracks <- vect(tracks,
               geom=c("x", "y"),
               crs=crs(depth)) #this ensures crs are the same as rasters

#extract
tracks$depth <- extract(depth, tracks, ID=F)

#remove rows where depth is NA - will be NA for every variable
plot(tracks, pch=".")
tracks <- tracks %>% drop_na(depth)
plot(tracks, pch=".")

###Slope###
slope <- rast("D:/Satellite_Data/static/slope/slope.nc")
tracks$slope <- extract(slope, tracks, ID=F)

###dShelf###
dshelf <- rast("D:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
tracks$dshelf <- extract(dshelf, tracks, ID=F)

#cleanup static
rm(depth, slope, dshelf)


#---------------
#Dynamic Variables

#dynamic_extract function from 05a script
source("code/05a_dynamic_extract_function.R")


###SST###
tracks <- dynamic_extract(predictor = "sst", tracks)

###MLD###
tracks <- dynamic_extract(predictor = "mld", tracks)

###SAL###
tracks <- dynamic_extract(predictor = "sal", tracks)

###SSH###
tracks <- dynamic_extract(predictor = "ssh", tracks)

###SIC###
tracks <- dynamic_extract(predictor = "sic", tracks)
tracks$sic[is.na(tracks$sic)] <- 0 #SIC values of 0 print as NA in GLORYS

###CURR###
tracks <- dynamic_extract(predictor = "uo", tracks) #eastward velocity
tracks <- dynamic_extract(predictor = "vo", tracks) #northward velocity
tracks$curr <- sqrt((tracks$uo^2) + (tracks$vo^2)) #current speed

###EKE###
tracks$eke <- 0.5 * ((tracks$uo^2) + (tracks$vo^2))

###CHL### 
#Satellite data unavailable before 04-09-1997 and needs adjusted function for other dates in 1997

cutoff_97 <- as_date("1997-09-04")
cutoff_98 <- as_date("1998-01-01")

if(min_date < cutoff_98){
tracks_pre97 <- tracks %>% filter(date < cutoff_97)
tracks_97 <- tracks %>% filter(date >= cutoff_97 & date < cutoff_98)
tracks <- tracks %>% filter(date >= cutoff_98)
}

source("code/05b_dynamic_chlorophyll_function.R") #unique function for different file structure
tracks <- dynamic_chlorophyll(predictor = "chl", tracks)

if(min_date < cutoff_98){
source("code/05d_dynamic_chlorophyll_1997_function.R")
try(tracks_97 <- dynamic_chlorophyll_1997(predictor="chl", tracks=tracks_97))

tracks <- bind_spat_rows(tracks, tracks_97, tracks_pre97)
min(tracks$date)
rm(cutoff_97, cutoff_98, tracks_pre97, tracks_97)
}

###WIND###
#Satellite data unavailable before 01-08-1999 and needs adjusted function for other dates in 1999

cutoff_99 <- as_date("1999-08-01")
cutoff_00 <- as_date("2000-01-01")

if(min_date < cutoff_00){
tracks_pre99 <- tracks %>% filter(date < cutoff_99)
tracks_99 <- tracks %>% filter(date >= cutoff_99 & date < cutoff_00)
tracks <- tracks %>% filter(date >= cutoff_00)
}

source("code/05c_dynamic_wind_function.R")
tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "east")
tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "north")
tracks$wind <- sqrt(tracks$wind_east^2 + tracks$wind_north^2)

if(min_date < cutoff_00){
source("code/05e_dynamic_wind_1999_function.R")
try(tracks_99 <- dynamic_wind_1999(predictor="wind", tracks=tracks_99, direction = "east"))
try(tracks_99 <- dynamic_wind_1999(predictor="wind", tracks=tracks_99, direction = "north"))
tracks_99$wind <- sqrt(tracks_99$wind_east^2 + tracks_99$wind_north^2)

tracks <- bind_spat_rows(tracks, tracks_99, tracks_pre99)
min(tracks$date)
rm(cutoff_99, cutoff_00, tracks_pre99, tracks_99)
}

#---------------
#Export
plot(tracks, pch=".")
tracks <- as.data.frame(tracks, geom="XY")

write.csv(tracks, 
          file=paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/presences.csv"))
