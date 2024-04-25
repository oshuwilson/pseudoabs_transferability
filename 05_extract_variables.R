#script for extracting environmental variables to tracks and pseudoabsences
#RUN SOES WAP POST-MOULT
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(terra)
  library(dplyr)
  library(lubridate)
  library(tidyterra)
  library(ggplot2)
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
min(tracks$date) #important for wind and chlorophyll

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
ggplot(tracks, aes(x=depth)) + geom_density()

#remove rows where depth is NA - will be NA for every variable
plot(tracks, pch=".")
tracks <- tracks %>% drop_na(depth)
plot(tracks, pch=".")


###Slope###
slope <- rast("D:/Satellite_Data/static/slope/slope.nc")
tracks$slope <- extract(slope, tracks, ID=F)
ggplot(tracks, aes(x=slope)) + geom_density()

###dShelf###
dshelf <- rast("D:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
tracks$dshelf <- extract(dshelf, tracks, ID=F)
ggplot(tracks, aes(x=dshelf)) + geom_density()

#cleanup static
rm(depth, slope, dshelf)


#---------------
#Dynamic Variables

#dynamic_extract function from 05a script
source("code/05a_dynamic_extract_function.R")


###SST###
tracks <- dynamic_extract(predictor = "sst", tracks)
ggplot(tracks, aes(x=sst)) + geom_density()

###MLD###
tracks <- dynamic_extract(predictor = "mld", tracks)
ggplot(tracks, aes(x=mld)) + geom_density()

###SAL###
tracks <- dynamic_extract(predictor = "sal", tracks)
ggplot(tracks, aes(x=sal)) + geom_density()

###SSH###
tracks <- dynamic_extract(predictor = "ssh", tracks)
ggplot(tracks, aes(x=ssh)) + geom_density()

###SIC###
tracks <- dynamic_extract(predictor = "sic", tracks)
tracks$sic[is.na(tracks$sic)] <- 0 #SIC values of 0 print as NA in GLORYS
ggplot(tracks, aes(x=sic)) + geom_density()

###CURR###
tracks <- dynamic_extract(predictor = "uo", tracks) #eastward velocity
tracks <- dynamic_extract(predictor = "vo", tracks) #northward velocity
tracks$curr <- sqrt((tracks$uo^2) + (tracks$vo^2)) #current speed
ggplot(tracks, aes(x=curr)) + geom_density()

###EKE###
tracks$eke <- 0.5 * ((tracks$uo^2) + (tracks$vo^2))
ggplot(tracks, aes(x=eke)) + geom_density()

###CHL### 
#Does not work for data before 04-09-1997 and needs adjusted function for other dates in 1997
#only use hashtagged code if data starts in 1997 or earlier

cutoff_97 <- as_date("1997-09-04")
cutoff_98 <- as_date("1998-01-01")
tracks_pre97 <- tracks %>% filter(date < cutoff_97)
plot(tracks_pre97)
tracks_97 <- tracks %>% filter(date >= cutoff_97 & date < cutoff_98)
plot(tracks_97)
tracks <- tracks %>% filter(date >= cutoff_98)
plot(tracks, pch=".")


source("code/05b_dynamic_chlorophyll_function.R") #unique function for different file structure

tracks <- dynamic_chlorophyll(predictor = "chl", tracks)
ggplot(tracks, aes(x=chl)) + geom_density()

source("code/05d_dynamic_chlorophyll_1997_function.R")
tracks_97 <- dynamic_chlorophyll_1997(predictor="chl", tracks=tracks_97)
ggplot(tracks_97, aes(x=chl)) + geom_density()

tracks <- bind_spat_rows(tracks, tracks_97, tracks_pre97)
min(tracks$date)
rm(cutoff_97, cutoff_98, tracks_pre97, tracks_97)


###WIND###
#Does not work for data before 01-08-1999 and needs adjusted function for other dates in 1999
#only use hashtagged code if data starts in 1999 or earlier

cutoff_99 <- as_date("1999-08-01")
cutoff_00 <- as_date("2000-01-01")
tracks_pre99 <- tracks %>% filter(date < cutoff_99)
plot(tracks_pre99)
tracks_99 <- tracks %>% filter(date >= cutoff_99 & date < cutoff_00)
plot(tracks_99)
tracks <- tracks %>% filter(date >= cutoff_00)
plot(tracks, pch=".")

source("code/05c_dynamic_wind_function.R")
tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "east")
tracks <- dynamic_wind(predictor = "wind", tracks = tracks, direction = "north")
tracks$wind <- sqrt(tracks$wind_east^2 + tracks$wind_north^2)
ggplot(tracks, aes(x=wind)) + geom_density()

# source("code/05e_dynamic_wind_1999_function.R")
# tracks_99 <- dynamic_wind_1999(predictor="wind", tracks=tracks_99, direction = "east")
# tracks_99 <- dynamic_wind_1999(predictor="wind", tracks=tracks_99, direction = "north")
# tracks_99$wind <- sqrt(tracks_99$wind_east^2 + tracks_99$wind_north^2)
# ggplot(tracks_99, aes(x=wind)) + geom_density()

tracks <- bind_spat_rows(tracks, tracks_99, tracks_pre99)
min(tracks$date)
rm(cutoff_99, cutoff_00, tracks_pre99, tracks_99)

#---------------
#Export
plot(tracks, pch=".")
tracks <- as.data.frame(tracks, geom="XY")

write.csv(tracks, 
          file=paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/presences.csv"))
