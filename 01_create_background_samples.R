#create background pseudo-absences from tracking data

#clear workspace and set working directory
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

#load required packages
{
  library(dynamicSDM)
  library(sf)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(CCAMLRGIS)
  library(ggplot2)
}

#import coast file from CCAMLRGIS - used to mask out land points in background sampling
coast <- load_Coastline()
coast_v <- vect(coast)

#refresh from here for each new breeding stage
rm(list=setdiff(ls(), c("coast", "coast_v")))

#change species, site, breeding stage, and season variable
this.species <- "ADPE"
this.site <- "Pointe_Geologie"
this.stage <- "chick-rearing"
season <- TRUE #TRUE if breeding stage overlaps the new year, FALSE if not

# 1. Format data for dynamicSDM
#import and format tracks
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
tracks <- tracks %>% select(individual_id, date, x, y)

#isolate a day month and year column for dynamicSDM package
tracks$day <- as.numeric(day(tracks$date))
tracks$month <- month(tracks$date)
tracks$year <- year(tracks$date)

#filter to remove points with NA or invalid coordinates and dates
tracks <- spatiotemp_check(occ.data = tracks,
                           na.handle = "exclude",
                           date.handle = "exclude",
                           date.res = "day",
                           coord.handle = "exclude",
                           duplicate.handle = "exclude")


# 2. Run background sampling for each test run
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
  
  #plot tracks
  tracks_terra <- vect(train_tracks,
                       geom = c("x", "y"),
                       crs = "epsg:4326")
  tracks_terra <- project(tracks_terra, "EPSG:6932")
  plot(tracks_terra, pch = ".")
  
  
  # 3. Fit minimum convex hull
  #convert to sf - terra convHull incompatible with dynamicSDM at times
  tracks_sf <- st_as_sf(tracks_terra)
  
  #create minimum convex hull
  mch_sf <- st_convex_hull(st_union(tracks_sf))
  plot(mch_sf, add=T)
  
  #crop coast to land within minimum convex hull
  e <- ext(tracks_terra)
  crop_coast <- crop(coast_v,e)
  crop_coast <- st_as_sf(crop_coast)
  plot(crop_coast, add=T)
  
  #gets rid of self-intersections - eliminates bugs in dynamicSDM
  coast_buff <- st_buffer(crop_coast,0)
  
  #erase coastline from minimum convex hull
  mch_masked <- st_difference(mch_sf, st_union(coast_buff))
  
  #check final hull and tracks
  plot(mch_masked)
  plot(tracks_sf, pch=".", add=T, col = "black")
  
  #cleanup
  rm(mch_sf, tracks, tracks_terra, crop_coast, coast_buff)
  
  
  # 4. Create background points
  #use dynamicSDM function to create background points
  suppressMessages( #suppressing message that prints when temporal.buffer = 0
    background <- spatiotemp_pseudoabs(occ.data=train_tracks, 
                                       spatial.method = "random",
                                       spatial.ext = mch_masked,
                                       temporal.method = "buffer",
                                       temporal.buffer = 0,
                                       n.pseudoabs = nrow(train_tracks),
                                       prj = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))
  
  
  #plot mask, background points, and tracks
  plot(mch_masked)
  terra::plot(terra::vect(background[, c("x", "y")],
                          geom = c("x", "y"),
                          crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"),
              pch = ".", col = "red", add=T)
  plot(tracks_sf, add=T, pch=".", col="black")
  
  
  # 5. Export
  #format dataframe
  background$date <- as.Date(with(background, paste(year, month, day, sep="-")), "%Y-%m-%d")
  background <- select(background, x, y, date)
  
  #project to EPSG 4326 for later extraction
  background_vect <- terra::vect(background[, c("x", "y", "date")],
                                 geom = c("x", "y"),
                                 crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
  background_vect <- project(background_vect, "EPSG:4326")
  plot(background_vect, pch=".")
  
  #reconvert to dataframe
  background <- as.data.frame(background_vect, geom="XY")
  
  #export
  write.csv(background, paste0("output/background/", this.species, "/", this.site, "/", this.stage, "_", i, ".csv")) 
}


# 6. Repeat the process for all tracks for spatial transfer
train_tracks <- tracks

#plot tracks
tracks_terra <- vect(train_tracks,
                     geom = c("x", "y"),
                     crs = "epsg:4326")
tracks_terra <- project(tracks_terra, "EPSG:6932")
plot(tracks_terra, pch = ".")

#convert to sf - terra convHull incompatible with dynamicSDM at times
tracks_sf <- st_as_sf(tracks_terra)

#create minimum convex hull
mch_sf <- st_convex_hull(st_union(tracks_sf))
plot(mch_sf, add=T)

#crop coast to land within minimum convex hull
e <- ext(tracks_terra)
crop_coast <- crop(coast_v,e)
crop_coast <- st_as_sf(crop_coast)
plot(crop_coast, add=T)

#gets rid of self-intersections - eliminates bugs in dynamicSDM
coast_buff <- st_buffer(crop_coast,0)

#erase coastline from minimum convex hull
mch_masked <- st_difference(mch_sf, st_union(coast_buff))

#check final hull and tracks
plot(mch_masked)
plot(tracks_sf, pch=".", add=T, col = "black")

#cleanup
rm(mch_sf, tracks, tracks_terra, crop_coast, coast_buff)

#use dynamicSDM function to create background points
suppressMessages( #suppressing message that prints when temporal.buffer = 0
  background <- spatiotemp_pseudoabs(occ.data=train_tracks, 
                                     spatial.method = "random",
                                     spatial.ext = mch_masked,
                                     temporal.method = "buffer",
                                     temporal.buffer = 0,
                                     n.pseudoabs = nrow(train_tracks),
                                     prj = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))


#plot mask, background points, and tracks
plot(mch_masked)
terra::plot(terra::vect(background[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"),
            pch = ".", col = "red", add=T)
plot(tracks_sf, add=T, pch=".", col="black")

#format dataframe
background$date <- as.Date(with(background, paste(year, month, day, sep="-")), "%Y-%m-%d")
background <- select(background, x, y, date)

#project to EPSG 4326 for later extraction
background_vect <- terra::vect(background[, c("x", "y", "date")],
                               geom = c("x", "y"),
                               crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
background_vect <- project(background_vect, "EPSG:4326")
plot(background_vect, pch=".")

#reconvert to dataframe
background <- as.data.frame(background_vect, geom="XY")

#only when happy export
write.csv(background, paste0("output/background/", this.species, "/", this.site, "/", this.stage, ".csv")) 