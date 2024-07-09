#Create buffer pseudo-absences from tracking data

#clear workspace and set working directory
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

#load required passages
{
  library(dynamicSDM)
  library(sf)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(CCAMLRGIS)
}

#read in coast file for masking
coast <- load_Coastline()
coast_v <- vect(coast)

#refresh from here
rm(list=setdiff(ls(), c("coast", "coast_v")))

#read in species/site/stage - Need to figure out MAPE late chick-rearing, SOES Marion post-breeding and ANFS post-moult
this.species <- "ADPE"
this.site <- "Pointe_Geologie"
this.stage <- "chick-rearing"
buff.value <- 40708 #75th percentile step length for the species

# 1. Format data for dynamicSDM
#read in tracks
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
tracks$individual_id <- as.factor(tracks$individual_id)
tracks <- tracks %>% select(date, x, y)

#isolate a day month and year column
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

#plot tracks
tracks_terra <- vect(tracks,
                     geom = c("x", "y"),
                     crs = "epsg:4326")
tracks_terra <- project(tracks_terra, "EPSG:6932")
plot(tracks_terra, pch=".")

# 2. Exclude land points from tracks
#create extent of tracks for smaller coast file
e <- ext(tracks_terra)
crop_coast <- crop(coast_v, e)
plot(crop_coast, add=T)

#crop tracks to exclude land-based points
tracks_terra <- mask(tracks_terra, crop_coast, inverse=T)

#check changes
plot(tracks_terra, pch=".")
plot(crop_coast, add=T, col="gold")

#convert cropped tracks to a dataframe for dynamicSDM
train_tracks <- as.data.frame(tracks_terra, geom="XY")

# 3. Create mask file
#convert to sf
tracks_sf <- st_as_sf(tracks_terra)

#create minimum convex hull
mch_sf <- st_convex_hull(st_union(tracks_sf))
plot(mch_sf)

#buffer to allow for full range of buffer samples
mch_buff <- st_buffer(mch_sf, buff.value)

#create new crop_coast file
mch_buff_v <- vect(mch_buff)
e2 <- ext(mch_buff_v)
crop_coast2 <- crop(coast_v, e2)
crop_coast2 <- st_as_sf(crop_coast2)

#plot both
plot(mch_buff)
plot(crop_coast2, add=T)

#buffer of 0 removes self-intersections
coast_buff <- st_buffer(crop_coast2,0)

#create MCH with coast masked out
mch_masked <- st_difference(mch_buff, st_union(coast_buff))

#plot mask and tracks
plot(mch_masked)
plot(tracks_sf, pch=".", add=T, col = "black")

# 4. Sample buffer pseudo-absences 

#run buffer creation script
suppressMessages( #suppress message when temporal.buffer=0
  buffers <- spatiotemp_pseudoabs(occ.data = train_tracks,
                    spatial.method = "buffer",
                    spatial.ext = mch_masked, 
                    temporal.method = "buffer",
                    spatial.buffer = c(9000,buff.value), #9000 ensures that pseudoabs falls in a different cell
                    temporal.buffer = 0,
                    n.pseudoabs = nrow(train_tracks),
                    prj = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))


#plot
terra::plot(terra::vect(buffers[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"),
            pch = ".", col = "red", add=T) 


# 5. Export 
#format dataframe
buffers$date <- as.Date(with(buffers, paste(year, month, day, sep="-")), "%Y-%m-%d")
buffers <- select(buffers, x, y, date)

#project to epsg 4326 for extracting covariates
buffers_vect <- terra::vect(buffers[, c("x", "y", "date")],
                               geom = c("x", "y"),
                               crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
buffers_vect <- project(buffers_vect, "EPSG:4326")
plot(buffers_vect, pch=".")

#reconvert to dataframe
buffers <- as.data.frame(buffers_vect, geom="XY")

#export
write.csv(buffers, paste0("output/buffers/", this.species, "/", this.site, "/", this.stage, ".csv"))
