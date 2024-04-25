#create background pseudoabsences from tracking data

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(dynamicSDM)
  library(sf)
  library(terra)
  library(dplyr)
  library(lubridate)
  library(CCAMLRGIS)
  library(ggplot2)
  library(parallel)
  library(parallelMap)
}

#sf coast file 
coast <- load_Coastline()
coast_v <- vect(coast)

#refresh from here 
rm(list=setdiff(ls(), c("coast", "coast_v")))

#change species, site, and stage
this.species <- "SOES"
this.site <- "South_Georgia"
this.stage <- "post-moult"

# 1. Format data for dynamicSDM
#change spreadsheet name for each species
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
tracks$individual_id <- as.factor(tracks$individual_id)
tracks <- tracks %>% select(-X, -year)

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


# 2. Isolate training data - do this by study for the real runs
train_tracks <- tracks

#plot tracks
tracks_terra <- vect(train_tracks,
                     geom = c("x", "y"),
                     crs = "epsg:4326")
tracks_terra <- project(tracks_terra, "EPSG:6932")
plot(tracks_terra, pch = ".")


# 3. Fit minimum convex hull -- ADD BUFFER????
#convert to sf
tracks_sf <- st_as_sf(tracks_terra)

#MCH
mch_sf <- st_convex_hull(st_union(tracks_sf))

plot(mch_sf, add=T)

e <- ext(tracks_terra)
crop_coast <- crop(coast_v,e)
crop_coast <- st_as_sf(crop_coast)
plot(crop_coast, add=T)

#gets rid of self-intersections
coast_buff <- st_buffer(crop_coast,0)

mch_masked <- st_difference(mch_sf, st_union(coast_buff))

plot(mch_masked)
plot(tracks_sf, pch=".", add=T, col = "black")

#cleanup
rm(mch_sf, tracks, tracks_terra, crop_coast, coast_buff)


# 4. Create background points - can take a while 
parallelStartSocket(cpus = 8) #parallelize to speed up

suppressMessages( #suppressing message that prints when temporal.buffer = 0
  background <- spatiotemp_pseudoabs(occ.data=train_tracks, 
                     spatial.method = "random",
                     spatial.ext = mch_masked,
                     temporal.method = "buffer",
                     temporal.buffer = 0,
                     n.pseudoabs = nrow(train_tracks),
                     prj = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"))

parallelStop()

plot(mch_masked)
terra::plot(terra::vect(background[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"),
            pch = ".", col = "red", add=T)
plot(tracks_sf, add=T, pch=".", col="black")


# 5. Format for export 
background$date <- as.Date(with(background, paste(year, month, day, sep="-")), "%Y-%m-%d")
background <- select(background, x, y, date)

#project to lat/lon
background_vect <- terra::vect(background[, c("x", "y", "date")],
                               geom = c("x", "y"),
                               crs = "+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
background_vect <- project(background_vect, "EPSG:4326")
plot(background_vect, pch=".")

background <- as.data.frame(background_vect, geom="XY")

#only when happy export
write.csv(background, paste0("output/background/", this.species, "/", this.site, "/", this.stage, ".csv")) 
