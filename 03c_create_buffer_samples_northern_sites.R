#creating buffer samples for when tracks exceed ESPG 6932 projection limits
#used for Marion SOES, Marion ANFS, Marion SUFS and South_Georgia MAPE
#shorter script, but takes spatiotemp_pseudoabs longer

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
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(rnaturalearthhires)
}

#read in coast file for masking
coast <- load_Coastline()
coast_v <- vect(coast)

#download oceans file for masking
oceans <- ne_download(scale = "medium", category = "physical", type = "ocean", returnclass = "sf")

#refresh from here
rm(list=setdiff(ls(), c("oceans","coast", "coast_v")))

#read in species/site/stage - Need to do SOES Marion post-moult
this.species <- "ANFS"
this.site <- "Marion"
this.stage <- "breeding"
buff.value <- 58968

# 1. Format data for dynamicSDM
#change spreadsheet name for each species
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

#project back to EPSG 4326
tracks_terra <- project(tracks_terra, "EPSG:4326")
plot(tracks_terra, pch=".")

#convert cropped tracks to a dataframe for dynamicSDM
train_tracks <- as.data.frame(tracks_terra, geom="XY")

# 3. Create mask file
#use oceans file for these examples

# 4. Sample buffer pseudo-absences 
#run buffer creation script
suppressMessages( #suppress message when temporal.buffer=0
  buffers <- spatiotemp_pseudoabs(occ.data = train_tracks,
                                  spatial.method = "buffer",
                                  spatial.ext = oceans, 
                                  temporal.method = "buffer",
                                  spatial.buffer = c(12000,buff.value), #9000 ensures that pseudoabs falls in a different cell
                                  temporal.buffer = 0,
                                  n.pseudoabs = nrow(train_tracks),
                                  prj = "+proj=longlat +datum=WGS84"))


#plot
terra::plot(terra::vect(buffers[, c("x", "y")],
                        geom = c("x", "y"),
                        crs = "+proj=longlat +datum=WGS84"),
            pch = ".", col = "red") 


# 5. Export buffers
#format dataframe
buffers$date <- as.Date(with(buffers, paste(year, month, day, sep="-")), "%Y-%m-%d")
buffers <- select(buffers, x, y, date)

#export
write.csv(buffers, paste0("output/buffers/", this.species, "/", this.site, "/", this.stage, ".csv"))
