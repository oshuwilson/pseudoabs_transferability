#code to calculate step lengths
#clear workspace and set working directory
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")
rm(list=ls())

#load required packages
{
  library(lubridate)
  library(dplyr)
  library(amt)
  library(ggplot2)
}

#set species
this.species <- "ADPE"

# 1. read in tracks and format variables
#change spreadsheet for diff species
tracks <- read.csv(paste0("data/tracks_no_gls/", this.species, ".csv"))
tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
tracks$individual_id <- as.factor(tracks$individual_id)

# 2. format tracks
#ensure all observations are complete
if(all(complete.cases(tracks)) == FALSE) {
  tracks <- tracks[complete.cases(tracks),]
}

#make tracks
trks <- make_track(tracks, decimal_longitude, decimal_latitude, date, id = individual_id,
                   crs=4326)

#project coordinates to meters from lat lon
trks <- transform_coords(trks, 6932)

#nest by ID
trks <- trks |> nest(data = -"id")

#resample to daily rate and run steps_by_burst for all individuals
trks2 <- trks |>
  mutate(steps = map(data, function(x)
    x |> track_resample(rate = days(1), tolerance = hours(1)) |> steps_by_burst()))
 

# 3. Calculate daily step length metrics
#unnest for stats 
trks3 <- trks2 |> amt::select(id, steps) |> unnest(cols = steps)

#visualise
ggplot(trks3, aes(x=sl_)) + geom_density()

#75th percentile step length
quantile(trks3$sl_, probs=0.75)
