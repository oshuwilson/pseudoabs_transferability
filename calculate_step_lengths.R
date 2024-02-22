#code to calculate step lengths
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")
rm(list=ls())

{
  library(lubridate)
  library(dplyr)
  library(amt)
  library(ggplot2)
}

# 1. read in tracks and format variables
#change spreadsheet for diff species
coords <- read.csv("data/tracks_no_gls/GHAL.csv")
coords$date <- as.POSIXct(coords$date, format = "%Y-%m-%d %H:%M:%S")
coords$individual_id <- as.factor(coords$individual_id)

# 2. format tracks
#ensure all observations are complete
if(all(complete.cases(coords)) == FALSE) {
  coords <- coords[complete.cases(coords),]
}

#make tracks
trks <- make_track(coords, decimal_longitude, decimal_latitude, date, id = individual_id,
                   crs=4326)

#transform coordinates to meters from lat lon
trks <- transform_coords(trks, 6932)

#nest by ID
trks <- trks |> nest(data = -"id")
trks

#resample to daily rate and run steps_by_burst for all individuals
trks2 <- trks |>
  mutate(steps = map(data, function(x)
    x |> track_resample(rate = days(1), tolerance = hours(1)) |> steps_by_burst()))
trks2
 

# 3. Calculate daily step length metrics
#unnest for stats 
trks3 <- trks2 |> amt::select(id, steps) |> unnest(cols = steps)

#visualise
ggplot(trks3, aes(x=sl_)) + geom_density()

#stats
mean(trks3$sl_) 
median(trks3$sl_)       
max(trks3$sl_)
quantile(trks3$sl_, probs=0.75)


# 4. Try 2 hour sampling
trks4 <- trks |>
  mutate(steps = map(data, function(x)
    x |> track_resample(rate = hours(2), tolerance = minutes(30)) |> steps_by_burst()))
trks4

trks5 <- trks4 |> amt::select(id, steps) |> unnest(cols = steps)
ggplot(trks5, aes(x=sl_)) + geom_density()

#multiply stats by 12 for theoretical daily travel
mean(trks5$sl_) * 12
median(trks5$sl_) * 12