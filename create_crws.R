#create CRW simulation from tracks

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(aniMotum)
  library(dplyr)
  library(lubridate)
  library(terra)
  library(sf)
  library(tidyr)
}


#refresh from here
rm(list=ls())

#read in species/site/stage - not sure why EMPE Mawson is too small - redo ADPE???
this.species <- "SOES"
this.site <- "Marion"
this.stage <- "post-moult"
#central place forager?
CPF <- TRUE
#timestep of state-space-modeled tracks for this species (in hours)
timestep <- 2


# 1. Read in and format tracking data
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, "_at_sea.csv"))
tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
tracks$individual_id <- as.factor(tracks$individual_id)

#create lc column to tell aniMotum that errors are included
tracks$lc <- "GL"

#rename and select columns for aniMotum model fitting
tracks <- tracks %>% rename(id = individual_id, lon = x, lat = y, x.sd = longitude_se, y.sd = latitude_se)
tracks <- tracks %>% select(id, date, lc, lon, lat, x.sd, y.sd)

#remove individuals with only one point 
IDs <- tracks %>% group_by(id) %>% summarise(n = n())
IDs <- filter(IDs, n < 2)
IDs$id <- as.character(IDs$id)
IDs$id <- as.factor(IDs$id)
solo <- levels(IDs$id)

tracks <- tracks %>% filter(!id %in% solo)

#2. Split each ID into trips (only needed for CPFs??)
#use temporal filter? - should work with no land points but be vigilant in checking
#will this work when animals rest on ice?

tracks_by_ID <- tracks %>% group_by(id) %>%  #group_by_ID
  mutate(date0 = lag(date)) %>% #create column for date at point t-1
  mutate(timegap = interval(date0, date)/hours(1)) %>% #create column for the diff between t and t-1
  mutate(timegap = if_else(row_number() == 1, 1, timegap)) %>% #makes sure that first value is 1 not NA
  mutate(trip = as.factor(cumsum(c(1, diff(timegap) > timestep)))) %>% #creates trip number column, allows for one missed location fix
  mutate(trip_ID = as.factor(paste(id, trip, sep="_"))) #create trip ID columns

#identify trips too short for aniMotum (n<5)
trip_lengths <- tracks_by_ID %>% group_by(trip_ID) %>% #group by ID and trip
  summarise(n=n()) %>% #count row numbers for each trip
  filter(n < 5) #isolate short trips

#extract short trip IDs
trip_lengths$trip_ID <- as.factor(as.character(trip_lengths$trip_ID)) #removes deleted levels
solo_trips <- levels(trip_lengths$trip_ID)

#filter out short trips
tracks_by_ID <- tracks_by_ID %>% filter(!trip_ID %in% solo_trips)

#only keep relevant columns
tracks_by_ID <- tracks_by_ID %>% ungroup() %>% 
  select(trip_ID, date, lc, lon, lat, x.sd, y.sd) %>%
  rename(id = trip_ID)

tracks <- tracks_by_ID

# 3. Create CRWs

#plot tracks
tracks_terra <- vect(tracks,
                     geom = c("lon", "lat"),
                     crs = "epsg:4326")
tracks_terra <- project(tracks_terra, "EPSG:6932")
plot(tracks_terra, pch=".")

#run CRW model
fit <- fit_ssm(tracks, model="crw", control=ssm_control(verbose=0), spdf = FALSE)
st <- sim_fit(fit, what="fitted", rep=20, cpf=CPF) #simulates 20 tracks per ID
plot(st[102,])

#keep the most representative track (limits number of tracks that were over land)
st_filter <- sim_filter(st, keep=0.05) #use min(lat) and then null filter for antarctic CPFs
plot(st_filter[103:105,])

#reroute away from land
st_routed <- route_path(st_filter)
plot(st_routed[106:109,])

#3 ID erroneous tracks - can be unrealistic (e.g. circumpolar) or not CPFs among CPF data
plot(st_routed[114,])

#change by study
resample <- c(3, 22, 51:52, 54:56, 71:73, 78, 90:91, 93, 95:101, 103:105, 106:109, 111:112)
nonCPF <- c(4, 8, 36, 44, 45, 70)
discard <- c(2:3, 6:7, 9, 22, 25:28, 30:32, 40:42, 51:52, 54:56, 58, 61:62, 71:73,
             78, 81, 83, 88, 90, 91, 94)


#isolate error tracks and refit CRWs
errors <- filter(st_routed, row_number() %in% nonCPF)
errors$id <- as.factor(errors$id)
error_IDs <- levels(errors$id)

#isolate error trips 
error_tracks <- tracks %>% filter(id %in% error_IDs)

#run CRWs
error_fit <- fit_ssm(error_tracks, model="crw", control=ssm_control(verbose=0), spdf = FALSE)
error_st <- sim_fit(error_fit, what="fitted", rep=10, cpf=CPF) #simulates multiple tracks per ID
plot(error_st[3,])
error_st_filter <- sim_filter(error_st, keep = 0.1)
plot(error_st_filter[3,])
error_st_routed <- route_path(error_st_filter, centroids = TRUE)

#check all have been fixed
plot(error_st_routed[6,])

st_nonCPF <- error_st_routed

st_pure <- st_routed %>% filter(!id %in% error_IDs)
fixed_st <- rbind(st_pure, st_resampled, st_nonCPF)

#4. extract data
CRW <- fixed_st %>% unnest(cols = c(sims)) %>% 
  filter(rep!=0) %>%
  select(id, date, lon, lat)

#plot CRW and tracks
CRW_terra <- vect(CRW,
                  geom = c("lon", "lat"),
                  crs = "epsg:4326")
CRW_terra <- project(CRW_terra, "EPSG:6932")

plot(CRW_terra, col="red", pch=".")
plot(tracks_terra, add=T, pch=".", col="black")

# 3. Format and Export
CRW <- CRW %>% rename(x=lon, y=lat)
tracks <- tracks %>% select(id, date, lon, lat) %>% rename(x=lon, y=lat)

write.csv(CRW, paste0("output/CRWs/", this.species, "/", this.site, "/", this.stage, ".csv"))
write.csv(tracks, paste0("output/CRWs/", this.species, "/", this.site, "/", this.stage, "_presences.csv"))

