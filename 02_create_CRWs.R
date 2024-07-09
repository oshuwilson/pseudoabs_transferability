#create Correlated Random Walk simulations from tracks

#clear workspace and set working directory
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

#load required packages
{
  library(aniMotum)
  library(dplyr)
  library(lubridate)
  library(terra)
  library(sf)
  library(tidyr)
}


#refresh from here for each new breeding stage
rm(list=ls())

#read in species/site/stage
this.species <- "ADPE"
this.site <- "Pointe_Geologie"
this.stage <- "chick-rearing"
CPF <- TRUE #central place forager? TRUE for majority of tracks except for in HUWH and SOES post-moult


# 1. Read in and format tracking data for aniMotum
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))
tracks$date <- as.POSIXct(tracks$date, format = "%Y-%m-%d %H:%M:%S")
tracks$individual_id <- as.factor(tracks$individual_id)
tracks <- tracks %>% select(individual_id, date, x, y, longitude_se, latitude_se)

#extract timestep from state-spaced modeled tracks (in hours)
timestep <- as.numeric(tracks[2,]$date - tracks[1,]$date)

#create lc column to tell aniMotum that errors are included
tracks$lc <- "GL"

#rename and select columns for aniMotum model fitting
tracks <- tracks %>% rename(id = individual_id, lon = x, lat = y, x.sd = longitude_se, y.sd = latitude_se)


#2. Split each ID into trips - unenecessary for HUWH
tracks_by_ID <- tracks %>% group_by(id) %>%  #group_by_ID
  mutate(date0 = lag(date)) %>% #create column for date at point t-1
  mutate(timegap = interval(date0, date)/hours(1)) %>% #create column for the diff between t and t-1
  mutate(timegap = if_else(row_number() == 1, 1, timegap)) %>% #makes sure that first value is 1 not NA
  mutate(trip = as.factor(cumsum(c(1, diff(timegap) > timestep)))) %>% #creates trip number column, allows for one missed location fix
  mutate(trip_ID = as.factor(paste(id, trip, sep="_"))) #create trip ID columns

#identify trips too short for aniMotum (n<10)
trip_lengths <- tracks_by_ID %>% group_by(trip_ID) %>% #group by ID and trip
  summarise(n=n()) %>% #count row numbers for each trip
  filter(n < 10) #isolate short trips

#extract short trip IDs
trip_lengths$trip_ID <- as.factor(as.character(trip_lengths$trip_ID)) #removes deleted levels
short_trips <- levels(trip_lengths$trip_ID)

#filter out short trips
tracks_by_ID <- tracks_by_ID %>% filter(!trip_ID %in% short_trips)

#only keep relevant columns for aniMotum
tracks_by_ID <- tracks_by_ID %>% ungroup() %>% 
  select(trip_ID, date, lc, lon, lat, x.sd, y.sd) %>%
  rename(id = trip_ID)

#replace tracks dataset with tracks_by_ID
tracks <- tracks_by_ID


# 3. Create CRWs
#plot tracks
tracks_terra <- vect(tracks,
                     geom = c("lon", "lat"),
                     crs = "epsg:4326")
tracks_terra <- project(tracks_terra, "EPSG:6932")
plot(tracks_terra, pch=".")

#fit movement model 
fit <- fit_ssm(tracks, model="crw", control=ssm_control(verbose=0), spdf = FALSE)

#simulate 10 tracks per trip
st <- sim_fit(fit, what="fitted", rep=10, cpf=CPF)
plot(st[1,])

#keep the most representative track
#uses bearing and distance travelled as in Hazen et al. 2017 
#for Antarctic CPFs (e.g. EMPE), min latitude can be used first to filter out inland tracks
#this would be sim_filter(st, keep = 0.5, var = "lat", FUN = "min")
st_filter <- sim_filter(st, keep=0.1)
plot(st_filter[1,])

#reroute CRW tracks away from land
st_routed <- route_path(st_filter)
plot(st_routed[1,])


# 4. ID erroneous tracks (5 at a time)
#these can be unrealistic (e.g. circumpolar), where track ends before returning to colony, or condensed around coastlines
plot(st_routed[1:5, ])

#record error tracks
resample <- c(1, 3, 15) #unrealistic CRW distances or compressed around coast
nonCPF <- c(21, 43) #tracks that terminate without returning to/near the colony
discard <- c(22) #tracks that appear to be resting with little movement


# 5. Remove erroneous tracks 
#extract error trip_IDs
error <- c(resample, nonCPF, discard)
errors <- filter(st_routed, row_number() %in% error)
errors$id <- as.factor(errors$id)
error_IDs <- levels(errors$id)

#remove error IDs to keep good CRWs
st_pure <- st_routed %>% filter(!id %in% error_IDs)


# 6. Resample CRWs
#extract IDs to be resampled
resamples <- filter(st_routed, row_number() %in% resample)
resamples$id <- as.factor(resamples$id)
resample_IDs <- levels(resamples$id)

#isolate resampling trips 
resample_tracks <- tracks %>% filter(id %in% resample_IDs)

#rerun CRWs
resample_fit <- fit_ssm(resample_tracks, model="crw", control=ssm_control(verbose=0), spdf = FALSE)
resample_st <- sim_fit(resample_fit, what="fitted", rep=10, cpf=CPF) 
plot(resample_st[1,])

#can filter by lat/lon sd if prone to unrealistic distances
resample_st_filter <- sim_filter(resample_st, keep=0.1, var=c("lat"), FUN = "sd")
plot(resample_st_filter[1,])

#reroute around land
resample_st_routed <- route_path(resample_st_filter)
plot(resample_st_routed[1,])

#check all refit CRWs have been fixed
plot(resample_st_routed[1,])

#if yes, proceed
st_resampled <- resample_st_routed


# 7. Resample Non-CPF CRWs
#extract IDs to be resampled
nonCPFs <- filter(st_routed, row_number() %in% nonCPF)
nonCPFs$id <- as.factor(nonCPFs$id)
nonCPF_IDs <- levels(nonCPFs$id)

#isolate resampling trips 
nonCPF_tracks <- tracks %>% filter(id %in% nonCPF_IDs)

#rerun CRWs
nonCPF_fit <- fit_ssm(nonCPF_tracks, model="crw", control=ssm_control(verbose=0), spdf = FALSE)
nonCPF_st <- sim_fit(nonCPF_fit, what="fitted", rep=20, cpf=FALSE)
plot(nonCPF_st[1,])

#can filter by lat/lon sd if prone to unrealistic distances
nonCPF_st_filter <- sim_filter(nonCPF_st, keep = 0.5, var="lat", FUN="min")
nonCPF_st_filter <- sim_filter(nonCPF_st_filter, keep=0.1, var=c("lat", "lon"), FUN="sd")
plot(nonCPF_st_filter[1,])

#reroute around land
nonCPF_st_routed <- route_path(nonCPF_st_filter)
plot(nonCPF_st_routed[1,])

#check all refit CRWs have been fixed
plot(nonCPF_st_routed[1,])

#if yes, proceed
st_nonCPF <- nonCPF_st_routed

#combine all good and resampled tracks
fixed_st <- rbind(st_pure, st_resampled, st_nonCPF)


# 8. Visualise CRWs
#extract key info from CRW models
CRW <- fixed_st %>% unnest(cols = c(sims)) %>% 
  filter(rep!=0) %>%
  select(id, date, lon, lat)

#plot CRWs and tracks
CRW_terra <- vect(CRW,
                  geom = c("lon", "lat"),
                  crs = "epsg:4326")
CRW_terra <- project(CRW_terra, "EPSG:6932")
plot(CRW_terra, col="red", pch=".")
plot(tracks_terra, add=T, pch=".", col="black")

# 9. Format and Export
CRW <- CRW %>% rename(x=lon, y=lat)
write.csv(CRW, paste0("output/CRWs/", this.species, "/", this.site, "/", this.stage, ".csv"))
