#-------------------------------------------------------------------------------
# Extract environmental covariates 
#-------------------------------------------------------------------------------

#clear workspace and set working directory
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

#load required packages
{
  library(terra)
  library(dplyr)
  library(lubridate)
  library(tidyterra)
}

# read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

for(z in 1:21){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "z")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  #read in thinned tracks
  tracks <- readRDS(paste0("output/thinned_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_thinned.rds"))
  
  # read in background samples
  back <- read.csv(paste0("output/background/", this.species, "_", this.site, "_", this.stage, ".csv")) %>%
    select(date, individual_id, x, y) %>%
    mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"))
  
  # read in buffer samples
  buffers <- read.csv(paste0("output/buffers/", this.species, "_", this.site, "_", this.stage, ".csv")) %>%
    select(date, individual_id, x, y) %>%
    mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"))
  
  # read in correlated random walks
  crws <- read.csv(paste0("output/CRWs/", this.species, "_", this.site, "_", this.stage, ".csv")) %>%
    rename(individual_id = id) %>%
    select(date, individual_id, x, y) %>%
    mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"))
  
  # list all background yearly sampling files
  pattern <- paste0(this.species, "_", this.site, "_", this.stage)
  back_files <- list.files("output/background/", full.names = T,
                           pattern = pattern)
  
  # for each file
  for(back_file in back_files[2:length(back_files)]){
    
    # read in the file
    bg <- read.csv(back_file) %>%
      select(date, individual_id, x, y) %>%
      mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"))
    
    # get year name from filename
    bg_year <- gsub(paste0("output/background/", pattern, "_"), "", back_file)
    bg_year <- gsub(".csv", "", bg_year)
    gsub(paste0("output/background/", pattern, "_"), "", bg_year)
    
    # create year column
    bg <- bg %>%
      mutate(test_year = bg_year,
             pb = "background") 
    
    # join to all
    if(back_file == back_files[2]){
      back_yrs <- bg
    } else {
      back_yrs <- rbind(back_yrs, bg) 
    }
  }
  
  # create additional column for presence/pseudo-absence type
  tracks <- tracks %>%
    mutate(pb = "presence")
  buffers <- buffers %>% 
    mutate(pb = "buffer")
  back <- back %>%
    mutate(pb = "background")
  crws <- crws %>%
    mutate(pb = "crw")
  
  # combine all together
  data <- bind_rows(tracks, buffers, back, crws)
  
  # create test_year column to bind with yearly_background samples
  data <- data %>%
    mutate(test_year = NA)
  
  # bind together
  data <- bind_rows(data, back_yrs)
  
  #get minimum date for chlorophyll extraction later on
  min_date <- min(data$date) 
  
  
  #-------------------------------------------------------------------------------
  # Static Variables
  #-------------------------------------------------------------------------------
  
  ###Depth###
  depth <- rast("E:/Satellite_Data/static/depth/depth.nc")
  
  #create SpatVector for data
  data <- vect(data,
               geom=c("x", "y"),
               crs=crs(depth)) #this ensures crs are the same as rasters
  
  #extract
  data$depth <- extract(depth, data, ID=F)
  
  #remove rows where depth is NA - will be NA for every variable
  plot(data, pch=".")
  data <- data %>% drop_na(depth)
  plot(data, pch=".")
  
  ###Slope###
  slope <- rast("E:/Satellite_Data/static/slope/slope.nc")
  data$slope <- extract(slope, data, ID=F)
  
  ###dShelf###
  dshelf <- rast("E:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
  data$dshelf <- extract(dshelf, data, ID=F)
  
  #cleanup static
  rm(depth, slope, dshelf)
  
  
  #-------------------------------------------------------------------------------
  # Dynamic Variables
  #-------------------------------------------------------------------------------
  
  #load functions for dynamic extraction
  source("code/R/extraction_functions.R")
  
  ###SST###
  data <- dynamic_extract(predictor = "sst", data, crop = F)
  print("sst")
  
  ###MLD###
  data <- dynamic_extract(predictor = "mld", data, crop = F)
  print("mld")
  
  ###SAL###
  data <- dynamic_extract(predictor = "sal", data, crop = F)
  print("sal")
  
  ###SSH###
  data <- dynamic_extract(predictor = "ssh", data, crop = F)
  print("ssh")
  
  ###SIC###
  data <- dynamic_extract(predictor = "sic", data, crop = F)
  data$sic[is.na(data$sic)] <- 0 #SIC values of 0 print as NA in GLORYS
  print("sic")
  
  ###CURR###
  data <- dynamic_extract(predictor = "uo", data, crop = F) #eastward velocity
  data <- dynamic_extract(predictor = "vo", data, crop = F) #northward velocity
  data$curr <- sqrt((data$uo^2) + (data$vo^2)) #current speed
  print("curr/eke")
  
  ###EKE###
  data$eke <- 0.5 * ((data$uo^2) + (data$vo^2))
  
  ###CHL### 
  #Satellite data unavailable before 04-09-1997 and needs adjusted function for other dates in 1997
  
  #cutoff dates for data from before chl data available and in 1997 needing unique function
  cutoff_97 <- as_date("1997-09-04")
  cutoff_98 <- as_date("1998-01-01")
  
  #if data predate 1998, split into 3 objects for extraction
  if(min_date < cutoff_98){
    data_pre97 <- data %>% filter(date < cutoff_97)
    data_97 <- data %>% filter(date >= cutoff_97 & date < cutoff_98)
    data <- data %>% filter(date >= cutoff_98)
  }
  
  #extract data post-1997
  data <- dynamic_chlorophyll(predictor = "chl", data)
  
  #extract data in 1997 and rejoin data together
  if(min_date < cutoff_98){
    try(data_97 <- dynamic_chlorophyll_1997(predictor="chl", data=data_97))
    
    data <- bind_spat_rows(data, data_97, data_pre97)
    min(data$date)
    rm(cutoff_97, cutoff_98, data_pre97, data_97)
  }
  print("chl")
  
  #-------------------------------------------------------------------------------
  # Export
  #-------------------------------------------------------------------------------
  
  #convert to dataframe
  plot(data, pch=".")
  data <- as.data.frame(data, geom="XY")
  
  #export
  saveRDS(data, 
          file=paste0("output/extraction/", this.species, "_", this.site, "_", this.stage, "_extracted.rds"))
  
  #show species has finished
  print(paste0(this.species, " ", this.site, " ", this.stage, " success"))
  
}


#-------------------------------------------------------------------------------
# Repeat the Process for the Spatial Testing Populations
#-------------------------------------------------------------------------------

# cleanup
rm(list = ls())

# read in table with info for each species, site and stage
meta <- read.csv("data/spatial_site_metadata.csv")

for(z in 1:25){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "z")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  
  #read in thinned tracks
  tracks <- readRDS(paste0("output/thinned_tracks/spatial/", this.species, "_", this.site, "_", this.stage, "_tracks_thinned.rds")) %>%
    mutate(pb = "presence")
  
  # read in background samples
  back <- read.csv(paste0("output/background/spatial/", this.species, "_", this.site, "_", this.stage, ".csv")) %>%
    select(date, individual_id, x, y) %>%
    mutate(date = as_datetime(date, format = "%Y-%m-%d %H:%M:%S"),
           pb = "background")
  
  # combine all together
  data <- bind_rows(tracks, back)
  
  #get minimum date for chlorophyll extraction later on
  min_date <- min(data$date, na.rm = T) 
  
  
  #-------------------------------------------------------------------------------
  # Static Variables
  #-------------------------------------------------------------------------------
  
  ###Depth###
  depth <- rast("E:/Satellite_Data/static/depth/depth.nc")
  
  #create SpatVector for data
  data <- vect(data,
               geom=c("x", "y"),
               crs=crs(depth)) #this ensures crs are the same as rasters
  
  #extract
  data$depth <- extract(depth, data, ID=F)
  
  #remove rows where depth is NA - will be NA for every variable
  plot(data, pch=".")
  data <- data %>% drop_na(depth)
  plot(data, pch=".")
  
  ###Slope###
  slope <- rast("E:/Satellite_Data/static/slope/slope.nc")
  data$slope <- extract(slope, data, ID=F)
  
  ###dShelf###
  dshelf <- rast("E:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
  data$dshelf <- extract(dshelf, data, ID=F)
  
  #cleanup static
  rm(depth, slope, dshelf)
  
  
  #-------------------------------------------------------------------------------
  # Dynamic Variables
  #-------------------------------------------------------------------------------
  
  #load functions for dynamic extraction
  source("code/R/extraction_functions.R")
  
  ###SST###
  data <- dynamic_extract(predictor = "sst", data, crop = F)
  print("sst")
  
  ###MLD###
  data <- dynamic_extract(predictor = "mld", data, crop = F)
  print("mld")
  
  ###SAL###
  data <- dynamic_extract(predictor = "sal", data, crop = F)
  print("sal")
  
  ###SSH###
  data <- dynamic_extract(predictor = "ssh", data, crop = F)
  print("ssh")
  
  ###SIC###
  data <- dynamic_extract(predictor = "sic", data, crop = F)
  data$sic[is.na(data$sic)] <- 0 #SIC values of 0 print as NA in GLORYS
  print("sic")
  
  ###CURR###
  data <- dynamic_extract(predictor = "uo", data, crop = F) #eastward velocity
  data <- dynamic_extract(predictor = "vo", data, crop = F) #northward velocity
  data$curr <- sqrt((data$uo^2) + (data$vo^2)) #current speed
  print("curr/eke")
  
  ###EKE###
  data$eke <- 0.5 * ((data$uo^2) + (data$vo^2))
  
  ###CHL### 
  #Satellite data unavailable before 04-09-1997 and needs adjusted function for other dates in 1997
  
  #cutoff dates for data from before chl data available and in 1997 needing unique function
  cutoff_97 <- as_date("1997-09-04")
  cutoff_98 <- as_date("1998-01-01")
  
  #if data predate 1998, split into 3 objects for extraction
  if(min_date < cutoff_98){
    data_pre97 <- data %>% filter(date < cutoff_97)
    data_97 <- data %>% filter(date >= cutoff_97 & date < cutoff_98)
    data <- data %>% filter(date >= cutoff_98)
  }
  
  #extract data post-1997
  data <- dynamic_chlorophyll(predictor = "chl", data)
  
  #extract data in 1997 and rejoin data together
  if(min_date < cutoff_98){
    try(data_97 <- dynamic_chlorophyll_1997(predictor="chl", data=data_97))
    
    data <- bind_spat_rows(data, data_97, data_pre97)
    min(data$date)
    rm(cutoff_97, cutoff_98, data_pre97, data_97)
  }
  print("chl")
  
  #-------------------------------------------------------------------------------
  # Export
  #-------------------------------------------------------------------------------
  
  #convert to dataframe
  plot(data, pch=".")
  data <- as.data.frame(data, geom="XY")
  
  #export
  saveRDS(data, 
          file=paste0("output/extraction/spatial/", this.species, "_", this.site, "_", this.stage, "_extracted.rds"))
  
  #show species has finished
  print(paste0(this.species, " ", this.site, " ", this.stage, " success"))
  
}
