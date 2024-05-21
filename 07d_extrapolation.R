#calculate extrapolation for each temporal test
rm(list=ls())
setwd("/mainfs/home/jcw2g17/Chapter 01/")

{
  library(dplyr)
  library(lubridate)
  library(flexsdm)
  library(ggplot2)
  library(foreach)
  library(doParallel)
}

#read in metadata
meta <- read.csv("data/species_site_stage_metadata_BRT.csv")

#define predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")

#parallel
registerDoParallel(cores = 21)

foreach(z = 1:21) %dopar% {
  try({
    
    #define parameters in loop
    rm(list=setdiff(ls(), c("meta", "predictors", "z")))
    this.species <- meta[z, 1]
    this.site <- meta[z, 2]
    this.stage <- meta[z, 3]
    season <- meta[z, 4]

# 1. Formatting
#read in presences and pseudo-absences
tracks <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/presences.csv"))
buff <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/buffers.csv"))
back <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/background.csv"))
crw <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/CRWs.csv"))

#create presence/absence column
tracks$pa <- 1
buff$pa <- 0
back$pa <- 0
crw$pa <- 0

#remove extra columns to allow rbind
columns <- c("date", "depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope", "x", "y", "pa")
tracks <- tracks %>% select(all_of(columns))
buff <- buff %>% select(all_of(columns))
back <- back %>% select(all_of(columns))
crw <- crw %>% select(all_of(columns))

#keep background points for testing
back_test <- back

#rbind for models
buff <- rbind(tracks, buff)
back <- rbind(tracks, back)
crw <- rbind(tracks, crw)

#setup date column
buff$date <- as_date(buff$date)
back$date <- as_date(back$date)
crw$date <- as_date(crw$date)
tracks$date <- as_date(tracks$date)


# 2. Isolate training and testing data

#if season = FALSE, separate by year
if(season == FALSE){
  buff$season <- year(buff$date)
  back$season <- year(back$date)
  crw$season <- year(crw$date)
  tracks$season <- year(tracks$date)
}

#if season = TRUE, separate by season
if(season == TRUE){
  buff$season <- year(round_date(buff$date, unit="year"))
  back$season <- year(round_date(back$date, unit="year"))
  crw$season <- year(round_date(crw$date, unit="year"))
  tracks$season <- year(round_date(tracks$date, unit="year"))
}

#extract seasons for loop
seasons <- levels(as.factor(tracks$season))

#create empty table
shape_values <- NULL

#define year
for(i in seasons){
  this.test <- i
  
  #extract training data
  buff_train <- buff %>% filter(season != this.test) %>% select(all_of(predictors), pa)
  back_train <- back %>% filter(season != this.test) %>% select(all_of(predictors), pa)
  crw_train <- crw %>% filter(season != this.test) %>% select(all_of(predictors), pa)
  tracks_train <- tracks %>% filter(season != this.test) %>% select(all_of(predictors), pa)
  
  #extract testing data
  back_test <- back %>% filter(season == this.test) %>% select(all_of(predictors))
  
  #fix sea ice if all values = 0
  if(sum(back_test$sic) == 0){
    back_test <- back_test %>% select(-sic)
  }
  
  #calculate shape
  buff_shape <- extra_eval(training_data = buff_train, pr_ab = "pa", projection_data = back_test)
  back_shape <- extra_eval(training_data = back_train, pr_ab = "pa", projection_data = back_test)
  crw_shape <- extra_eval(training_data = crw_train, pr_ab = "pa", projection_data = back_test)
  
  #store scores
  shape_scores <- cbind(buff_shape[,1], back_shape[,1], crw_shape[,1])
  names(shape_scores) <- c("buff", "back", "crw")
  
  buffers <- quantile(shape_scores$buff, probs=c(0.25, 0.5, 0.75))
  background <- quantile(shape_scores$back, probs=c(0.25, 0.5, 0.75))
  CRWs <- quantile(shape_scores$crw, probs=c(0.25, 0.5, 0.75))
  
  shape_IQR <- as.data.frame(rbind(buffers, background, CRWs))
  names(shape_IQR) <- c("LQ", "Median", "UQ")
  shape_IQR$season <- this.test
  shape_IQR$pseudo <- c("buffers", "background", "CRWs")
  
  shape_values <- rbind(shape_values, shape_IQR)
  
}

saveRDS(shape_values,
        file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/shape_values.RDS"))

})
  
}