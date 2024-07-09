#automated script to test spatial extrapolation

#clear workspace and set working directory
rm(list=ls())
#setwd("/mainfs/home/jcw2g17/Chapter 01/")
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

#load required packages
{
  library(lubridate)
  library(flexsdm)
  library(dplyr)
}

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")
meta2 <- read.csv("output/spatial/spatial_site_metadata.csv")

#define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")

#remove ANPE, EMPE, SUFS, and MAPE incubation/post-breeding (no spatial transfer)
meta <- meta %>% filter(Species != "ANPE" & Species != "EMPE" & Species != "SUFS")
meta <- meta %>% filter(Species != "MAPE" | 
                          (Species == "MAPE" & Stage != "incubation" & Stage != "post-breeding"))

#loop over every species, site, and stage
for(z in 1:nrow(meta)){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "predictors", "z", "meta2")))
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
  
  #rbind for models
  buff <- rbind(tracks, buff)
  back <- rbind(tracks, back)
  crw <- rbind(tracks, crw)
  
  #setup date column
  buff$date <- as_date(buff$date)
  back$date <- as_date(back$date)
  crw$date <- as_date(crw$date)
  tracks$date <- as_date(tracks$date)
  
  # 2. Load in test data
  #filter spatial metadata to this species and stage
  meta3 <- meta2 %>% filter(Species == this.species & Stage == this.stage)
  
  #extract list of sites for this species and stage
  meta3$Site <- as.factor(meta3$Site)
  sites <- levels(meta3$Site)
  
  #create empty table
  shape_values <- NULL
  shape_all <- NULL
  
  #loop over each season
  for(i in sites){
    test.site <- i
    
    #load in test data
    back_test <- read.csv(paste0("output/spatial/", this.species, "/", this.stage, "/extraction/", test.site, "_background.csv"))
    tracks_test <- read.csv(paste0("output/spatial/", this.species, "/", this.stage, "/extraction/", test.site, "_presences.csv"))
    
    #remove chl and/or wind from predictors list if relevant
    meta4 <- meta3 %>% filter(Site == test.site)
    if(meta4$Missing[1] == "windchl"){
      shape_predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "slope")
    }
    
    if(meta4$Missing[1] == "chl"){
      shape_predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "wind", "eke", "slope")
    }
    
    if(meta4$Missing[1] == ""){
      shape_predictors <- predictors
    }
    
    #fix sea ice if all values = 0
    if(sum(back_test$sic) == 0 & sum(back$sic == 0)){
      shape_predictors <- shape_predictors[shape_predictors != "sic"]
    }
    
    #only select predictors for training and testing
    back_test <- back_test %>% select(all_of(shape_predictors))
    tracks_test <- tracks_test %>% select(all_of(shape_predictors))
    
    buff_train <- buff %>% select(all_of(shape_predictors), pa)
    back_train <- back %>% select(all_of(shape_predictors), pa)
    crw_train <- crw %>% select(all_of(shape_predictors), pa)
    
    #merge testing data with PA column
    tracks_test$pa <- 1
    back_test$pa <- 0
    test_data <- rbind(tracks_test, back_test)
    
    #BUFFERS
    #calculate shape
    buff_shape <- extra_eval(training_data = buff_train, pr_ab = "pa", projection_data = test_data, n_cores = 7)
    
    #BACKGROUND
    #calculate shape
    back_shape <- extra_eval(training_data = back_train, pr_ab = "pa", projection_data = test_data, n_cores = 7)

    #CRW
    #calculate shape
    crw_shape <- extra_eval(training_data = crw_train, pr_ab = "pa", projection_data = test_data, n_cores = 7)
    
    #store scores
    shape_scores <- cbind(buff_shape[,1], back_shape[,1], crw_shape[,1])
    names(shape_scores) <- c("buff", "back", "crw")
    
    #extract LQ, median, and UQ for each set of scores
    buffers <- quantile(shape_scores$buff, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
    background <- quantile(shape_scores$back, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
    CRWs <- quantile(shape_scores$crw, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
    
    #format data table with all LQ, median, and UQ values
    shape_IQR <- as.data.frame(rbind(buffers, background, CRWs))
    names(shape_IQR) <- c("LQ", "Median", "UQ")
    shape_IQR$region <- test.site
    shape_IQR$pseudo <- c("buffers", "background", "CRWs")
    
    #bind to dataset with information for every season
    shape_values <- rbind(shape_values, shape_IQR)
    
    shape_scores$region <- test.site
    shape_all <- rbind(shape_all, shape_scores)
    
  }
  
  #export dataset
  saveRDS(shape_values,
          file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/shape_medians.RDS"))
  saveRDS(shape_all,
          file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/shape_scores.RDS"))
  
  #print to show that this species has completed
  print(paste0(this.species, " ", this.site, " ", this.stage, " completed"))
  
  #print dataset
  print(shape_values)
  
}
