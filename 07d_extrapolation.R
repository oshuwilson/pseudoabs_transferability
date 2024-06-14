#calculate extrapolation for each temporal test
#SOES Marion post-moult incomplete. Need to do SOES WAP and SG ones too.

rm(list=ls())
#setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")
setwd("/mainfs/home/jcw2g17/Chapter 01/")

{
  library(lubridate)
  library(flexsdm)
  library(dplyr)
}

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

#define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")
meta <- meta[-c(1:15, 20:21),]

#loop over every species, site, and stage
for(z in 1:21){
  
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
  shape_all <- NULL
  
  #loop over each season
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
    
    #if test data is missing over 10% of a predictor, remove predictor
    pred_check <- back_test
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    predictors2 <- names(pred_check)
    
    #BUFFERS
    #remove predictors if over 10% of column is NA
    pred_check <- buff_train %>% select(all_of(predictors2))
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    buff_predictors <- names(pred_check)
    
    #make testing data the same columns
    buff_testing <- back_test %>% select(all_of(buff_predictors))
    
    #calculate shape
    buff_shape <- extra_eval(training_data = buff_train, pr_ab = "pa", projection_data = buff_testing, n_cores = 40)
    
    
    #BACKGROUND
    #remove predictors if over 10% of column is NA
    pred_check <- back_train %>% select(all_of(predictors2))
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    back_predictors <- names(pred_check)
    
    #make testing data the same columns
    back_testing <- back_test %>% select(all_of(back_predictors))
    
    #calculate shape
    back_shape <- extra_eval(training_data = back_train, pr_ab = "pa", projection_data = back_testing, n_cores = 40)
    
    
    #CRW
    #remove predictors if over 10% of column is NA
    pred_check <- crw_train %>% select(all_of(predictors2))
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    crw_predictors <- names(pred_check)
    
    #make testing data the same columns
    crw_testing <- back_test %>% select(all_of(crw_predictors))
    
    #calculate shape
    crw_shape <- extra_eval(training_data = crw_train, pr_ab = "pa", projection_data = crw_testing, n_cores = 40)

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
    shape_IQR$season <- this.test
    shape_IQR$pseudo <- c("buffers", "background", "CRWs")
    
    #bind to dataset with information for every season
    shape_values <- rbind(shape_values, shape_IQR)
    
    #do the same for scores of entire dataset
    shape_scores$season <- this.test
    shape_all <- rbind(shape_all, shape_scores)
    
  }
  
  #export dataset
  saveRDS(shape_values,
          file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/shape_medians.RDS"))
  saveRDS(shape_all,
          file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/shape_scores.RDS"))
  
  #print to show that this species has completed
  print(paste0(this.species, " ", this.site, " ", this.stage, " completed"))
  
  #print dataset
  print(shape_values)
  
}
