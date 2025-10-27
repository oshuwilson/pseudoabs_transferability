#-------------------------------------------------------------------------------
# Temporal Environmental Extrapolation
#-------------------------------------------------------------------------------

#clear workspace and set working directory
rm(list=ls())
setwd("/mainfs/scratch/jcw2g17/dsdm/")

#load required packages
{
  library(lubridate)
  library(flexsdm)
  library(dplyr)
}

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

#define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")
#meta <- meta[-c(1:15, 20:21),]

#loop over every species, site, and stage
for(z in 19:21){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "predictors", "z")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  #-------------------------------------------------------------------------------
  # 1. Formatting
  #-------------------------------------------------------------------------------
  
  #read in presences and pseudo-absences with environmental data
  data <- readRDS(paste0("output/extraction/", this.species, "_", this.site, "_", this.stage, "_extracted.rds"))
  
  #if season = FALSE, separate by year
  if(season == FALSE){
    data$season <- year(data$date)
  }
  
  #if season = TRUE, separate by season
  if(season == TRUE){
    data$season <- year(round_date(data$date, unit="year"))
  }
  
  #extract seasons for loop
  seasons <- unique(data$season)
  
  #create empty tables
  shape_values <- NULL
  shape_all <- NULL
  
  #loop over each season
  for(i in seasons){
    this.test <- i
    
    #extract training data
    train <- data %>%
      filter(season != this.test)
    
    #create test dataset
    test <- data %>%
      filter(season == this.test & is.na(test_year) & pb %in% c("presence", "background")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background"))) %>%
      select(all_of(predictors), pb)
    
    #fix sea ice if all values = 0
    if(sum(test$sic) == 0){
      test <- test %>% select(-sic)
      sicpreds <- predictors[predictors != "sic"]
    } else {
      sicpreds <- predictors
    }
    
    #if test data is missing over 10% of a predictor, remove predictor
    pred_check <- test %>% select(all_of(sicpreds))
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    predictors2 <- names(pred_check)
    
    #BUFFERS
    # get buffer training data
    buff_train <- train %>% 
      filter(pb %in% c("buffer", "presence")) %>%
      select(all_of(predictors2), pb) %>%
      mutate(pb = ifelse(pb == "buffer", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #remove predictors if over 10% of column is NA
    pred_check <- buff_train %>% select(all_of(predictors2))
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    buff_predictors <- names(pred_check)
    
    #make testing data the same columns
    buff_testing <- test %>% select(all_of(buff_predictors))
    
    #calculate shape
    buff_shape <- extra_eval(training_data = buff_train, pr_ab = "pb", projection_data = buff_testing, n_cores = 40)
    
    #BACKGROUND
    # get background training data
    back_train <- train %>% 
      filter(pb == "background" & test_year == this.test | pb == "presence") %>%
      select(all_of(predictors2), pb) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #remove predictors if over 10% of column is NA
    pred_check <- back_train %>% select(all_of(predictors2))
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    back_predictors <- names(pred_check)
    
    #make testing data the same columns
    back_testing <- test %>% select(all_of(back_predictors))
    
    #calculate shape
    back_shape <- extra_eval(training_data = back_train, pr_ab = "pb", projection_data = back_testing, n_cores = 40)
    
    
    #CRW
    # remove trip IDs
    crw_data <- train %>%
      filter(pb == "crw")
    crw_data$individual_id <- stringr::str_remove(crw_data$individual_id, "_\\d+$")
    
    # get crw training data
    crw_train <- train %>% 
      filter(pb == "presence") %>%
      bind_rows(crw_data) %>%
      select(all_of(predictors2), pb, individual_id) %>%
      mutate(pb = ifelse(pb == "crw", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #remove predictors if over 10% of column is NA
    pred_check <- crw_train %>% select(all_of(predictors2))
    if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    crw_predictors <- names(pred_check)
    
    #make testing data the same columns
    crw_testing <- test %>% select(all_of(crw_predictors))
    
    #calculate shape
    crw_shape <- extra_eval(training_data = crw_train, pr_ab = "pb", projection_data = crw_testing, n_cores = 40)
    
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
    
    # print season
    print(i)
  }
  
  #export dataset
  saveRDS(shape_values,
          file = paste0("output/temporal/extrapolation/", this.species, "_", this.site, "_", this.stage, "_shape_medians.RDS"))
  saveRDS(shape_all,
          file = paste0("output/temporal/extrapolation/", this.species, "_", this.site, "_", this.stage, "_shape_scores.RDS"))
  
  #print to show that this species has completed
  print(paste0(this.species, " ", this.site, " ", this.stage, " completed"))
  
}