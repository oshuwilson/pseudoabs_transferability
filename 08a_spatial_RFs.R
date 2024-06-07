#automated script to run RFs for spatial transfer validation
rm(list=ls())
#setwd("/mainfs/home/jcw2g17/Chapter 01/")
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(dplyr)
  library(caret)
  library(ranger)
  library(enmSdmX)
  library(lubridate)
  library(miceRanger)
}

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")
meta2 <- read.csv("output/spatial/spatial_site_metadata.csv")

#define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")

#isolate subsets where all predictors are present in test data - do the same for those missing chl and/or wind and change predictors
meta <- meta[c(1, 3, 9, 10, 11, 13, 15:19),]
meta2 <- meta2[c(4, 8, 9, 16:21, 25),]

#loop to run through each species, stage, and site iteratively
for(z in 1:11) {
  try({
    
    #define parameters in loop
    rm(list=setdiff(ls(), c("meta", "meta2", "predictors", "z")))
    this.species <- meta[z, 1]
    this.site <- meta[z, 2]
    this.stage <- meta[z, 3]
    
    # 1. Formatting 
    #read in presences and pseudo-absences
    tracks <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/presences.csv"))
    buff <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/buffers.csv"))
    back <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/background.csv"))
    crw <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/CRWs.csv"))
    
    #create presence/absence column
    tracks$pa <- "presence"
    buff$pa <- "absence"
    back$pa <- "absence"
    crw$pa <- "absence"
    
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
    
    # 2. Create RFs
    
    #create parameter grid to vary mtry between 2, 3, and 4
    param_grid <- expand.grid(mtry=2:4, splitrule = "gini", min.node.size=1)
    
    #setup 10-fold cross-validation
    cv_scheme <- trainControl(method = "cv", number = 10, verboseIter = FALSE,
                              summaryFunction = twoClassSummary, classProbs = TRUE)
    
    #BUFFER
    #remove non-predictor columns
    buff_sel <- buff %>% select(all_of(predictors), pa)
    
    #check for NA - less than 10% of training data okay for imputing
    if(sum(is.na(buff_sel)) < 0.1*nrow(buff_sel) & sum(is.na(buff_sel)) > 0){
      buff_mice <- miceRanger(buff_sel, m=1)
      buff_sel <- completeData(buff_mice)[[1]]
    }
    
    #remove columns where missing data is over 10% of rows then impute
    if(sum(is.na(buff_sel)) > 0.1*nrow(buff_sel)){
      buff_sel <- buff_sel[colSums(is.na(buff_sel)) < 0.1*nrow(buff_sel)]
      buff_mice <- miceRanger(buff_sel, m=1)
      buff_sel <- completeData(buff_mice)[[1]]
    }
    
    #perform tuning search
    X <- buff_sel %>% select(-pa)
    Y <- as.factor(buff_sel$pa)
    buff_rf <- train(x = X, y = Y, method = "ranger", metric = "ROC", trControl = cv_scheme, 
                     tuneGrid = param_grid, num.trees = 1000, importance = "impurity")
    
    #save mtry results
    buff_mtry <- buff_rf$results[,c(1,4)]
    buff_mtry$pseudo <- "buff"
    
    #save model
    saveRDS(buff_rf, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/buff_rf.RDS"))
    
    #remove unnecessary parameters to continue
    rm(buff_sel, buff_mice, X, Y)
    
    
    #BACKGROUND
    #remove non-predictor columns
    back_sel <- back %>% select(all_of(predictors), pa)
    
    #check for NA - less than 10% of training data okay for imputing
    if(sum(is.na(back_sel)) < 0.1*nrow(back_sel) & sum(is.na(back_sel)) > 0){
      back_mice <- miceRanger(back_sel, m=1)
      back_sel <- completeData(back_mice)[[1]]
    }
    
    #remove columns where missing data is over 10% of rows then impute
    if(sum(is.na(back_sel)) > 0.1*nrow(back_sel)){
      back_sel <- back_sel[colSums(is.na(back_sel)) < 0.1*nrow(back_sel)]
      back_mice <- miceRanger(back_sel, m=1)
      back_sel <- completeData(back_mice)[[1]]
    }
    
    #perform tuning search
    X <- back_sel %>% select(-pa)
    Y <- as.factor(back_sel$pa)
    back_rf <- train(x = X, y = Y, method = "ranger", metric = "ROC", trControl = cv_scheme, 
                     tuneGrid = param_grid, num.trees = 1000, importance = "impurity")
    
    #save mtry results
    back_mtry <- back_rf$results[,c(1,4)]
    back_mtry$pseudo <- "back"
    
    #save model
    saveRDS(back_rf, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/back_rf.RDS"))
    
    #remove unnecessary parameters to continue
    rm(back_sel, back_mice, X, Y)
    
    
    #CRWs
    #remove non-predictor columns
    crw_sel <- crw %>% select(all_of(predictors), pa)
    
    #check for NA - less than 10% of training data okay for imputing
    if(sum(is.na(crw_sel)) < 0.1*nrow(crw_sel) & sum(is.na(crw_sel)) > 0){
      crw_mice <- miceRanger(crw_sel, m=1)
      crw_sel <- completeData(crw_mice)[[1]]
    }
    
    #remove columns where missing data is over 10% of rows then impute
    if(sum(is.na(crw_sel)) > 0.1*nrow(crw_sel)){
      crw_sel <- crw_sel[colSums(is.na(crw_sel)) < 0.1*nrow(crw_sel)]
      crw_mice <- miceRanger(crw_sel, m=1)
      crw_sel <- completeData(crw_mice)[[1]]
    }
    
    #perform tuning search
    X <- crw_sel %>% select(-pa)
    Y <- as.factor(crw_sel$pa)
    crw_rf <- train(x = X, y = Y, method = "ranger", metric = "ROC", trControl = cv_scheme, 
                     tuneGrid = param_grid, num.trees = 1000, importance = "impurity")
    
    #save mtry results
    crw_mtry <- crw_rf$results[,c(1,4)]
    crw_mtry$pseudo <- "crw"
    
    #save model
    saveRDS(crw_rf, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/crw_rf.RDS"))
    
    #remove unnecessary parameters to continue
    rm(crw_sel, crw_mice, X, Y)
    
    
    #EXPORT
    #mtry values
    mtry_values <- rbind(buff_mtry, back_mtry, crw_mtry)
    saveRDS(mtry_values, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/mtry_values.RDS"))
    
    
    # 3. Test RFs
    
    #filter spatial metadata to this species and stage
    meta3 <- meta2 %>% filter(Species == this.species & Stage == this.stage)
    
    #extract list of sites for this species and stage
    meta3$Site <- as.factor(meta3$Site)
    sites <- levels(meta3$Site)
    
    #null table for output
    rf_boyce_final <- NULL
    
    #run for loop to test each site
    for(i in sites){
      test.site <- i
      
      #load in test data
      back_test <- read.csv(paste0("output/spatial/", this.species, "/", this.stage, "/extraction/", test.site, "_background.csv"))
      tracks_test <- read.csv(paste0("output/spatial/", this.species, "/", this.stage, "/extraction/", test.site, "_presences.csv"))
      
      #only select predictors for testing
      back_test <- back_test %>% select(all_of(predictors))
      tracks_test <- tracks_test %>% select(all_of(predictors))
      
      #check for NA - less than 10% of training data okay for imputing
      if(sum(is.na(back_test)) < 0.1*nrow(back_test) & sum(is.na(back_test)) > 0){
        back_test_mice <- miceRanger(back_test, m=1)
        back_test <- completeData(back_test_mice)[[1]]
      }

      if(sum(is.na(tracks_test)) < 0.1*nrow(tracks_test) & sum(is.na(tracks_test)) > 0){
        tracks_test_mice <- miceRanger(tracks_test, m=1)
        tracks_test <- completeData(tracks_test_mice)[[1]]
      }
      
      #predict and evaluate buffers
      p1 <- predict(buff_rf, tracks_test, type = "prob")[,2]
      p2 <- predict(buff_rf, back_test, type = "prob")[,2]
      buff_rf_boyce <- evalContBoyce(p1, p2)
      
      #predict and evaluate background
      p1 <- predict(back_rf, tracks_test, type = "prob")[,2]
      p2 <- predict(back_rf, back_test, type = "prob")[,2]
      back_rf_boyce <- evalContBoyce(p1, p2)
      
      #predict and evaluate crws
      p1 <- predict(crw_rf, tracks_test, type = "prob")[,2]
      p2 <- predict(crw_rf, back_test, type = "prob")[,2]
      crw_rf_boyce <- evalContBoyce(p1, p2)
      
      #FINAL DATA
      #boyce scores
      rf_boyce <- expand.grid(buff = buff_rf_boyce, back = back_rf_boyce, crw = crw_rf_boyce)
      rf_boyce$site <- i
      rf_boyce_final <- rbind(rf_boyce_final, rf_boyce)
    }
    
    #export boyce scores
    saveRDS(rf_boyce_final, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/boyce_scores_rf.RDS"))
    
  })
  
}