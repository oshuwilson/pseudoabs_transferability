#automated script to run Boosted Regression Trees for leave-year-out validation
rm(list=ls())
#setwd("/mainfs/home/jcw2g17/Chapter 01/")
setwd("/mainfs/home/jcw2g17/Chapter 01/")

{
  library(dplyr)
  library(caret)
  library(gbm)
  library(enmSdmX)
  library(lubridate)
  library(foreach)
  library(doParallel)
}

meta <- read.csv("data/species_site_stage_metadata_BRT.csv")
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")

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
    
    #null tables for output
    gbm_hypers <- NULL
    gbm_boyce_final <- NULL
    gbm_hyper_meta <- NULL
    
    #loop starts here
    for(i in seasons){
      this.test <- i
      
      #extract training data
      buff_train <- buff %>% filter(season != this.test)
      back_train <- back %>% filter(season != this.test)
      crw_train <- crw %>% filter(season != this.test)
      tracks_train <- tracks %>% filter(season != this.test)
      
      #extract testing data
      back_test <- back %>% filter(season == this.test, pa == 0)
      tracks_test <- tracks %>% filter(season == this.test)
      
      #if test data is missing over 10% of a predictor, remove predictor from models
      pred_check <- back %>% filter(season == this.test) %>% select(all_of(predictors))
      
      if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
        pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
      }
      
      predictors <- names(pred_check)
      
      
      # 3. gbm Predictions - Tune parameters
      
      #create parameter grid to vary hyperparameters
      param_grid <- expand.grid(interaction.depth = c(1, 3, 5),
                                n.trees = seq(1000, 10000, 1000),
                                shrinkage = c(0.005, 0.01, 0.5),
                                n.minobsinnode = 20)
      
      #setup 10-fold cross-validation
      cv_scheme <- trainControl(method = "cv", number = 10, verboseIter = FALSE,
                                summaryFunction = twoClassSummary, classProbs = TRUE)
      
      #remove NAs from test data
      back_test <- back_test %>% select(all_of(predictors)) %>% na.omit()
      tracks_test <- tracks_test %>% select(all_of(predictors)) %>% na.omit()
      
      
      #BUFFERS
      #remove non-predictor columns
      buff_sel <- buff_train %>% select(all_of(predictors), pa)
      
      #make presence-absence a character name
      buff_sel$pa <- if_else(buff_sel$pa == 1, "presence", "absence")
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(buff_sel)) > 0.1*nrow(buff_sel)){
        buff_sel <- buff_sel[colSums(is.na(buff_sel)) < 0.1*nrow(buff_sel)]
      }
      
      #perform tuning search
      X <- buff_sel %>% select(-pa)
      Y <- as.factor(buff_sel$pa)
      buff_gbm <- train(x = X, y = Y, method = "gbm", metric = "ROC", trControl = cv_scheme, 
                        tuneGrid = param_grid)
      
      #save parameter results
      buff_params <- buff_gbm$bestTune
      buff_params$pseudo <- "buffer"
      buff_params$season <- i
      
      #predict and evaluate
      p1 <- predict(buff_gbm, tracks_test, type = "prob")[,2]
      p2 <- predict(buff_gbm, back_test, type = "prob")[,2]
      buff_gbm_boyce <- evalContBoyce(p1, p2)
      
      #save model
      saveRDS(buff_gbm, 
              file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/buff_gbm_", this.test, ".RDS"))
      
      
      #remove unnecessary parameters to continue
      rm(buff_gbm,buff_sel, X, Y, p1, p2)
      
      
      #BACKGROUND
      #remove non-predictor columns
      back_sel <- back_train %>% select(all_of(predictors), pa)
      
      #make presence-absence a character name
      back_sel$pa <- if_else(back_sel$pa == 1, "presence", "absence")
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(back_sel)) > 0.1*nrow(back_sel)){
        back_sel <- back_sel[colSums(is.na(back_sel)) < 0.1*nrow(back_sel)]
      }
      
      #perform tuning search
      X <- back_sel %>% select(-pa)
      Y <- as.factor(back_sel$pa)
      back_gbm <- train(x = X, y = Y, method = "gbm", metric = "ROC", trControl = cv_scheme, 
                        tuneGrid = param_grid)
      
      #save parameter results
      back_params <- back_gbm$bestTune
      back_params$pseudo <- "backer"
      back_params$season <- i
      
      #predict and evaluate
      p1 <- predict(back_gbm, tracks_test, type = "prob")[,2]
      p2 <- predict(back_gbm, back_test, type = "prob")[,2]
      back_gbm_boyce <- evalContBoyce(p1, p2)
      
      #save model
      saveRDS(back_gbm, 
              file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/back_gbm_", this.test, ".RDS"))
      
      
      #remove unnecessary parameters to continue
      rm(back_gbm,back_sel, X, Y, p1, p2)
      
      
      #CRWs
      #remove non-predictor columns
      crw_sel <- crw_train %>% select(all_of(predictors), pa)
      
      #make presence-absence a character name
      crw_sel$pa <- if_else(crw_sel$pa == 1, "presence", "absence")
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(crw_sel)) > 0.1*nrow(crw_sel)){
        crw_sel <- crw_sel[colSums(is.na(crw_sel)) < 0.1*nrow(crw_sel)]
      }
      
      #perform tuning search
      X <- crw_sel %>% select(-pa)
      Y <- as.factor(crw_sel$pa)
      crw_gbm <- train(x = X, y = Y, method = "gbm", metric = "ROC", trControl = cv_scheme, 
                       tuneGrid = param_grid)
      
      #save parameter results
      crw_params <- crw_gbm$bestTune
      crw_params$pseudo <- "crwer"
      crw_params$season <- i
      
      #predict and evaluate
      p1 <- predict(crw_gbm, tracks_test, type = "prob")[,2]
      p2 <- predict(crw_gbm, back_test, type = "prob")[,2]
      crw_gbm_boyce <- evalContBoyce(p1, p2)
      
      #save model
      saveRDS(crw_gbm, 
              file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/crw_gbm_", this.test, ".RDS"))
      
      
      #remove unnecessary parameters to continue
      rm(crw_gbm,crw_sel, X, Y, p1, p2)
      
      
      #FINAL DATA
      gbm_boyce <- expand.grid(buff = buff_gbm_boyce, back = back_gbm_boyce, crw = crw_gbm_boyce)
      gbm_boyce$season <- i
      gbm_boyce_final <- rbind(gbm_boyce_final, gbm_boyce)
      
      hyper_values <- rbind(buff_params, back_params, crw_params)
      gbm_hyper_meta <- rbind(gbm_hyper_meta, hyper_values)
      
    }
    
    
    # 4. Export Boyce, Mtry, and Metadata
    saveRDS(gbm_boyce_final, 
            file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/boyce_scores_gbm.RDS"))
    saveRDS(gbm_hyper_meta, 
            file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/brt_hyperparameter_values.RDS"))
    
    print(paste0(this.species, " ", this.site, " ", this.stage, " success"))
    
  }) 
  
}