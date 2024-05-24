#automated script to run GAMs for leave-year-out validation
rm(list=ls())
setwd("/mainfs/home/jcw2g17/Chapter 01/")

{
  library(dplyr)
  library(mgcv)
  library(enmSdmX)
  library(lubridate)
  library(foreach)
  library(doParallel)
}

#GAM function to only use predictors that remain after later steps remove some
pred_gam <- function(df){
  mgcv::gam(
    as.formula(
    paste0(
      "pa ~ s(", 
      setdiff(names(df), "pa") %>% paste0(collapse = ", bs = 'ts', k=5) + s("),
      ", bs = 'ts', k=5)"
    )), 
    family=binomial, data=df)
}

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata_GAMs.csv")

#define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")

#setup parallel programming
registerDoParallel(cores = 21)

#loop to run through each species, stage, and site iteratively
foreach(z=1:21) %dopar% {
  try({
    
    #define parameters in loop
    rm(list=setdiff(ls(), c("meta", "predictors", "z", "pred_gam")))
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
    
    #null tables for output
    gam_boyce_final <- NULL
    
    #loop over each season
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
      
      
      # 3. Run GAMs
      
      #BUFFER
      #remove non-predictor columns
      buff_sel <- buff_train %>% select(all_of(predictors), pa)
      
      #remove SIC if less than 5 distinct values
      if(n_distinct(buff_sel$sic) < 5){
        buff_sel <- buff_sel %>% select(-sic)
      }
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(buff_sel)) > 0.1*nrow(buff_sel)){
        buff_sel <- buff_sel[colSums(is.na(buff_sel)) < 0.1*nrow(buff_sel)]
      }
      
      #run gam
      buff_gam <- pred_gam(buff_sel)
      
      #predict and evaluate
      p1 <- predict.gam(buff_gam, tracks_test, type = "response")
      p2 <- predict.gam(buff_gam, back_test, type = "response")
      buff_gam_boyce <- evalContBoyce(p1, p2, na.rm=TRUE)
      
      #save model
      saveRDS(buff_gam, 
              file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/buff_gam_", this.test, ".RDS"))
      
      #remove unnecessary parameters to continue
      rm(buff_gam, buff_sel, p1, p2)
      
      
      #BACKGROUND
      #remove non-predictor columns
      back_sel <- back_train %>% select(all_of(predictors), pa)
      
      #remove SIC if less than 5 distinct values
      if(n_distinct(back_sel$sic) < 5){
        back_sel <- back_sel %>% select(-sic)
      }
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(back_sel)) > 0.1*nrow(back_sel)){
        back_sel <- back_sel[colSums(is.na(back_sel)) < 0.1*nrow(back_sel)]
      }
      
      #run gam
      back_gam <- pred_gam(back_sel)
      
      #predict and evaluate
      p1 <- predict.gam(back_gam, tracks_test, type = "response")
      p2 <- predict.gam(back_gam, back_test, type = "response")
      back_gam_boyce <- evalContBoyce(p1, p2, na.rm=TRUE)
      
      #save model
      saveRDS(back_gam, 
              file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/back_gam_", this.test, ".RDS"))
      
      #remove unnecessary parameters to continue
      rm(back_gam, back_sel, p1, p2)
      
      
      #CRWs
      #remove non-predictor columns
      crw_sel <- crw_train %>% select(all_of(predictors), pa)
      
      #remove SIC if less than 5 distinct values
      if(n_distinct(crw_sel$sic) < 5){
        crw_sel <- crw_sel %>% select(-sic)
      }
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(crw_sel)) > 0.1*nrow(crw_sel)){
        crw_sel <- crw_sel[colSums(is.na(crw_sel)) < 0.1*nrow(crw_sel)]
      }
      
      #run gam
      crw_gam <- pred_gam(crw_sel)
      
      #predict and evaluate
      p1 <- predict.gam(crw_gam, tracks_test, type = "response")
      p2 <- predict.gam(crw_gam, back_test, type = "response")
      crw_gam_boyce <- evalContBoyce(p1, p2, na.rm=TRUE)
      
      #save model
      saveRDS(crw_gam, 
              file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/crw_gam_", this.test, ".RDS"))
      
      #remove unnecessary parameters to continue
      rm(crw_gam, crw_sel, p1, p2)
      
      #FINAL DATA
      gam_boyce <- expand.grid(buff = buff_gam_boyce, back = back_gam_boyce, crw = crw_gam_boyce)
      gam_boyce$season <- i
      gam_boyce_final <- rbind(gam_boyce_final, gam_boyce)
    }
    
    saveRDS(gam_boyce_final, 
            file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/boyce_scores_gam.RDS"))
    
    
  })
  
}