#-------------------------------------------------------------------------------
# Temporal marginalised space-time point process model
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(INLA)
library(elevatr)
library(terra)
library(spdep)
library(spatialreg)
library(lubridate)
library(dplyr)
library(tidysdm)
library(ctmm)
library(terra)
library(tidyterra)
library(miceRanger)

# read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

# run over each row of metadata
for(z in 13){
  
  # define initial predictors
  predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")
  
  # define parameters in loop
  rm(list=setdiff(ls(), c("meta", "predictors", "z", "cores")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  # read in pseudo-absences with environmental data
  pseudo <- readRDS(paste0("output/extraction/", this.species, "_", this.site, "_", this.stage, "_extracted.rds")) %>%
    filter(pb != "presence")
  
  # read in presences with environmental data
  pres <- read.csv(paste0("output_old/extraction/", this.species, "/", this.site, "/", this.stage, "/presences.csv"))
  
  # combine presences and pesudoabsences
  data <- bind_rows(
    pseudo %>% select(pb, all_of(predictors), individual_id, date, x, y, test_year),
    pres %>% mutate(pb = "presence", date = as_datetime(date)) %>% 
      select(pb, all_of(predictors), individual_id, date, x, y)
  )
  
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
  
  # if ANFS post-moult, remove 2008 + 2009 (too few data for modelling to converge)
  if(this.species == "ANFS" & this.stage == "post-moult"){
    seasons <- seasons[!seasons %in% c(2008, 2009)]
  }
  
  #-------------------------------------------------------------------------------
  # 1. Creating kernel density for STPP
  #-------------------------------------------------------------------------------
  
  # get presences from the training data
  tracks <- pres %>%
      select(individual_id, date, x, y) %>%
      mutate(date = as_datetime(date))
  
  # rename columns to movebank format for ctmm
  tracks <- tracks %>%
    rename(individual.local.identifier = individual_id,
           timestamp = date,
           location.long = x,
           location.lat = y)
  
  # arrange by ID and date
  tracks <- tracks %>%
    arrange(individual.local.identifier, timestamp) %>%
    distinct()
  
  # make sure over 5 points per individual
  track_counts <- tracks %>%
    group_by(individual.local.identifier) %>%
    summarise(n = n()) %>%
    ungroup() %>% 
    filter(n < 5) %>%
    pull(individual.local.identifier)
  tracks <- tracks %>%
    filter(!individual.local.identifier %in% track_counts)
  
  # convert to telemetry object
  tel <- as.telemetry(tracks)
  
  # fit KDEs to each individual
  fits <- lapply(tel, function(x) {
    m.iid <- ctmm.fit(x)
  })
  
  # fit UDs - weighted to account for irregular sampling
  UDs <- akde(tel, fits, weights = T)
  
  #-------------------------------------------------------------------------------
  # 2. Run temporal leave-year-out tests
  #-------------------------------------------------------------------------------
  
  # null tables for output
  stpp_boyce_final <- NULL
  all_metrics <- NULL
  
  # loop over seasons
  for(i in seasons){
    this.test <- i
    print(paste0("Season started: ", this.test))
    
    # define initial predictors
    predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")
    
    # if SUFS, remove sea ice
    if(this.species == "SUFS"){
      predictors <- predictors[!predictors %in% c("sic")]
    }
    
    #extract training data
    train <- data %>%
      filter(season != this.test)
    
    # isolate self-testing dataset from within training data
    self_test_inds <- train %>%
      filter(pb == "presence") %>%
      pull(individual_id) %>%
      unique() %>%
      sample(size = length(.)*0.2)
    self_test <- train %>% 
      filter(pb == "background" & test_year == this.test | pb == "presence") %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background"))) %>%
      filter(individual_id %in% self_test_inds)
    
    # impute if any missing values
    if(sum(is.na(self_test)) > 0){
      self_mice <- miceRanger(self_test, m=1)
      self_test <- completeData(self_mice)[[1]]
    }
    
    # remove self-test data from training data
    train <- train %>%
      filter(!individual_id %in% self_test_inds)
    
    #create test dataset
    test <- data %>%
      filter(season == this.test & is.na(test_year) & pb %in% c("presence", "background")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #if test data is missing over 10% of a predictor, remove predictor from models
    pred_check <- test %>% 
      select(all_of(predictors))
    if(sum(is.na(pred_check)) > 0.01*nrow(pred_check)){
      pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
    }
    predictors <- names(pred_check)
    
    #remove any NAs from test data to allow prediction
    test <- test %>% 
      select(all_of(predictors), pb, individual_id) %>%
      drop_na()
    
    # if no test data skip
    if(nrow(test) == 0){
      next
    }
    
    # get a list of all individuals in the training data
    train_inds <- train %>%
      filter(pb == "presence") %>%
      pull(individual_id) %>%
      unique()
    
    # isolate the UDs of this season's training individuals
    UDs_season <- UDs[train_inds]
    
    # compute population UD
    popUD <- mean(UDs_season)
    
    # convert to SpatRaster
    popUD <- popUD %>% raster::raster() %>% rast()
    
    # project to epsg:4326
    popUD <- popUD %>% project("epsg:4326")
    
    # extract the values of the kernel density at the locations of the training data
    train_terra <- train %>%
      vect(geom = c("x", "y"), crs = "epsg:4326")
    train_terra$mov.kern <- terra::extract(popUD, train_terra, ID = F)
    
    # reconvert to dataframe
    train <- as.data.frame(train_terra, geom = "XY")
    
    # invert mov.kern values
    train$mov.kern <- 1-train$mov.kern
    
    # if NA, revalue to 0 (i.e. outside of the kernel)
    train$mov.kern[is.na(train$mov.kern)] <- 0
    
    
    #-------------------------------------------------------------------------------
    # 3. Fit STPP using INLA - following Eisaguirre et al. 2025
    #-------------------------------------------------------------------------------
    
    # get minimum number of psuedo-absences
    min_n <- train %>% filter(pb != "presence") %>%
      group_by(pb) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      pull(n) %>%
      min()
    
    # subsample all pseudoabsences
    pres <- train %>% filter(pb == "presence")
    train <- train %>%  filter(pb != "presence") %>%
      group_by(pb) %>%
      sample_n(min_n) %>%
      ungroup() %>%
      bind_rows(pres)
    
    # 1. Background Sampling
    bg <- train %>%
      filter(pb %in% c("background", "presence"))
    
    # convert pb to 1-0
    bg$use <- ifelse(bg$pb == "presence", 1, 0)
    
    # base model of covariates - use remaining predictors only
    b.f <- paste0("use ~ (", paste(predictors, collapse = " + "), ")")
    
    # get residuals for movement kernel
    mov.mod <- glm(update.formula(b.f,log(mov.kern+.000001)~.),data=bg,na.action='na.exclude')
    bg$mov.resid <- scale(residuals(mov.mod))[,1]
    bg$mov.resid.raw <- residuals(mov.mod)
    
    # add movement kernel
    f.inla <- update.formula(b.f,~.+mov.resid)
    
    # fit INLA model
    bg.sp <- inla(f.inla,    ## should take <1 min to run on most machines
                  control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                     mean=0,prec=1/1),
                  family='binomial',
                  data=bg,
                  control.compute=list(config=FALSE))
    
    # create mov.kern and mov.resid values of 0 on test data
    test$mov.kern <- 0
    test$mov.resid <- 0
    
    # reformat test data to have NA response variable
    test$use <- NA
    
    # append to bg training dataset
    df_bg <- bind_rows(
      bg %>% select(all_of(predictors), use, pb, mov.kern, mov.resid),
      test %>% select(all_of(predictors), use, pb, mov.kern, mov.resid)
    )
    
    # refit INLA model
    bg.sp2 <- inla(f.inla,    ## should take <1 min to run on most machines
                   control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                      mean=0,prec=1/1),
                   family='binomial',
                   data=df_bg,
                   control.compute=list(config=FALSE))
    
    # get predictions for new rows
    n_train <- nrow(bg)
    pred <- bg.sp2$summary.fitted.values$mean[(n_train+1):nrow(df_bg)]
    
    # back-transform predictions
    pred <- plogis(pred)
    
    # append predictions to test data
    test$pred <- pred
    
    # calculate boyce index
    back_boyce <- boyce_cont(test, pb, pred) %>%
      pull(.estimate)
    
    # repeat the process for the self-test individuals
    self_test$mov.kern <- 0
    self_test$mov.resid <- 0
    self_test$use <- NA
    
    df_self <- bind_rows(
      bg %>% select(all_of(predictors), use, pb, mov.kern, mov.resid),
      self_test %>% select(all_of(predictors), use, pb, mov.kern, mov.resid)
    )
    
    bg.sp3 <- inla(f.inla,    ## should take <1 min to run on most machines
                   control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                      mean=0,prec=1/1),
                   family='binomial',
                   data=df_self,
                   control.compute=list(config=FALSE))
    
    n_train_self <- nrow(bg)
    pred_self <- bg.sp3$summary.fitted.values$mean[(n_train_self+1):nrow(df_self)]
    pred_self <- plogis(pred_self)
    self_test$pred <- pred_self
    
    back_boyce0 <- boyce_cont(self_test, pb, pred) %>%
      pull(.estimate)
    
    
    # 2. Buffer Sampling
    buff <- train %>%
      filter(pb %in% c("buffer", "presence"))
    
    # convert pb to 1-0
    buff$use <- ifelse(buff$pb == "presence", 1, 0)
    
    # base model of covariates - use remaining predictors only
    b.f <- paste0("use ~ (", paste(predictors, collapse = " + "), ")")
    
    # get residuals for movement kernel
    mov.mod <- glm(update.formula(b.f,log(mov.kern+.000001)~.),data=buff,na.action='na.exclude')
    buff$mov.resid <- scale(residuals(mov.mod))[,1]
    buff$mov.resid.raw <- residuals(mov.mod)
    
    # add movement kernel
    f.inla <- update.formula(b.f,~.+mov.resid)
    
    # fit INLA model
    buff.sp <- inla(f.inla,    ## should take <1 min to run on most machines
                    control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                       mean=0,prec=1/1),
                    family='binomial',
                    data=buff,
                    control.compute=list(config=FALSE))
    
    # append to buff training dataset
    df_buff <- bind_rows(
      buff %>% select(all_of(predictors), use, pb, mov.kern, mov.resid),
      test %>% select(all_of(predictors), use, pb, mov.kern, mov.resid)
    )
    
    # refit INLA model
    buff.sp2 <- inla(f.inla,    ## should take <1 min to run on most machines
                     control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                        mean=0,prec=1/1),
                     family='binomial',
                     data=df_buff,
                     control.compute=list(config=FALSE))
    
    # get predictions for new rows
    n_train_buff <- nrow(buff)
    pred_buff <- buff.sp2$summary.fitted.values$mean[(n_train_buff+1):nrow(df_buff)]
    pred_buff <- plogis(pred_buff)
    test$pred_buff <- pred_buff
    
    # calculate boyce index
    buff_boyce <- boyce_cont(test, pb, pred_buff) %>%
      pull(.estimate)
    
    # repeat the process for the self-test individuals
    df_self <- bind_rows(
      buff %>% select(all_of(predictors), use, pb, mov.kern, mov.resid),
      self_test %>% select(all_of(predictors), use, pb, mov.kern, mov.resid)
    )
    
    buff.sp3 <- inla(f.inla,    ## should take <1 min to run on most machines
                     control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                        mean=0,prec=1/1),
                     family='binomial',
                     data=df_self,
                     control.compute=list(config=FALSE))
    
    n_train_self <- nrow(buff)
    pred_self_buff <- buff.sp3$summary.fitted.values$mean[(n_train_self+1):nrow(df_self)]
    pred_self_buff <- plogis(pred_self_buff)
    self_test$pred_buff <- pred_self_buff
    
    buff_boyce0 <- boyce_cont(self_test, pb, pred_buff) %>%
      pull(.estimate)
    
    
    # 3. CRW
    crw <- train %>%
      filter(pb %in% c("crw", "presence"))
    
    # convert pb to 1-0
    crw$use <- ifelse(crw$pb == "presence", 1, 0)
    
    # base model of covariates - use remaining predictors only
    b.f <- paste0("use ~ (", paste(predictors, collapse = " + "), ")")
    
    # get residuals for movement kernel
    mov.mod <- glm(update.formula(b.f,log(mov.kern+.000001)~.),data=crw,na.action='na.exclude')
    crw$mov.resid <- scale(residuals(mov.mod))[,1]
    crw$mov.resid.raw <- residuals(mov.mod)
    
    # add movement kernel
    f.inla <- update.formula(b.f,~.+mov.resid)
    
    # fit INLA model
    crw.sp <- inla(f.inla,    ## should take <1 min to run on most machines
                   control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                      mean=0,prec=1/1),
                   family='binomial',
                   data=crw,
                   control.compute=list(config=FALSE))
    
    # append to crw training dataset
    df_crw <- bind_rows(
      crw %>% select(all_of(predictors), use, pb, mov.kern, mov.resid),
      test %>% select(all_of(predictors), use, pb, mov.kern, mov.resid)
    )
    
    # refit INLA model
    crw.sp2 <- inla(f.inla,    ## should take <1 min to run on most machines
                    control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                       mean=0,prec=1/1),
                    family='binomial',
                    data=df_crw,
                    control.compute=list(config=FALSE))
    
    # get predictions for new rows
    n_train_crw <- nrow(crw)
    pred_crw <- crw.sp2$summary.fitted.values$mean[(n_train_crw+1):nrow(df_crw)]
    pred_crw <- plogis(pred_crw)
    test$pred_crw <- pred_crw
    
    # calculate boyce index
    crw_boyce <- boyce_cont(test, pb, pred_crw) %>%
      pull(.estimate)
    
    # repeat the process for the self-test individuals
    df_self <- bind_rows(
      crw %>% select(all_of(predictors), use, pb, mov.kern, mov.resid),
      self_test %>% select(all_of(predictors), use, pb, mov.kern, mov.resid)
    )
    
    crw.sp3 <- inla(f.inla,    ## should take <1 min to run on most machines
                    control.fixed=list(mean.intercept=0,prec.intercept=1/1000^2,
                                       mean=0,prec=1/1),
                    family='binomial',
                    data=df_self,
                    control.compute=list(config=FALSE))
    
    n_train_self <- nrow(crw)
    pred_self_crw <- crw.sp3$summary.fitted.values$mean[(n_train_self+1):nrow(df_self)]
    pred_self_crw <- plogis(pred_self_crw)
    self_test$pred_crw <- pred_self_crw
    
    crw_boyce0 <- boyce_cont(self_test, pb, pred_crw) %>%
      pull(.estimate)
    
    # collate self test scores
    self <- expand.grid(buff = buff_boyce0, back = back_boyce0, crw = crw_boyce0) %>%
      mutate(season = i)
    
    # join to other self test scores
    if(i == seasons[1]){
      all_self <- self
    } else {
      all_self <- rbind(all_self, self)
    }
    
    #combine boyce scores
    stpp_boyce <- expand.grid(buff = buff_boyce, back = back_boyce, crw = crw_boyce)
    stpp_boyce$season <- i
    stpp_boyce_final <- rbind(stpp_boyce_final, stpp_boyce)
    
    #print season completion
    print(i)
  }
  
  # 4. Export
  
  # export boyce scores
  saveRDS(stpp_boyce_final, 
          file = paste0("output/temporal/stpp/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_stpp.RDS"))
  
  # export self test 
  saveRDS(all_self,
          file = paste0("output/temporal/stpp/", this.species, "_", this.site, "_", this.stage, "_self_test_stpp.RDS"))
  
  #show species has finished
  print(paste0(this.species, " ", this.site, " ", this.stage, " success"))
  
}
