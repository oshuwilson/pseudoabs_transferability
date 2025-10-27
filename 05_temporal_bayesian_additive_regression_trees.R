#-------------------------------------------------------------------------------
# Temporal Bayesian Additive Regression Trees
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/dsdm/")

{
  library(dplyr)
  library(lubridate)
  library(tidymodels)
  library(themis)
  library(tidysdm)
  library(miceRanger)
  library(future)
}

# set number of cores for parallelisation
cores <- 78

# read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

# remove completed colonies
#meta <- meta[-c(1:7, 9, 10, 12, 15),]

# run over each row of metadata
for(z in 3:4){
  
  # define initial predictors
  predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "predictors", "z", "cores")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  #-------------------------------------------------------------------------------
  # 1. Formatting
  #-------------------------------------------------------------------------------
  
  #read in presences and pseudo-absences with environmental data
  data <- readRDS(paste0("output/extraction/", this.species, "_", this.site, "_", this.stage, "_extracted.rds"))
  
  #format data
  data <- data %>%
    mutate(date = as_date(date),
           pb = as.factor(pb)) 
  
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
  
  # if ANFS, remove 2008 + 2009 (too few data for modelling to converge)
  if(this.species == "ANFS"){
    seasons <- seasons[!seasons %in% c(2008, 2009)]
  }
  
  #null tables for output
  bart_boyce_final <- NULL
  all_metrics <- NULL
  
  # loop over seasons
  for(i in seasons){
    this.test <- i
    
    # define initial predictors
    predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")
    
    # if subantarctic fur seals, remove sea ice
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
    
    #-------------------------------------------------------------------------------
    # 2. Hyperparameter Tuning
    #-------------------------------------------------------------------------------
    
    # set ideal number of cross-validation folds to 5
    v <- 5
    
    #if number of individuals is fewer than 5, change v to n_ind
    n_ind <- train %>% filter(pb == "presence") %>% pull(individual_id) %>% unique() %>% length()
    if(n_ind < 5){
      v <- n_ind
    }
    
    #define bayesian additive regression tree settings
    bart_mod <- parsnip::bart() %>%
      set_mode("classification") %>%
      set_engine("dbarts") %>%
      set_args(trees = tune()) #tune trees
    
    #create workflow
    bart_wf <- workflow() %>%
      add_model(bart_mod)
    
    #define tree values to vary over 
    trees <- c(50, 100, 200, 300)
    grid <- expand_grid(trees = trees)
    
    #BUFFERS
    #isolate buffer data
    buff_sel <- train %>% 
      filter(pb %in% c("buffer", "presence")) %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ifelse(pb == "buffer", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #check for NA - less than 10% of training data okay for imputing
    if(sum(is.na(buff_sel)) < 0.1*nrow(buff_sel) & sum(is.na(buff_sel)) > 0){
      buff_mice <- miceRanger(buff_sel, m=1)
      buff_sel <- completeData(buff_mice)[[1]]
    }
    
    #remove columns where missing data is over 10% of rows then impute
    if(sum(is.na(buff_sel)) > 0.1*nrow(buff_sel)){
      buff_sel <- buff_sel[colSums(is.na(buff_sel)) < 0.1*nrow(buff_sel)]
      if(sum(is.na(buff_sel)) > 0){
        buff_mice <- miceRanger(buff_sel, m=1)
        buff_sel <- completeData(buff_mice)[[1]]
      }
    }
    
    #create cross-validation folds
    buff_folds <- group_vfold_cv(data = buff_sel, 
                                 group = individual_id, #split training/testing data by individual
                                 v = v, #number of folds
                                 balance = "observations" #roughly the same number of points in each fold
    )
    
    #define formula for modelling
    buff_rec <- recipe(pb ~ ., data = buff_sel) %>%
      update_role(individual_id, new_role = "ID")  %>% #let model know that id is not a predictor
      step_downsample(pb)
    
    #update workflow
    buff_wf <- bart_wf %>%
      add_recipe(buff_rec)
    
    # enable parallelisation
    plan(multisession, workers = cores)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    buff_tun <- tune_grid(buff_wf,
                          resamples = buff_folds,
                          grid = grid,
                          metrics = sdm_metric_set()) #includes boyce index as a tuning parameter
    
    # end parallelisation
    plan(sequential)
    
    #get metric scores
    buff_metrics <- collect_metrics(buff_tun, summarize = T) %>%
      mutate(pb = "buffer")
    
    #extract best model
    buff_best <- select_best(buff_tun, metric = "boyce_cont")
    
    #set up model
    buff_best_mod <- bart() %>%
      set_engine(engine = "dbarts") %>%
      set_mode("classification") %>%
      set_args(trees = buff_best$trees[1])
    
    #update workflow
    buff_best_wf <- buff_wf %>%
      update_model(buff_best_mod)
    
    #run best model on all data
    buff_fit <- buff_best_wf %>%
      fit(buff_sel)
    
    #predict habitat suitability to training data
    test$buff_suitability <- predict(buff_fit, test, type = "prob") %>%
      pull(.pred_presence)
    
    #calculate boyce index
    buff_boyce <- boyce_cont(test, pb, buff_suitability) %>%
      pull(.estimate)
    
    #remove unnecessary parameters to continue
    rm(buff_sel, buff_folds, buff_rec, buff_wf, buff_grid, buff_tun, buff_best, buff_best_mod, buff_best_wf)
    
    
    #BACKGROUND
    #remove non-predictor columns
    back_sel <- train %>% 
      filter(pb == "background" & test_year == this.test | pb == "presence") %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #check for NA - less than 10% of training data okay for imputing
    if(sum(is.na(back_sel)) < 0.1*nrow(back_sel) & sum(is.na(back_sel)) > 0){
      back_mice <- miceRanger(back_sel, m=1)
      back_sel <- completeData(back_mice)[[1]]
    }
    
    #remove columns where missing data is over 10% of rows then impute
    if(sum(is.na(back_sel)) > 0.1*nrow(back_sel)){
      back_sel <- back_sel[colSums(is.na(back_sel)) < 0.1*nrow(back_sel)]
      if(sum(is.na(back_sel)) > 0){
        back_mice <- miceRanger(back_sel, m=1)
        back_sel <- completeData(back_mice)[[1]]
      }
    }
    
    #create cross-validation folds
    back_folds <- group_vfold_cv(data = back_sel, 
                                 group = individual_id, #split training/testing data by individual
                                 v = v, #number of folds
                                 balance = "observations" #roughly the same number of points in each fold
    )
    
    #define formula for modelling
    back_rec <- recipe(pb ~ ., data = back_sel) %>%
      update_role(individual_id, new_role = "ID") %>% #let model know that id is not a predictor
      step_downsample(pb)
    
    #update workflow
    back_wf <- bart_wf %>%
      add_recipe(back_rec)
    
    # enable parallelisation
    plan(multisession, workers = cores)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    back_tun <- tune_grid(back_wf,
                          resamples = back_folds,
                          grid = grid,
                          metrics = sdm_metric_set()) #includes boyce index as a tuning parameter
    
    # end parallelisation
    plan(sequential)
    
    #get metric scores
    back_metrics <- collect_metrics(back_tun, summarize = T) %>%
      mutate(pb = "background")
    
    #extract best model
    back_best <- select_best(back_tun, metric = "boyce_cont")
    
    #set up model
    back_best_mod <- bart() %>%
      set_engine(engine = "dbarts") %>%
      set_mode("classification") %>%
      set_args(trees = back_best$trees[1])
    
    #update workflow
    back_best_wf <- back_wf %>%
      update_model(back_best_mod)
    
    #run best model on all data
    back_fit <- back_best_wf %>%
      fit(back_sel)
    
    #predict habitat suitability to training data
    test$back_suitability <- predict(back_fit, test, type = "prob") %>%
      pull(.pred_presence)
    
    #calculate boyce index
    back_boyce <- boyce_cont(test, pb, back_suitability) %>%
      pull(.estimate)
    
    #remove unnecessary parameters to continue
    rm(back_folds, back_rec, back_wf, back_tun, back_best, back_best_mod, back_best_wf)
    
    
    #CRWs
    #get CRW data and remove trip numbers
    crw_data <- train %>%
      filter(pb == "crw")
    crw_data$individual_id <- stringr::str_remove(crw_data$individual_id, "_\\d+$")
    
    # isolate CRW data
    crw_sel <- train %>% 
      filter(pb == "presence") %>%
      bind_rows(crw_data) %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ifelse(pb == "crw", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #check for NA - less than 10% of training data okay for imputing
    if(sum(is.na(crw_sel)) < 0.1*nrow(crw_sel) & sum(is.na(crw_sel)) > 0){
      crw_mice <- miceRanger(crw_sel, m=1)
      crw_sel <- completeData(crw_mice)[[1]]
    }
    
    #remove columns where missing data is over 10% of rows then impute
    if(sum(is.na(crw_sel)) > 0.1*nrow(crw_sel)){
      crw_sel <- crw_sel[colSums(is.na(crw_sel)) < 0.1*nrow(crw_sel)]
      if(sum(is.na(crw_sel)) > 0){
        crw_mice <- miceRanger(crw_sel, m=1)
        crw_sel <- completeData(crw_mice)[[1]]
      }
    }
    
    #create cross-validation folds
    crw_folds <- group_vfold_cv(data = crw_sel, 
                                group = individual_id, #split training/testing data by individual
                                v = v, #number of folds
                                balance = "observations" #roughly the same number of points in each fold
    )
    
    #define formula for modelling
    crw_rec <- recipe(pb ~ ., data = crw_sel) %>%
      update_role(individual_id, new_role = "ID") %>% #let model know that id is not a predictor
      step_downsample(pb)
    
    #update workflow
    crw_wf <- bart_wf %>%
      add_recipe(crw_rec)
    
    # enable parallelisation
    plan(multisession, workers = cores)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    crw_tun <- tune_grid(crw_wf,
                         resamples = crw_folds,
                         grid = grid,
                         metrics = sdm_metric_set()) #includes boyce index as a tuning parameter
    
    # end parallelisation
    plan(sequential)
    
    #get metric scores
    crw_metrics <- collect_metrics(crw_tun, summarize = T) %>%
      mutate(pb = "crw")
    
    #extract best model
    crw_best <- select_best(crw_tun, metric = "boyce_cont")
    
    #set up model
    crw_best_mod <- bart() %>%
      set_engine(engine = "dbarts") %>%
      set_mode("classification") %>%
      set_args(trees = crw_best$trees[1])
    
    #update workflow
    crw_best_wf <- crw_wf %>%
      update_model(crw_best_mod)
    
    #run best model on all data
    crw_fit <- crw_best_wf %>%
      fit(crw_sel)
    
    #predict habitat suitability to training data
    test$crw_suitability <- predict(crw_fit, test, type = "prob") %>%
      pull(.pred_presence)
    
    #calculate boyce index
    crw_boyce <- boyce_cont(test, pb, crw_suitability) %>%
      pull(.estimate)
    
    #combine boyce scores
    bart_boyce <- expand.grid(buff = buff_boyce, back = back_boyce, crw = crw_boyce)
    bart_boyce$season <- i
    bart_boyce_final <- rbind(bart_boyce_final, bart_boyce)
    
    # self test models
    self_test$buff_suitability <- predict(buff_fit, self_test, type = "prob") %>%
      pull(.pred_presence)
    self_test$back_suitability <- predict(back_fit, self_test, type = "prob") %>%
      pull(.pred_presence)
    self_test$crw_suitability <- predict(crw_fit, self_test, type = "prob") %>%
      pull(.pred_presence)
    
    buff_boyce0 <- boyce_cont(self_test, pb, buff_suitability) %>%
      pull(.estimate)
    back_boyce0 <- boyce_cont(self_test, pb, back_suitability) %>%
      pull(.estimate)
    crw_boyce0 <- boyce_cont(self_test, pb, crw_suitability) %>%
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
    
    # combine metrics
    metrics <- bind_rows(buff_metrics, back_metrics, crw_metrics) %>%
      mutate(season = i)
    all_metrics <- bind_rows(all_metrics, metrics)
  }
  
  #-------------------------------------------------------------------------------
  # 3. Export
  #-------------------------------------------------------------------------------
  
  # export boyce scores
  saveRDS(bart_boyce_final, 
          file = paste0("output/temporal/bart/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_bart.RDS"))
  
  # export metrics
  saveRDS(all_metrics,
       file = paste0("output/temporal/bart/", this.species, "_", this.site, "_", this.stage, "_metrics_bart.RDS"))
  
  # export self test 
  saveRDS(all_self,
          file = paste0("output/temporal/bart/", this.species, "_", this.site, "_", this.stage, "_self_test_bart.RDS"))
  
  #show species has finished
  print(paste0(this.species, " ", this.site, " ", this.stage, " success"))
}
