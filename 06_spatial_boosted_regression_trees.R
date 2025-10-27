#-------------------------------------------------------------------------------
# Spatial Boosted Regression Trees
#-------------------------------------------------------------------------------

#clear workspace and set working directory
rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/dsdm/")

#load required packages
{
  library(dplyr)
  library(tidymodels)
  library(tidysdm)
  library(themis)
  library(bonsai)
  library(lubridate)
  library(future)
  library(miceRanger)
}

#cores for parallelisation
cores <- 78

#set seed
set.seed(777)

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")
meta2 <- read.csv("data/spatial_site_metadata.csv")

#remove ANPE, EMPE, SUFS, and MAPE incubation/post-breeding (no spatial transfer)
meta <- meta %>% filter(Species != "ANPE" & Species != "EMPE" & Species != "SUFS")
meta <- meta %>% filter(Species != "MAPE" | 
                          (Species == "MAPE" & Stage != "incubation" & Stage != "post-breeding"))

for(z in 1:nrow(meta)){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "meta2", "predictors", "z", "cores", "missing")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  
  # filter spatial meta to test sites
  spatial_tests <- meta2 %>%
    filter(Species == this.species, Stage == this.stage)
  
  # define missing var options
  missing_options <- spatial_tests %>%
    pull(Missing) %>%
    unique()
  
  # for each option
  for(missing in missing_options){
    
    #define initial predictors
    predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")
    
    #remove predictors if missing
    if(missing == "chl"){
      predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "slope")
    }
    
    #-------------------------------------------------------------------------------
    # 1. Formatting
    #-------------------------------------------------------------------------------
    
    #read in presences and pseudo-absences with environmental data
    data <- readRDS(paste0("output/extraction/", this.species, "_", this.site, "_", this.stage, "_extracted.rds"))
    
    #format data
    data <- data %>%
      mutate(date = as_date(date),
             pb = as.factor(pb),
             sic = as.numeric(sic)) 
    
    # isolate self-testing dataset from within training data
    self_test_inds <- data %>%
      filter(pb == "presence") %>%
      pull(individual_id) %>%
      unique() %>%
      sample(size = length(.)*0.2)
    self_test <- data %>% 
      filter(pb == "background" & is.na(test_year) | pb == "presence") %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background"))) %>%
      filter(individual_id %in% self_test_inds)
    
    # impute if any missing values
    if(sum(is.na(self_test)) > 0){
      self_mice <- miceRanger(self_test, m=1)
      self_test <- completeData(self_mice)[[1]]
    }
    
    # remove self-test data from training data
    data <- data %>%
      filter(!individual_id %in% self_test_inds)
    
    #-------------------------------------------------------------------------------
    # 2. Hyperparameter Tuning
    #-------------------------------------------------------------------------------
    
    # set ideal number of cross-validation folds to 10
    v <- 10
    
    #if number of individuals is fewer than 10, change v to n_ind
    n_ind <- length(unique(data$individual_id))
    if(n_ind < 10){
      v <- n_ind
    }
    
    #define random forest settings
    brt_mod <- boost_tree() %>%
      set_mode("classification") %>%
      set_engine("lightgbm" #use lightgbm package
      ) %>%
      set_args(trees = tune(),
               tree_depth = tune(), 
               learn_rate = tune(), 
               min_n = 20) 
    
    #create workflow
    brt_wf <- workflow() %>%
      add_model(brt_mod)
    
    #create tuning grid
    learn.rate <- c(0.005, 0.01, 0.5)
    tree.depth <- c(1, 3, 5)
    trees <- seq(2000, 10000, 2000)
    grid <- expand_grid(learn_rate = learn.rate, tree_depth = tree.depth, trees = trees)
    
    #BUFFER
    #isolate buffer data
    buff <- data %>% 
      filter(pb %in% c("buffer", "presence")) %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ifelse(pb == "buffer", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #create cross-validation folds
    buff_folds <- group_vfold_cv(data = buff, 
                                 group = individual_id, #split training/testing data by individual
                                 v = v, #number of folds
                                 balance = "observations" #roughly the same number of points in each fold
    )
    
    #define formula for modelling
    buff_rec <- recipe(pb ~ ., data = buff) %>%
      update_role(individual_id, new_role = "ID") %>% #let model know that id is not a predictor
      step_downsample(pb)
    
    #update workflow
    buff_wf <- brt_wf %>%
      add_recipe(buff_rec)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    buff_tun <- tune_grid(buff_wf,
                          resamples = buff_folds,
                          grid = grid,
                          metrics = sdm_metric_set(), #includes boyce index as a tuning parameter
                          control = control_grid(allow_par = F, verbose = T))
    
    #get metric scores for each tuning value
    buff_metrics <- collect_metrics(buff_tun, summarize = T) %>%
      mutate(pb = "buffer")
    
    #extract best model
    buff_best <- select_best(buff_tun, metric = "boyce_cont")
    
    #set up model
    buff_best_mod <- boost_tree() %>%
      set_engine(engine = "lightgbm") %>%
      set_mode("classification") %>%
      set_args(min_n = 20, trees = buff_best$trees[1], 
               tree_depth = buff_best$tree_depth[1], 
               learn_rate = buff_best$learn_rate[1])
    
    #update workflow
    buff_best_wf <- buff_wf %>%
      update_model(buff_best_mod)
    
    #run best model on all data
    buff_fit <- buff_best_wf %>%
      fit(buff)
    
    #print buffer fit
    print("buffer")
    
    #BACKGROUND
    #remove non-predictor columns
    back <- data %>% 
      filter(pb == "background" & is.na(test_year) | pb == "presence") %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #create cross-validation folds
    back_folds <- group_vfold_cv(data = back, 
                                 group = individual_id, #split training/testing data by individual
                                 v = v, #number of folds
                                 balance = "observations" #roughly the same number of points in each fold
    )
    
    #define formula for modelling
    back_rec <- recipe(pb ~ ., data = back) %>%
      update_role(individual_id, new_role = "ID") %>% #let model know that id is not a predictor
      step_downsample(pb)
    
    #update workflow
    back_wf <- brt_wf %>%
      add_recipe(back_rec)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    back_tun <- tune_grid(back_wf,
                          resamples = back_folds,
                          grid = grid,
                          metrics = sdm_metric_set(), #includes boyce index as a tuning parameter
                          control = control_grid(allow_par = F))
    
    #get metric scores for each tuning value
    back_metrics <- collect_metrics(back_tun, summarize = T) %>%
      mutate(pb = "background")
    
    #extract best model
    back_best <- select_best(back_tun, metric = "boyce_cont")
    
    #set up model
    back_best_mod <- boost_tree() %>%
      set_engine(engine = "lightgbm") %>%
      set_mode("classification") %>%
      set_args(min_n = 20, trees = back_best$trees[1], 
               tree_depth = back_best$tree_depth[1], 
               learn_rate = back_best$learn_rate[1])
    
    #update workflow
    back_best_wf <- back_wf %>%
      update_model(back_best_mod)
    
    #run best model on all data
    back_fit <- back_best_wf %>%
      fit(back)
    
    #print background fit
    print("background")
    
    
    #CRWs
    #get CRW data and remove trip numbers
    crw_data <- data %>%
      filter(pb == "crw")
    crw_data$individual_id <- stringr::str_remove(crw_data$individual_id, "_\\d+$")
    
    # isolate CRW data
    crw <- data %>% 
      filter(pb == "presence") %>%
      bind_rows(crw_data) %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ifelse(pb == "crw", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #create cross-validation folds
    crw_folds <- group_vfold_cv(data = crw, 
                                group = individual_id, #split training/testing data by individual
                                v = v, #number of folds
                                balance = "observations" #roughly the same number of points in each fold
    )
    
    #define formula for modelling
    crw_rec <- recipe(pb ~ ., data = crw) %>%
      update_role(individual_id, new_role = "ID") %>% #let model know that id is not a predictor
      step_downsample(pb)
    
    #update workflow
    crw_wf <- brt_wf %>%
      add_recipe(crw_rec)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    crw_tun <- tune_grid(crw_wf,
                         resamples = crw_folds,
                         grid = grid,
                         metrics = sdm_metric_set(), #includes boyce index as a tuning parameter
                         control = control_grid(allow_par = F))
    
    #get metric scores for each tuning value
    crw_metrics <- collect_metrics(crw_tun, summarize = T) %>%
      mutate(pb = "crw")
    
    #extract best model
    crw_best <- select_best(crw_tun, metric = "boyce_cont")
    
    #set up model
    crw_best_mod <- boost_tree() %>%
      set_engine(engine = "lightgbm") %>%
      set_mode("classification") %>%
      set_args(min_n = 20, trees = crw_best$trees[1], 
               tree_depth = crw_best$tree_depth[1], 
               learn_rate = crw_best$learn_rate[1])
    
    #update workflow
    crw_best_wf <- crw_wf %>%
      update_model(crw_best_mod)
    
    #run best model on all data
    crw_fit <- crw_best_wf %>%
      fit(crw)
    
    #print crw fit
    print("crw")
    
    
    
    #-------------------------------------------------------------------------------
    # 3. Test Models
    #-------------------------------------------------------------------------------
    
    #filter spatial metadata to this species and stage
    meta3 <- meta2 %>% filter(Species == this.species & Stage == this.stage & Missing == missing)
    
    #extract list of sites for this species and stage
    sites <- unique(meta3$Site)
    
    #null table for output
    brt_boyce_final <- NULL
    
    #loop to test each site
    for(i in sites){
      test.site <- i
      
      #load in test data
      test <- readRDS(paste0("output/extraction/spatial/", this.species, "_", test.site, "_", this.stage, "_extracted.rds"))
      
      #only select predictors for testing
      test <- test %>%
        select(all_of(predictors), pb) %>%
        mutate(individual_id = "test")
      
      #remove NA values from test data
      test <- test %>% drop_na()
      
      #order PA
      test$pb <- ordered(as.factor(test$pb), levels = c("presence", "background"))
      
      # predict
      test$buff_suitability <- predict(buff_fit, test, type = "prob") %>%
        pull(.pred_presence)
      test$back_suitability <- predict(back_fit, test, type = "prob") %>%
        pull(.pred_presence)
      test$crw_suitability <- predict(crw_fit, test, type = "prob") %>%
        pull(.pred_presence)
      
      #evaluate buffers
      buff_boyce <- boyce_cont(test, pb, buff_suitability) %>%
        pull(.estimate)
      
      #evaluate background
      back_boyce <- boyce_cont(test, pb, back_suitability) %>%
        pull(.estimate)
      
      #evaluate crws
      crw_boyce <- boyce_cont(test, pb, crw_suitability) %>%
        pull(.estimate)
      
      #boyce scores
      brt_boyce <- expand.grid(buff = buff_boyce, back = back_boyce, crw = crw_boyce)
      brt_boyce$test_site <- i
      brt_boyce_final <- bind_rows(brt_boyce_final, brt_boyce)
    }
    
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
      mutate(missing_covar = missing)
    
    # join to other self test scores
    if(missing == missing_options[1]){
      all_self <- self
    } else {
      all_self <- rbind(all_self, self)
    }
    
    # join together boyce scores
    if(missing == missing_options[1]){
      all_boyce <- brt_boyce_final
    } else {
      all_boyce <- rbind(all_boyce, brt_boyce_final)
    }
    
    # collate model metrics
    metrics <- bind_rows(buff_metrics, back_metrics, crw_metrics) %>%
      filter(.metric == "boyce_cont") %>%
      mutate(missing_covar = missing)
    
    # join to other metrics
    if(missing == missing_options[1]){
      all_metrics <- metrics
    } else {
      all_metrics <- rbind(all_metrics, metrics)
    }
    
  }
  
  #export boyce scores
  saveRDS(all_boyce, 
          file = paste0("output/spatial/brt/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_brt.RDS"))
  
  # export metrics
  saveRDS(all_metrics,
          file = paste0("output/spatial/brt/", this.species, "_", this.site, "_", this.stage, "_metrics_brt.RDS"))
  
  # export self test
  saveRDS(all_self,
          file = paste0("output/spatial/brt/", this.species, "_", this.site, "_", this.stage, "_self_test_brt.RDS"))
  
  #print completion
  print(paste0("Completed ", this.species, " ", this.stage))
  
}
