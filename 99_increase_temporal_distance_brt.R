#-------------------------------------------------------------------------------
# Test different temporal intervals for Southern Elephant Seals
#-------------------------------------------------------------------------------
# isolate the final 3 years of data and train models on preceding years
# iteratively increase the temporal distance between training and testing data

rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/dsdm/")

{
  library(dplyr)
  library(lubridate)
  library(tidymodels)
  library(tidysdm)
  library(parsnip)
  library(lightgbm)
  library(bonsai)
  library(themis)
  library(miceRanger)
  library(future)
}

# set seed
set.seed(777)

# set number of cores for parallelisation
cores <- 10

# read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

# limit to rows with enough data
meta <- meta %>%
  filter(Species == "SOES" & Site %in% c("Marion", "WAP"))

# define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "slope")

# run over each row of metadata
for(z in 1:4){
  print(z)
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "predictors", "z", "cores", "final_results")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  # 1. Formatting 
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
  
  # limit to background and presence info
  data <- data %>%
    filter(pb %in% c("presence", "background"))
  
  # find individuals with > 15 presences following thinning (> 15 tracked days)
  data_rich_inds <- data %>%
    filter(pb == "presence") %>% 
    group_by(individual_id) %>% 
    summarise(n = n()) %>% 
    ungroup() %>% 
    filter(n > 15) %>% 
    pull(individual_id)
  
  # limit data to these inds
  data <- data %>% filter(individual_id %in% data_rich_inds)
  
  # find final season 
  final_season <- max(data$season)
  if(z == 4){
    final_season <- 2010
  }
  
  # isolate test data within 2 years of final season
  test <- data %>% filter(season >= (final_season - 2) & is.na(test_year)) %>%
    drop_na(all_of(predictors)) %>%
    mutate(pb = factor(pb)) %>%
    mutate(pb = ordered(pb, levels = c("presence", "background")))
  
  # isolate training data from earlier seasons
  train <- data %>% filter(season < (final_season - 2))
  
  # list each season in training data in reverse order
  seasons <- train %>% arrange(-season) %>% pull(season) %>% unique()
  
  # for each season
  for(this_season in seasons){
    print(this_season)
    
    # if(z == 2 & this_season <= 2006){
    #   next
    # }
    # if(z == 3 & this_season <= 1998){
    #   next
    # }
    # if(z == 4 & this_season == 2006){
    #   next
    # }
    
    # isolate training data from this season and earlier
    train_season <- train %>% filter(season <= this_season)
    
    # get unique number of individuals
    n_inds <- length(unique(train_season$individual_id))
    
    # if n_inds under 5 move on to next loop
    if(n_inds < 5){
      next
    }
    
    # if fewer than 100 presences, move on to next loop
    n_presences <- train_season %>% filter(pb == "presence") %>% nrow()
    if(n_presences < 100){
      next
    }
    
    # randomly sample 100 presence locations
    presences <- train_season %>% filter(pb == "presence") %>% sample_n(100)
    
    # get temporally corresponding background points
    background <- train_season %>% filter(pb == "background" & 
                                            as_date(date) %in% as_date(presences$date))
    
    # join together
    train_season <- bind_rows(presences, background)
    
    # drop NAs
    train_season <- train_season %>% drop_na(all_of(predictors))
    
    #remove non-predictor columns
    train_season <- train_season %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    # define random forest settings
    brt_mod <- boost_tree() %>%
      set_mode("classification") %>%
      set_engine("lightgbm" #use lightgbm package
      ) %>%
      set_args(trees = tune(),
               tree_depth = tune(), 
               learn_rate = tune(), 
               min_n = 20) 
    
    # create workflow
    brt_wf <- workflow() %>%
      add_model(brt_mod)
    
    #create tuning grid
    learn.rate <- c(0.005, 0.01, 0.5)
    tree.depth <- c(1, 3, 5)
    trees <- seq(2000, 10000, 2000)
    grid <- expand_grid(learn_rate = learn.rate, tree_depth = tree.depth, trees = trees)
    
    #create cross-validation folds
    folds <- group_vfold_cv(data = train_season, 
                            group = individual_id, #split training/testing data by individual
                            v = 5, #number of folds
                            balance = "observations" #roughly the same number of points in each fold
    )
    
    #define formula for modelling
    rec <- recipe(pb ~ ., data = train_season) %>%
      update_role(individual_id, new_role = "ID") %>% #let model know that id is not a predictor
      step_downsample(pb)
    
    #update workflow
    brt_wf <- brt_wf %>%
      add_recipe(rec)
    
    # run in sequence
    plan(sequential)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    tun <- tune_grid(brt_wf,
                     resamples = folds,
                     grid = grid,
                     metrics = sdm_metric_set())
    
    #get metric scores for each mtry value
    metrics <- collect_metrics(tun, summarize = T) 
    
    #extract best model
    best <- select_best(tun, metric = "boyce_cont")
    
    #set up model
    best_mod <- boost_tree() %>%
      set_engine(engine = "lightgbm") %>%
      set_mode("classification") %>%
      set_args(min_n = 20, trees = best$trees[1], 
               tree_depth = best$tree_depth[1], 
               learn_rate = best$learn_rate[1])
    
    #update workflow
    best_wf <- brt_wf %>%
      update_model(best_mod)
    
    #run best model on all data
    fit <- best_wf %>%
      fit(train_season)
    
    #predict habitat suitability to training data
    test$suitability <- predict(fit, test, type = "prob") %>%
      pull(.pred_presence)
    
    #calculate boyce index
    boyce <- boyce_cont(test, pb, suitability) %>%
      pull(.estimate)
    
    # calculate gap between training and testing data
    test_seasons <- test$season %>% unique() %>% sort()
    gap <- min(test_seasons) - this_season
    
    # create table with gap and boyce index
    results <- data.frame(
      species = this.species,
      site = this.site,
      stage = this.stage,
      gap = gap,
      boyce = boyce,
      algo = "BRT"
    )
    
    # join to other results
    if(!exists("all_results")){
      all_results <- results
    } else {
      all_results <- bind_rows(all_results, results)
    }

  }
  
  # join to all case studies
  if(z == 1){
    final_results <- all_results
  } else {
    final_results <- bind_rows(final_results, all_results)
  }
}

# export
saveRDS(final_results, "output/temporal/distance/brt.rds")
