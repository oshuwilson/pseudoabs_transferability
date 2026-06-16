#-------------------------------------------------------------------------------
# Test different temporal intervals for Southern Elephant Seals
#-------------------------------------------------------------------------------
# isolate the final 3 years of data and train models on preceding years
# iteratively increase the temporal distance between training and testing data

rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/resubmission/")

{
  library(dplyr)
  library(lubridate)
  library(tidymodels)
  library(tidysdm)
  library(ranger)
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
for(z in 4){
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
    #drop_na(all_of(predictors)) %>%
    mutate(pb = factor(pb)) %>%
    mutate(pb = ordered(pb, levels = c("presence", "background")))
  
  # isolate training data from earlier seasons
  train <- data %>% filter(season < (final_season - 2))
  
  # list each season in training data in reverse order
  seasons <- train %>% arrange(-season) %>% pull(season) %>% unique()
  
  # for each season
  for(this_season in seasons){
    print(this_season)
    
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
                                            date %in% presences$date)
    
    # join together
    train_season <- bind_rows(presences, background)
    
    # drop NAs
    train_season <- train_season %>% drop_na(all_of(predictors))
    
    #remove non-predictor columns
    train_season <- train_season %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    # define random forest settings
    rf_mod <- rand_forest() %>%
      set_mode("classification") %>%
      set_engine("ranger", #use ranger package
                 importance = "impurity" #gini index for importance
      ) %>%
      set_args(trees = 1000, #1000 trees
               mtry = tune(), #tune mtry
               min_n = 1) #minimum number of samples in a node
    
    # create workflow
    rf_wf <- workflow() %>%
      add_model(rf_mod)
    
    #create tuning grid
    mtry <- 2:4
    grid <- expand_grid(mtry = mtry)
    
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
    rf_wf <- rf_wf %>%
      add_recipe(rec)
    
    # enable parallelisation
    plan(multisession, workers = cores)
    
    # set seed
    set.seed(777)
    
    #run models with tuning
    tun <- tune_grid(rf_wf,
                     resamples = folds,
                     grid = grid,
                     metrics = sdm_metric_set()) #includes boyce index as a tuning parameter
    
    # end parallelisation
    plan(sequential)
    
    #get metric scores for each mtry value
    metrics <- collect_metrics(tun, summarize = T) 
    
    #extract best model
    best <- select_best(tun, metric = "boyce_cont")
    
    #set up model
    best_mod <- rand_forest() %>%
      set_engine(engine = "ranger", importance = "impurity") %>%
      set_mode("classification") %>%
      set_args(trees = 1000, mtry = best$mtry[1], min_n = 1)
    
    #update workflow
    best_wf <- rf_wf %>%
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
      algo = "RF"
    )
    
    # join to other results
    if(this_season == seasons[1]){
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

ggplot(final_results, aes(x = gap, y = boyce, col = paste0(species, site, stage))) +
  geom_point() +
  geom_smooth(aes(group = paste0(species, site, stage)), method = "lm", se = F) +
  theme_bw() +
  labs(x = "Temporal gap (years)", y = "Boyce index", col = "Case study") +
  theme(legend.position = "bottom")

# export
saveRDS(final_results, "output/temporal/distance/rf.rds")
