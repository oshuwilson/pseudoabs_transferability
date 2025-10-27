#automated script to run GLMMs for leave-year-out validation

#clear workspace and set working directory
rm(list=ls())
setwd("/iridisfs/scratch/jcw2g17/dsdm/")

#load required packages
{
  library(dplyr)
  library(lme4)
  library(tidysdm)
  library(lubridate)
  library(foreach)
  library(doParallel)
}
select <- dplyr::select

# GLM function to only use predictors that remain after later steps remove some
pred_glm <- function(df){
  lme4::glmer(
    as.formula(
      paste0(
        "pb ~ ", 
        setdiff(names(df), c("pb", "individual_id")) %>% paste0(collapse = " + "),
        " + (1|individual_id)")),
    family = binomial, data = df
  )
}

#read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")

#loop to run through each species, stage, and site iteratively
for(z in 21) {
  try({
    
    #define initial predictors
    predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")
    
    #define parameters in loop
    rm(list=setdiff(ls(), c("meta", "predictors", "z", "pred_glm", "cores")))
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
    glm_boyce_final <- NULL
    
    #setup parallel programming
    #registerDoParallel(cores = length(seasons))
    
    #loop over each season
    #glm_boyce_final <- foreach(i = seasons, .combine = rbind, .packages = c("dplyr", "lme4", "tidysdm", "lubridate")) %dopar% {
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
        mutate(pb = ifelse(pb == "background", 0, 1)) %>%
        mutate(pb = ordered(pb, levels = c("1", "0"))) %>%
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
      
      
      # 3. Run GLMMs
      
      #BUFFER
      #isolate buffer data
      buff_sel <- train %>% 
        filter(pb %in% c("buffer", "presence")) %>%
        select(all_of(predictors), pb, individual_id) %>%
        mutate(pb = ifelse(pb == "buffer", 0, 1))
      
      #remove SIC if less than 5 distinct values
      if(n_distinct(buff_sel$sic) < 5 & this.species != "SUFS"){
        buff_sel <- buff_sel %>% select(-sic)
      }
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(buff_sel)) > 0.1*nrow(buff_sel)){
        buff_sel <- buff_sel[colSums(is.na(buff_sel)) < 0.1*nrow(buff_sel)]
      }
      
      #run glm
      buff_glm <- pred_glm(buff_sel)
      
      #predict and evaluate
      self_test$buff_suitability <- predict(buff_glm, self_test, re.form = NA, type = "response")
      test$buff_pred <- predict(buff_glm, test, re.form = NA, type = "response")
      buff_glm_boyce <- boyce_cont(test, pb, buff_pred) %>%
        pull(.estimate)
      buff_glm_boyce0 <- boyce_cont(self_test, pb, buff_suitability) %>%
        pull(.estimate)
      
      #remove unnecessary parameters to continue
      rm(buff_glm, buff_sel)
      
      
      #BACKGROUND
      #remove non-predictor columns
      back_sel <- train %>% 
        filter(pb == "background" & test_year == this.test | pb == "presence") %>%
        select(all_of(predictors), pb, individual_id) %>%
        mutate(pb = ifelse(pb == "background", 0, 1))
      
      #remove SIC if less than 5 distinct values
      if(n_distinct(back_sel$sic) < 5 & this.species != "SUFS"){
        back_sel <- back_sel %>% select(-sic)
      }
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(back_sel)) > 0.1*nrow(back_sel)){
        back_sel <- back_sel[colSums(is.na(back_sel)) < 0.1*nrow(back_sel)]
      }
      
      #run glm
      back_glm <- pred_glm(back_sel)
      
      #predict and evaluate
      self_test$back_suitability <- predict(back_glm, self_test, re.form = NA, type = "response")
      test$back_pred <- predict(back_glm, test, re.form = NA, type = "response")
      back_glm_boyce <- boyce_cont(test, pb, back_pred) %>%
        pull(.estimate)
      back_glm_boyce0 <- boyce_cont(self_test, pb, back_suitability) %>%
        pull(.estimate)
      
      #remove unnecessary parameters to continue
      rm(back_glm, back_sel)
      
      
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
        mutate(pb = ifelse(pb == "crw", 0, 1))
      
      #remove SIC if less than 5 distinct values
      if(n_distinct(crw_sel$sic) < 5 & this.species != "SUFS"){
        crw_sel <- crw_sel %>% select(-sic)
      }
      
      #remove columns where missing data is over 10% of rows
      if(sum(is.na(crw_sel)) > 0.1*nrow(crw_sel)){
        crw_sel <- crw_sel[colSums(is.na(crw_sel)) < 0.1*nrow(crw_sel)]
      }
      
      #run glm
      crw_glm <- pred_glm(crw_sel)
      
      #predict and evaluate
      self_test$crw_suitability <- predict(crw_glm, self_test, re.form = NA, type = "response")
      test$crw_pred <- predict(crw_glm, test, re.form = NA, type = "response")
      crw_glm_boyce <- boyce_cont(test, pb, crw_pred) %>%
        pull(.estimate)
      crw_glm_boyce0 <- boyce_cont(self_test, pb, crw_suitability) %>%
        pull(.estimate)
      
      #remove unnecessary parameters to continue
      rm(crw_glm, crw_sel)
      
      # collate self test scores
      self <- expand.grid(buff = buff_glm_boyce0, back = back_glm_boyce0, crw = crw_glm_boyce0) %>%
        mutate(season = i)
      
      # join to other self test scores
      if(i == seasons[1]){
        all_self <- self
      } else {
        all_self <- rbind(all_self, self)
      }
      
      #FINAL DATA
      glm_boyce <- expand.grid(buff = buff_glm_boyce, back = back_glm_boyce, crw = crw_glm_boyce)
      glm_boyce$season <- i
      glm_boyce_final <- rbind(glm_boyce_final, glm_boyce)
      # glm_boyce
    }
    
    # export CBI scores
    saveRDS(glm_boyce_final, 
            file = paste0("output/temporal/glmm/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_glmm.RDS"))
    
    # export self-test scores
    saveRDS(all_self,
            file = paste0("output/temporal/glmm/", this.species, "_", this.site, "_", this.stage, "_self_test_glmm.RDS"))
    
    #print completion
    print(paste0("Completed ", this.species, " ", this.site, " ", this.stage))
    
  })
  
}
