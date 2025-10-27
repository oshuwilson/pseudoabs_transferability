#-------------------------------------------------------------------------------
# Spatial Generalised Linear Mixed Models
#-------------------------------------------------------------------------------

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
  library(miceRanger)
}

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
meta2 <- read.csv("data/spatial_site_metadata.csv")

#remove ANPE, EMPE, SUFS, and MAPE incubation/post-breeding (no spatial transfer)
meta <- meta %>% filter(Species != "ANPE" & Species != "EMPE" & Species != "SUFS")
meta <- meta %>% filter(Species != "MAPE" | 
                          (Species == "MAPE" & Stage != "incubation" & Stage != "post-breeding"))

# for each population
for(z in 5:nrow(meta)){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "meta2", "predictors", "z", "missing", "pred_glm")))
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
    
    #setup parallel programming
    #registerDoParallel(cores = 21)
    
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
      mutate(pb = ifelse(pb == "background", 0, 1)) %>%
      mutate(pb = as.factor(pb)) %>%
      mutate(pb = ordered(pb, levels = c("1", "0"))) %>%
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
    # 2. Create Models
    #-------------------------------------------------------------------------------
    
    #BUFFER
    #isolate buffer data
    buff <- data %>% 
      filter(pb %in% c("buffer", "presence")) %>%
      select(all_of(predictors), pb, individual_id) %>%
      mutate(pb = ifelse(pb == "buffer", 0, 1))
    
    #remove SIC if less than 5 distinct values
    if(n_distinct(buff$sic) < 5){
      buff <- buff %>% select(-sic)
    }
    
    #remove columns where missing data is over 10% of rows
    if(sum(is.na(buff)) > 0.1*nrow(buff)){
      buff <- buff[colSums(is.na(buff)) < 0.1*nrow(buff)]
    }
    
    #convert individual_id to factor
    buff$individual_id <- as.factor(buff$individual_id)
    
    #run glm
    buff_glm <- pred_glm(buff)
    
    #print buffer completion
    print("buff")
    
    #BACKGROUND
    #remove non-predictor columns
    back <- data %>% 
      filter(pb == "background" & is.na(test_year) | pb == "presence") %>%
      select(all_of(predictors), pb, individual_id)  %>%
      mutate(pb = ifelse(pb == "background", 0, 1))
    
    #remove SIC if less than 5 distinct values
    if(n_distinct(back$sic) < 5){
      back <- back %>% select(-sic)
    }
    
    #remove columns where missing data is over 10% of rows
    if(sum(is.na(back)) > 0.1*nrow(back)){
      back <- back[colSums(is.na(back)) < 0.1*nrow(back)]
    }
    
    #convert individual_id to factor
    back$individual_id <- as.factor(back$individual_id)
    
    #run glm
    back_glm <- pred_glm(back)
    
    #print background completion
    print("back")
    
    
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
      mutate(pb = ifelse(pb == "crw", 0, 1))
    
    #remove SIC if less than 5 distinct values
    if(n_distinct(crw$sic) < 5){
      crw <- crw %>% select(-sic)
    }
    
    #remove columns where missing data is over 10% of rows
    if(sum(is.na(crw)) > 0.1*nrow(crw)){
      crw <- crw[colSums(is.na(crw)) < 0.1*nrow(crw)]
    }
    
    #convert individual_id to factor
    crw$individual_id <- as.factor(crw$individual_id)
    
    #run glm
    crw_glm <- pred_glm(crw)
    
    #print crw completion
    print("crw")
    
    #-------------------------------------------------------------------------------
    # 3. Test Models
    #-------------------------------------------------------------------------------
    
    #filter spatial metadata to this species and stage
    meta3 <- meta2 %>% filter(Species == this.species & Stage == this.stage & Missing == missing)
    
    #extract list of sites for this species and stage
    sites <- unique(meta3$Site)
    
    #null table for output
    glmm_boyce_final <- NULL
    
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
      test$buff_suitability <- predict(buff_glm, test, re.form = NA, type = "response")
      test$back_suitability <- predict(back_glm, test, re.form = NA, type = "response")
      test$crw_suitability <- predict(crw_glm, test, re.form = NA, type = "response")
      
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
      glmm_boyce <- expand.grid(buff = buff_boyce, back = back_boyce, crw = crw_boyce)
      glmm_boyce$test_site <- i
      glmm_boyce_final <- bind_rows(glmm_boyce_final, glmm_boyce)
    }
    
    # join together boyce scores
    if(missing == missing_options[1]){
      all_boyce <- glmm_boyce_final
    } else {
      all_boyce <- rbind(all_boyce, glmm_boyce_final)
    }
    
    # self test models
    self_test$buff_suitability <- predict(buff_glm, self_test, re.form = NA, type = "response")
    self_test$back_suitability <- predict(back_glm, self_test, re.form = NA, type = "response")
    self_test$crw_suitability <- predict(crw_glm, self_test, re.form = NA, type = "response")
    
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
    
  }
  
  #export boyce scores
  saveRDS(all_boyce, 
          file = paste0("output/spatial/glmm/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_glmm.RDS"))
  
  # export self test
  saveRDS(all_self,
          file = paste0("output/spatial/glmm/", this.species, "_", this.site, "_", this.stage, "_self_test_glmm.RDS"))
  
  # print completion
  print(paste0(this.species, " ", this.stage, " ", this.site, " completed"))
  
}
