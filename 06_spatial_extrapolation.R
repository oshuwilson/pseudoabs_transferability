#-------------------------------------------------------------------------------
# Spatial Environmental Extrapolation
#-------------------------------------------------------------------------------

#clear workspace and set working directory
rm(list=ls())
setwd("/mainfs/scratch/jcw2g17/dsdm/")

#load required packages
{
  library(lubridate)
  library(flexsdm)
  library(dplyr)
  library(tidyr)
}

#read in table with info for each species, site and stage
meta <- readRDS("data/species_site_stage_metadata.RDS")
meta2 <- readRDS("data/spatial_site_metadata.RDS")

#define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "slope")

#remove ANPE, EMPE, SUFS, and MAPE incubation/post-breeding (no spatial transfer)
meta <- meta %>% filter(Species != "ANPE" & Species != "EMPE" & Species != "SUFS")
meta <- meta %>% filter(Species != "MAPE" | 
                          (Species == "MAPE" & Stage != "incubation" & Stage != "post-breeding"))

#loop over every species, site, and stage
for(z in 5:nrow(meta)){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "predictors", "z", "meta2")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  # filter spatial meta to test sites
  spatial_tests <- meta2 %>%
    filter(Species == this.species, Stage == this.stage)
  
  #-------------------------------------------------------------------------------
  # 1. Formatting
  #-------------------------------------------------------------------------------
  
  #read in presences and pseudo-absences with environmental data
  data <- readRDS(paste0("output/extraction/", this.species, "_", this.site, "_", this.stage, "_extracted.rds"))
  
  #format data
  data <- data %>%
    mutate(pb = as.factor(pb),
           sic = as.numeric(sic)) %>%
    select(all_of(predictors), pb, test_year)
  
  #extract list of sites for this species and stage
  sites <- unique(spatial_tests$Site)
  
  #create empty table
  shape_values <- NULL
  shape_all <- NULL
  
  #loop over each season
  for(i in sites){
    test.site <- i
    
    #load in test data
    test <- readRDS(paste0("output/extraction/spatial/", this.species, "_", test.site, "_", this.stage, "_extracted.rds"))
    
    #remove chl from predictors if missing
    missing_chl <- spatial_tests %>%
      filter(Site == test.site) %>%
      pull(Missing)
    if(missing_chl == "chl"){
      shape_predictors <- predictors[predictors != "chl"]
    } else {
      shape_predictors <- predictors
    }
    
    #fix sea ice if all values = 0
    if(sum(test$sic) == 0 & sum(data$sic == 0)){
      shape_predictors <- shape_predictors[shape_predictors != "sic"]
    }
    
    #only select predictors for testing data
    test <- test %>%
      select(all_of(shape_predictors))
    
    # isolate buffer training data
    buff_train <- data %>% 
      filter(pb == "presence" | pb == "buffer") %>%
      select(all_of(shape_predictors), pb) %>%
      mutate(pb = ifelse(pb == "buffer", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    # isolate background training data
    back_train <- data %>% 
      filter(pb == "background" & is.na(test_year) | pb == "presence") %>%
      select(all_of(shape_predictors), pb) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    # isolate crw training data
    crw_train <- data %>% 
      filter(pb == "presence" | pb == "crw") %>%
      select(all_of(shape_predictors), pb) %>%
      mutate(pb = ifelse(pb == "crw", "background", "presence")) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background")))
    
    #BUFFERS
    #calculate shape
    buff_shape <- extra_eval(training_data = buff_train, pr_ab = "pb", projection_data = test, n_cores = 40)
    
    #BACKGROUND
    #calculate shape
    back_shape <- extra_eval(training_data = back_train, pr_ab = "pb", projection_data = test, n_cores = 40)

    #CRW
    #calculate shape
    crw_shape <- extra_eval(training_data = crw_train, pr_ab = "pb", projection_data = test, n_cores = 40)
    
    #store scores
    shape_scores <- cbind(buff_shape[,1], back_shape[,1], crw_shape[,1])
    names(shape_scores) <- c("buff", "back", "crw")
    
    #extract LQ, median, and UQ for each set of scores
    buffers <- quantile(shape_scores$buff, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
    background <- quantile(shape_scores$back, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
    CRWs <- quantile(shape_scores$crw, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
    
    #format data table with all LQ, median, and UQ values
    shape_IQR <- as.data.frame(rbind(buffers, background, CRWs))
    names(shape_IQR) <- c("LQ", "Median", "UQ")
    shape_IQR$region <- test.site
    shape_IQR$pseudo <- c("buffers", "background", "CRWs")
    
    #bind to dataset with information for every season
    shape_values <- rbind(shape_values, shape_IQR)
    
    shape_scores$region <- test.site
    shape_all <- rbind(shape_all, shape_scores)
    
  }
  
  #export dataset
  saveRDS(shape_values,
          file = paste0("output/spatial/extrapolation/", this.species, "_", this.site, "_", this.stage, "_shape_medians.RDS"))
  saveRDS(shape_all,
          file = paste0("output/spatial/extrapolation/", this.species, "_", this.site, "_", this.stage, "_shape_scores.RDS"))
  
  #print to show that this species has completed
  print(paste0(this.species, " ", this.site, " ", this.stage, " completed"))
  
}
