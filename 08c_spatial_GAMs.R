#automated script to run GAMs for spatial transfer validation
rm(list=ls())
#setwd("/mainfs/home/jcw2g17/Chapter 01/")
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(dplyr)
  library(mgcv)
  library(enmSdmX)
  library(lubridate)
  #library(foreach)
  #library(doParallel)
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
meta <- read.csv("data/species_site_stage_metadata.csv")
meta2 <- read.csv("output/spatial/spatial_site_metadata.csv")

#isolate subsets where all predictors are present in test data - change for those missing >10% of chl and/or wind
meta2 <- meta2 %>% filter(Missing == "") #possible options "", "chl", or "windchl"
meta <- meta %>% filter(Species %in% meta2$Species & Stage %in% meta2$Stage)

#define initial predictors
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")

#remove predictors if missing
missing <- meta2$Missing[1]

if(missing == "chl"){
  predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "wind", "eke", "slope")
}

if(missing == "windchl"){
  predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "slope")
}

#setup parallel programming
#registerDoParallel(cores = 21)

#loop to run through each species, stage, and site
for(z in 1:nrow(meta)) { #change to foreach to run in parallel
  try({
    
    #define parameters in loop
    rm(list=setdiff(ls(), c("meta", "meta2", "predictors", "z", "pred_gam")))
    this.species <- meta[z, 1]
    this.site <- meta[z, 2]
    this.stage <- meta[z, 3]

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
    
    # 2. Create GAMs
    
    #BUFFER
    #remove non-predictor columns
    buff_sel <- buff %>% select(all_of(predictors), pa)
    
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
    
    #save model
    saveRDS(buff_gam, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/buff_gam.RDS"))
    
    #remove unnecessary parameters to continue
    rm(buff_sel)
    
    
    #BACKGROUND
    #remove non-predictor columns
    back_sel <- back %>% select(all_of(predictors), pa)
    
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
    
    #save model
    saveRDS(back_gam, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/back_gam.RDS"))
    
    #remove unnecessary parameters to continue
    rm(back_sel)
    
    
    #CRWs
    #remove non-predictor columns
    crw_sel <- crw %>% select(all_of(predictors), pa)
    
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
    
    #save model
    saveRDS(crw_gam, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/crw_gam.RDS"))
    
    #remove unnecessary parameters to continue
    rm(crw_sel)
    
    
    
    # 3. Test GAMs
    
    #filter spatial metadata to this species and stage
    meta3 <- meta2 %>% filter(Species == this.species & Stage == this.stage)
    
    #extract list of sites for this species and stage
    meta3$Site <- as.factor(meta3$Site)
    sites <- levels(meta3$Site)
    
    #null table for loop
    gam_boyce_final <- NULL
    
    #run for loop to test each site
    for(i in sites){
      test.site <- i
      
      #load in test data
      back_test <- read.csv(paste0("output/spatial/", this.species, "/", this.stage, "/extraction/", test.site, "_background.csv"))
      tracks_test <- read.csv(paste0("output/spatial/", this.species, "/", this.stage, "/extraction/", test.site, "_presences.csv"))
      
      #predict and evaluate buffers
      p1 <- predict.gam(buff_gam, tracks_test, type = "response")
      p2 <- predict.gam(buff_gam, back_test, type = "response")
      buff_gam_boyce <- evalContBoyce(p1, p2, na.rm=TRUE)
      
      #predict and evaluate background
      p1 <- predict.gam(back_gam, tracks_test, type = "response")
      p2 <- predict.gam(back_gam, back_test, type = "response")
      back_gam_boyce <- evalContBoyce(p1, p2, na.rm=TRUE)
      
      #predict and evaluate crws
      p1 <- predict.gam(crw_gam, tracks_test, type = "response")
      p2 <- predict.gam(crw_gam, back_test, type = "response")
      crw_gam_boyce <- evalContBoyce(p1, p2, na.rm=TRUE)
      
      #FINAL DATA
      #boyce scores
      gam_boyce <- expand.grid(buff = buff_gam_boyce, back = back_gam_boyce, crw = crw_gam_boyce)
      gam_boyce$site <- i
      gam_boyce_final <- rbind(gam_boyce_final, gam_boyce)
    }
    
    #export boyce scores
    saveRDS(gam_boyce_final, 
            file = paste0("output/spatial/", this.species, "/", this.site, "/", this.stage, "/boyce_scores_gam.RDS"))
    
    print(paste0(this.species, " ", this.stage, " ", this.site, " completed"))
    
  })
  
}
