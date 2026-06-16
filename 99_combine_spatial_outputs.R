rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(tidyverse)

# read in metadata
meta <- read.csv("data/species_site_stage_metadata.csv")
meta2 <- read.csv("data/spatial_site_metadata.csv") %>%
  rename(species = Species, test_site = Site, stage = Stage, missing_covar = Missing)

# remove case studies without spatial tests available
meta <- meta %>% filter(Species != "ANPE" & Species != "EMPE" & Species != "SUFS")
meta <- meta %>% filter(Species != "MAPE" | 
                          (Species == "MAPE" & Stage != "incubation" & Stage != "post-breeding"))

# read in sample sizes
sample_sizes <- readRDS("output/spatial/sample_sizes/spatial_sample_sizes.RDS")

# for each case study
for(i in 1:nrow(meta)){
  
  # define species, site, and stage
  this.species <- meta$Species[i]
  this.site <- meta$Site[i]
  this.stage <- meta$Stage[i]
  
  # limit sample sizes to relevant species, site, and stage
  ss <- sample_sizes %>%
    filter(species == this.species,
           site == this.site,
           stage == this.stage)
  
  # reorder to put species, site, and stage first in dataframe
  ss <- ss %>%
    select(species, site, stage, everything())
  
  # read in extrapolation scores
  extra <- readRDS(paste0("output/spatial/extrapolation/", this.species, "_", this.site, "_", this.stage, "_shape_medians.RDS"))
  
  # join sample sizes and extrapolation scores
  ss <- ss %>%
    left_join(extra %>% mutate(species = this.species))
  
  # reorder to move pseudo up
  ss <- ss %>%
    select(species, site, stage, region, pseudo, everything()) %>%
    rename(test_site = region)
  
  # add missing covars from meta2
  ss <- ss %>%
    left_join(meta2 %>%
                select(species, test_site, stage, missing_covar))
  
  # read in random forests
  rf <- readRDS(paste0("output/spatial/rf/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_rf.RDS"))
  rf_self <- readRDS(paste0("output/spatial/rf/", this.species, "_", this.site, "_", this.stage, "_self_test_rf.RDS"))
  
  # read in boosted regression trees
  brt <- readRDS(paste0("output/spatial/brt/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_brt.RDS"))
  brt_self <- readRDS(paste0("output/spatial/brt/", this.species, "_", this.site, "_", this.stage, "_self_test_brt.RDS"))
  
  # read in bayesian additive regression trees
  bart <- readRDS(paste0("output/spatial/bart/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_bart.RDS"))
  bart_self <- readRDS(paste0("output/spatial/bart/", this.species, "_", this.site, "_", this.stage, "_self_test_bart.RDS"))
  
  # read in generalised additive mixed models
  gamm <- readRDS(paste0("output/spatial/gamm/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_gamm.RDS"))
  gamm_self <- readRDS(paste0("output/spatial/gamm/", this.species, "_", this.site, "_", this.stage, "_self_test_gamm.RDS"))
  
  # read in generalised linear mixed models
  glmm <- readRDS(paste0("output/spatial/glmm/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_glmm.RDS"))
  glmm_self <- readRDS(paste0("output/spatial/glmm/", this.species, "_", this.site, "_", this.stage, "_self_test_glmm.RDS"))
  
  # read in maxent
  maxent <- readRDS(paste0("output/spatial/maxent/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_maxent.RDS"))
  maxent_self <- readRDS(paste0("output/spatial/maxent/", this.species, "_", this.site, "_", this.stage, "_self_test_maxent.RDS"))
  
  # read in stpp
  stpp <- readRDS(paste0("output/spatial/stpp/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_stpp.RDS"))
  stpp_self <- readRDS(paste0("output/spatial/stpp/", this.species, "_", this.site, "_", this.stage, "_self_test_stpp.RDS"))
  
  # read in sdpe
  sdpe <- readRDS(paste0("output/spatial/sdpe/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_sdpe.RDS"))
  sdpe_self <- readRDS(paste0("output/spatial/sdpe/", this.species, "_", this.site, "_", this.stage, "_self_test_sdpe.RDS"))
  
  # assign the algorithm to each dataframe
  rf$algorithm <- "RF"
  rf_self$algorithm <- "RF"
  brt$algorithm <- "BRT"
  brt_self$algorithm <- "BRT"
  bart$algorithm <- "BART"
  bart_self$algorithm <- "BART"
  gamm$algorithm <- "GAMM"
  gamm_self$algorithm <- "GAMM"
  glmm$algorithm <- "GLMM"
  glmm_self$algorithm <- "GLMM"
  maxent$algorithm <- "MaxEnt"
  maxent_self$algorithm <- "MaxEnt"
  stpp$algorithm <- "mSTPP"
  stpp_self$algorithm <- "mSTPP"
  sdpe$algorithm <- "INLA-SPDE"
  sdpe_self$algorithm <- "INLA-SPDE"
  
  # combine all algorithm dataframes
  algos <- rbind(rf, brt, bart, gamm, glmm, maxent, stpp, sdpe)
  algos_self <- rbind(rf_self, brt_self, bart_self, gamm_self, glmm_self, maxent_self, stpp_self, sdpe_self)
  
  # pivot longer 
  algos <- algos %>%
    pivot_longer(cols = 1:3,
                 names_to = "pseudo",
                 values_to = "transferability") %>%
    mutate(pseudo = case_when(pseudo == "buff" ~ "buffers",
                              pseudo == "back" ~ "background",
                              pseudo == "crw" ~ "CRWs"))
  algos_self <- algos_self %>%
    pivot_longer(cols = 1:3,
                 names_to = "pseudo",
                 values_to = "self_test") %>%
    mutate(pseudo = case_when(pseudo == "buff" ~ "buffers",
                              pseudo == "back" ~ "background",
                              pseudo == "crw" ~ "CRWs"))
  
  # add species, site, and stage
  algos$species <- this.species
  algos$site <- this.site
  algos$stage <- this.stage
  algos_self$species <- this.species
  algos_self$site <- this.site
  algos_self$stage <- this.stage
  
  # join to sample sizes and extrapolation scores
  ss <- ss %>%
    left_join(algos, by = c("species", "site", "stage", "test_site", "pseudo"))
  ss <- ss %>%
    left_join(algos_self, by = c("species", "site", "stage", "pseudo", "algorithm", "missing_covar"))
  
  # join to all other case studies
  if(i == 1){
    all_cases <- ss
  } else {
    all_cases <- rbind(all_cases, ss)
  }
}

# change pseudo-absence naming
all_cases <- all_cases %>%
  mutate(pseudo = case_when(pseudo == "background" ~ "Background",
                              pseudo == "buffers" ~ "Buffer",
                              pseudo == "CRWs" ~ "CRWs")) 

# create pseudo-algorithm column
all_cases <- all_cases %>%
  mutate(algopseudo = paste(algorithm, pseudo, sep = " "))

# remove self-tests
all_cases <- all_cases %>%
  filter(site != test_site)

# export
saveRDS(all_cases, "output/spatial/transferability_scores.RDS")
