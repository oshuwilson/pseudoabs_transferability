rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(tidyverse)

# read in metadata
meta <- read.csv("data/species_site_stage_metadata.csv")

# read in sample sizes
sample_sizes <- readRDS("output/temporal/sample_sizes/temporal_sample_sizes.RDS")

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
  extra <- readRDS(paste0("output/temporal/extrapolation/", this.species, "_", this.site, "_", this.stage, "_shape_medians.RDS"))
  
  # join sample sizes and extrapolation scores
  ss <- ss %>%
    left_join(extra, by = "season")
  
  # reorder to move pseudo up
  ss <- ss %>%
    select(species, site, stage, season, pseudo, everything())
  
  # read in random forests
  rf <- readRDS(paste0("output/temporal/rf/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_rf.RDS"))
  rf_self <- readRDS(paste0("output/temporal/rf/", this.species, "_", this.site, "_", this.stage, "_self_test_rf.RDS"))
  
  # read in boosted regression trees
  brt <- readRDS(paste0("output/temporal/brt/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_brt.RDS"))
  brt_self <- readRDS(paste0("output/temporal/brt/", this.species, "_", this.site, "_", this.stage, "_self_test_brt.RDS"))
  
  # read in bayesian additive regression trees
  bart <- readRDS(paste0("output/temporal/bart/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_bart.RDS"))
  bart_self <- readRDS(paste0("output/temporal/bart/", this.species, "_", this.site, "_", this.stage, "_self_test_bart.RDS"))
  
  # read in generalised additive mixed models
  gamm <- readRDS(paste0("output/temporal/gamm/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_gamm.RDS"))
  gamm_self <- readRDS(paste0("output/temporal/gamm/", this.species, "_", this.site, "_", this.stage, "_self_test_gamm.RDS"))
  
  # read in generalised linear mixed models
  glmm <- readRDS(paste0("output/temporal/glmm/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_glmm.RDS"))
  glmm_self <- readRDS(paste0("output/temporal/glmm/", this.species, "_", this.site, "_", this.stage, "_self_test_glmm.RDS"))
  
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
  
  # combine all algorithm dataframes
  algos <- rbind(rf, brt, bart, gamm, glmm)
  algos_self <- rbind(rf_self, brt_self, bart_self, gamm_self, glmm_self)
  
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
    left_join(algos, by = c("species", "site", "stage", "season", "pseudo"))
  ss <- ss %>%
    left_join(algos_self, by = c("species", "site", "stage", "season", "pseudo", "algorithm"))
  
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

# remove case studies with too small a sample size
cases <- all_cases %>%
  filter(test_n_ind > 2, n > 50)

# export
saveRDS(cases, "output/temporal/transferability_scores.RDS")
saveRDS(all_cases, "output/temporal/transferability_all_cases.RDS")
