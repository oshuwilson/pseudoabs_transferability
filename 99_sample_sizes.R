#-------------------------------------------------------------------------------
# Calculate Sample Sizes
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(tidyverse)


#-------------------------------------------------------------------------------
# Temporal
#-------------------------------------------------------------------------------

# read in metadata
meta <- read.csv("data/species_site_stage_metadata.csv")
for(z in 1:nrow(meta)){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "counts_final", "z")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  # read in thinned tracks
  tracks <- readRDS(paste0("output/thinned_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_thinned.RDS"))
  
  #if season = FALSE, separate by year
  if(season == FALSE){
    tracks$season <- year(tracks$date)
  }
  
  #if season = TRUE, separate by season
  if(season == TRUE){
    tracks$season <- year(round_date(tracks$date, unit="year"))
  }
  
  # total count
  ntot <- nrow(tracks)
  
  # total number of individuals
  ntot_ind <- n_distinct(tracks$individual_id)
  
  # count number of presences per season
  counts <- tracks %>%
    group_by(season) %>%
    summarise(n = ntot - n(),
              test_n = n(),
              n_ind = n_distinct(individual_id),
              train_n_ind = ntot_ind - n_distinct(individual_id),
              test_n_ind = n_distinct(individual_id)) %>%
    mutate(species = this.species,
           site = this.site,
           stage = this.stage)
  
  # combine to all other case studies
  if(z == 1){
    counts_final <- counts
  } else {
    counts_final <- bind_rows(counts_final, counts)
  }
}

# write out
saveRDS(counts_final, "output/temporal/sample_sizes/temporal_sample_sizes.RDS")


#-------------------------------------------------------------------------------
# Spatial
#-------------------------------------------------------------------------------

rm(list=ls())

# read in metadata
meta <- read.csv("data/species_site_stage_metadata.csv")
z <- 1
for(z in 1:21){
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "counts_final", "z")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  season <- meta[z, 4]
  
  # read in thinned tracks
  tracks <- readRDS(paste0("output/thinned_tracks/", this.species, "_", this.site, "_", this.stage, "_tracks_thinned.RDS"))
  
  # get sample size
  ntot <- nrow(tracks)
  
  # get number of individuals
  ntot_ind <- n_distinct(tracks$individual_id)
  
  # create data frame
  df <- data.frame(species = this.species,
                   site = this.site,
                   stage = this.stage,
                   n = ntot,
                   n_ind = ntot_ind)
  
  # combine to all other case studies
  if(z == 1){
    counts_final <- df
  } else {
    counts_final <- bind_rows(counts_final, df)
  }
}

# export
saveRDS(counts_final, "output/spatial/sample_sizes/spatial_sample_sizes.RDS")
