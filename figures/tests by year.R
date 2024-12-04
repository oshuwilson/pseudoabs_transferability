#tile plot of the years per stage
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01/output/leave-year-out")

library(tidyverse)

#read in the data
boyce <- readRDS("boyce_final.RDS")

#remove ANFS GLS years
boyce <- boyce %>% filter(species != "ANFS" | stage != "post-moult" | !season %in% c("2008", "2009"))

#summarise seasons recorded for each species/site/stage
seasons <- boyce %>% 
  select(species, site, stage, season) %>%
  distinct()

#add number of individuals per year
n_ind <- readRDS("n_ind.RDS")
seasons <- left_join(seasons, n_ind, by = c("species", "site", "stage", "season"))

#rename levels of species
seasons <- seasons %>%
  mutate(species = as.factor(species))
levels(seasons$species) <- c("Adélie Penguin", "Antarctic Fur Seal", "Antarctic Petrel", 
                             "Crabeater Seal", "Emperor Penguin", "Grey-headed Albatross", 
                             "Humpback Whale", "Macaroni Penguin", "Southern Elephant Seal",
                             "Subantarctic Fur Seal")

#rename levels of site
seasons <- seasons %>%
  mutate(site = as.factor(site))
levels(seasons$site) <- c("Marion Island", "Mawson Coast", "Pointe Geologie", "Svarthamaren",
                          "South Georgia", "Antarctic Peninsula")

#rename levels of stage
seasons <- seasons %>%
  mutate(stage = as.factor(stage))
levels(seasons$stage) <- c("Breeding", "Chick-Rearing", "Early Chick-Rearing", 
                           "Non-Central Place Foraging", "Incubation", "Late Chick-Rearing",
                           "No Stage", "Post-Breeding", "Post-Moult")

#change crabeater stage to non-central place foraging
seasons <- seasons %>%
  mutate(stage = ifelse(species == "Crabeater Seal", "Non-Central Place Foraging", as.character(stage)))

#add a column to combine species, site and stage
seasons <- seasons %>%
  mutate(test = as.factor(paste(species, site, stage, sep = " - ")),
         season = as.factor(season))



#plot
ggplot(seasons, aes(x = season, y = fct_rev(test), fill = n_ind)) +
  geom_tile(aes(height = 0.95, width = 0.96)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Assigned Year", y = "Study Species, Site, and Stage",
       fill = "Number of Tracked Individuals") +
  scale_fill_viridis_c(trans = "log",
                       breaks = c(1, 5, 15, 50)) + 
  theme(axis.text.x=element_text(angle = 0, hjust=0.5),
        axis.text.y=element_text(size=12),
        axis.title = element_text(size = 14),
        legend.position = "bottom") 
ggsave(filename = "~/OneDrive - University of Southampton/Documents/Chapter 01/text/figures/tests by year.png",
       width = 14, height = 10, dpi = 300)

###----------------------###
# Repeat for Spatial Tests #
###----------------------###

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01/output/spatial")

#read in the stage info
meta <- read.csv("~/OneDrive - University of Southampton/Documents/Chapter 01/data/species_site_stage_metadata.csv")
meta <- meta %>% 
  select(-Site) %>%
  distinct() %>%
  rename(species = Species, stage = Stage, season = Season)

#read in spatial test info
boyce <- readRDS("boyce_final.RDS")
boyce <- boyce %>%
  select(species, site, stage) %>%
  distinct()

#append stage info
boyce <- boyce %>%
  left_join(meta, by = c("species", "stage"))

#remove SOES stages where tracks were taken from other temporal test populations
soes_sites <- c("WAP", "South_Georgia", "Marion")
boyce <- boyce %>%
  filter(species != "SOES" | !site %in% soes_sites)

#null variable for loop
alltracks <- NULL

#for each species/site/stage, list each season where tracks were present
for(i in 1:nrow(boyce)){
  
  #define this species, site, stage, and season
  this.species <- boyce$species[i]
  this.site <- boyce$site[i]
  this.stage <- boyce$stage[i]
  this.season <- boyce$season[i]
  
  #read in tracks
  tracks <- read.csv(paste0(this.species, "/", this.stage, "/", this.site, ".csv"))
  
  #if season is FALSE, summarise the number of years
  if(this.season == FALSE) {
    tracks <- tracks %>% 
      group_by(year) %>% 
      summarise(n_ind = n_distinct(individual_id)) %>%
      rename(season = year)
  }
  
  #if season is TRUE, summarise the number of years when the dates are rounded to the nearest new year
  if(this.season == TRUE){
    tracks <- tracks %>% 
      mutate(season = round_date(as.Date(date), "year")) %>%
      mutate(season = year(season)) %>%
      group_by(season) %>%
      summarise(n_ind = n_distinct(individual_id))
  }
  
  #format table to include species, site, and stage
  tracks <- tracks %>%
    mutate(species = this.species,
           site = this.site,
           stage = this.stage)
  
  #append to all tracks
  alltracks <- bind_rows(alltracks, tracks)
}

#summarise seasons recorded for each species/site/stage
seasons <- alltracks %>%
  filter(!is.na(season))

#rename levels of species
seasons <- seasons %>%
  mutate(species = as.factor(species))
levels(seasons$species) <- c("Adélie Penguin", "Antarctic Fur Seal", "Crabeater Seal", "Grey-headed Albatross", 
                             "Humpback Whale", "Macaroni Penguin", "Southern Elephant Seal")

#rename levels of site
seasons <- seasons %>%
  mutate(site = as.factor(site))
levels(seasons$site) <- c("Amundsen Sea", "East Antarctica", "Bouvet Island", "Campbell Island", "East Antarctica",
                          "Heard Island", "Kerguelen", "Macquarie Island", "Marion", "Mawson Coast", "Raoul Island", 
                          "Ross Sea", "South Orkney", "South Georgia", "South Shetland", "Vestfold Hills", 
                          "Antarctic Peninsula", "Weddell Sea")

#rename levels of stage
seasons <- seasons %>%
  mutate(stage = as.factor(stage))
levels(seasons$stage) <- c("Breeding", "Chick-Rearing", "Early Chick-Rearing",
                           "Incubation", "Late Chick-Rearing", "No Stage", 
                           "Post-Breeding", "Post-Moult")

#change crabeater seal stages to non-central place foraging
seasons <- seasons %>%
  mutate(stage = ifelse(species == "Crabeater Seal", "Non-Central Place Foraging", as.character(stage)))

#add a column to combine species, site and stage
seasons <- seasons %>%
  mutate(test = as.factor(paste(species, site, stage, sep = " - ")),
         season = as.factor(season))

ggplot(seasons, aes(x = season, y = fct_rev(test), fill = n_ind)) +
  geom_tile(aes(height = 0.95, width = 0.96)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Assigned Year", y = "Study Species, Site, and Stage",
       fill = "Number of Tracked Individuals") +
  scale_fill_viridis_c(trans = "log",
                       breaks = c(1, 5, 15, 50)) + 
  theme(axis.text.x=element_text(angle = 0, hjust=0.5),
        axis.text.y=element_text(size=12),
        axis.title = element_text(size = 14),
        legend.position = "bottom")
ggsave(filename = "~/OneDrive - University of Southampton/Documents/Chapter 01/text/figures/spatial tests by year.png",
       width = 14, height = 12, dpi = 300)
