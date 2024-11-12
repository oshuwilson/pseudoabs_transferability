#visualisations for schematic
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(terra)
  library(ggplot2)
  library(tidyterra)
  library(CCAMLRGIS)
  library(lubridate)
  library(dplyr)
}

#coast file
coast <- load_Coastline()
coast <- vect(coast)

#model species is ANFS Marion breeding
this.species <- "MAPE"
this.site <- "South_Georgia"
this.stage <- "incubation"

#read in tracks
tracks <- read.csv(paste0("data/tracks_by_stage/", this.species, "/", this.site, "/", this.stage, ".csv"))

#create terra object
tracks <- vect(tracks, geom=c("x", "y"), crs = "epsg:4326")
tracks <- project(tracks, "EPSG:6932")
plot(tracks, pch=".")

#read in pseudoabs
buffer <- read.csv(paste0("output/buffers/", this.species, "/", this.site, "/", this.stage, ".csv"))
background <- read.csv(paste0("output/background/", this.species, "/", this.site, "/", this.stage, ".csv"))
crw <- read.csv(paste0("output/CRWs/", this.species, "/", this.site, "/", this.stage, ".csv"))

#create terra objects
buffer <- vect(buffer,
               geom = c("x", "y"),
               crs = "epsg:4326")
buffer <- project(buffer, "EPSG:6932")

background <- vect(background,
                   geom = c("x", "y"),
                   crs = "epsg:4326")
background <- project(background, "EPSG:6932")

crw <- vect(crw,
            geom = c("x", "y"),
            crs = "epsg:4326")
crw <- project(crw, "EPSG:6932")

#crop coast to fit
e <- ext(crw) 
crop_coast <- crop(coast, e)

#plots
ggplot() +
  geom_spatvector(data = buffer, size=0.1, col = "cadetblue4") + 
  geom_spatvector(data = tracks, size=0.1, col = "black") +
  geom_spatvector(data=crop_coast, fill="grey") + 
  theme_void() + 
  xlim(-3066418.83637174, -2183192.10410952) +
  ylim(2858741.12954235, 3718285.10207138)

ggsave(filename = "text/figures/pseudo-abs/buffer.png",
       width = 8, height = 6)

ggplot() +
  geom_spatvector(data = background, size=0.1, col = "cadetblue") + 
  geom_spatvector(data = tracks, size=0.1, col = "black") +
  geom_spatvector(data=crop_coast, fill="grey") + 
  theme_void() +
  xlim(-3066418.83637174, -2183192.10410952) +
  ylim(2858741.12954235, 3718285.10207138)

ggsave(filename = "text/figures/pseudo-abs/background.png",
       width = 8, height = 6)

ggplot() +
  geom_spatvector(data = crw, size=0.1, col = "cadetblue") + 
  geom_spatvector(data = tracks, size=0.1, col = "black") +
  geom_spatvector(data=crop_coast, fill="grey") + 
  theme_void() +
  xlim(-3066418.83637174, -2183192.10410952) +
  ylim(2858741.12954235, 3718285.10207138)

ggsave(filename = "text/figures/pseudo-abs/crw.png",
       width = 8, height = 6)

#South Georgia S/O Map
bbox <- ext(-3066418.83637174, -2183192.10410952, 2858741.12954235, 3718285.10207138)
bbox <- as.polygons(bbox, crs = crs(tracks))

ggplot() + 
  geom_spatvector(data = coast, fill="white") + theme_void() +
  geom_spatvector(data = bbox, fill=NA, col="red3") + 
  xlim(-5180778.6221, 5180778.6221) +
  ylim(-5180778.6221, 5180778.6221) + 
  theme(plot.background = element_rect(fill = "lightblue2"))

ggsave(filename = "text/figures/pseudo-abs/inset.svg",
       width = 8, height = 6)

ext(coast)
