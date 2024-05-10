#script for extracting environmental variables to pseudo and pseudoabsences
#Finished ADPE, GHAL, EMPE, MAPE, ANPE, SUFS, ANFS, CRAS, SOES marion + south_georgia, and HUWH
#SOES WAP POST-MOULT BACKGROUND STILL NEEDS CHL AND WIND
#DOES SLOPE NEED TO BE RESAMPLED??? (i.e. from GEBCO and resampled to 0.083deg)
#only SOES WAP left - really slow for some reason
#rerun slope for all but soes wap crws post-moult

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(terra)
  library(dplyr)
  library(lubridate)
  library(tidyterra)
  library(ggplot2)
}

#setup
rm(list=ls())
this.species <- "SOES"
this.site <- "WAP"
this.stage <- "post-moult"
pseudo.type <- "background" #either background, buffers, or CRWs

#read in pseudo
pseudo <- read.csv(paste0("output/", pseudo.type, "/", this.species, "/", this.site, "/", this.stage, ".csv"))

#only keep relevant columns
pseudo <- pseudo %>% select(-X)

#format date (for temporal extraction later)
pseudo$date <- as_datetime(pseudo$date)
min(pseudo$date) #important for wind and chlorophyll - check dates if before 2000

#---------------
#Static Variables

###Depth###
depth <- rast("D:/Satellite_Data/static/depth/depth.nc")

#create SpatVector for pseudo
pseudo <- vect(pseudo,
               geom=c("x", "y"),
               crs=crs(depth)) #this ensures crs are the same as rasters

#extract
pseudo$depth <- extract(depth, pseudo, ID=F)
ggplot(pseudo, aes(x=depth)) + geom_density()

#remove rows where depth is NA - will be NA for every variable
plot(pseudo, pch=".")
pseudo <- pseudo %>% drop_na(depth)
plot(pseudo, pch=".")


###Slope###
slope <- rast("D:/Satellite_Data/static/slope/slope.nc")
pseudo$slope <- extract(slope, pseudo, ID=F)
ggplot(pseudo, aes(x=slope)) + geom_density()

###dShelf###
dshelf <- rast("D:/Satellite_Data/static/dshelf/dshelf_resampled.nc")
pseudo$dshelf <- extract(dshelf, pseudo, ID=F)
ggplot(pseudo, aes(x=dshelf)) + geom_density()

#cleanup static
rm(depth, slope, dshelf)


#---------------
#Dynamic Variables

#dynamic_extract function from 05a script
source("code/05a_dynamic_extract_function.R")


###SST###
pseudo <- dynamic_extract(predictor = "sst", pseudo)
ggplot(pseudo, aes(x=sst)) + geom_density()

###MLD###
pseudo <- dynamic_extract(predictor = "mld", pseudo)
ggplot(pseudo, aes(x=mld)) + geom_density()

###SAL###
pseudo <- dynamic_extract(predictor = "sal", pseudo)
ggplot(pseudo, aes(x=sal)) + geom_density()

###SSH###
pseudo <- dynamic_extract(predictor = "ssh", pseudo)
ggplot(pseudo, aes(x=ssh)) + geom_density()

###SIC###
pseudo <- dynamic_extract(predictor = "sic", pseudo)
pseudo$sic[is.na(pseudo$sic)] <- 0 #SIC values of 0 print as NA in GLORYS
ggplot(pseudo, aes(x=sic)) + geom_density()

###CURR###
pseudo <- dynamic_extract(predictor = "uo", pseudo) #eastward velocity
pseudo <- dynamic_extract(predictor = "vo", pseudo) #northward velocity
pseudo$curr <- sqrt((pseudo$uo^2) + (pseudo$vo^2)) #current speed
ggplot(pseudo, aes(x=curr)) + geom_density()

###EKE###
pseudo$eke <- 0.5 * ((pseudo$uo^2) + (pseudo$vo^2))
ggplot(pseudo, aes(x=eke)) + geom_density()

###CHL### 
#Does not work for data before 04-09-1997 and needs adjusted function for other dates in 1997
#only use hashtagged code if data starts in 1997 or earlier
# 
cutoff_97 <- as_date("1997-09-04")
cutoff_98 <- as_date("1998-01-01")
pseudo_pre97 <- pseudo %>% filter(date < cutoff_97)
plot(pseudo_pre97)
pseudo_97 <- pseudo %>% filter(date >= cutoff_97 & date < cutoff_98)
plot(pseudo_97)
pseudo <- pseudo %>% filter(date >= cutoff_98)
plot(pseudo, pch=".")


source("code/05b_dynamic_chlorophyll_function.R") #unique function for different file structure

pseudo <- dynamic_chlorophyll(predictor = "chl", pseudo)
ggplot(pseudo, aes(x=chl)) + geom_density()

source("code/05d_dynamic_chlorophyll_1997_function.R")
pseudo_97 <- dynamic_chlorophyll_1997(predictor="chl", pseudo_97)
ggplot(pseudo_97, aes(x=chl)) + geom_density()

pseudo <- bind_spat_rows(pseudo, pseudo_97, pseudo_pre97)
min(pseudo$date)
rm(cutoff_97, cutoff_98, pseudo_pre97, pseudo_97)


###WIND###
#Does not work for data before 01-08-1999 and needs adjusted function for other dates in 1999
#only use hashtagged code if data starts in 1999 or earlier

cutoff_99 <- as_date("1999-08-01")
cutoff_00 <- as_date("2000-01-01")
pseudo_pre99 <- pseudo %>% filter(date < cutoff_99)
plot(pseudo_pre99)
pseudo_99 <- pseudo %>% filter(date >= cutoff_99 & date < cutoff_00)
#plot(pseudo_99)
pseudo <- pseudo %>% filter(date >= cutoff_00)
plot(pseudo, pch=".")

source("code/05c_dynamic_wind_function.R")
pseudo <- dynamic_wind(predictor = "wind", pseudo, direction = "east")
pseudo <- dynamic_wind(predictor = "wind", pseudo, direction = "north")
pseudo$wind <- sqrt(pseudo$wind_east^2 + pseudo$wind_north^2)
ggplot(pseudo, aes(x=wind)) + geom_density()

# source("code/05e_dynamic_wind_1999_function.R")
# pseudo_99 <- dynamic_wind_1999(predictor="wind", pseudo_99, direction = "east")
# pseudo_99 <- dynamic_wind_1999(predictor="wind", pseudo_99, direction = "north")
# pseudo_99$wind <- sqrt(pseudo_99$wind_east^2 + pseudo_99$wind_north^2)
# ggplot(pseudo_99, aes(x=wind)) + geom_density()
# 
pseudo <- bind_spat_rows(pseudo, pseudo_99, pseudo_pre99)
min(pseudo$date)
rm(cutoff_99, cutoff_00, pseudo_pre99, pseudo_99)

#---------------
#Export
plot(pseudo, pch=".")
pseudo <- as.data.frame(pseudo, geom="XY")

write.csv(pseudo, 
          file=paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/", pseudo.type, ".csv"))

