#functions used in dynamic extraction

#--------------------------------------------
# 1. dynamic_extract - for general use

#extracts variables to points based on the day and year of that point
#requires terra, tidyterra, dplyr, and lubridate
#tracks must have a datetime column called date and be a SpatVector

dynamic_extract <- function(predictor, tracks){
  
  #first create a list of the years within the tracks
  tracks$year <- as.factor(year(tracks$date))
  years <- levels(tracks$year)
  
  #empty variable for loop to feed into
  tracks_extracted <- NULL
  
  #for loop by year - bind = TRUE
  for(z in years){
    trax <- filter(tracks, year==z)
    pred <- rast(paste0("D:/Satellite_Data/daily/", predictor, "/", predictor, "_", z, ".nc"))
    
    e <- ext(trax) + c(0.5,0.5,0.5,0.5) #create SpatExtent for cropping raster
    pred_crop <- crop(pred, e) #crop to increase speed
    
    trax$yday <- as.factor(yday(trax$date)) #extract all yday numbers from data
    ydays <- levels(trax$yday) #different levels of ydays
    
    xtractions <- NULL #create empty list for next loop to feed into
    
    #for loop by yday
    for(i in ydays){
      points <- filter(trax, yday==i) #subsets by yday
      yday.date <- as_date(first(points$date)) #extracts date for yday
      slice <- pred_crop[[time(pred_crop) == yday.date]] #slices raster by yday
      xtracted <- extract(slice, points, ID=F, bind=T) #extract values from slice
      xtracted_df <- as.data.frame(xtracted, geom="XY") #create dataframe for binding
      names(xtracted_df)[length(names(xtracted_df))-2] <- predictor #rename column to predictor name
      xtractions <- rbind(xtractions, xtracted_df) #bind with previous ydays
    }
    
    tracks_extracted <- rbind(tracks_extracted, xtractions) #bind together all years
    
  }
  
  #remove yday column for next predictor to work
  tracks_extracted <- dplyr::select(tracks_extracted, -yday, -year)
  
  #reformat into SpatVector
  tracks_extracted <- vect(tracks_extracted,
                           geom=c("x", "y"),
                           crs=crs(tracks))
  
  return(tracks_extracted)
  
}



#--------------------------------------------
# 2. dynamic_chlorophyll - for unique file structure of chlorophyll

#extracts variables to points based on the day and year of that point
#requires terra, tidyterra, dplyr, and lubridate
#tracks must have a datetime or date column called date and be a SpatVector

dynamic_chlorophyll <- function(predictor, tracks){
  
  #first create a list of the years within the tracks
  tracks$year <- as.factor(year(tracks$date))
  years <- levels(tracks$year)
  
  #empty variable for loop to feed into
  tracks_extracted <- NULL
  
  #for loop by year - bind = TRUE
  for(z in years){
    trax <- filter(tracks, year==z)
    pred <- rast(paste0("D:/Satellite_Data/daily/chl/resampled/chl_", z, "_resampled.nc"))
    
    e <- ext(trax) + c(0.5,0.5,0.5,0.5) #create SpatExtent for cropping raster
    pred_crop <- crop(pred, e) #crop to increase speed
    
    trax$yday <- as.factor(yday(trax$date)) #extract all yday numbers from data
    ydays <- levels(trax$yday) #different levels of ydays
    
    xtractions <- NULL #create empty list for next loop to feed into
    
    #for loop by yday
    for(i in ydays){
      points <- filter(trax, yday==i) #subsets by yday
      yday.date <- as_date(first(points$date)) #extracts date for yday
      slice <- pred_crop[[time(pred_crop) == yday.date]] #slices raster by yday
      xtracted <- extract(slice, points, ID=F, bind=T) #extract values from slice
      xtracted_df <- as.data.frame(xtracted, geom="XY") #create dataframe for binding
      names(xtracted_df)[length(names(xtracted_df))-2] <- predictor #rename column to predictor name
      xtractions <- rbind(xtractions, xtracted_df) #bind with previous ydays
    }
    
    tracks_extracted <- rbind(tracks_extracted, xtractions) #bind together all years
    
  }
  
  #remove yday column for next predictor to work
  tracks_extracted <- dplyr::select(tracks_extracted, -yday, -year)
  
  #reformat into SpatVector
  tracks_extracted <- vect(tracks_extracted,
                           geom=c("x", "y"),
                           crs=crs(tracks))
  
  return(tracks_extracted)
  
}



#--------------------------------------------
# 3. dynamic_chlorophyll_1997 - for unique file structure of chlorophyll in 1997 (starts in September)

#extracts variables to points based on the day and year of that point
#requires terra, tidyterra, dplyr, and lubridate
#tracks must have a datetime column called date and be a SpatVector

dynamic_chlorophyll_1997 <- function(predictor, tracks){
  
  trax <- tracks
  pred <- rast(paste0("D:/Satellite_Data/daily/chl/resampled/chl_1997_resampled.nc"))
  
  e <- ext(trax) + c(0.5,0.5,0.5,0.5) #create SpatExtent for cropping raster
  pred_crop <- crop(pred, e) #crop to increase speed
  
  trax$yday <- as.factor(yday(trax$date)) #extract all yday numbers from data
  ydays <- levels(trax$yday) #different levels of ydays
  
  xtractions <- NULL #create empty list for next loop to feed into
  
  #for loop by yday
  for(i in ydays){
    points <- filter(trax, yday==i) #subsets by yday
    slice <- pred_crop[[as.numeric(i)-246]] #slices raster by yday offset for 1999-09-01 start day
    xtracted <- extract(slice, points, ID=F, bind=T) #extract values from slice
    xtracted_df <- as.data.frame(xtracted, geom="XY") #create dataframe for binding
    names(xtracted_df)[length(names(xtracted_df))-2] <- predictor #rename column to predictor name
    xtractions <- rbind(xtractions, xtracted_df) #bind with previous ydays
  }
  
  tracks_extracted <- xtractions #bind together all years
  
  #remove yday column for next predictor to work
  tracks_extracted <- dplyr::select(tracks_extracted, -yday)
  
  #reformat into SpatVector
  tracks_extracted <- vect(tracks_extracted,
                           geom=c("x", "y"),
                           crs=crs(tracks))
  
  return(tracks_extracted)
  
}




#--------------------------------------------
# 4. dynamic_wind - for monthly resolution of wind

#extracts variables to points based on the month and year of that point
#requires terra, tidyterra, dplyr, and lubridate
#tracks must have a datetime column called date and be a SpatVector

dynamic_wind <- function(predictor, tracks, direction){
  
  #first create a list of the years within the tracks
  tracks$year <- as.factor(year(tracks$date))
  years <- levels(tracks$year)
  
  #empty variable for loop to feed into
  tracks_extracted <- NULL
  
  #for loop by year - bind = TRUE
  for(z in years){
    trax <- filter(tracks, year==z)
    pred <- rast(paste0("D:/Satellite_Data/monthly/wind/", direction, "/", direction, "_resampled_", z, ".nc"))
    
    e <- ext(trax) + c(0.5,0.5,0.5,0.5) #create SpatExtent for cropping raster
    pred_crop <- crop(pred, e) #crop to increase speed
    
    trax$month <- as.factor(month(trax$date)) #extract all month numbers from data
    months <- levels(trax$month) #different levels of months
    
    xtractions <- NULL #create empty list for next loop to feed into
    
    #for loop by month
    for(i in months){
      points <- filter(trax, month==i) #subsets by month
      slice <- pred_crop[[month(time(pred_crop)) == as.numeric(i)]] #slices raster by month
      xtracted <- extract(slice, points, ID=F, bind=T) #extract values from slice
      xtracted_df <- as.data.frame(xtracted, geom="XY") #create dataframe for binding
      names(xtracted_df)[length(names(xtracted_df))-2] <- paste0(predictor, "_", direction) #rename column to predictor name
      xtractions <- rbind(xtractions, xtracted_df) #bind with previous months
    }
    
    tracks_extracted <- rbind(tracks_extracted, xtractions) #bind together all years
    
  }
  
  #remove month column for next predictor to work
  tracks_extracted <- dplyr::select(tracks_extracted, -month, -year)
  
  #reformat into SpatVector
  tracks_extracted <- vect(tracks_extracted,
                           geom=c("x", "y"),
                           crs=crs(tracks))
  
  return(tracks_extracted)
  
}




#--------------------------------------------
# 5. dynamic_wind_1999 - for monthly resolution of wind and later start date in 1999

#extracts variables to points based on the day and year of that point
#requires terra, tidyterra, dplyr, and lubridate
#tracks must have a datetime column called date and be a SpatVector from 1999 data only
#direction can be either east or north

dynamic_wind_1999 <- function(predictor, tracks, direction){
  
  trax <- tracks
  pred <- rast(paste0("D:/Satellite_Data/monthly/wind/", direction, "/", direction, "_resampled_1999.nc"))
  
  e <- ext(trax) + c(0.5,0.5,0.5,0.5) #create SpatExtent for cropping raster
  pred_crop <- crop(pred, e) #crop to increase speed
  
  trax$month <- as.factor(month(trax$date)) #extract all month numbers from data
  months <- levels(trax$month) #different levels of months
  
  xtractions <- NULL #create empty list for next loop to feed into
  
  #for loop by month
  for(i in months){
    points <- filter(trax, month==i) #subsets by month
    slice <- pred_crop[[as.numeric(i)-7]] #slices raster by month offset for Aug start
    xtracted <- extract(slice, points, ID=F, bind=T) #extract values from slice
    xtracted_df <- as.data.frame(xtracted, geom="XY") #create dataframe for binding
    names(xtracted_df)[length(names(xtracted_df))-2] <- paste0(predictor, "_", direction) #rename column to predictor name
    xtractions <- rbind(xtractions, xtracted_df) #bind with previous months
  }
  
  tracks_extracted <- xtractions  
  
  #remove month column for next predictor to work
  tracks_extracted <- dplyr::select(tracks_extracted, -month)
  
  #reformat into SpatVector
  tracks_extracted <- vect(tracks_extracted,
                           geom=c("x", "y"),
                           crs=crs(tracks))
  
  return(tracks_extracted)
  
}