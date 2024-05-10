#dynamic_extract function for chlorophyll file structure
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
      slice <- pred_crop[[as.numeric(i)]] #slices raster by yday
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