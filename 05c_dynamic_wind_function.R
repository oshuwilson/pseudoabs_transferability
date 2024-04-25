#dynamic_extract function for wind file structure
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
      slice <- pred_crop[[as.numeric(i)]] #slices raster by month
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