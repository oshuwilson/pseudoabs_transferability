#dynamic_extract function for wind - 1 AUG to 31 DEC 1999 ONLY
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