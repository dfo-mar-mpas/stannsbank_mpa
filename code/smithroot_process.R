## Function to process SmithRoot water sampler

smithroot_process <- function(x,longform=FALSE){
  
  # x = file path to the .csv output from smithroot
  # longform = logical specifying simplifed (DEFAULT - FALSE) or complete (TRUE) formatted output for the output
  
  #required libraries
  require(dplyr)
  require(sf)
  require(lubridate)
  
  #Assumed projection
  latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  
 #Data read-in
  data_dim <- read.csv(x,row.names = NULL) #full data read in
  
  data_ind <- which(grepl("Date",data_dim[,1])) #whats the first column of the data file (no headers)
  
  data_df <- read.csv(x,skip=data_ind) #read in just the longform data
  names(data_df) <- gsub("\\..*","",names(data_df)) #clean up the headers
  
  #check for date errors if the machine hasn't been turned on in a while it can catch up with an erroneous old date followed by the current dates
  date_check <- !is.null(data_dim[grepl("Start Time",data_dim[,1]),3])
  
  #location of the sample taken as the centriod of the recorded coordinates
  location <- data_df%>%
              st_as_sf(coords=c("Longitude","Latitude"),crs=latlong)%>%
              summarise(st_union(geometry))%>%
              st_centroid()%>%
              mutate(Longitude = st_coordinates(.)[,1],
                     Latitude = st_coordinates(.)[,2])%>%
              data.frame()%>%
              dplyr::select(Longitude,Latitude)
  
  #extract date
  if(!date_check){date_val <- mean(as.POSIXct(data_df$Date),na.rm=TRUE)} #taken as the mean time of sampling
  
  if(date_check){
    
    time_series <- data_df%>%
                   mutate(Date = as.POSIXct(Date),
                          timestep = difftime(Date,lag(Date,1)))
    
    if(length(which(time_series$timestep>10))>1){ #each step should be 2 seconds buffer with 10 seconds. looking for big breaks
      
      message(paste0("There is an issue with the date in this file (",x,"), date time may be innacurate, please check."))
      
      date_val = NA
      
    }else{
      
      message(paste0("There is an issue with the date in this file (",x,"), date time may be innacurate, please check."))
      
      date_val <- mean(as.POSIXct(time_series[-(which.max(time_series$timestep)-1),"Date"]),na.rm=TRUE) #gets rid of that one erronous out of sequence value at the beginning
      
    }
    
  }
  
  
  #simple output
  if(!longform){
    
    out <- data.frame(date=date_val,
                      longitude=location$Longitude,
                      latitude=location$Latitude,
                      duration = data_dim[grepl("Duration",data_dim[,1]),2],
                      volume = data_dim[grepl("Total",data_dim[,1]),2],
                      distance = data_dim[grepl("Distance",data_dim[,1]),2],
                      peak_pressure = data_dim[grepl("Peak",data_dim[,1]),2],
                      flow = data_dim[grepl("Avg Flow",data_dim[,1]),2],
                      flow_rate = data_dim[grepl("Avg Rate",data_dim[,1]),2],
                      speed = data_dim[grepl("Avg Speed",data_dim[,1]),2])
    
  }
  
  if(longform){
    
    out <- data_df #this could have erroneous dates
    
  }
  
  #return the specified output
  return(out)
  
}
