#Functions to read BNAM - matlab data and generate rasters of a specified resolution. 

bnam_raster_gen <- function(x,region=NULL,res_val=0.09,crs_val){
  
  #function to pull data from the matlab file, extract at a specified resolution (0.09 degrees), and write a raster
  
  #x is the file path for the matlab file. 
  #region is a shape files (sf) that is used to get the extent. if NULL no trimming is conducted (default). 
  #res_val is the resolution of the raster grid (DEFAULT = 0.09). 
  #crs_val is the decimal degree projection of BNAM
  
  #load libraries
  require(R.matlab)
  require(raster)
  require(dplyr)
  require(sp)
  require(sf)
  
  #progress message
  message(paste0("Reading in matlab file ",x))
  
  #load the matlab data
  BNAM_MAT <- readMat(x)
  
  #pull coordinates (based on the fact that long and lat are the first elements in these matlab containers)
  lon <- BNAM_MAT[[1]] %>% `dim<-`(.,NULL)  #flatten longitude array 
  lat <- BNAM_MAT[[2]] %>% `dim<-`(.,NULL)  #flatten latitude array 
  
  array_dim <- length(lon)
  #same as the dimensions of the matrix here - array_dim == dim(BNAM_MAT[[10]])[2]*dim(BNAM_MAT[[10]])[3]
  
  #for this function we will extract the annual averages but this step can be modified
  select_names <- names(BNAM_MAT)[grepl("ann",names(BNAM_MAT))]
  
  #progress message
  message(paste("Extracting",paste(select_names,collapse = " and "),sep=" "))
  
  #extract data and write a geotiff for each layer (layer loop)
  for(i in 1:length(select_names)){
    
    #create an xyz data.frame
    extract_df <- data.frame(lon=lon,lat=lat,z=as.vector(BNAM_MAT[[select_names[[i]]]]))
    
    if(!is.null(region)){
      
      focus_extent <- region%>%
        st_bbox()%>%
        as.numeric()
      
      #trim the data down to the focal region
      extract_df <- extract_df%>%
        filter(between(lon,focus_extent[1],focus_extent[3]))%>%
        filter(between(lat,focus_extent[2],focus_extent[4]))
      
    } #end of region clip
    
    #create a grid to do the extraction
    grid <- extract_df%>%st_as_sf(coords=c("lon","lat"),crs=crs_val)%>%raster::extent(.)%>%raster::raster(.)
    raster::res(grid) <- res_val
    
    #now do the grid extraction
    out <- extract_df%>% #this takes a second. fasterize' is an option but I find it buggy
      st_as_sf(coords=c("lon","lat"),crs=crs_val)%>%
      st_coordinates()%>%
      raster::rasterize(.,grid,field=extract_df[,3],fun=mean)%>%
      raster::rasterToPoints(.)%>%data.frame()%>%
      rename(lon = x, lat = y, MAP = layer)
    
    #convert to a raster
    raster_out=rasterFromXYZ(out,crs=crs_val)
    
    #directory check
    if(!dir.exists("output/bnam_rasters/")){dir.create("output/bnam_rasters/")}
    
    #file name for the geotiff
    filename <- paste0("output/bnam_rasters/",paste0(gsub("[.]","_",select_names[i]),".tif"))
    
    #progress message
    message(paste0("Writing ",filename))
    
    #save raster
    writeRaster(raster_out,filename,overwrite=TRUE)
    
  } #end of the layer loop
  
} #end of function

sum_rasters <- function(file1,file2,outputfile){
  
  require(dplyr)
  require(raster)
  
  #file1/file2 - full file paths for rasters to be added
  
  #read in the rasters
  ras1 <- raster(file1)
  ras2 <- raster(file2)
  
  #stack and combine
  rs <- stack(list(ras1,ras2))%>%calc(.,sum)
  
  writeRaster(rs,outputfile)
  
}

bnam_raster_gen_rcp_season <- function(x,region=NULL,res_val=0.09,crs_val,rasterfun = min,funname="min",output_dir){
  
  #function to pull data from the matlab file, extract at a specified resolution (0.09 degrees), and write a raster
  
  #x is the file path for the matlab file. 
  #region is a shape files (sf) that is used to get the extent. if NULL no trimming is conducted (default). 
  #res_val is the resolution of the raster grid (DEFAULT = 0.09). 
  #crs_val is the decimal degree projection of BNAM
  #rasterfun is a function that is used for the grid for each season. Must be mean, min, max, or med. -- something that can be applied to derive one value from a vector
  #funname is the name of the function applied (character) ** I can't figure a way to convert 'mean' into a character 
  #output_dir is the directory to store the rasters
  
  #load libraries
  require(R.matlab)
  require(raster)
  require(dplyr)
  require(sp)
  require(sf)
  
  #progress message
  message(paste0("Reading in matlab file ",x))
  
  #load the matlab data
  BNAM_MAT <- readMat(x)
  
  #pull coordinates (based on the fact that long and lat are the first elements in these matlab containers)
  lon <- BNAM_MAT[[1]] %>% `dim<-`(.,NULL)  #flatten longitude array 
  lat <- BNAM_MAT[[2]] %>% `dim<-`(.,NULL)  #flatten latitude array 
  
  array_dim <- length(lon)
  #same as the dimensions of the matrix here - array_dim == dim(BNAM_MAT[[10]])[2]*dim(BNAM_MAT[[10]])[3]
  
  #for this function we will extract the annual averages but this step can be modified
  select_names <- names(BNAM_MAT)[grepl("rcp",tolower(names(BNAM_MAT)))]
  select_names <- select_names[!grepl("ann",select_names)] #get rid of the seasonal data
  
  #progress message
  message(paste("Extracting",paste(select_names,collapse = " and "),sep=" "))
  
  #extract data and write a geotiff for each layer (layer loop)
  for(i in 1:length(select_names)){
    
    month_extract <- BNAM_MAT[[select_names[[i]]]]
    
    month_extract_flat <- NULL
    for(j in 1:12){ #extract monthly data and put it into a long, flattened data.frame
      
      #create an xyz data.frame
      month_extract_flat <- rbind(month_extract_flat,
                                  data.frame(lon=lon,lat=lat,mon=j,z=as.vector(month_extract[j,,])))
      
    } #end month loop
    
    if(!is.null(region)){
      
      focus_extent <- region%>%
        st_bbox()%>%
        as.numeric()
      
      #trim the data down to the focal region
      month_extract_flat <- month_extract_flat%>%
        filter(between(lon,focus_extent[1],focus_extent[3]))%>%
        filter(between(lat,focus_extent[2],focus_extent[4]))
      
    } #end of region clip
    
    month_extract_flat <- month_extract_flat%>%
      mutate(season = case_when(mon %in% 1:3 ~ "winter",
                                mon %in% 4:6 ~ "spring",
                                mon %in% 7:9 ~ "summer",
                                mon %in% 10:12 ~"fall"))
    
    month_extract_flat <- rbind(month_extract_flat,
                                month_extract_flat%>%mutate(season="annual"))# this is so the annual average can be taken from the same flattened data without needing logical code in the raster gen part below. 
    
    #create a grid to do the extraction
    grid <- month_extract_flat%>%
      filter(mon==1)%>% #only need the grid based on one month since all months are the same
      st_as_sf(coords=c("lon","lat"),crs=crs_val)%>%
      raster::extent(.)%>%
      raster::raster(.)
    
    #set the resolution of the grid to the function specified variable 'res_val'
    raster::res(grid) <- res_val
    
    for(j in c("winter","spring","summer","fall","annual")){
      
      #the file you are working with
      fn <- strsplit(x,"./")%>%unlist()
      fn <- fn[length(fn)]
      
      #variable you are working with
      varnam <- strsplit(select_names[i],"[.]")%>%unlist()
      varnam <- gsub("d","",varnam[length(varnam)])
      
      #the RCP scenario you are working with
      rcp <- strsplit(x,"/")%>%
        unlist()%>%
        data.frame(names=.)%>%
        filter(grepl("rcp",tolower(names)))%>%
        mutate(names=gsub(" ","",names),
               names=tolower(names),
               names=gsub("\\.","",names))%>%
        pull(names)
      
      
      message(paste0("Working on ",j," ",varnam,":  ",gsub(".mat","",fn)))
      
      #now do the grid extraction
      out <- month_extract_flat%>%#this takes a second. fasterize' is an option but I find it buggy
        filter(season == j)%>%
        st_as_sf(coords=c("lon","lat"),crs=crs_val)%>%
        st_coordinates()%>%
        raster::rasterize(.,grid,field=month_extract_flat%>%filter(season == j)%>%pull(z),fun=rasterfun,na.omit=TRUE)%>%
        raster::rasterToPoints(.)%>%
        data.frame()%>%
        rename(lon = x, lat = y, MAP = layer)%>%
        suppressWarnings()# will return warnings for the NaN 'land' values
      
      #convert to a raster
      raster_out=rasterFromXYZ(out,crs=crs_val)
      
      #directory check
      if(!dir.exists(output_dir)){dir.create(output_dir)}
      
      #file name for the geotiff
      filename <- paste0(output_dir,j,"_",funname,"_",varnam,"_",rcp,".tif")
      
      #progress message
      message(paste0("Writing ",filename))
      
      #save raster
      writeRaster(raster_out,filename,overwrite=TRUE)
      
    }#end raster season loop
    
  } #end of the layer loop
  
} #end of function

bnam_raster_gen_season <- function(x,region=NULL,res_val=0.09,crs_val,rasterfun = min,funname="min",output_dir){
  
  #function to pull data from the matlab file, extract at a specified resolution (0.09 degrees), and write a raster
  
  #x is the file path for the matlab file. 
  #region is a shape files (sf) that is used to get the extent. if NULL no trimming is conducted (default). 
  #res_val is the resolution of the raster grid (DEFAULT = 0.09). 
  #crs_val is the decimal degree projection of BNAM
  #rasterfun is a function that is used for the grid for each season. Must be mean, min, max, or med. -- something that can be applied to derive one value from a vector
  #funname is the name of the function applied (character) ** I can't figure a way to convert 'mean' into a character 
  #output_dir is the directory to store the rasters
  
  #load libraries
  require(R.matlab)
  require(raster)
  require(dplyr)
  require(sp)
  require(sf)
  
  #progress message
  message(paste0("Reading in matlab file ",x))
  
  #load the matlab data
  BNAM_MAT <- readMat(x)
  
  #pull coordinates (based on the fact that long and lat are the first elements in these matlab containers)
  lon <- BNAM_MAT[[1]] %>% `dim<-`(.,NULL)  #flatten longitude array 
  lat <- BNAM_MAT[[2]] %>% `dim<-`(.,NULL)  #flatten latitude array 
  
  array_dim <- length(lon)
  #same as the dimensions of the matrix here - array_dim == dim(BNAM_MAT[[10]])[2]*dim(BNAM_MAT[[10]])[3]
  
  #for this function we will extract the annual averages but this step can be modified
  select_names <- names(BNAM_MAT)[grepl("mon",names(BNAM_MAT))]
  
  #progress message
  message(paste("Extracting",paste(select_names,collapse = " and "),sep=" "))
  
  #extract data and write a geotiff for each layer (layer loop)
  for(i in 1:length(select_names)){
    
    month_extract <- BNAM_MAT[[select_names[[i]]]]
    
    month_extract_flat <- NULL
    for(j in 1:12){ #extract monthly data and put it into a long, flattened data.frame
      
      #create an xyz data.frame
      month_extract_flat <- rbind(month_extract_flat,
                                  data.frame(lon=lon,lat=lat,mon=j,z=as.vector(month_extract[j,,])))
      
    } #end month loop
    
    if(!is.null(region)){
      
      focus_extent <- region%>%
        st_bbox()%>%
        as.numeric()
      
      #trim the data down to the focal region
      month_extract_flat <- month_extract_flat%>%
        filter(between(lon,focus_extent[1],focus_extent[3]))%>%
        filter(between(lat,focus_extent[2],focus_extent[4]))
      
    } #end of region clip
    
    month_extract_flat <- month_extract_flat%>%
      mutate(season = case_when(mon %in% 1:3 ~ "winter",
                                mon %in% 4:6 ~ "spring",
                                mon %in% 7:9 ~ "summer",
                                mon %in% 10:12 ~"fall"))
    
    month_extract_flat <- rbind(month_extract_flat,
                                month_extract_flat%>%mutate(season="annual"))# this is so the annual average can be taken from the same flattened data without needing logical code in the raster gen part below. 
    
    
    #create a grid to do the extraction
    grid <- month_extract_flat%>%
      filter(mon==1)%>% #only need the grid based on one month since all months are the same
      st_as_sf(coords=c("lon","lat"),crs=crs_val)%>%
      raster::extent(.)%>%
      raster::raster(.)
    
    
    raster::res(grid) <- res_val
    
    for(j in c("winter","spring","summer","fall","annual")){
      
      message(paste0("Working on ",j," ",gsub(".mat","",strsplit(x,"./")%>%unlist()%>%.[4])))
      
      #now do the grid extraction
      out <- month_extract_flat%>%#this takes a second. fasterize' is an option but I find it buggy
        filter(season == j)%>%
        st_as_sf(coords=c("lon","lat"),crs=crs_val)%>%
        st_coordinates()%>%
        raster::rasterize(.,grid,field=month_extract_flat%>%filter(season == j)%>%pull(z),fun=rasterfun,na.omit=TRUE)%>%
        raster::rasterToPoints(.)%>%
        data.frame()%>%
        rename(lon = x, lat = y, MAP = layer)%>%
        suppressWarnings()# will return warnings for the NaN 'land' values
      
      #convert to a raster
      raster_out=rasterFromXYZ(out,crs=crs_val)
      
      #directory check
      if(!dir.exists(output_dir)){dir.create(output_dir)}
      
      #file name for the geotiff
      filename <- paste0(output_dir,j,"_",funname,"_",paste0(gsub(".mat","",strsplit(x,"./")%>%unlist()%>%.[4]),".tif"))
      
      #progress message
      message(paste0("Writing ",filename))
      
      #save raster
      writeRaster(raster_out,filename,overwrite=TRUE)
      
    }#end raster season loop
    
  } #end of the layer loop
  
} #end of function