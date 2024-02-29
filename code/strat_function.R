stratFun <- function(dem,min_depth=NULL,max_depth=NULL){
  
  require(stars)
  require(sf)
  require(raster)
  
  x <- dem
  
  if(is.null(max_depth) & !is.null(min_depth)){
    
    x[!x>min_depth] <- NA
    ind <- !is.na(x)
    x[!ind] <- NA
    x[ind] <- 1
    x <- x > -Inf
    
    x <- st_as_sf(stars::st_as_stars(x),as_points = FALSE, merge = TRUE) 
    
  }
  
  if(!is.null(max_depth) & !is.null(min_depth)){
    
    x[x>=min_depth] <- NA
    x[x<=max_depth] <- NA
    ind <- !is.na(x)
    x[!ind] <- NA
    x[ind] <- 1
    x <- x > -Inf
    
    x <- st_as_sf(stars::st_as_stars(x),as_points = FALSE, merge = TRUE) 
    
  }
  
  if(!is.null(max_depth) & is.null(min_depth)){
    
    x[x>max_depth] <- NA
    ind <- !is.na(x)
    x[!ind] <- NA
    x[ind] <- 1
    x <- x > -Inf
    
    x <- st_as_sf(stars::st_as_stars(x),as_points = FALSE, merge = TRUE) 
    
  }
  
  return(x)
  
}