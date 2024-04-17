
correlated_random_walk <- function(initial_position, fish_bear, speed, prey_noise_sd, num_steps) {
  
  # Convert initial bearing from degrees to radians
  initial_bearing <- ## Load required libraries
library(leaflet)
library(lwgeom)
library(readr)
library(gganimate)
library(viridis)
library(utils)
library(tidyterra)
library(sf)
library(nngeo)
library(tidyverse)
library(terra)
library(animation)
library(sfheaders)
library(rnaturalearth)
library(geosphere)

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load polygons
sab <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
        st_transform(utm)%>%
        st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
        st_union()%>% #gets rid of the zones
        st_as_sf()

sab_banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%st_transform(utm)

#load coastline and make basemap ------
coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")

basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
           filter(woe_name !="Nova Scotia")%>% #this is in HR with coast_hr
           st_as_sf()%>%
           st_union()%>%
           st_transform(utm)%>%
           st_as_sf()

plot_lims <- sab%>%
             st_buffer(150)%>%
             st_bbox()

#now get coordinates for a grid
sab_grid <- sab%>%
            st_make_grid(cellsize=10,what="centers")%>%
            st_as_sf()%>%
            st_intersection(sab) * pi / 180
  
  # Extract x and y coordinates from the initial position sf object
  start_x <- st_coordinates(initial_position)[, 1]
  start_y <- st_coordinates(initial_position)[, 2]
  
  # Initialize variables
  bearings <- numeric(num_steps)
  x_coords <- numeric(num_steps)
  y_coords <- numeric(num_steps)
  
  # Initialize starting position
  x_coords[1] <- start_x
  y_coords[1] <- start_y
  
  # Generate correlated random walk
  for (i in 2:num_steps) {
    # Generate noise for the current step
    noise <- rnorm(1, mean = 0, sd = prey_noise_sd)
    
    # Compute the correlated bearing for the next step
    correlated_bearing <- (initial_bearing + noise) %% (2 * pi)
    
    # Update positions based on the correlated bearing and speed
    x_coords[i] <- x_coords[i - 1] + speed * cos(correlated_bearing)
    y_coords[i] <- y_coords[i - 1] + speed * sin(correlated_bearing)
    
    # Store the correlated bearing
    bearings[i] <- correlated_bearing
  }
  
  # Return the result as a list of x and y coordinates and bearings
  return(list(x_coords = x_coords, y_coords = y_coords, bearings = bearings))
}





rand_step <- function(bear,start,speed,bound){
  
  require(dplyr)
  require(sf)
  require(geosphere)
  
  #function to identify the next location based on a correlated random walk
  
  #bear - bearing in degrees that the simulated animal will move
  #start - location that the animal is starting from - this is an sf data.frame
  #speed - swim speed
  #bound - this is the bounding polygon to ensure that the position doesn't jump to land
  
  #set up a correlated random bearing 
  
  correlated_bearing <- (bearing  + rnorm(1, mean = 0, sd = prey.noise_sd))%% (2 * pi)
  
  
  
  
}