## Create a survey design for the offshore wind area (Sydney Bite)

#load libraries -----
library(tidyverse)
library(sf)
library(ggspatial)
library(marmap)
library(terra)

#function for making a grid
make_centered_waypoints <- function(x, #assumes a m based projection
                                    n = 20,
                                    ncol = NULL,
                                    nrow = NULL,
                                    width = 20*1000,
                                    height = 20*1000,
                                    return_sf = TRUE,
                                    clip = FALSE) {
  
  library(sf)
  library(dplyr)
  
  x <- st_union(x) %>% st_make_valid()
  
  # centroid
  cent <- st_centroid(x)
  cx <- st_coordinates(cent)[1,1]
  cy <- st_coordinates(cent)[1,2]
  
  # determine grid shape
  if (is.null(ncol) & is.null(nrow)) {
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)
  }
  
  if (is.null(nrow)) nrow <- ceiling(n / ncol)
  if (is.null(ncol)) ncol <- ceiling(n / nrow)
  
  # recompute actual number of stations
  n_total <- ncol * nrow
  
  # spacing based on fixed footprint
  cell_w <- width / ncol
  cell_h <- height / nrow
  
  # bounding box centered on centroid
  bbox <- x%>%
    st_bbox()
  
  bbox[1] <- cx - width/2
  bbox[3] <- cx + width/2
  bbox[2] <- cy - height/2
  bbox[4] <- cy + height/2
  
  bbox_poly <- st_as_sfc(bbox)
  
  # grid
  grid <- st_make_grid(
    bbox_poly,
    n = c(ncol, nrow),
    square = TRUE
  )
  
  pts <- st_centroid(grid)
  
  if (clip) {
    pts <- pts[st_within(pts, x, sparse = FALSE), ]
  }
  
  if (!return_sf) {
    pts <- st_coordinates(pts) %>%
      as.data.frame()
  }
  
  attr(pts, "grid_info") <- list(
    n_requested = n,
    n_actual = n_total,
    ncol = ncol,
    nrow = nrow,
    cell_w = cell_w,
    cell_h = cell_h
  )
  
  return(pts)
}

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- 32620 

#load the Offshore wind polygons -----
osw <- read_sf("C:/Users/stanleyr/Documents/GitHub/offfshore_wind/data/shapefiles/Designated_WEAs_25_07_29.shp")%>%
  st_transform(latlong)

#Create a geographic stratification of the MPA  ----------------

sb_osw <- osw %>%
  filter(WEA == "Sydney Bight") %>%
  st_transform(utm)%>%
  st_make_valid() 

wps <- make_centered_waypoints(
  sb_osw,
  ncol = 6,
  nrow = 6,
  width = 44*1000,
  height = 44*1000,
  clip=TRUE
)


#get the depth

# 1. transform to lon/lat
wps_ll <- st_transform(wps, 4326)
coords <- st_coordinates(wps_ll)

# 2. get bathy
bathy <- getNOAA.bathy(
  lon1 = min(coords[,1]) - 1,
  lon2 = max(coords[,1]) + 1,
  lat1 = min(coords[,2]) - 1,
  lat2 = max(coords[,2]) + 1,
  resolution = 1
)

# 3. convert to raster
bathy_r <- marmap::as.raster(bathy)

# 4. extract depth values (NO interactive mode)
depths <- terra::extract(
  terra::rast(bathy_r),
  cbind(coords[,1], coords[,2])
)[,1]

#now assemble the outputs
wps_num <- wps %>%
  st_as_sf()%>%
  mutate(
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2]
  ) %>%
  mutate(
    Xr = round(X, 3),
    Yr = round(Y, 3),
    depth=depths
  ) %>%
  arrange(desc(Xr), Yr) %>%   # EAST → WEST, then SOUTH → NORTH
  mutate(
    station = paste0("sb_osw_", row_number())
  ) %>%
  st_transform(latlong)%>%
  mutate(lon=st_coordinates(.)[,1],
         lat=st_coordinates(.)[,2])%>%
  data.frame()%>%
  dplyr::select(station,depth,lon,lat)


#save outputs
write.csv(wps_num,file = "data/2026 Mission/osw_stations.csv",row.names = FALSE)



