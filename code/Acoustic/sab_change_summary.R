#code to get the average warming conditions on SAB

#Figure 1 for the Acoustic Report

#load libraries ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(stars)

s2_as_sf = FALSE

#source functions
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/trim_img_ws.R")

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load slope change object from Remi
sab_df <- readRDS('data/Acoustic/CJFAS paper/SABslopes.RDS')

sst_df <- rast(sab_df["sst_slope"])
sst_p <- rast(sab_df["sst_p"])
bt_df <- rast(sab_df["bottom_slope"])
bt_p <- rast(sab_df["bottom_p"])

#just a check since there are no insignificant slopes
sst_df[values(sst_p>=0.05)] <- NA
bt_df[values(bt_p)>=0.05] <- NA

#load the sab polygon and banks -- -
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

sab_banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%st_transform(latlong)

#load bathymetry data --- 

sab_bathy <- rast("data/Bathymetry/TotalSABbathy/wholeBathy.tif")

sab_bathy <- mask(sab_bathy,sab%>%st_transform(st_crs(sab_bathy)))


#trim out bottom temperatures that are deeper than 150m 
shallow_mask <- sab_bathy >= -150

shallow_poly <- as.polygons(shallow_mask, dissolve = TRUE)

shallow_poly <- shallow_poly[shallow_poly$wholeBathy == 1, ]

# Convert to sf
shallow_sf <- st_as_sf(shallow_poly)%>%
              st_make_valid()%>%
              st_transform(st_crs(sst_df))

#convert to terra formats
shallow_vect <- shallow_sf%>%vect()

sab_vect <- sab%>%
            st_transform(st_crs(sst_df))%>%
            vect()

#extract summary stats

#sea surface temp
sst_vals <- extract(
  sst_df,
  sab_vect,
  weights = TRUE,
  cells = TRUE
)

sst_sd <- sst_vals %>%
  group_by(ID) %>%
  summarise(
    wm = sum(sst_slope * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE),
    wsd = sqrt(
      sum(weight * (sst_slope - wm)^2, na.rm = TRUE) /
        (sum(weight, na.rm = TRUE) - 
           sum(weight^2, na.rm = TRUE) / sum(weight, na.rm = TRUE))
    )
  )

#Bottom temp (within threshold)
bt_vals <- extract(
  bt_df,
  sab_vect,
  weights = TRUE,
  cells = TRUE
)

bt_sd <- bt_vals %>%
  group_by(ID) %>%
  summarise(
    wm = sum(bottom_slope * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE),
    wsd = sqrt(
      sum(weight * (bottom_slope - wm)^2, na.rm = TRUE) /
        (sum(weight, na.rm = TRUE) - 
           sum(weight^2, na.rm = TRUE) / sum(weight, na.rm = TRUE))
    )
  )




