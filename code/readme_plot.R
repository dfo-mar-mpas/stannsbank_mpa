## Code to make plots for the first St. Anns Bank MPA Working Group Workshop - November 16-17 2022 Membertou
#Make plot for the github readme

#Load libraries ----
library(ggplot2)
library(raster)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(stars)
library(tidyr)
library(ggspatial)

sf_use_s2=FALSE

#Projections ------------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm_mar <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load the St. Anns Bank MPA Polygon
sab <- read_sf("data/Shapefiles/StAnnsBank/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab_nozones <- sab%>%
  st_transform(utm_mar)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

bioregion <- read_sf("data/Shapefiles/MaritimesPlanningArea.shp")%>%st_transform(latlong)

bioregion_box <- bioregion%>%
  st_transform(utm_mar)%>%
  st_buffer(5)%>%
  st_transform(latlong)%>%
  st_bbox()%>%
  st_as_sfc()

plotlims <- c(-60.2,45.7,-58.25,46.5)

#basemap
basemap_atlantic <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="Canada"),
                          ne_states(country = "United States of America",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="USA"))%>%
  st_intersection(.,bioregion_box)

basemap_hr <- read_sf("data/Shapefiles/SAB_hr_coast.shp")%>%
              st_transform(latlong)

#load high resolution bathymetry (50m) from multibeam 
sab_bathy <- raster("data/Bathymetry/bathy50/w001001.adf", RAT=FALSE)%>%
             projectRaster(.,crs=latlong)

sab_benthoscape <- read_sf("data/Shapefiles/StAnns_Benthoscape/benthoscape/benthoscape.shp")%>%st_transform(latlong)

#load the eDNA sampling points
sample_coords <- read.csv("data/sampling_2022_coords.csv")%>%
                 filter(grepl("SAB",station))%>%
                 dplyr::select(station,lon,lat)%>%
                 st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                 mutate(type="eDNA")

#load reciever locations
reciever_coords <- read.csv("data/OTN_redesign_coords.csv")%>%
                   rename(lon=long)%>%
                   st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                   mutate(type="Acoustic recievers")


#make the map
p1 <- ggplot()+
      geom_sf(data=basemap_hr)+
      geom_sf(data=sab_benthoscape,aes(fill=Assigned_c),show.legend = FALSE)+
      geom_sf(data=sab,fill=NA)+
      geom_sf(data=reciever_coords,shape=20)+
      geom_sf(data=sample_coords,shape=19)+
      geom_sf(data=sample_coords%>%filter(station == "SAB_trans_1"),shape=21,fill="white",size=1.5)+
      coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
      theme_bw()+
      theme(legend.position = "none")+
      annotation_scale()

ggsave("output/sab_readme_plot.png",p1,width=6,height=6,dpi=300,units="in")      
