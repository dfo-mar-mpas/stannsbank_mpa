## map of imagery locations -- 2009 Campod survey 

#load libraries
library(ggplot2)
library(raster)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(stars)
library(tidyr)
library(ggspatial)
library(xlsx)
library(ggnewscale)

sf_use_s2=FALSE

#Projections ------------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm_mar <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load the St. Anns Bank MPA Polygon
sab <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab_nozones <- sab%>%
  st_transform(utm_mar)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()%>%
  mutate(name="St. Anns Bank MPA")

bounding_area <- sab%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_sf()%>%
  st_buffer(0.5)%>%# 1 degree buffer
  st_bbox()%>%
  st_as_sfc()%>%
  st_sf()%>%
  suppressWarnings()%>%
  suppressMessages()

plotlims <- c(-60.2,45.7,-58.25,46.5)

bioregion <- read_sf("data/Shapefiles/MaritimesPlanningArea.shp")%>%st_transform(latlong)

bioregion_box <- bioregion%>%
  st_transform(utm_mar)%>%
  st_buffer(5)%>%
  st_transform(latlong)%>%
  st_bbox()%>%
  st_as_sfc()

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

#load high resolution bathymetry (50m) from multibeam 
sab_bathy <- raster("data/Bathymetry/bathy50/w001001.adf", RAT=FALSE)%>%
  projectRaster(.,crs=latlong)

sab_benthoscape <- read_sf("data/Shapefiles/benthoscape.shp")%>%
  st_transform(latlong)%>%
  mutate(class = Assigned_c,
         class = gsub("A - ","",class), #clean up the classification labels
         class = gsub("Asp - ","",class),
         class = gsub("B - ","",class),
         class = gsub("C - ","",class),
         class = gsub("D - ","",class),
         class = gsub("E - ","",class),
         class = gsub("F - ","",class))

#load bathymetry 
ras <- raster("data/Bathymetry/sab_dem.tif")
rasproj <- proj4string(ras)

#load coordinates - downloaded form https://data.mendeley.com/datasets/vhj7cg9y68/2
photo_coords <- read_sf("data/shapefiles/Kenchington2023_SAB_MPA_Photo_Archive_Metadata.shp")%>%
                   st_transform(latlong)%>%
                   mutate(type="photo")

photo_centroid <- photo_coords%>%
                  group_by(Station)%>%
                  summarise(geometry = st_union(geometry)) %>%
                  st_centroid()%>%
                  st_as_sf() # get the centroid of point clusters

video_coords <- read_sf("data/shapefiles/Kenchington2023_SAB_MPA_Video_Transect_Metadata.shp")%>%
                 st_transform(latlong)%>%
                 mutate(type="video")

imagery_coords <- rbind(photo_coords%>%dplyr::select(Station,type,geometry),video_coords%>%dplyr::select(Station,type,geometry))

sab_imagery_plot <- ggplot()+
                    geom_sf(data=basemap_atlantic)+
                    geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
                    geom_sf(data=sab,fill=NA)+
                    # new_scale_fill() +
                    # geom_sf(data=imagery_coords,aes(fill=type),pch=21,size=3,col="black")+
                    geom_sf(data=photo_centroid,pch=21,fill="white",col="black",size=3)+
                    coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
                    theme_bw()+
                    theme(legend.position = "bottom",
                          axis.text = element_blank())+
                    labs(fill="")+
                    guides(fill = guide_legend(nrow = 3))+
                    annotation_scale();sab_imagery_plot

ggsave("output/sab_imagery_plot.png",sab_imagery_plot,height=6,width=7,units="in",dpi=300)

