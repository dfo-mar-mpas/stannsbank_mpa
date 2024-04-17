## Load required libraries
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(patchwork)

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load polygons
sab <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab_banks <- read_sf("data/Shapefiles/sab_banks.shp")

#load coastline and make basemap ------
coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")%>%
            st_transform(latlong)

basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
  filter(woe_name !="Nova Scotia")%>% #this is in HR with coast_hr
  st_as_sf()%>%
  st_union()%>%
  st_transform(utm)%>%
  st_as_sf()%>%
  st_transform(latlong)


#load acoustic stations
sab_array <- read.csv("data/Acoustic/OTN_redesign_coords.csv")%>%
             st_as_sf(coords=c("long","lat"),crs=latlong)

sab_new <- read.csv("data/Acoustic/2024_Additional_Station_Design.csv")%>%
           st_as_sf(coords=c("lon","lat"),crs=latlong)


plot_lims <- sab_array%>%
  st_transform(utm)%>%
  st_buffer(10)%>%
  st_transform(latlong)%>%
  st_bbox()

plot_lims_big <- sab%>%
                 st_transform(utm)%>%
                 st_buffer(50)%>%
                 st_transform(latlong)%>%
                 st_bbox()

#now get coordinates for a grid
p1 <- ggplot()+
  geom_sf(data=sab,fill=NA)+
  geom_sf(data=sab%>%filter(Zone==1),fill="grey90")+
  geom_sf(data=sab_banks,fill="grey20")+
  geom_sf(data=sab_array)+
  geom_sf(data=sab_new,col="red",size=2)+
  theme_bw()+
  coord_sf(expand=0,xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)])

p2 <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=coast_hr)+
  geom_sf(data=sab,fill=NA)+
  geom_sf(data=sab%>%filter(Zone==1),fill="grey90")+
  geom_sf(data=sab_banks,fill="grey20")+
  geom_sf(data=sab_array)+
  geom_sf(data=sab_new,col="red",size=2)+
  geom_sf(data=plot_lims%>%st_as_sfc(),fill=NA,lty=2)+
  theme_bw()+
  coord_sf(expand=0,xlim=plot_lims_big[c(1,3)],ylim=plot_lims_big[c(2,4)])+
  annotation_scale()
  
p3 <- p1 + p2 + plot_layout(ncol=2)

ggsave("output/Acoustic/2024_newstations.png",p3,width=10,height=7,units="in",dpi=300)
