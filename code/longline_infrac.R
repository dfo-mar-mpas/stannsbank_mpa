#load libraries ----
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)
library(ggplot2)
library(raster)
library(stars)
library(ggspatial)
library(viridis)

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

##SAB
sab_zones <- read_sf("Data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

sab_bounds <- st_bbox(sab)

#SAB Benthoscape
sab_benthoscape <- read_sf("Data/Shapefiles/benthoscape.shp")%>%
  st_transform(latlong)%>%
  rename(class=Assigned_c)%>%
  dplyr::select(class,geometry)%>%
  st_make_valid()%>%
  group_by(class)%>%
  summarise()

#infraction coordiantes
46 03.49 N, 58 51.01 W

infrac_loc <- data.frame(lat=46+(3.49/60),lon = (58 + (51.01/60))*-1)%>%
              st_as_sf(coords=c("lon","lat"),crs=latlong)

%>%
              st_intersection(.,sab_benthoscape%>%dplyr::select(class,geometry))

wolffish <- read_sf("M:/OCMD/02_Projects/MPA/NetworkPlanning/ConservationPriorities/Shapefiles/AtlanticWolffish.shp")%>%
            st_transform(latlong)

wolffish2 <- read_sf("M:/OCMD/02_Projects/MPA/NetworkPlanning/ConservationPriorities/Shapefiles/atwolf7884f_fall_poly.shp")%>%
  st_transform(latlong)

wolffish3 <- read_sf("M:/OCMD/02_Projects/MPA/NetworkPlanning/ConservationPriorities/Shapefiles/atwolf7884s_spring_poly.shp")%>%
  st_transform(latlong)

cod <- read_sf("M:/OCMD/02_Projects/MPA/NetworkPlanning/ConservationPriorities/TechReportLayers/Ecological/FineFilter/DepletedSpecies/AtlanticCod4Vn.shp")%>%
  st_transform(latlong)

tt <- read_sf("M:/OCMD/02_Projects/MPA/NetworkPlanning/ConservationPriorities/TechReportLayers/Ecological/FineFilter/DepletedSpecies/AtlanticWolffish.shp")%>%st_transform(latlong)

ggplot()+
  geom_sf(data=sab_zones,fill=NA)+
  geom_sf(data=tt,fill="red")+
  geom_sf(data=wolffish2,fill="blue",alpha=0.5)+
  geom_sf(data=wolffish3,fill="green",alpha=0.5)+
  geom_sf(data=cod,fill="yellow")+
  geom_sf(data=infrac_loc)+
  theme_bw()+
  coord_sf(xlim=sab_bounds[c(1,3)],ylim=sab_bounds[c(2,4)])
