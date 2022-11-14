#plot snow crab stations 

#libraries ---
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthhires)
library(marmap)

#map projection
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load polygons
sab <- rgdal::readOGR("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/SAB_boundary_zones_2017.shp",verbose = FALSE)%>%
  maptools::unionSpatialPolygons(.,IDs=rep(1,nrow(sp::coordinates(.))))%>%
  st_as_sf()%>%
  st_transform(latlong)%>%
  mutate(name="St Anns Bank")

gully <- st_read("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Gully_Boundary_Redone2014.shp")%>%
  st_transform(latlong)%>%
  mutate(name="Gully")

NovaScotia <- st_read("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/NS_coastline_project_Erase1.shp")%>%
  st_transform(latlong) #better resolution

MPAs <- rbind(sab%>%select(name,geometry),
              gully%>%select(name,geometry))

#load diet data stations
sdat <- read.csv("R:/Science/CESD/HES_MPAGroup/Projects/St Anns Bank/Monitoring/Diets/SnCb_Stomach_Data.csv")%>%
  select(MISSION,SETNO,FISHSET_ID,YEAR,LAT2,LONG2)%>%
  mutate(latvar = LAT2)%>%
  rename(Longitude = LONG2, Latitude = LAT2)%>%
  filter(Latitude<46.5)%>%
  mutate(settrip = paste(Longitude,Latitude,YEAR))%>%
  distinct(settrip,.keep_all = T)%>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=latlong)%>%
  st_join(.,MPAs,join=st_intersects)%>%
  mutate(name=ifelse(is.na(name),ifelse(latvar<44.5,"Outside Gully","Outside SAB"),name))

#summary by year and inside outside
table(sdat$name,sdat$YEAR)

#Map extent (Nova Scotia)
rextent <- sdat%>%
  raster::extent(.)

#set limits
Long.lim  <-  c(rextent[1]-2, rextent[2]+1)
Lat.lim <-  c(rextent[3]-0.5, rextent[4]+0.5)

#download bathy
  curdir <- getwd()
  setwd("data/")
  bathy <- getNOAA.bathy(Long.lim[1],Long.lim[2],Lat.lim[1],Lat.lim[2],res=1,keep=T)
  setwd(curdir)

p1 <- ggplot()+
  geom_sf(data=NovaScotia,fill="white")+
  geom_sf(data=sab,alpha=0.5,fill="cornflowerblue")+
  geom_sf(data=gully,alpha=0.5,fill="cornflowerblue")+
  geom_sf(data=sdat)+
  geom_contour(data=bathy,aes(x=x,y=y,z=z),breaks=c(-200),lwd=0.05,colour="grey50")+
  coord_sf(xlim = Long.lim,  ylim = Lat.lim, expand=FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x=expression(paste("Longitude ",degree,"W",sep="")),
       y=expression(paste("Latitude ",degree,"N",sep="")))

ggsave("CrabSurveyStations.png",p1,dpi=600,height=6,width=6,units="in")
