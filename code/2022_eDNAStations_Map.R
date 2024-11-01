coords<-read.csv("data/eDNA/PER-2022-473/Perley2022_eDNA_Coords.csv", header = T)
head(coords)

#coords<-coords %>% mutate(Latitude=paste(Lat,X,sep=" "), Longitude=paste(Long,X.1, sep=" "))

#make coords into something sf and ggplot2 can use
coords_
#angle2dec <- function(angle) {
  angle <- as.character(angle)
  x <- do.call(rbind, strsplit(angle, split=' '))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    y[1] + y[2]/60 + y[3]/3600
  })
  return(x)
}
#Projections ------------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm_mar <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"



#load MPA polygons
sabshape <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)%>%
  mutate(name="St Anns Bank")

#load high resolution bathymetry (50m) from multibeam 
sab_bathy <- raster::raster("data/Bathymetry/bathy50/w001001.adf", RAT=FALSE)%>%
  raster::projectRaster(.,crs=latlong)

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
sab_benthoscape$class<-gsub("cobblles","cobbles",sab_benthoscape$class)
#load bathymetry 
ras <- raster::raster("data/Bathymetry/sab_dem.tif") %>%
  raster::projectRaster(.,crs=latlong)
rasproj <- sp::proj4string(ras)
#basemap of Nova Scotia
novsco <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Coastline/NS_coastline_project_Erase1.shp")%>%st_transform(latlong)%>%
  mutate(name="Nova Scotia")%>%
  dplyr::select(name,geometry)

sab_depth_plot <- ggplot()+
  geom_stars(data=stars_dem)+
  geom_sf(data=novsco)+
  geom_sf(data=sab,fill=NA,linewidth=1.25,colour="blue")+
  #coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.6,46.5), expand=F)+
  geom_point(data=coords, aes(x=DD_Long, y=DD_Lat))+
  theme_bw()+
  scale_fill_viridis(option="A",direction = -1,na.value="transparent")+
  labs(x="",y="",fill="Depth (m)")+
  theme(#axis.text = element_blank(),
    legend.position = "bottom", text=element_text(size=15), legend.key.width = unit(1.5, "cm"))
sab_depth_plot

ggsave("2022_eDNAStations.png",plot=sab_depth_plot, device = "png", path = "output/", width = 8, height=6, units = "in", dpi=320)
