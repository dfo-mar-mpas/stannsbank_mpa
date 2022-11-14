#load libraries 
library(sf)
library(ggplot2)
library(dplyr)
library(raster)
library(rnaturalearth)
library(rnaturalearthhires)
library(stars)
library(tidyr)
library(ggspatial)

sf_use_s2=FALSE

#coefficient for converting from m to nm 
nm_coeff <- 0.539956803/1000

#Projections ------------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm_mar <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load the St. Anns Bank MPA Polygon
sab <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab_nozones <- sab%>%
  st_transform(utm_mar)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

bioregion <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/MaritimesPlanningArea.shp")%>%st_transform(latlong)

#Load bathymetry
    # load("r:/Science/CESD/HES_MPAGroup/Data/Bathymetry/PredSpaceData.RData")
    # bathydat <- predSpace%>%dplyr::select("plon","plat","z")
    # colnames(bathydat) <- c("long","lat","depth")
    # coordinates(bathydat) <- ~long+lat
    # proj4string(bathydat) <- utm_mar
    # gridded(bathydat) <- TRUE
    # dem <- crop(raster(bathydat),bounding_area%>%st_transform(utm_mar)%>%as_Spatial())
    # writeRaster(dem,"r:/Science/CESD/HES_MPAGroup/Projects/eDNA Genome Quebec/data/sab_dem.tif")
dem <- raster("r:/Science/CESD/HES_MPAGroup/Projects/eDNA Genome Quebec/data/sab_dem.tif")

#Load the NSCC Benthoscape layers
sab_bs <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/StAnns_Benthoscape/benthoscape/benthoscape.shp")%>%
          st_transform(latlong)%>%
  group_by(Assigned_c)%>%
  summarise(st_combine(geometry))

#create area to search for stations
bounding_area <- sab%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_sf()%>%
  st_buffer(1)%>%# 1 degree buffer
  st_bbox()%>%
  st_as_sfc()%>%
  st_sf()%>%
  suppressWarnings()%>%
  suppressMessages()

#create a buffer for locations that the Perley could reach within one 12 hour day based on a steam time of 12 hours. 55 nms will be the furthest one can be away
louisbourg <- data.frame(lon=c(-59.971415,-59.970650,-59.955728),lat=c(45.916957,45.904229,45.901368))%>%
              st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)

#distance in nautical miles from the dock to the area where a bearing could be made to a station. 
louisbourg_approach <- louisbourg%>%summarize(do_union=TRUE) %>% st_cast("LINESTRING")%>%st_length()%>%as.numeric(.)*nm_coeff

#Make a bounding distance (km) from Louisbourg for stations that could be too far (55 nm - approach)

maxdist <- (55 - louisbourg_approach)/nm_coeff/1000

louisbourg_bounding <- louisbourg[3,]%>%
                       st_transform(utm_mar)%>%
                       st_buffer(maxdist)%>%
                       st_transform(latlong)

#load the 

#make a new polygon to select stations from
sab_bounded <- st_intersection(sab,louisbourg_bounding)

sab_bs_bounded <- st_intersection(sab_bs,louisbourg_bounding)%>%
                  st_intersection(.,sab)%>%
                  suppressMessages()%>%
                  suppressWarnings()

#select a series of 24 stations at random from this area. These can serve as placeholders. 
#Note I can't get group_by to work with st_sample thus the for loop

#Function to shrink polygons so that the values aren't selected at the edges, which for some reason they seem to be by default

  # https://gis.stackexchange.com/questions/392505/can-i-use-r-to-do-a-buffer-inside-polygons-shrink-polygons-negative-buffer
    # NB size is given as a positive value
    shrinkIfPossible <- function(sf, size) {
      # compute inward buffer
      sg <- st_buffer(st_geometry(sf), -size)
      
      # update geometry only if polygon is not degenerate
      st_geometry(sf)[!st_is_empty(sg)] = sg[!st_is_empty(sg)]
      
      # return updated dataset
      return(sf)
    }


station_design <- NULL
for(i in sab_bs$Assigned_c){
  
  temp <- sab_bs_bounded%>%
          filter(Assigned_c == i)%>%
          st_transform(utm_mar)%>%
          shrinkIfPossible(.,0.3)%>% #this just is to try to get the points away from the boundaries
          st_transform(latlong)%>%
          st_sample(size=4)%>%
          suppressMessages()%>%
          st_as_sf()%>%
          mutate(Assigned_c = i,
                 lb_dist = as.numeric(st_distance(.,louisbourg[3,])*nm_coeff),
                 lon=st_coordinates(.)[,1],
                 lat=st_coordinates(.)[,2],
                 stat_num = 1:4)%>%
          rename(geometry=x)%>%
          st_transform(proj4string(dem))%>%
          mutate(depth=raster::extract(dem,as_Spatial(.)))%>%
          st_transform(latlong)%>%
          dplyr::select(lon,lat,Assigned_c,lb_dist,depth,stat_num,geometry)
  
  station_design = rbind(station_design,temp)
  
  
}

station_design$station <- paste0("SAB_eDNA_",1:nrow(station_design))

p1 <- ggplot()+
  geom_sf(data=sab_bs%>%rename(benthoscape = Assigned_c),aes(fill=benthoscape))+
  geom_sf(data=sab,fill=NA)+
  #geom_sf(data=sab_bs_bounded,aes(fill=Assigned_c))+
  geom_sf(data=station_design)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid=element_blank())+
  annotation_north_arrow(location="tl",height = unit(0.3,"in"),width=unit(0.3,"in"))+
  annotation_scale(location="br");p1

ggsave("r:/Science/CESD/HES_MPAGroup/Projects/St Anns Bank/figures/2022_design.png",p1,width=6,height=6)

#save the output

output <- station_design%>%
          mutate(station_num = 1:n(),
                 Station_ID = paste0("SAB_eDNA_",station_num),
                 Longitude=round(st_coordinates(.)[,1],2),
                 Latitude=round(st_coordinates(.)[,2],2))%>%
          st_transform(proj4string(dem))%>%
          mutate(Nominal_Depth=raster::extract(dem,as_Spatial(.))%>%round(2))%>%
          data.frame()%>%
          dplyr::select(Station_ID,Longitude,Latitude,Nominal_Depth)

write.csv(output,
          "r:/Science/CESD/HES_MPAGroup/Projects/St Anns Bank/data/StAnns_eDNA_2022.csv",row.names=F)



#load proposed stations for St. Anns Bank. THese came from the 2020 ship time application 
sab_edna <- read.csv("r:/Science/CESD/HES_MPAGroup/Projects/St Anns Bank/data/eDNA_stations_SAB.csv")%>%
            dplyr::select("Station.Name","Station.Type","Latitude","Longitude")%>% #there is some wierd excel column things going on that I am too lazy to fix. 
            st_as_sf(coords=c("Longitude","Latitude"),crs=latlong,remove=FALSE)
            
    
#load basemap 
NScoast <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/NS_coastline_project_Erase1.shp")

basemap <- NScoast%>%
  st_intersection(bounding_area%>%st_transform(st_crs(NScoast)))%>%
  st_transform(latlong)%>%
  suppressWarnings()

#Plot stations
ggplot()+
  geom_sf(data=louisbourg_bounding,fill=NA)+
  geom_sf(data=basemap)+
  geom_sf(data=sab,fill=NA)+
  geom_sf(data=sab_edna)
  theme_bw()
  
  
### spatially balanced sampling
  #remotes::install_github("paul-vdb/DFO-master-sample")
  library(BASMasterSample)
  
  NS_master <- buildMS(shp = bioregion,rotate=TRUE)
  NS_master_rot <- rotate.shp(NS_master,NS_master,back=TRUE)
  
  small_bound <- sab%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()%>%
                buildMS(.,rotate=TRUE)%>%
                rotate.shp(.,.,back=TRUE)
                
    
   
  
  ggplot()+
    geom_sf(data=small_bound%>%st_transform(latlong))+
    geom_sf(data=sab)
    
  
  test=masterSample(sab_bs,n=4)
  test=masterSample(sab%>%st_transform(st_crs(small_bound)),n=20)
