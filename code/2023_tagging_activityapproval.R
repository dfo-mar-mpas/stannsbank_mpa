#load libraries ----
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)
library(ggplot2)
library(raster)
library(stars)
library(geosphere)
library(ggspatial)
library(ggnewscale)
library(viridis)
library(patchwork)
library(tmaptools)
library(patchwork)

s2_as_sf = FALSE

source("r:/Science/CESD/HES_MPAGroup/R/Functions/box_grid.R")

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load bathymetry ------
dem_sab <- raster("data/Bathymetry/sab_dem.tif") # SAB based on the 'predSpace' data
load("data/Bathymetry/250m_stars_object.RData") # STARS object based on the 250m depth contour extracted (converted) from GEBCO

#load St Anns Bank MPA shapefile 
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

#load coastline and make basemap ------
coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")

bounding_area <- sab%>%
  st_transform(utm)%>%
  st_buffer(dist = 100)%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_as_sf()%>%
  st_transform(latlong)%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_as_sf()%>%
  suppressWarnings()%>%
  suppressMessages()

basemap <- coast_hr%>%
  st_intersection(bounding_area%>%st_transform(st_crs(coast_hr)))%>%
  st_transform(latlong)%>%
  suppressWarnings()

plot_boundaries <- c(c(-61,-58,45.6,46.6))

#tagging locations -----

#2021 activity approval start positions
halibut_2021 <- read.csv("data/Acoustic/2021_HalibutStations.csv")%>%
                st_as_sf(coords=c("Longitude","Latitude"),crs=latlong)%>%
                mutate(type="Halibut Survey")%>%
                dplyr::select(type,geometry)

halibut_2022 <- read.csv("data/Acoustic/Halibut tagging sheet Nov 2022 OTN.csv")%>%
                filter(!is.na(TAG_SERIAL_NUMBER))%>% #tailing blank spaces in the Excel sheet - this will be inclusive of cod and halibut where sets were deployed
                mutate(lon.deg=as.numeric(substr(RELEASE_LONGITUDE, 1, 2)),
                       lon.min=RELEASE_LONGITUDE - (lon.deg*100),
                       lon = (lon.deg + lon.min/60)*-1,
                       lat.deg=as.numeric(substr(RELEASE_LATITUDE, 1, 2)),
                       lat.min=RELEASE_LATITUDE - (lat.deg*100),
                       lat = lat.deg + lat.min/60,
                       type="2022 tag locations")%>%
                st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                dplyr::select(type,geometry)

surveysets <- rbind(halibut_2021,halibut_2022)%>%
              st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
              mutate(lon=round(st_coordinates(.)[,1],3),
                     lat=round(st_coordinates(.)[,2],3))%>%
              arrange(type,Zone)

#save coordinates for the activity approval applcation - to simplify the included table, the 2022 tag locations will be presented as centriods of the three clusters
mod_sets <- surveysets%>%
            filter(type=="2022 tag locations")%>%
            mutate(cluster = case_when(lon> -59.26 ~ 1, #rougly 4 clusters
                                       lon> -59.3 & lon < (-59.26) ~ 2,
                                       lon< (-59.5) ~ 3,
                                       TRUE ~ 4))%>%
            group_by(cluster)%>%
            summarise(geometry = st_union(geometry))%>%
            st_centroid()%>%
            ungroup()%>%
            st_as_sf()%>%
            mutate(lon=round(st_coordinates(.)[,1],3),
                   lat=round(st_coordinates(.)[,2],3),
                   type="2022 tag locations")%>%
            st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
            suppressWarnings()

surveysets_df <- rbind(surveysets%>%
                         filter(type=="Halibut Survey")%>%
                         data.frame()%>%
                         dplyr::select(type,Zone,lon,lat),
                       mod_sets%>%
                         data.frame()%>%
                         dplyr::select(type,Zone,lon,lat))%>%
                  arrange(type,Zone)

#save coordinates
write.csv(surveysets_df,file="output/Acoustic/2023_ActivityApproval_taglocations.csv",row.names = FALSE)

#make plot
p1 <- ggplot()+
  geom_sf(data=shelfbreak,fill=NA,lwd=0.2)+
  geom_sf(data=basemap)+
  geom_sf(data=sab_zones,fill="cornflowerblue",alpha=0.30,col="grey10")+
  geom_sf(data=surveysets,aes(fill=type),col="black",pch=21,size=2)+
  geom_sf(data=mod_sets,pch=21,fill="white",col="black")+
  scale_fill_manual(values=c("grey30","grey80"))+
  coord_sf(xlim=plot_boundaries[c(1,2)],ylim=plot_boundaries[c(3,4)],expand=0)+
  labs(fill="")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position="bottom")

ggsave("output/Acoustic/2023_ActivityApproval_taglocations.png",p1,height=5,width=8,units="in",dpi=300)
