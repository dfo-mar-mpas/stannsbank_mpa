##Load libraries --- 

library(tidyverse)
library(sf)
library(terra)
library(ggspatial)
library(patchwork)
library(ggnewscale)

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load the St. Anns Bank MPA shapefiles
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

#load the reciever locations
sab_acoustic_new <- read.csv("data/Acoustic/SABMPAnew_June2025.csv")%>%
                    st_as_sf(coords=c("DEPLOY_LONG_DD","DEPLOY_LAT_DD"),crs=latlong)%>%
                    mutate(type="New array")

sab_acoustic <- read.csv("data/Acoustic/OTN_redesign_coords.csv")%>%
                st_as_sf(coords=c("long","lat"),crs=latlong)%>%
                mutate(type="Core array")

sab_acoustic_coords <- rbind(sab_acoustic%>%dplyr::select(type),
                             sab_acoustic_new%>%dplyr::select(type))

#load the banks polygon
sab_banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%
             st_transform(latlong)


#load bathymetric data
sab_dem <- rast("data/Bathymetry/sab_dem.tif")%>%project(latlong)

contour_sab <- as.contour(sab_dem, levels = 150)%>% # could also use 237 from 'deep' but just round out to 250
  st_as_sf()%>%
  st_transform(latlong)

#Basemap ----
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

#eDNA stations
edna_coords <- read.csv("output/2024_mission/station_coords.csv")%>%
  dplyr::select(Station,Depth,Longitude,Latitude)%>%
  filter(!grepl("Cam",Station))%>%
  rename(name=Station,depth=Depth,long=Longitude,lat=Latitude)%>%
  st_as_sf(coords=c("long","lat"),crs=latlong)%>%
  mutate(sample="eDNA")

#Map boundaries --- 
curdo_plot_bound <- sab_acoustic_new%>%
  filter(grepl("CB",STATION_NO))%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_transform(utm)%>%
  st_buffer(0.5)%>%
  st_transform(latlong)%>%
  st_bbox()

curdo_box <- curdo_plot_bound%>%st_as_sfc()

p1_curdo <- ggplot()+
  geom_sf(data=sab_banks%>%filter(name=="Curdo Bank"),fill="#F8766D",col="black")+
  geom_sf(data=edna_coords,shape=3,size=0.75,stroke=0.2)+
  geom_sf(data=sab_acoustic_coords,aes(fill=type),shape=21,size=2.5)+
  theme_bw()+
  theme(axis.text=element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size=10,vjust=-2),
        legend.position = "none",
        panel.grid = element_blank(),
        title = element_blank())+
  coord_sf(xlim=curdo_plot_bound[c(1,3)],ylim=curdo_plot_bound[c(2,4)],expand=0)+
  #labs(title="Curdo Bank",fill="")+
  labs(fill="")+
  annotation_scale(location="br")+
  scale_fill_manual(values=c("New array" = "royalblue4",
                             "Core array" = "darkorange"))

#Sactarie Bank

scatarie_plot_bound <- sab_acoustic_new%>%
  filter(!STATION_NO %in% c("SAB_SB2516","SAB_CB2501","SAB_CB2502","SAB_CB2503"))%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_transform(utm)%>%
  st_buffer(4)%>%
  st_transform(latlong)%>%
  st_bbox()

scatarie_box <- scatarie_plot_bound%>%st_as_sfc()

p1_scatarie <- ggplot()+
  geom_sf(data=sab_banks%>%filter(name=="Scatarie Bank"),fill="#00BFC4",col="black")+
  geom_sf(data=edna_coords,shape=3,size=0.75,stroke=0.2)+
  geom_sf(data=sab_acoustic_coords,aes(fill=type),shape=21,size=2.5)+
  theme_bw()+
  theme(axis.text=element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size=10,vjust=-2),
        legend.position = "none",
        panel.grid = element_blank(),
        title = element_blank())+
  coord_sf(xlim=scatarie_plot_bound [c(1,3)],ylim=scatarie_plot_bound [c(2,4)],expand=0)+
  #labs(title="Scatarie Bank",fill="")+
  labs(fill="")+
  annotation_scale(location="br")+
  scale_fill_manual(values=c("New array" = "royalblue4",
                             "Core array" = "darkorange"))


## whole array

sab_lims <- sab_acoustic_coords%>%
            st_bbox()%>%
            st_as_sfc()%>%
            st_transform(utm)%>%
            st_buffer(10)%>%
            st_transform(latlong)%>%
            st_bbox()

p1 <- ggplot()+
  geom_sf(data=basemap)+
  #geom_sf(data=contour_sab,linewidth=0.25,linetype=1,col="grey")+
  geom_sf(data=curdo_box,fill=NA)+
  geom_sf(data=scatarie_box,fill=NA)+
  geom_sf(data=sab_zones,fill=NA)+
  geom_sf(data=sab_banks,fill=c("#F8766D","#00BFC4"))+
  geom_sf(data=edna_coords,shape=3,size=1.1,stroke=0.2)+
  new_scale_fill()+
  geom_sf(data=sab_acoustic_coords,aes(fill=type),shape=21,size=1.1,stroke=0.2)+
  theme_bw()+
  #coord_sf(xlim=sab_lims[c(1,3)],ylim=sab_lims[c(2,4)],expand=0)+
  coord_sf(xlim=c(-59.8,-58.3),ylim=c(45.75,46.5),expand=0)+
  theme(plot.margin = margin(0, 0, 0, 0),
        legend.position="inside",
        legend.position.inside = c(0.1,0.9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        axis.text = element_blank())+
  annotation_scale(location="br")+
  scale_fill_manual(values=c("New array" = "royalblue4",
                             "Core array" = "darkorange"))+
  guides(fill = guide_legend(override.aes = list(size = 2.75)))

mosaic_map <- p1 + (p1_scatarie/p1_curdo) + plot_layout(ncol=2,widths=c(3.5,1.4))

ggsave("output/acoustic/sab_acoustic_2025.png",mosaic_map,width=8,height=4,units="in",dpi=600)
