### Activity report for PER-2023-767 mission in St. Anns Bank

#load libraries ----
library(dplyr)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(ggnewscale)
library(patchwork)

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#Sample stations ----
stations_2023 <- read.csv("data/eDNA/PER-2023-767/SampleDataSheet.csv")%>%#Stations sampled in 2023 - this is the table from form-c
                 filter(Station != "BLANK",
                        Location == "St Anns Bank")%>%
                 mutate(lat=as.numeric(substring(Lat,1,2)) + as.numeric(substring(Lat,4))/60,
                        lon=(as.numeric(substring(Long,1,2)) + as.numeric(substring(Long,4))/60)*-1,
                        type=ifelse(type !="eDNA","Paired",type),
                        type=factor(type,levels=c("Paired","eDNA")))%>%
                  st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                 distinct(Long,.keep_all=TRUE)

#Benthoscape classification
benthoscape_classes <- data.frame(Assigned_c=c("A - Mud",
                                               "Asp - Mud and seapens",
                                               "B - Gravelly sand/mud <50% cobblles/gravel",
                                               "C - Till >50% cobbles/gravel",
                                               "D - Till with coraline algae",
                                               "E - Gravel with crinoids",
                                               "F - Sand with Sand dollars"),
                                  classn=c("Mud", #simplified names
                                           "Mud/Seapens",
                                           "Gravelly sand/mud",
                                           "Till >50% cob/grav",
                                           "Till/coraline algae",
                                           "Gravel/crinoids",
                                           "Sand/Sand dollars"))

#load benthoscape classification
sab_benthoscape <- read_sf("data/Shapefiles/benthoscape.shp")%>%
  st_transform(latlong)%>%
  st_make_valid()%>%
  left_join(.,benthoscape_classes)


#load St Anns Bank MPA shapefile 
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

#sab map coast
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

plot_boundaries <- c(c(-60,-58.3,45.75,46.5))

plot_lim2 <- stations_2023%>%
             st_bbox()%>%
             st_as_sfc()%>%
             st_transform(utm)%>%
             st_buffer(2)%>%
             st_transform(latlong)%>%
             st_bbox()


p1 <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=sab_zones,fill=NA)+
  geom_sf(data=stations_2023,aes(fill=type),pch=21)+
  geom_sf(data=plot_lim2%>%st_as_sfc(),fill=NA,lty=2)+
  coord_sf(expand=0,xlim=plot_boundaries[c(1,2)],ylim=plot_boundaries[c(3,4)])+
  scale_fill_manual(values=c("black","white"))+
  theme_bw()+
  labs(fill="")+
  theme(legend.position="bottom")

p2 <- ggplot()+
  geom_sf(data=sab_benthoscape,aes(fill=classn))+
  geom_sf(data=sab_zones,fill=NA)+
  new_scale_fill()+
  geom_sf(data=stations_2023,aes(fill=type),pch=21,size=2.5,show.legend = FALSE)+
  scale_fill_manual(values=c("black","white"))+
  coord_sf(expand=0,xlim=plot_lim2[c(1,3)],ylim=plot_lim2[c(2,4)])+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position="bottom")+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))

p3 <- p1+p2+plot_layout(ncol=2)

ggsave("output/2023_mission/activity_report_plot.png",height=5,width=10,p3,dpi=300)

## formated table for report

out <- stations_2023%>%
       data.frame()%>%
       dplyr::select(Station,Date,type,lat,lon)

write.csv(out,"output/2023_mission/survey_stations_activity_report.csv",row.names=F)

#full list

read.csv("data/eDNA/PER-2023-767/SampleDataSheet.csv")%>%#Stations sampled in 2023 - this is the table from form-c
  filter(Station != "BLANK")%>%
  mutate(lat=as.numeric(substring(Lat,1,2)) + as.numeric(substring(Lat,4))/60,
         lon=(as.numeric(substring(Long,1,2)) + as.numeric(substring(Long,4))/60)*-1,
         type=ifelse(type !="eDNA","Paired",type),
         type=factor(type,levels=c("Paired","eDNA")))%>%
  st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
  distinct(Long,.keep_all=TRUE)%>%
  data.frame()%>%
  dplyr::select(Station,Date,type,lat,lon)%>%
  write.csv(.,"output/2023_mission/survey_stations_formc.csv",row.names=F)
