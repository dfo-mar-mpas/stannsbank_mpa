#Get the fish species code from ANDES - these get updated from time to time so may need to replace this csv eventually
fishcodes <- read.csv("data/CrabSurvey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv")


catchdat <- read.csv("data/CrabSurvey/SABMPA2023export.csv", header = T) 

arsw <- read.csv("C:/Users/JEFFERYN/Documents/GitHub/stannsbank_mpa/data/CrabSurvey/sabmpa_area_swept.csv", header = T)
arsw <- arsw %>% unite(TRIP_STATION, c("TRIP", "STATION"), remove = F)
arsw %>% summarise(mean=mean(AREA_SWEPT,na.rm = T), 
                   median=median(AREA_SWEPT, na.rm=T), 
                   stdev=sd(AREA_SWEPT, na.rm=T),
                   se=(sd(AREA_SWEPT)/sqrt(length(AREA_SWEPT))))
#mean        median        stdev          se 
#0.00435671 0.004279407 0.0008242672    6.730114-05

#merge catch data with species codes
catchdat <- merge(catchdat, fishcodes, by.x="SPECCD_ID", by.y="SPECCD_ID", all.x=T) #Amy Glass did a more recent pull in December 2023 so comparing this csv file to the MPATOWS RData object
catchdata <- left_join(arsw, catchdat, by = "TRIP_STATION", relationship = "many-to-many", keep = F)
#remove redundant columns
catchdata <- catchdata[,-c(8, 11)]

catch<-catchdata %>% mutate(Year=format(as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y"))


sf_use_s2 = FALSE
latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"


#basemap of Nova Scotia
novsco <- read_sf("data/Shapefiles/NS_coastline_project_Erase1.shp")%>%st_transform(latlong)%>%
  mutate(name="Nova Scotia")%>%
  dplyr::select(name,geometry)

ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1.25)+
  #geom_sf(data=gully2, colour="red", fill=NA)+
  coord_sf(xlim=c(-61, -58), ylim=c(45.25,47.1), expand=F)+
  geom_point(data=catch, aes(x=LONGITUDE*-1, y=LATITUDE),size=3)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="lightblue1"), text=element_text(size=20))

#get the most common species



#Plot Wolffish
m1 <-ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25, alpha=0.75)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=0.75)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="STRIPED ATLANTIC WOLFFISH"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m1

ggsave(filename = "output/CrabSurvey/Wolffish_Catch_byYear.png", plot = m1, device = "png", width = 10, height =6, units = "in", dpi=500)


#Plot Snow Crab
m2 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="SNOW CRAB  QUEEN; INCL CHIONOECETES SP"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m2

ggsave(filename = "output/CrabSurvey/Snowcrab_Catch_byYear.png", plot = m2, device = "png", width = 10, height =6, units = "in", dpi=500)

#Plot Silver Hake
m3 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="SILVER HAKE"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m3

ggsave(filename = "output/CrabSurvey/SilverHake_Catch_byYear.png", plot = m3, device = "png", width = 10, height =6, units = "in", dpi=500)

#Plot cod
m4 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="COD ATLANTIC"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m4

ggsave(filename = "output/CrabSurvey/Cod_Catch_byYear.png", plot = m4, device = "png", width = 10, height =6, units = "in", dpi=500)

#American Plaice
m5 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="AMERICAN PLAICE"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m5

ggsave(filename = "output/CrabSurvey/AmericanPlaice_Catch_byYear.png", plot = m5, device = "png", width = 10, height =6, units = "in", dpi=500)

#Witch Flounder
m6 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="WITCH FLOUNDER"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m6

ggsave(filename = "output/CrabSurvey/WitchFlounder_Catch_byYear.png", plot = m6, device = "png", width = 10, height =6, units = "in", dpi=500)

#Redfish
m7 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="REDFISH UNSEPARATED"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m7

ggsave(filename = "output/CrabSurvey/Redfish_Catch_byYear.png", plot = m7, device = "png", width = 10, height =6, units = "in", dpi=500)

#Sea Pens
m8 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="SEA PEN; PENNATULACEA"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m8

ggsave(filename = "output/CrabSurvey/SeaPens_Catch_byYear.png", plot = m8, device = "png", width = 10, height =6, units = "in", dpi=500)

#Sponges
m9 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="SPONGES"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m9

ggsave(filename = "output/CrabSurvey/Sponges_Catch_byYear.png", plot = m9, device = "png", width = 10, height =6, units = "in", dpi=500)

#Thorny Skate
m10 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="THORNY SKATE"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m10

ggsave(filename = "output/CrabSurvey/ThornySkate_Catch_byYear.png", plot = m10, device = "png", width = 10, height =6, units = "in", dpi=500)


#White Hake
m11 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(COMM=="WHITE HAKE"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m11

ggsave(filename = "output/CrabSurvey/WhiteHake_Catch_byYear.png", plot = m11, device = "png", width = 10, height =6, units = "in", dpi=500)

#LEPTASTERIAS POLARIS
m12 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(SPEC=="LEPTASTERIAS (HEXASTERIAS) POLARIS"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m12

ggsave(filename = "output/CrabSurvey/LepPolaris_Catch_byYear.png", plot = m12, device = "png", width = 10, height =6, units = "in", dpi=500)

#Pennatula
m13 <- ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.5,46.5), expand=F)+
  geom_point(data=catch %>% filter(SPEC=="PENNATULOIDEA"), aes(x=LONGITUDE*-1, y=LATITUDE, size=EST_NUM_CAUGHT))+
  facet_wrap(vars(Year),nrow=2)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(legend.position="bottom",
        panel.background = element_rect(fill="white"), text=element_text(size=12))+
  guides(size=guide_legend(title = "Number per trawl"), fill=guide_legend(title = ""))
m13

ggsave(filename = "output/CrabSurvey/Pennatuloidea_Catch_byYear.png", plot = m13, device = "png", width = 10, height =6, units = "in", dpi=500)
