## Code to make plots for the first St. Anns Bank MPA Working Group Workshop - November 16-17 2022 Membertou

#Load libraries ----
library(ggplot2)
library(raster)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(stars)
library(tidyr)
library(ggspatial)
library(wdpar)
library(viridis)
library(fasterize)
library(stars)
library(patchwork)
library(lme4)
library(foster)

sf_use_s2=FALSE

#source functions
source("code/MCT_table.R") #from CanadaMPAs analysis
source("code/raster_gen.R") #from Eelgrass Poolseq

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

mar_network <- read_sf("c:/Users/stanleyr/Documents/Github/CAN_MCN_Vulnerability/data/shapefiles/networksites_proposed_OEM_MPA_20220617.shp")%>%
  st_transform(latlong)%>%
  mutate(TYPE = ifelse(TYPE == "TBD","Draft",TYPE))

mar_plotlims <- st_bbox(bioregion)

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

#load bathymetry 
ras <- raster::raster("data/Bathymetry/sab_dem.tif") %>%
  raster::projectRaster(.,crs=latlong)
rasproj <- proj4string(ras)

#bathymetry comparison
x <- sab_nozones%>%st_transform(rasproj)
cr <- x%>%extent()
b1 <- crop(ras,cr,snap="out") 

ind <- b1%>%
  fasterize(x,.)%>% #this is reverse engineering 
  values(.)%>%
  as.data.frame()%>%
  rename(val=1)%>%
  mutate(val=!is.na(val))%>%
  pull(val)

sab_bathy_values <- data.frame(depth=values(b1)[ind]*-1,
                               site="St Anns Bank")
#create an output raster --
values(b1)[!ind] <- NA #masks out the non polygon raster values

##make it a plot
stars_dem <- ras%>%
             projectRaster(.,crs=latlong)%>%
             crop(.,extent(sab_nozones))%>%
             mask(.,sab_nozones%>%as_Spatial())%>%
             st_as_stars()

sab_depth_plot <- ggplot()+
  geom_stars(data=stars_dem)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=sab,fill=NA)+
  coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  theme_bw()+
  scale_fill_viridis(option="A",direction = -1,na.value="transparent")+
  labs(x="",y="",fill="Depth (m)")+
  theme(axis.text = element_blank(),
        legend.position = "bottom")

ggsave("output/sab_bathymetry_plot.png",sab_depth_plot,width=8,height=7,units="in",dpi=600)

sab_depth_ridge_plot <- ggplot(sab_bathy_values,aes(x=depth,y=1,fill=stat(x)))+
  geom_density_ridges_gradient(scale=3,size=0.3,rel_min_height=0.001)+
  scale_fill_viridis_c(option="C")+
  theme_bw()+
  scale_x_continuous()+
  coord_flip(expand = 0)+
  labs(x=expression(paste(log[10]," Depth (m)",sep="")),y="",fill="")+
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank())

ggsave("output/sab_bathymetry_ridges.png",sab_depth_ridge_plot,width=8,height=8,units="in",dpi=600)

##SAB vs MPAs analysis

global_pas <- wdpa_fetch("global")


#Canadian MPAs
MCT <- MCT_table(update=TRUE)%>%
  rename(waterbody = 1,
         type=4,
         agency=5,
         area=6,
         precent_coverage=7)%>%
  mutate(area=as.numeric(gsub(",","",area)),
         type_closure = gsub("Oceans Act MPA","MPA",type),
         type_closure = ifelse(type_closure %in% c("Marine Refuge","MPA"),type_closure,"Other"))

sum_mct <- MCT%>%group_by(type_closure)%>%summarise(area=sum(area))%>%ungroup()%>%data.frame()

#empirical frequency distribution
can_cas_ecdf <- ecdf(MCT$area)
sab_area <- MCT%>%filter(grepl("Anns",Name))%>%pull(area)
sab_ecdf <- can_cas_ecdf(sab_area)

sab_ecdf_plot <- ggplot(data=MCT,aes(area))+
  geom_point(aes(x=sab_area,y=sab_ecdf),col="red",size=4)+
  geom_segment(aes(x=0,xend=sab_area,y=sab_ecdf,yend=sab_ecdf),lty=2)+
  geom_segment(aes(x=sab_area,xend=sab_area,y=-Inf,yend=sab_ecdf),lty=2)+
  stat_ecdf(geom="point")+
  scale_x_log10()+
  scale_y_continuous(labels=percent)+
  theme_bw()+
  labs(y="Cumulative Frequency Distribution",x=expression ("Area "~km^2))+
  annotation_logticks(sides="b")

ggsave("output/sab_ecdf_plot.png",sab_ecdf_plot,height=6,width=6,units="in",dpi=600)

##Benthoscape plot

sab_benthoscape_plot <- ggplot()+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sab,fill=NA)+
  coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text = element_blank())+
  labs(fill="")+
  guides(fill = guide_legend(nrow = 3))

ggsave("output/sab_benthoscape_plot.png",width=8,height=7,units="in",dpi=300)

##Climate Change plots - this uses eelgrass Poolseq data 

#2075
# bnam_raster_gen(x="C:/Users/stanleyr/Desktop/Eelgrass_Poolseq/Data/bnam/2075/RCP 8.5/dTbtm_dSbtm_F_R85_2066-2085.mat",
#                 region=bioregion_box,
#                 crs_val = latlong)
# 
# bnam_raster_gen(x="C:/Users/stanleyr/Desktop/Eelgrass_Poolseq/Data/bnam/2075/RCP 4.5/dTbtm_dSbtm_F_R45_2066-2085.mat",
#                 region=bioregion_box,
#                 crs_val = latlong)
# 
# #2055
# bnam_raster_gen(x="C:/Users/stanleyr/Desktop/Eelgrass_Poolseq/Data/bnam/2055/RCP 8.5/dTbtm_dSbtm_F_R85_2046-2065.mat",
#                 region=bioregion_box,
#                 crs_val = latlong)
# 
# bnam_raster_gen(x="C:/Users/stanleyr/Desktop/Eelgrass_Poolseq/Data/bnam/2055/RCP 4.5/dTbtm_dSbtm_F_R45_2046-2065.mat",
#                 region=bioregion_box,
#                 crs_val = latlong)

rcp45_55 <- raster("output/bnam_rasters/RCP45_2055_dTbtm_ann.tif")%>%
  st_as_stars()%>%
  st_transform(latlong)

rcp85_55 <- raster("output/bnam_rasters/RCP85_2055_dTbtm_ann.tif")%>%
  st_as_stars()%>%
  st_transform(latlong)

rcp45_75 <- raster("output/bnam_rasters/RCP45_2075_dTbtm_ann.tif")%>%
  st_as_stars()%>%
  st_transform(latlong)

rcp85_75 <- raster("output/bnam_rasters/RCP85_2075_dTbtm_ann.tif")%>%
  st_as_stars()%>%
  st_transform(latlong)

#range of temperatures for plotting
plotrange2 <- c(-0.1,2.7)

rcp_45_2055 <- ggplot()+
  geom_stars(data=rcp45_55)+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),lwd=0.25)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="grey70",lwd=0.4)+
  geom_sf(data=basemap_atlantic,lwd=0.01)+
  geom_sf(data=bioregion,lwd=0.02,fill=NA)+
  coord_sf(expand=0,xlim=mar_plotlims[c(1,3)],ylim=mar_plotlims[c(2,4)])+
  theme_bw()+
  scale_fill_viridis(option = "B",limits=plotrange2)+
  labs(title="RCP 4.5 2055",
       fill=expression("Bottom temperature anomaly " ( degree*C)))

rcp_85_2055 <- ggplot()+
  geom_stars(data=rcp85_55)+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),lwd=0.25)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="grey70",lwd=0.4)+
  geom_sf(data=basemap_atlantic,lwd=0.01)+
  geom_sf(data=bioregion,lwd=0.02,fill=NA)+
  coord_sf(expand=0,xlim=mar_plotlims[c(1,3)],ylim=mar_plotlims[c(2,4)])+
  theme_bw()+
  scale_fill_viridis(option = "B",limits=plotrange2)+
  labs(title="RCP 8.5 2055",
       fill=expression("Bottom temperature anomaly " ( degree*C)))

rcp_45_2075 <- ggplot()+
  geom_stars(data=rcp45_75)+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),lwd=0.25)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="grey70",lwd=0.4)+
  geom_sf(data=basemap_atlantic,lwd=0.01)+
  geom_sf(data=bioregion,lwd=0.02,fill=NA)+
  coord_sf(expand=0,xlim=mar_plotlims[c(1,3)],ylim=mar_plotlims[c(2,4)])+
  theme_bw()+
  scale_fill_viridis(option = "B",limits=plotrange2)+
  labs(title="RCP 4.5 2075",
       fill=expression("Bottom temperature anomaly " ( degree*C)))

rcp_85_2075 <- ggplot()+
  geom_stars(data=rcp85_75)+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),lwd=0.25)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="grey70",lwd=0.4)+
  geom_sf(data=basemap_atlantic,lwd=0.01)+
  geom_sf(data=bioregion,lwd=0.02,fill=NA)+
  coord_sf(expand=0,xlim=mar_plotlims[c(1,3)],ylim=mar_plotlims[c(2,4)])+
  theme_bw()+
  scale_fill_viridis(option = "B",limits=plotrange2)+
  labs(title="RCP 8.5 2075",
       fill=expression("Bottom temperature anomaly " ( degree*C)))

#bring them together with patchwork
output_plot <- rcp_45_2055 + rcp_45_2075 + rcp_85_2055 + rcp_85_2075 + plot_layout(nrow=2,guides="collect") & theme(legend.position = "bottom",axis.text=element_blank())

ggsave("output/SAB_Climate_Scenarios.png",output_plot,width=8,height=8,units="in",dpi=300)

#zoomed in plot using RCP_85
sab_rcp_85_2075 <- ggplot()+
  geom_stars(data=rcp85_75)+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),lwd=0.25)+
  geom_sf(data=sab,fill=NA,col="grey70",lwd=0.4)+
  geom_sf(data=basemap_atlantic,lwd=0.01)+
  geom_sf(data=bioregion,lwd=0.02,fill=NA)+
  coord_sf(expand=0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  theme_bw()+
  scale_fill_viridis(option = "B")+
  labs(fill=expression("Bottom temperature anomaly " ( degree*C)))+
  theme(legend.position="bottom",
        axis.text=element_blank())

sab_rcp_85_2075_full <- ggplot()+
  geom_stars(data=rcp85_75)+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),lwd=0.25)+
  geom_sf(data=sab,fill=NA,col="grey70",lwd=0.4)+
  geom_sf(data=basemap_atlantic,lwd=0.01)+
  geom_sf(data=bioregion,lwd=0.02,fill=NA)+
  coord_sf(expand=0,xlim=mar_plotlims[c(1,3)],ylim=mar_plotlims[c(2,4)])+
  theme_bw()+
  scale_fill_viridis(option = "B",limits=plotrange2)+
  labs(fill=expression("Bottom temperature anomaly " ( degree*C)))+
  theme(legend.position="bottom",
        axis.text=element_blank())

ggsave("output/SAB_RCP85_2075.png",sab_rcp_85_2075,width=8,height=8,units="in",dpi=300)
ggsave("output/SAB_RCP85_2075_full.png",sab_rcp_85_2075_full,width=8,height=8,units="in",dpi=300)

sab_overlay <- rbind(data.frame(rcp=45,year=2055,lab="RCP 4.5",lab_full="RCP 4.5 2055",vals=rcp45_55$RCP45_2055_dTbtm_ann[st_intersects(rcp45_55,sab,sparse =FALSE)]),
                     data.frame(rcp=45,year=2075,lab="RCP 4.5",lab_full="RCP 4.5 2075",vals=rcp45_75$RCP45_2075_dTbtm_ann[st_intersects(rcp45_75,sab,sparse =FALSE)]),
                     data.frame(rcp=85,year=2055,lab="RCP 8.5",lab_full="RCP 8.5 2055",vals=rcp85_55$RCP85_2055_dTbtm_ann[st_intersects(rcp85_55,sab,sparse =FALSE)]),
                     data.frame(rcp=85,year=2075,lab="RCP 8.5",lab_full="RCP 8.5 2075",vals=rcp85_75$RCP85_2075_dTbtm_ann[st_intersects(rcp85_75,sab,sparse =FALSE)]))%>%
  mutate(lab_full = factor(lab_full,levels=c("RCP 4.5 2055","RCP 8.5 2055","RCP 4.5 2075","RCP 8.5 2075")))

sab_overlay%>%group_by(lab_full)%>%summarise(min=min(vals,na.rm=T),max=max(vals,na.rm=T),mean=mean(vals,na.rm=T))%>%ungroup()%>%data.frame()

#Look at the Lewis et al. results

sab_emerge <- read.csv("data/ENSEMBLE_species_lost.csv")%>% #from Lewis et al 
  filter(abbreviation == "SABMPA")

#Bioclassification overlay analysis 

climate.sensitive <- st_read('r:/Science/CESD/HES_MPAGroup/Projects/SPERA/Biological Classification/John_OB/bioclassification/Output/Mar_ClimateSensitive.shp')%>%
  st_transform(latlong)

mar.classes <- st_read("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Bioclassification/MaritimesBioclassificationPolygons.shp")%>%st_transform(latlong)

sab_focus <- sab%>%
  st_transform(utm_mar)%>%
  st_buffer(200)%>%
  st_transform(latlong)%>%
  st_bbox()

climate.sensitive%>%st_make_valid()%>%st_union()%>%st_intersection(.,sab_nozones)%>%st_area()/st_area(sab_nozones)*100

mc_1 <- ggplot()+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=mar.classes%>%filter(name %in% c("ESS","ESS: Banks","Laurentian Channel/Shelf Break")),aes(fill=name))+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),col="grey100",lwd=0.25)+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=climate.sensitive,fill=alpha("grey",0.08),col="black",lwd=1.15)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="black",lwd=0.5)+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text=element_blank())+
  coord_sf(expand=0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  labs(fill="")

mc_2 <- ggplot()+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=mar.classes%>%filter(name %in% c("ESS","ESS: Banks","Laurentian Channel/Shelf Break")),aes(fill=name))+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),col="grey100",lwd=0.25)+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="black",lwd=0.5)+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.text=element_blank())+
  coord_sf(expand=0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  labs(fill="")

mc_3 <- ggplot()+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=mar.classes,aes(fill=name))+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),col="grey100",lwd=0.25)+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="black",lwd=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")+
  coord_sf(expand=0,xlim=mar_plotlims[c(1,3)],ylim=mar_plotlims[c(2,4)])+
  labs(fill="")

#save the figures
ggsave("output/mar_class_vuln_sab.png",mc_1,height=6,width=6,units = "in",dpi=300)
ggsave("output/mar_class_sab.png",mc_2,height=6,width=6,units = "in",dpi=300)
ggsave("output/mar_class_network.png",mc_3,height=6,width=6,units = "in",dpi=300)

#Crab survey data 

#load crab survey data
load("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/2022-10-31_SAB_crabsurvey_mpastations.RData")

crabdf <- crabdat%>%
  st_as_sf(coords=c("LONGITUDE","LATITUDE"),crs=latlong,remove=FALSE)%>%
  st_intersection(bounding_area)%>%
  st_join(.,sab_nozones,left=TRUE)%>%
  mutate(location = ifelse(is.na(name),"Outside","Inside"))

crab_year <- data.frame()
for(i in unique(crabdf$YEAR)){
  
  inside <- crabdf%>%
    data.frame()%>%
    filter(YEAR == i,
           location=="Inside")%>%
    distinct(SET_NO)%>%
    nrow()
  
  outside <- crabdf%>%
    data.frame()%>%
    filter(YEAR == i,
           location=="Outside")%>%
    distinct(SET_NO)%>%
    nrow()
  
  crab_year <- rbind(crab_year,data.frame(year=i,inside=inside,outside=outside))
  
}

#year plot
year_df <- crab_year%>%gather(location,value,2:3)

year_plot <- ggplot(year_df,aes(x=year,y=value,fill=location))+
  geom_bar(position="stack",stat="identity",col="black")+
  theme_bw()+
  scale_x_continuous(breaks=2015:2021)+
  scale_y_continuous(breaks=seq(1,14,2),expand=c(0,0.1))+
  labs(fill="",y="Number of sets",x="")+
  theme(legend.position = "bottom")

ggsave("output/SAB_crabsurvey_by_year.png",year_plot,height = 8,width=6,units="in",dpi=300)

#richness
rich_df <- crabdf%>%
  data.frame()%>%
  group_by(YEAR,STATION)%>%
  summarise(nspecies = length(unique(SCIENTIFIC)))%>%
  ungroup()%>%
  left_join(.,crabdf%>%
              data.frame()%>%
              mutate(id=paste(YEAR,STATION,sep="_"))%>%
              distinct(id,.keep_all=T)%>%
              dplyr::select(STATION,DATE_TIME2,DEP2,LONGITUDE,LATITUDE,YEAR,location))%>%
  mutate(station=paste("ST",STATION,sep="_"),
         strata=ifelse(DEP2>150,">150m","<150m"),
         year_strat=ifelse(YEAR>2018,"2019-2021","2015-2018"))


sab_rich_depth <- ggplot(rich_df%>%filter(nspecies>9),aes(x=DEP2,y=nspecies,col=location,group=location))+
  geom_point(size=2)+
  stat_smooth(method="lm")+
  theme_bw()+
  labs(x="Depth (m)",y="Richness",col="")+
  theme(legend.position = "none")

ggsave("output/SAB_crabsurvey_richness_by_depth.png",sab_rich_depth,height = 8,width=6,units="in",dpi=300)

#model 

rich_count <- rich_df%>%
  group_by(STATION)%>%
  summarise(nyears=n(),
            minyear=min(YEAR),
            maxyear=max(YEAR))

rich_mod <- lmer(nspecies ~ location + (1 + location | station),data=rich_df) #doesn't work

ggplot(rich_df,aes(x=YEAR,y=nspecies,group=station))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+
  labs(x="Year",y="Richness",col="")+
  facet_grid(~location)

#very-very preliminary - no actual stats.   
richness_plot <- ggplot(rich_df,aes(x=YEAR,y=nspecies))+
  geom_line(aes(col=station,group=station))+
  geom_point(size=2,aes(col=station))+
  stat_smooth(method="lm",aes(group=location))+
  theme_bw()+
  labs(x="Year",y="Richness",col="")+
  facet_grid(~location)+
  theme(legend.position="none")

ggsave("output/SAB_crabsurvey_richness_by_year.png",richness_plot,height = 8,width=6,units="in",dpi=300)


sab_rich_depth_strata <- ggplot(rich_df,aes(x=location,y=nspecies,fill=location))+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(~strata)+
  theme_bw()+
  labs(x="",y="Richness")+
  theme(legend.position="none")

ggsave("output/SAB_crabsurvey_richness_by_depthstrata.png",sab_rich_depth_strata ,height = 8,width=6,units="in",dpi=300)


inside_outside <- rich_df%>%
  group_by(YEAR,location,strata)%>%
  summarise(mean=mean(nspecies,na.rm=T))%>%
  ungroup()%>%
  spread(location,mean)%>%
  mutate(diff=Inside-Outside)


inside_outside_plot <- ggplot(inside_outside,aes(x=YEAR,y=diff,group=strata,col=strata))+
  geom_line()+
  geom_point(size=2)+
  geom_hline(yintercept = 0,lty=2)+
  theme_bw()+
  labs(x="",y="Average species richness difference Inside-Outside",col="Depth strata")+
  scale_x_continuous(breaks=2015:2021)

ggsave("output/SAB_crabsurvey_richness_inside_outside.png",inside_outside_plot ,height = 8,width=6,units="in",dpi=300)

## RVdata --

load("r:/Science/CESD/HES_MPAGroup/Data/RVdata/RVDATA_extraction.RData")

rv_data <- rvdat%>%
  filter(!is.na(LATITUDE),!is.na(LONGITUDE))%>%
  dplyr::select(LATITUDE,LONGITUDE)%>%
  mutate(id=paste(LONGITUDE,LATITUDE,sep="_"))%>%
  distinct(id,.keep_all = TRUE)%>%
  st_as_sf(coords=c("LONGITUDE","LATITUDE"),crs=latlong,remove=F)%>%
  st_intersection(bioregion)

rv_hull <- rv_data%>%st_union() %>% 
  st_buffer(0.5, endCapStyle = "SQUARE") %>% 
  st_buffer(-0.5, endCapStyle = "SQUARE")

sab_rv_map <- ggplot()+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=rv_data,size=0.1)+
  geom_sf(data=mar_network,fill=alpha("grey",0.07),col="grey100",lwd=0.25)+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=sab,fill=alpha("grey",0.15),col="black",lwd=0.5)+
  theme_bw()+
  coord_sf(expand=0)

ggsave("output/SAB_rv_map.png",sab_rv_map ,height = 8,width=8,units="in",dpi=300)

#Acoustic map

sab_dets <- read.csv("data/SABMPAfishdets_stations.csv")%>%
             st_as_sf(coords=c("deploy_long","deploy_lat"),crs=latlong)%>%
             st_transform(CanProj)

dets_bound <- sab_dets%>%
              st_bbox()%>%
              st_as_sfc()%>%
              st_as_sf()%>%
              st_buffer(1000*1000)

basemap_dets <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
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
  st_transform(CanProj)%>%
  st_intersection(.,dets_bound)

basemap_dets <- ne_states(country = c("Commonwealth of the Bahamas","Belize","Canada" ,                    
                                      "Republic of Cuba","Dominican Republic","Republic of Honduras","Republic of Haiti",          
                                      "Jamaica","United Mexican States","Commonwealth of Puerto Rico","United States of America"),returnclass = "sf")%>%
  st_transform(CanProj)%>%
  st_intersection(.,dets_bound)

ggplot()+
  geom_sf(data=basemap_dets)+
  geom_sf(data=sab_dets)

#cod residence time

cod_res <- read.csv("R:/Science/CESD/HES_MPAGroup/Presentations/data/SABMPA detection summaries Nov 2022/rik_sum_Acod.csv")%>%
           st_as_sf(coords=c("mean_longitude","mean_latitude"),crs=latlong,remove=FALSE)%>%
           filter(mean_latitude>40)

cod_domain <- st_bbox(cod_res)

cod_domain[1] <- -82
cod_domain[2] <- -75
cod_domain[3] <- 47
cod_domain[4] <- 56

cod_domain <- cod_domain%>%
              st_as_sfc()%>%
              st_as_sf()%>%
              st_transform(CanProj)%>%
              st_bbox()

ggplot()+
  geom_sf(data=basemap_dets%>%st_transform(latlong))+
  geom_sf(data=cod_res%>%st_transform(latlong),aes(size=ri_mean))+
  geom_sf(data=sab%>%st_transform(latlong),fill=NA)+
  coord_sf(xlim=c(-75,-82),ylim=c(47,56),expand=0)

