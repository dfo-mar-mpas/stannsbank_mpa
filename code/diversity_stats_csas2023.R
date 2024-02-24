##2024 SAB CSAS Snow Crab survey
#Diversity stats
library(dplyr)
library(ggplot2)
library(pwr2)
library(BiodiversityR)
library(vegan)
library(iNEXT)

load("data/CrabSurvey/2023CrabSurveyDat.RData")


#stole some code from Ryan "2022_Workshop_Plots.R" script
#use data1 dataframe from datareview script for richness
head(data1)
data1$LONGITUDE <- data1$LONGITUDE*-1 #check if this needs to be done depending on what data image you're loading
#Add inside vs outside
tt<-MPATOWS %>% filter(LATITUDE>45.1) #this gets us just the SAB tows that are 'enhanced' stations 
data2 <- data1 %>% mutate(enhanced=as.logical(data1$STATION.x %in% unique(as.integer(tt$STATION))), Year=format(as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y")) %>% replace(is.na(.),FALSE)
#data2$LONGITUDE <- data1$LONGITUDE
#data2$LATITUDE <- data1$LATITUDE
data2<-merge(data2, stns, by.x="STATION.x", by.y="STATION", all.x=T)
centroids <- st_centroid(sabshape)
centroids$X <- st_coordinates(centroids)[,1]
centroids$Y <- st_coordinates(centroids)[,2]


ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  #geom_text(data=centroids, aes(x=X, y=Y, label=Id),nudge_x = -0.019, nudge_y = -0.008, size=5, fontface="bold")+
  #geom_sf(data=gully2, colour="red", fill=NA)+
  coord_sf(xlim=c(-60.5, -58.3), ylim=c(45.25,47.1), expand=F)+
  geom_point(data=data2, aes(x=LONGITUDE, y=LATITUDE, colour=Inside, shape=enhanced, size=enhanced))+
  scale_size_manual(values=c(1.75,2.5))+
  scale_color_manual(values=c("grey46","black"))+
  labs(y="Latitude",x="Longitude")+
  theme_bw()+
  theme(panel.background = element_rect(fill="lightcyan1"), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title = element_text(size=16), legend.text = element_text(size=10), legend.position = "right")
  #ggrepel::geom_label_repel(data=data2, mapping=aes(x=LONGITUDE, y=LATITUDE, label=STATION.x),max.overlaps=2000)

ggsave(filename = "SurveyStationsMap2023_NoZONES.png",plot = last_plot(), device = "png", path = "output/CrabSurvey/", width = 14, height=12, units="in", dpi=600)

#Get depths for each station from the dem 
data2$depth <- raster::extract(ras, data.frame(data2$LONGITUDE, data2$LATITUDE))
data2 <- data2 %>% mutate(location=ifelse(Inside,"Inside","Outside"))
#richness
rich_df <- data2%>%
  data.frame()%>%
  group_by(Year,STATION.x)%>%
  summarise(nspecies = length(unique(SPEC)))%>%
  ungroup()%>%
  left_join(.,data2%>%
              data.frame()%>%
              mutate(id=paste(Year,STATION.x,sep="_"))%>%
              distinct(id,.keep_all=T)%>%
              dplyr::select(STATION.x, LONGITUDE,LATITUDE, Year, location, depth))%>%
  mutate(station=paste("ST",STATION.x,sep="_"),
         strata=ifelse(depth>150,">150m","<150m"),
         year_strat=ifelse(Year>2018,"2019-2023","2015-2018"))


sab_rich_depth <- ggplot(rich_df%>%filter(nspecies>9),aes(x=depth, y=nspecies, col=location, group=location))+
  geom_point(size=3)+
  #facet_wrap(vars(year), nrow=2)+
  stat_smooth(method="lm")+
  #scale_colour_manual(labels=c("Outside","Inside"), values=c("red","blue"))+
  theme_bw()+
  labs(x="Depth (m)",y="Richness",col="")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "bottom", text=element_text(size=20))
sab_rich_depth
ggsave("output/CrabSurvey/richness_by_depth.png",sab_rich_depth,height = 8,width=10,units="in",dpi=500)


rich_count <- rich_df%>%
  group_by(STATION.x)%>%
  summarise(nyears=n(),
            minyear=min(Year),
            maxyear=max(Year))

rich_mod <- glmer(nspecies ~ location + (1 + location | STATION.x),data=rich_df, family=poisson(link="log"))
# Linear mixed model fit by REML ['lmerMod']
# Formula: nspecies ~ Inside + (1 + Inside | STATION.x)
# Data: rich_df
# REML criterion at convergence: 806.3076
# Random effects:
#   Groups    Name        Std.Dev. Corr 
# STATION.x (Intercept) 2.812         
# InsideTRUE  4.698    -0.97
# Residual              3.282         
# Number of obs: 150, groups:  STATION.x, 19
# Fixed Effects:
#   (Intercept)   InsideTRUE  
# 18.3430       0.4695  

ggplot(rich_df,aes(x=Year,y=nspecies,group=STATION.x))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+
  labs(x="Year",y="Richness",col="")+
  facet_grid(~location)

#very-very preliminary - no actual stats.   
richness_plot <- ggplot(rich_df,aes(x=Year,y=nspecies))+
  geom_line(aes(col=STATION.x,group=STATION.x))+
  geom_point(size=2,aes(col=STATION.x))+
  stat_smooth(method="lm",aes(group=location))+
  theme_bw()+
  labs(x="Year",y="Richness",col="")+
  facet_grid(~location)+
  theme(strip.background = element_rect(fill="white"),
         legend.position = "none")
richness_plot
ggsave("output/CrabSurvey/SAB_crabsurvey_richness_by_year.png",richness_plot,height = 8,width=10,units="in",dpi=400)

sab_rich_depth_strata <- ggplot(rich_df,aes(x=location,y=nspecies,fill=location))+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(~strata)+
  theme_bw()+
  labs(x="",y="Richness")+
  theme(legend.position="none")

ggsave("output/CrabSurvey/Richness_by_depthstrata.png",sab_rich_depth_strata ,height = 8,width=6,units="in",dpi=400)


###This part doesn't work yet with my TRUE-FALSE for inside-outside the MPA - fix this later
inside_outside <- rich_df%>%
  group_by(Year,location,strata)%>%
  summarise(mean=mean(nspecies,na.rm=T))%>%
  ungroup()%>%
  spread(location,mean)%>%
  mutate(diff=TRUE-FALSE)


inside_outside_plot <- ggplot(inside_outside,aes(x=Year,y=diff,group=strata,col=strata))+
  geom_line()+
  geom_point(size=2)+
  geom_hline(yintercept = 0,lty=2)+
  theme_bw()+
  labs(x="",y="Average species richness difference Inside-Outside",col="Depth strata")+
  scale_x_continuous(breaks=2015:2023)

ggsave("output/SAB_crabsurvey_richness_inside_outside.png",inside_outside_plot ,height = 8,width=6,units="in",dpi=300)


#Plot net temp per year and station
ggplot()+
  geom_point(data=data2 %>% filter(NET_TEMPERATURE>0), aes(x=as.factor(year), y=NET_TEMPERATURE, colour=depth),size=3)+
  #geom_smooth(method="lm")+
  scale_colour_continuous(breaks=c(75, 125, 175, 225))+
  xlab(label="Year")+
  ylab(label = expression( "Net Temperature ("*degree*'C)'))+
  scale_x_discrete(name="Year", labels=c("2015", "2016", "2017", "2018", "2019", "2021", "2022", "2023"))+
  theme_bw()+
  theme(axis.line=element_line(colour="black"), panel.border = element_blank(),
                                          text=element_text(size=18))+
  guides(col=guide_legend(title = "Depth (m)"))


ggsave("output/CrabSurvey/NetTemperatures2.png",last_plot() ,height = 8,width=12,units="in",dpi=400)

##############################
####Species accumulation curves and NMDS
data3 <- data2 %>% dplyr::select(c(STATION.x, year, SPEC, EST_NUM_CAUGHT, Inside))
data3 <- data3[!(is.na(data3$SPEC) | data3$SPEC==""),]
data4 <- data3 %>% group_by(SPEC) %>% 
  mutate(row=row_number()) %>% 
  pivot_wider(names_from = SPEC, values_from = EST_NUM_CAUGHT,values_fill = 0) %>%
  dplyr::select(-row) %>% as.data.frame()

data5<-as.matrix(data4[,4:length(colnames(data4))])
rownames(data5) <- paste(data4$STATION.x, data4$year, sep="_")

#sabnmds <- metaMDS(data5, distance = "bray", k=2)