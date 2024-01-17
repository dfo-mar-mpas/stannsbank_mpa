##2024 SAB CSAS Snow Crab survey
#Diversity stats
library(dplyr)
library(ggplot2)
library(pwr2)
library(BiodiversityR)
library(vegan)
library(iNEXT)

load("data/2023CrabSurveyDat.RData")


#stole some code from Ryan "2022_Workshop_Plots.R" script
#use data1 dataframe from datareview script for richness
head(data1)
data1$LONGITUDE <- data1$LONGITUDE*-1 #check if this needs to be done depending on what data image you're loading
#Add inside vs outside
tt<-MPATOWS %>% filter(LATITUDE>45.1) #this gets us just the SAB tows that are 'enhanced' stations 
data2 <- data1 %>% mutate(enhanced=as.logical(data1$STATION.x %in% unique(as.integer(tt$STATION)))) %>% replace(is.na(.),FALSE)
#data2$LONGITUDE <- data1$LONGITUDE
#data2$LATITUDE <- data1$LATITUDE

  
ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1)+
  #geom_sf(data=gully2, colour="red", fill=NA)+
  coord_sf(xlim=c(-60.5, -58.3), ylim=c(45.25,47.1), expand=F)+
  geom_point(data=data2, aes(x=LONGITUDE, y=LATITUDE, colour=Inside, shape=enhanced, size=enhanced))+
  scale_size_manual(values=c(1.75,2.5))+
  scale_color_manual(values=c("grey46","black"))+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="lightcyan1"), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title = element_text(size=16), legend.text = element_text(size=10))

ggsave(filename = "SurveyStationsMap2023.png",plot = last_plot(), device = "png", path = "output/CrabSurvey/", width = 16, height=10, units="in", dpi=600)

#Get depths for each station from the dem 
data2$depth <- raster::extract(ras, data.frame(data2$LONGITUDE, data2$LATITUDE))
#richness
rich_df <- data2%>%
  data.frame()%>%
  group_by(year,STATION.x)%>%
  summarise(nspecies = length(unique(SPEC)))%>%
  ungroup()%>%
  left_join(.,data2%>%
              data.frame()%>%
              mutate(id=paste(year,STATION.x,sep="_"))%>%
              distinct(id,.keep_all=T)%>%
              dplyr::select(STATION.x, LONGITUDE,LATITUDE, year, Inside, depth))%>%
  mutate(station=paste("ST",STATION.x,sep="_"),
         strata=ifelse(depth>150,">150m","<150m"),
         year_strat=ifelse(year>2018,"2019-2023","2015-2018"))


sab_rich_depth <- ggplot(rich_df%>%filter(nspecies>9),aes(x=depth, y=nspecies, col=Inside, group=Inside))+
  geom_point()+
  #facet_wrap(vars(year), nrow=2)+
  stat_smooth(method="lm")+
  theme_bw()+
  labs(x="Depth (m)",y="Richness",col="")+
  theme(legend.position = "none")
sab_rich_depth
ggsave("output/CrabSurvey/richness_by_depth.png",sab_rich_depth,height = 8,width=8,units="in",dpi=300)


rich_count <- rich_df%>%
  group_by(STATION.x)%>%
  summarise(nyears=n(),
            minyear=min(year),
            maxyear=max(year))

rich_mod <- lmer(nspecies ~ Inside + (1 + Inside | STATION.x),data=rich_df)
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

ggplot(rich_df,aes(x=year,y=nspecies,group=STATION.x))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+
  labs(x="Year",y="Richness",col="")+
  facet_grid(~Inside)

#very-very preliminary - no actual stats.   
richness_plot <- ggplot(rich_df,aes(x=year,y=nspecies))+
  geom_line(aes(col=STATION.x,group=STATION.x))+
  geom_point(size=2,aes(col=STATION.x))+
  stat_smooth(method="lm",aes(group=Inside))+
  theme_bw()+
  labs(x="Year",y="Richness",col="")+
  facet_grid(~Inside)+
  theme(legend.position="none")
richness_plot
ggsave("output/CrabSurvey/SAB_crabsurvey_richness_by_year.png",richness_plot,height = 8,width=6,units="in",dpi=300)

sab_rich_depth_strata <- ggplot(rich_df,aes(x=Inside,y=nspecies,fill=Inside))+
  geom_boxplot()+
  geom_jitter()+
  facet_grid(~strata)+
  theme_bw()+
  labs(x="",y="Richness")+
  theme(legend.position="none")

ggsave("output/CrabSurvey/Richness_by_depthstrata.png",sab_rich_depth_strata ,height = 8,width=6,units="in",dpi=300)


###This part doesn't work yet with my TRUE-FALSE for inside-outside the MPA - fix this later
inside_outside <- rich_df%>%
  group_by(year,Inside,strata)%>%
  summarise(mean=mean(nspecies,na.rm=T))%>%
  ungroup()%>%
  spread(Inside,mean)%>%
  mutate(diff=TRUE-FALSE)


inside_outside_plot <- ggplot(inside_outside,aes(x=year,y=diff,group=strata,col=strata))+
  geom_line()+
  geom_point(size=2)+
  geom_hline(yintercept = 0,lty=2)+
  theme_bw()+
  labs(x="",y="Average species richness difference Inside-Outside",col="Depth strata")+
  scale_x_continuous(breaks=2015:2023)

ggsave("output/SAB_crabsurvey_richness_inside_outside.png",inside_outside_plot ,height = 8,width=6,units="in",dpi=300)


#Plot net temp per year and station
ggplot()+
  geom_point(data=data2, aes(x=year, y=NET_TEMPERATURE, colour=as.factor(STATION.x)),size=3)+
  #geom_smooth(method="lm")+
  ylab(label =expression( "Net Temperature ("*degree*'C)'))+
  theme_bw()+
  theme(text=element_text(size=15))+
  guides(col=guide_legend(title = "Survey Station"))


ggsave("output/CrabSurvey/NetTemperatures.png",last_plot() ,height = 8,width=12,units="in",dpi=400)

##############################
####Species accumulation curves
data3 <- data2 %>% dplyr::select(c(STATION.x, year, SPEC, EST_NUM_CAUGHT, Inside))
data3 <- data3[!(is.na(data3$SPEC) | data3$SPEC==""),]
data4 <- data3 %>% group_by(SPEC) %>% 
  mutate(row=row_number()) %>% 
  pivot_wider(names_from = SPEC, values_from = EST_NUM_CAUGHT,values_fill = 0) %>%
  dplyr::select(-row) %>% as.data.frame()

divdat <- iNEXT(data4, datatype = "abundance")



ggiNEXT(sad_data, type = 1, facet.var = "Year")
# Type 1: Chao1 richness estimator
# Type 2: Bootstrap-based richness estimator
# Type 3: Sample-size-based rarefaction/extrapolation
# Type 4: Sample-size-based extrapolation for incomplete sampling


#Power analyses - looking at how many trawl sets it takes to achieve high statistical power to determine biodiversity/richness (and abundance or biomass of some species?)

