##2024 SAB CSAS Snow Crab survey
#Diversity stats
library(dplyr)
library(ggplot2)
library(pwr2)
library(BiodiversityR)
library(vegan)
library(iNEXT)

load("data/2023CrabSurveyDat.RData")


#stole some code from Ryan
#use data1 dataframe from datareview script for richness
head(data1)
data1$LONGITUDE <- data1$LONGITUDE*-1 #check if this needs to be done depending on what data image you're loading
#Add inside vs outside
data2 <- data1 %>% st_as_sf(coords = c("LONGITUDE","LATITUDE"),crs=latlong) %>% 
  mutate(Inside=as.logical(st_intersects(.,sab, sparse=TRUE)), Year=format(as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y"), enhanced=as.logical(STATION %in% unique(as.integer(SABTOWS3$STATION))))  %>% replace(is.na(.),FALSE)
data2$LONGITUDE <- data1$LONGITUDE
data2$LATITUDE <- data1$LATITUDE

  
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

ggsave(filename = "SurveyStationsMap2023.png",plot = last_plot(), device = "png", path = "output/", width = 16, height=10, units="in", dpi=600)


#richness
rich_df <- data2%>%
  data.frame()%>%
  group_by(Year,STATION)%>%
  summarise(nspecies = length(unique(SPEC)))%>%
  ungroup()%>%
  left_join(.,data2%>%
              data.frame()%>%
              mutate(id=paste(Year,STATION,sep="_"))%>%
              distinct(id,.keep_all=T)%>%
              dplyr::select(STATION,LONGITUDE,LATITUDE,Year,Inside))#%>%
  mutate(station=paste("ST",STATION,sep="_"),
         strata=ifelse(DEPTH>150,">150m","<150m"),
         year_strat=ifelse(Year>2018,"2019-2022","2015-2018"))


sab_rich_depth <- ggplot(rich_df%>%filter(nspecies>9),aes(x=Inside,y=nspecies))+
  geom_boxplot()+
  #stat_smooth(method="lm")+
  theme_bw()+
  labs(x="Year",y="Richness",col="")+
  theme(legend.position = "none")
sab_rich_depth

#Plot net temp per year and station
ggplot()+
  geom_point(data=data2, aes(x=Year, y=NET_TEMPERATURE, colour=as.factor(STATION)),size=3)+
  ylab(label =expression( "Net Temperature ("*degree*'C)'))+
  theme_bw()+
  theme(text=element_text(size=15))+
  guides(col=guide_legend(title = "Survey Station"))

#Species accumulation curves
sad_data <- data3 %>% as.data.frame() %>%
  group_by(Year, STATION) %>%
  summarize(species_counts = n()) 


ggiNEXT(sad_data, type = 1, facet.var = "Year")
# Type 1: Chao1 richness estimator
# Type 2: Bootstrap-based richness estimator
# Type 3: Sample-size-based rarefaction/extrapolation
# Type 4: Sample-size-based extrapolation for incomplete sampling


#Power analyses - looking at how many trawl sets it takes to achieve high statistical power to determine biodiversity/richness (and abundance or biomass of some species?)

