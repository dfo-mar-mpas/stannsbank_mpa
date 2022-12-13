#Code to make demonstration figures for the SAB OTN Tagging data for the 2020 SAB AC Meeting March 11th 2020

#Load libaries
library(ggplot2)
library(dplyr)
library(treemapify)
library(viridis)
library(lubridate)
library(sf)
library(data.table)
library(scales)
library(ggridges)
library(tidyr)
library(ggridges)

#Temperature data ------
tempdat <- read.csv("c:/Users/stanleyr/Desktop/SAB AC/SABTempData.csv",stringsAsFactors = F)%>%
           mutate(datetime = mdy(date),
                  yr=year(datetime),
                  line=factor(line,levels=c("Southern Line","Northern Line")))%>%
           filter(Temperature.C<9,
                  !is.na(line))

p1 <- ggplot(tempdat,aes(x=datetime,y=Temperature.C,col=line))+
  stat_smooth()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        legend.position = "bottom")+
  scale_x_date(date_labels = "%b-%y",date_breaks = "1 month")+
  labs(x="",
       y=expression(paste("Temperature ",degree,"C")),
       col="")+
  guides(color=guide_legend(override.aes=list(fill=NA)))

ggsave("c:/Users/stanleyr/Desktop/SAB AC/TemperatureSummary.png",p1,dpi=200,width=6,units="in")

MonthSum <- tempdat%>%
          group_by(year,month,line)%>%
          summarise(mn=mean(Temperature.C,na.rm=T),
                    sd=sd(Temperature.C,na.rm=T))%>%
          ungroup()%>%
          mutate(datetime=dmy(paste("01",month,year,sep="-")))%>%
          data.frame()

p2 <- ggplot(MonthSum,aes(x=datetime,y=mn,col=line,group=line))+
  # geom_errorbar(aes(ymin=mn-sd,ymax=mn+sd))+
  geom_line()+
  geom_point(size=1.3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1))+
  scale_x_date(date_labels = "%b-%y",date_breaks = "1 month")+
  labs(x="",
       y=expression(paste("Temperature ",degree,"C")),
       col="")+
  guides(color=guide_legend(override.aes=list(fill=NA)));p2

ggsave("c:/Users/stanleyr/Desktop/SAB AC/TemperatureSummary_means.png",p2,dpi=200,width=6,units="in")


#Seasonal detections ------------
seasondat <- read.csv("c:/Users/stanleyr/Desktop/SAB AC/DetectionTransposed.csv",stringsAsFactors = F)%>%
            mutate(Month=trimws(Month),
                   datetime=dmy(paste("01",Month,Year,sep="-")),
                   datetimestand=dmy(paste("01",Month,2019,sep="-")))%>%
            gather(.,key="species",value="count",White_Shark:American_Eel)%>%
            mutate(species=gsub("_"," ",species),)%>%
            filter(count>0)%>%
            uncount(count)

specieslevel <- table(seasondat$species)%>%data.frame()%>%arrange(Freq)%>%mutate(Var1=as.character(Var1))

seasondat$species <- factor(seasondat$species,levels=specieslevel%>%pull(Var1))

p3 <- ggplot(filter(seasondat,species %in% specieslevel[specieslevel$Freq>2,"Var1"]),
             aes(x=datetime,y=species,fill=species))+
  geom_density_ridges( jittered_points = TRUE,
                       position = position_points_jitter(width = 0.05, height = 0),
                       point_shape = 19, point_size = 1, point_alpha = 1, alpha = 0.7)+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(legend.position="none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(angle=45, hjust=1, vjust = 1))+
  scale_x_date(date_labels = "%b-%y",date_breaks = "3 month")+
  scale_fill_viridis(discrete=T)+
  labs(y="",x="");p3

ggsave("c:/Users/stanleyr/Desktop/SAB AC/TimeSeries_Ridges.png",p3,dpi=600,width=6,units="in")

p4 <- ggplot(data=filter(seasondat,species %in% specieslevel[specieslevel$Freq>2,"Var1"])%>%
               arrange(desc(datetimestand)),
             aes(x=datetimestand,y=species,fill=species))+
  geom_density_ridges(alpha = 0.7)+
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  theme(legend.position="none",
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(angle=45, hjust=1, vjust = 1))+
  scale_x_date(date_labels = "%b",date_breaks = "1 month")+
  scale_fill_viridis(discrete=T)+
  labs(y="",x="");p4

ggsave("c:/Users/stanleyr/Desktop/SAB AC/Monthly_Ridges.png",p4,dpi=600,width=6,units="in")

#Treemap distribution of detections  --------------
treedat <- read.csv("c:/Users/stanleyr/Desktop/SAB AC/uniquecount.csv",stringsAsFactors = F)%>%
  mutate(species=gsub("Atlantic Striped","Atlantic",species),
         countlab = paste0(species," (",count,")"),
         percent = count/sum(count)*100,
         percentlab = paste0(species," (",percent,")"))%>%
  arrange(-count)%>%
  mutate(species=factor(species,levels=species))

p5 <- ggplot(treedat, aes(area = count, fill = species,label=countlab)) +
  geom_treemap()+
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                                   grow = TRUE)+
  theme(legend.position = "none")+scale_fill_viridis(discrete=T);p5

ggsave("c:/Users/stanleyr/Desktop/SAB AC/DetectionsTreemap.png",p5,dpi=600,width=8,height=4,units="in")

#Donut plot -----------

donutdat <- read.csv("c:/Users/stanleyr/Desktop/SAB AC/uniquecount.csv",stringsAsFactors = F)%>%
  arrange(-count)%>%
  mutate(species=factor(species,levels=species),
         countlab = paste0(species," (",count,")"),
         fraction = count/sum(count),
         ymax=cumsum(fraction),
         percent=paste0(species,round(fraction,2)*100,"%"))

donutdat$ymin <- c(0, head(donutdat$ymax, n=-1))

# Make the plot
p6 <- ggplot(donutdat, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=species)) +
  geom_rect(colour="black") +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4))+
  theme_void()+
  scale_fill_viridis(discrete = T)+
  labs(fill="");p6

ggsave("c:/Users/stanleyr/Desktop/SAB AC/Detections_donut.png",p6,dpi=600,width=8,units="in")
