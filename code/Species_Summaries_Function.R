summary_plot_fun <- function(x,stns,catchdat_stand,fishmorph_df){
  
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  
  #x is a species common name -- e.g., "Snow crab", "Atlantic cod" - this will only work in the github repo with the right file

#note that wolfish is sampled just inside the MPA
wf_crab_df <- catchdat_stand%>%
              mutate(common = ifelse(common == "Witch","Witch flounder",common),#this is a 'witch' in one dataset
                     common = ifelse(species_filter == "SEBASTES SP.","Redfish sp",common))%>% #redfish has some issues and has to be fixed
              filter(common == x)%>%
              mutate(year=year(as.POSIXct(BOARD_DATE)),
                     dataset = "catch")%>%
              select(station,type,location,number,weight,year,dataset)%>%
              filter(location=="Inside")

wf_size <- fishmorph_df%>%
            mutate(station=as.integer(STATION),
                   species2 = ifelse(species2 == "Redfish unseparated","Redfish sp",species2))%>% #redfish naming issue. 
            filter(station %in% stns$STATION,
                   species2 == x)%>%
            left_join(.,stns%>%rename(station=STATION))%>%
            mutate(dataset="length",
                   location = ifelse(Inside,"Inside","Outside"),
                   type = ifelse(Enhanced,"Enhanced","Standard"))%>%
            filter(location == "Inside")

wf_stations <- wf_size%>%
                group_by(station,year)%>%
                summarise(med=median(FISH_LENGTH,na.rm=T),
                          mn=mean(FISH_LENGTH,na.rm=T),
                          sd=sd(FISH_LENGTH,na.rm=T),
                          bigfish = quantile(FISH_LENGTH,0.9))%>%
                ungroup()

#make the plots
wf_size_plot <- ggplot()+
  geom_vline(xintercept = 2017,lty=2)+
  geom_line(aes(x=year,y=mn,group=factor(station),col=factor(station)),data=wf_stations,lwd=0.25,alpha=0.25)+
  geom_linerange(aes(x=year,y=mn,group=factor(station),col=factor(station),ymin=mn-sd,ymax=mn+sd),data=wf_stations)+
  geom_point(aes(x=year,y=mn,group=factor(station),col=factor(station)),data=wf_stations,size=2)+
  theme_bw()+
  scale_x_continuous(limits=c(2015,2023))+
  theme(legend.position = "none")+
  labs(y=bquote(bar(x) ~ " size (cm) ± sd"),x="",title="a)")+
  stat_smooth(data=wf_stations,aes(x=year,y=mn),method="lm")

wf_size_big_plot <- ggplot()+
  geom_vline(xintercept = 2017,lty=2)+
  geom_line(aes(x=year,y=bigfish,group=factor(station),col=factor(station)),data=wf_stations,lwd=0.25,alpha=0.25)+
  geom_point(aes(x=year,y=bigfish,group=factor(station),col=factor(station)),data=wf_stations,size=2)+
  theme_bw()+
  scale_x_continuous(limits=c(2015,2023))+
  theme(legend.position = "none")+
  labs(y="90th percentile size (cm)",x="",title="b)")+
  stat_smooth(data=wf_stations,aes(x=year,y=bigfish),method="lm")

wf_catch_plot <- ggplot()+
  geom_vline(xintercept = 2017,lty=2)+
  geom_line(aes(x=year,y=number,group=factor(station),col=factor(station)),data=wf_crab_df,lwd=0.25,alpha=0.25)+
  geom_point(aes(x=year,y=number,group=factor(station),col=factor(station)),data=wf_crab_df,size=2)+
  theme_bw()+
  scale_x_continuous(limits=c(2015,2023))+
  theme(legend.position = "none",)+
  labs(y="Number",x="",title="c)")+
  stat_smooth(data=wf_crab_df,aes(x=year,y=number),method="lm")

wf_biomass_plot <- ggplot()+
  geom_vline(xintercept = 2017,lty=2)+
  geom_line(aes(x=year,y=weight,group=factor(station),col=factor(station)),data=wf_crab_df,lwd=0.25,alpha=0.25)+
  geom_point(aes(x=year,y=weight,group=factor(station),col=factor(station)),data=wf_crab_df,size=2)+
  theme_bw()+
  scale_x_continuous(limits=c(2015,2023))+
  theme(legend.position = "none")+
  labs(y="Biomass (kg)",x="",title="d)")+
  stat_smooth(data=wf_crab_df,aes(x=year,y=weight),method="lm")

wf_combo_plot <- (wf_size_plot+wf_size_big_plot)/(wf_catch_plot+wf_biomass_plot)

return(wf_combo_plot)

} # end function 