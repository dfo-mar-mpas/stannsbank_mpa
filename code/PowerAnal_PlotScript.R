#rbind all results together 
library(dplyr)
library(purrr)
library(ggplot2)


 crab_power <- list.files(pattern=".RDS") %>%
   map_dfr(readRDS)%>%
   mutate(PERIODwrongDirection=if_else(trend=="Increase",
                                       PERIODBefore>0,
                                       PERIODBefore<0))


 resultssummary <- crab_power %>%
   group_by(sp,sites,effect_size,trend) %>%
   reframe(Power=1-sum(p>0.05)/n(),
           PERIODwrongDirection=mean(PERIODwrongDirection)) %>%
   mutate(`Number of Samples`=factor(sites)) %>%
   rename(`Effect Size`=effect_size)



 ggplot(resultssummary,aes(x=`Effect Size`,y=Power,color=`Number of Samples`,group=`Number of Samples`))+
   geom_point()+
   geom_line()+
   geom_abline(intercept = 0.8,slope=0,linetype = "dotted")+
   scale_y_continuous(breaks=seq(0,1,0.2))+
   scale_color_brewer(palette = "Paired")+
   facet_grid(trend~sp)+
   theme_bw()

  file.remove(list.files(pattern=".RDS"))
 saveRDS(crab_power,"data/crab_power.RDS")

crab_power <- readRDS("data/CrabSurvey/crab_power.RDS")
head(crab_power)


# normal power
resultssummary <- crab_power %>%
  group_by(sp,sites,effect_size,trend) %>%
  reframe(Power=sum(p<0.05)/n()) %>%
  mutate(`Number of Samples`=factor(sites)) %>%
  rename(`Effect Size`=effect_size)

ggplot(resultssummary,aes(x=`Effect Size`,y=Power,color=`Number of Samples`,group=`Number of Samples`))+
  geom_point()+
  geom_line()+
  geom_abline(intercept = 0.8,slope=0,linetype = "dotted")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_color_brewer(palette = "Paired")+
  facet_grid(trend~sp)+
  theme_bw()

# false positives removed
resultssummary <- crab_power %>%
  group_by(sp,sites,effect_size,trend) %>%
  reframe(Power=(sum(p<0.05)-sum(p<0.05&PERIODwrongDirection))/n()) %>%
  mutate(`Number of Samples`=factor(sites)) %>%
  rename(`Effect Size`=effect_size)

ggplot(resultssummary,aes(x=`Effect Size`,y=Power,color=`Number of Samples`,group=`Number of Samples`))+
  geom_point()+
  geom_line()+
  geom_abline(intercept = 0.8,slope=0,linetype = "dotted")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_color_brewer(palette = "Paired")+
  facet_grid(trend~sp)+
  theme_bw()

# false positives removed x2
resultssummary <- crab_power %>%
  group_by(sp,sites,effect_size,trend) %>%
  reframe(Power=(sum(p<0.05)-2*sum(p<0.05&PERIODwrongDirection))/n()) %>%
  mutate(`Number of Samples`=factor(sites)) %>%
  rename(`Effect Size`=effect_size)

resultssummary$sp <- stringr::str_wrap(resultssummary$sp, width=8)
resultssummary$sp <- factor(resultssummary$sp, levels=c("Richness", "AMBLYRAJA\nRADIATA","ANARHICHAS\nLUPUS","GADUS\nMORHUA",
                                                        "GLYPTOCEPHALUS\nCYNOGLOSSUS","HIPPOGLOSSOIDES\nPLATESSOIDES","SEBASTES","CHIONOECETES\nOPILIO"))
ggplot(resultssummary,aes(x=`Effect Size`,y=Power,color=`Number of Samples`,group=`Number of Samples`))+
  geom_point()+
  geom_line()+
  geom_abline(intercept = 0.8,slope=0,linetype = "dotted")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_color_brewer(palette = "Paired")+
  facet_grid(trend~sp)+
  theme_bw()+
  theme(strip.background=element_rect(fill="white"), legend.position = "bottom",
        text=element_text(size=12))

ggsave("Final_PowerAnalysis_Plot_8species.png",plot=last_plot(),path = "output/CrabSurvey/",width = 12, height = 8, dpi = 320, bg = "white")

