#rbind all results together 
library(dplyr)
library(purrr)
library(ggplot2)

rr <- list.files(pattern=".RDS") %>% 
  map_dfr(readRDS)


resultssummary <- rr %>%
  group_by(sp,sites,effect_size,trend) %>%
  reframe(Power=1-sum(p>0.05)/n()) %>%
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
