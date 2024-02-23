library(sf)
library(vegan)
library(robis)
library(ggplot2)
library(arcpullr)
library(dplyr)
library(tidyr)

# get SAB
url <- "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Oceans_Act_Marine_Protected_Areas/MapServer/0"
sab <- get_spatial_layer(url,where="NAME_E='St. Anns Bank Marine Protected Area'")

#load crab data and clean it up for sp acc
crab <- read.csv("../stannsbank_mpa/data/CrabSurvey/SABMPA2023export.csv") %>%
  mutate(LONGITUDE=-LONGITUDE) %>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"),
           crs = 4326) %>%
  mutate(inSAB=as.vector(st_intersects(.,st_union(sab),sparse = FALSE))) %>%
  group_by(inSAB,TRIP,STATION,SPECCD_ID) %>%
  reframe(n=sum(EST_NUM_CAUGHT,na.rm = TRUE)) %>%
  pivot_wider(id_cols = c(inSAB,TRIP,STATION),names_from = SPECCD_ID,values_from = n,values_fill = 0) %>%
  dplyr::select(-c(TRIP,STATION))

# sp acc for crabs
#in SAB
spacc <- specaccum(crab %>% filter(inSAB) %>% dplyr::select(-inSAB),"random")

# prep data for ggplot
plotdata <- data.frame(Sites=spacc$sites,
                       Richness=spacc$richness,
                       SD=spacc$sd,
                       Source = "Snow Crab - inside SAB")


# out SAB
spacc <- specaccum(crab %>% filter(!inSAB) %>% dplyr::select(-inSAB),"random")
plotdata <- bind_rows(plotdata,
                      data.frame(Sites=spacc$sites,
                                 Richness=spacc$richness,
                                 SD=spacc$sd,
                                 Source="Snow Crab - outside SAB"))


# plot
ggplot(plotdata %>%
         mutate(Source=factor(Source,
                              levels=c("Snow Crab - outside SAB",
                                       "Snow Crab - inside SAB")))) +
  geom_line(aes(x=Sites, y=Richness,col=Source),linewidth=1) +
  geom_ribbon(aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD),fill=Source),alpha=0.2)+
  coord_cartesian(ylim = c(0, 150), xlim = c(0,50))+
  theme_classic()


