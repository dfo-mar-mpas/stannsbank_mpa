library(robis)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(arcpullr)
library(vegan)
library(marmap)
library(stars)
library(lme4)
library(MASS)
library(dplyr)
library(tidyr)

# get SAB
url <- "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Oceans_Act_Marine_Protected_Areas/MapServer/0"
sab <- get_spatial_layer(url,where="NAME_E='St. Anns Bank Marine Protected Area'")

# get crab
crab_swept <- read.csv("data/CrabSurvey/sabmpa_area_swept.csv")
spid <- read.csv("data/CrabSurvey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv") %>%
  rowwise() %>%
  mutate(sp_name = if_else(SPEC=="",
                           if_else(COMM=="",as.character(SPECCD_ID),COMM),
                           SPEC))
crab <- read.csv("data/CrabSurvey/SABMPA2023export.csv") %>%
  left_join(crab_swept,by = join_by(TRIP, "SET_NO"=="SET", STATION)) %>%
  mutate(CPUE = EST_NUM_CAUGHT/AREA_SWEPT,
         station_date=paste(STATION,BOARD_DATE,sep = ": "),
         LONGITUDE=-LONGITUDE) %>%
  left_join(spid,by = join_by(SPECCD_ID))


# get bathy
bbsab <- st_bbox(sab)
noaabathy <- getNOAA.bathy(bbsab[1]-1,bbsab[3]+1,bbsab[2]+1,bbsab[4]-1,resolution = 0.25,keep=T) %>%
  fortify.bathy() %>%
  st_as_stars() %>%
  st_set_crs(4326)

# get benthoscape
benthoscape <- st_read("data/Shapefiles/benthoscape.shp")%>%
  st_make_valid() %>%
  st_transform(4326)

# set up crab_stations and join with benthoscape and bathy

crabstations <- crab %>%
  group_by(STATION) %>%
  reframe(LONGITUDE=mean(LONGITUDE),
          LATITUDE=mean(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"),
           crs = 4326) %>%
  st_join(noaabathy %>% st_as_sf()) %>%
  st_join(benthoscape)  %>%
  mutate(Class_name=if_else(is.na(Class_name),
                            "Unclassified",
                            Class_name)) %>%
  dplyr::select(STATION,z,Class_name) %>%
  mutate(STATION = as.character(STATION),
         z_class = if_else(z>=(-100),
                           "<100 m",
                           if_else(z>=(-150),
                                   "100 to 150 m",
                                   ">150 m"))) %>%
  st_join(sab %>% select(ZONE_E)) %>%
  mutate(Boundary=if_else(is.na(ZONE_E),
                 "Outside",
                 "Inside"))

###################### NMDS #####################

comm_crab_ID <- crab %>%
  mutate(sp_name=if_else(is.na(sp_name),
                         as.character(SPECCD_ID),
                         sp_name),
         ID = paste(BOARD_DATE,STATION),
         CPUE = if_else(is.na(CPUE),
                        0,
                        CPUE)) %>%
  pivot_wider(id_cols = c(ID,BOARD_DATE,STATION),
              names_from = sp_name,
              values_from = CPUE,
              values_fill = 0)
ID <- comm_crab_ID %>%
  dplyr::select(ID,BOARD_DATE,STATION) %>%
  mutate(STATION=as.character(STATION))
comm_crab <- comm_crab_ID %>%
  tibble::column_to_rownames("ID")%>%
  #select(-BOARD_DATE,-STATION) %>%
  dplyr::select(-BOARD_DATE,-STATION) %>%
  decostand(method = "total")

zerocatch <- apply(comm_crab_ID %>%
                     tibble::column_to_rownames("ID")%>%
                     dplyr::select(-BOARD_DATE,-STATION),1,sum)>0

comm_crab.env <- ID %>%
  left_join(crabstations,by="STATION") %>%
  mutate(Year=as.numeric(substr(ID,1,4)))
# distmat_crab <- vegdist(comm_crab[zerocatch,],method="bray")


ord <- metaMDS(comm_crab[zerocatch,],k=5)

data.scores <- as.data.frame(scores(ord,"sites")) %>%
  mutate(ID=rownames(.)) %>%
  left_join(comm_crab.env,by="ID")

species.scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

hull.data <- data.scores %>%
  as.data.frame() %>%
  group_by(Class_name) %>%
  slice(chull(x=NMDS1,y=NMDS2))

nmdstheme <- theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# regular ggplot
p <- ggplot() +
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=Class_name),size=3) +
  scale_colour_brewer(palette = "Paired") +
  coord_equal()+
  scale_shape_manual(values=c(15,7,18,16,10,8,17))+
  geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(ord$stress,3),"k =",ord$ndim)))+
  nmdstheme
p



# plotly interaction
plotly::ggplotly(p,tooltip = "text")

# species scores
# selectedsp <- c("HIPPOGLOSSOIDES PLATESSOIDES",
#                 "GLYPTOCEPHALUS CYNOGLOSSUS",
#                 "SEBASTES",
#                 "AMBLYRAJA RADIATA",
#                 "GADUS MORHUA",
#                 "CHIONOECETES OPILIO")

selectedsp <- species.scores %>%
  mutate(vector=sqrt(NMDS1^2+NMDS2^2)) %>%
  arrange(desc(vector)) %>%
  head(10) %>%
  dplyr::select(species) %>%
  unlist()
sp_multfactor <- 0.7

p2 <- p +
  geom_segment(data=species.scores[species.scores$species %in% selectedsp,],
               aes(x=0, y=0, xend=NMDS1*sp_multfactor, yend=NMDS2*sp_multfactor),
               colour="red", arrow=arrow()) +
  geom_text(data=species.scores[species.scores$species %in% selectedsp,],
                  aes(x=NMDS1*sp_multfactor, y=NMDS2*sp_multfactor, label=species),
                  colour="red")
p2


# depth surface #TODO

# surf <- ordisurf(ord,)

# time as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=Year),size=3) +
  scale_colour_distiller(palette = "YlGnBu") +
  coord_equal()+
  theme_bw()
p

# MPA boundary as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=Boundary),size=3) +
  scale_colour_brewer(palette = "YlGnBu") +
  coord_equal()+
  theme_bw()
p

# z class boundary as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=z_class),size=3) +
  scale_colour_brewer(palette = "YlGnBu") +
  coord_equal()+
  theme_bw()
p

############ Summer RV  with inverts ###################
# tempfn <- tempfile()
# download.file("https://api-proxy.edh.azure.cloud.dfo-mpo.gc.ca/catalogue/records/1366e1f1-e2c8-4905-89ae-e10f1be0a164/attachments/SUMMER_csv.zip",tempfn)
# dir.create("data/RV")
# unzip(tempfn,exdir="data/RV/")
RV_info <- read.csv("data/RV/SUMMER_2023_GSINF.csv")
RV_species <- read.csv("data/RV/SUMMER_2023_GSSPECIES.csv")
RV_catch <- read.csv("data/RV/SUMMER_2023_GSCAT.csv")
RV_mission <- read.csv("data/RV/SUMMER_2023_GSMISSIONS.csv")

RV <- RV_catch %>%
  left_join(RV_info,by=c("MISSION","SETNO")) %>%
  left_join(RV_species,by=c("SPEC"="CODE")) %>%
  left_join(RV_mission,by="MISSION") %>%
  filter(YEAR>=2010) %>%
  st_as_sf(coords = c("SLONG","SLAT"),
           crs = 4326,
           remove = FALSE) %>%
  st_intersection(sab %>%
                    st_transform(4326)) %>%
  mutate(sp_name = if_else(SCI_NAME=="",
                           if_else(COMM=="",
                                   as.character(SPEC),
                                   COMM),
                           SCI_NAME),
         CPUE = TOTNO/DIST)


# set up crab_stations and join with benthoscape and bathy

RVstations <- RV %>%
  mutate(ID = paste(MISSION,SLAT,SLONG)) %>%
  group_by(ID) %>%
  mutate(geometry=st_union(geometry)) %>%
  st_join(noaabathy %>% st_as_sf()) %>%
  st_join(benthoscape)  %>%
  mutate(Class_name=if_else(is.na(Class_name),
                            "Unclassified",
                            Class_name)) %>%
  dplyr::select(ID,z,Class_name) %>%
  mutate(z_class = if_else(z>=(-100),
                           "<100 m",
                           if_else(z>=(-150),
                                   "100 to 150 m",
                                   ">150 m"))) %>%
  st_join(sab %>% select(ZONE_E)) %>%
  mutate(Boundary=if_else(is.na(ZONE_E),
                          "Outside",
                          "Inside"))


comm_RV_ID <- RV %>%
  as.data.frame() %>%
  mutate(sp_name=if_else(is.na(sp_name),
                         as.character(SPEC),
                         sp_name),
         ID = paste(MISSION,SLAT,SLONG),
         CPUE = if_else(is.na(CPUE),
                        0,
                        CPUE)) %>%
  pivot_wider(id_cols = c(ID,YEAR),
              names_from = sp_name,
              values_from = CPUE,
              values_fill = 0)
ID <- comm_RV_ID %>%
  dplyr::select(ID,YEAR)

comm_RV <- comm_RV_ID %>%
  tibble::column_to_rownames("ID")%>%
  dplyr::select(-YEAR) #%>%
  #decostand(method = "total")

zerocatch <- apply(comm_RV_ID %>%
                     tibble::column_to_rownames("ID")%>%
                     dplyr::select(-YEAR),1,sum)>0

comm_RV.env <- ID %>%
  left_join(RVstations,by="ID")
# distmat_crab <- vegdist(comm_crab[zerocatch,],method="bray")


ord <- metaMDS(comm_RV[zerocatch,],k=5)

data.scores <- as.data.frame(scores(ord,"sites")) %>%
  mutate(ID=rownames(.)) %>%
  left_join(comm_RV.env,by="ID")

species.scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

hull.data <- data.scores %>%
  as.data.frame() %>%
  group_by(Class_name) %>%
  slice(chull(x=NMDS1,y=NMDS2))


# regular ggplot
p <- ggplot() +
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=Class_name),size=3) +
  scale_colour_brewer(palette = "Paired") +
  coord_equal()+
  theme_bw()
p

# plotly interaction
plotly::ggplotly(p,tooltip = "text")

# species scores
# selectedsp <- c("HIPPOGLOSSOIDES PLATESSOIDES",
#                 "GLYPTOCEPHALUS CYNOGLOSSUS",
#                 "SEBASTES",
#                 "AMBLYRAJA RADIATA",
#                 "GADUS MORHUA",
#                 "CHIONOECETES OPILIO")

selectedsp <- species.scores %>%
  mutate(vector=sqrt(NMDS1^2+NMDS2^2)) %>%
  arrange(desc(vector)) %>%
  head(10) %>%
  dplyr::select(species) %>%
  unlist()
sp_multfactor <- 0.7

p2 <- p +
  geom_segment(data=species.scores[species.scores$species %in% selectedsp,],
               aes(x=0, y=0, xend=NMDS1*sp_multfactor, yend=NMDS2*sp_multfactor),
               colour="red", arrow=arrow()) +
  geom_text(data=species.scores[species.scores$species %in% selectedsp,],
            aes(x=NMDS1*sp_multfactor, y=NMDS2*sp_multfactor, label=species),
            colour="red")
p2


# depth surface #TODO

# surf <- ordisurf(ord,)

# time as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=YEAR),size=3) +
  scale_colour_distiller(palette = "YlGnBu") +
  coord_equal()+
  theme_bw()
p

# MPA boundary as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=Boundary),size=3) +
  scale_colour_brewer(palette = "YlGnBu") +
  coord_equal()+
  theme_bw()
p

# z class boundary as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=z_class),size=3) +
  scale_colour_brewer(palette = "YlGnBu") +
  coord_equal()+
  scale_shape_manual(values=c(15,7,18,16,10,8,17))+
  theme_bw()
p


###### RV "fish" only #####

fish1980 <- RV_catch %>%
  left_join(RV_info,by=c("MISSION","SETNO")) %>%
  left_join(RV_species,by=c("SPEC"="CODE")) %>%
  left_join(RV_mission,by="MISSION") %>%
  mutate(sp_name = if_else(SCI_NAME=="",
                           if_else(COMM=="",
                                   as.character(SPEC),
                                   COMM),
                           SCI_NAME),
         sp_name=if_else(is.na(sp_name),
                         as.character(SPEC),
                         sp_name)) %>%
  filter(YEAR<1980) %>%
  select(sp_name) %>%
  unlist() %>%
  unique()

RV <- RV_catch %>%
  left_join(RV_info,by=c("MISSION","SETNO")) %>%
  left_join(RV_species,by=c("SPEC"="CODE")) %>%
  left_join(RV_mission,by="MISSION") %>%
  st_as_sf(coords = c("SLONG","SLAT"),
           crs = 4326,
           remove = FALSE) %>%
  st_intersection(sab %>%
                    st_transform(4326)) %>%
  mutate(sp_name = if_else(SCI_NAME=="",
                           if_else(COMM=="",
                                   as.character(SPEC),
                                   COMM),
                           SCI_NAME),
         sp_name=if_else(is.na(sp_name),
                         as.character(SPEC),
                         sp_name),
         CPUE = TOTNO/DIST) %>%
  filter(sp_name %in% fish1980)


# set up crab_stations and join with benthoscape and bathy

RVstations <- RV %>%
  mutate(ID = paste(MISSION,SLAT,SLONG)) %>%
  group_by(ID) %>%
  mutate(geometry=st_union(geometry)) %>%
  st_join(noaabathy %>% st_as_sf()) %>%
  st_join(benthoscape)  %>%
  mutate(Class_name=if_else(is.na(Class_name),
                            "Unclassified",
                            Class_name)) %>%
  dplyr::select(ID,z,Class_name) %>%
  mutate(z_class = if_else(z>=(-100),
                           "<100 m",
                           if_else(z>=(-150),
                                   "100 to 150 m",
                                   ">150 m"))) %>%
  st_join(sab %>% select(ZONE_E)) %>%
  mutate(Boundary=if_else(is.na(ZONE_E),
                          "Outside",
                          "Inside"))


comm_RV_ID <- RV %>%
  as.data.frame() %>%
  mutate(sp_name=if_else(is.na(sp_name),
                         as.character(SPEC),
                         sp_name),
         ID = paste(MISSION,SLAT,SLONG),
         CPUE = if_else(is.na(CPUE),
                        0,
                        CPUE)) %>%
  pivot_wider(id_cols = c(ID,YEAR),
              names_from = sp_name,
              values_from = CPUE,
              values_fill = 0)
ID <- comm_RV_ID %>%
  dplyr::select(ID,YEAR)

comm_RV <- comm_RV_ID %>%
  tibble::column_to_rownames("ID")%>%
  dplyr::select(-YEAR) #%>%
#decostand(method = "total")

zerocatch <- apply(comm_RV_ID %>%
                     tibble::column_to_rownames("ID")%>%
                     dplyr::select(-YEAR),1,sum)>0

comm_RV.env <- ID %>%
  left_join(RVstations,by="ID")
# distmat_crab <- vegdist(comm_crab[zerocatch,],method="bray")


ord <- metaMDS(comm_RV[zerocatch,],k=5)

data.scores <- as.data.frame(scores(ord,"sites")) %>%
  mutate(ID=rownames(.)) %>%
  left_join(comm_RV.env,by="ID")

species.scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

hull.data <- data.scores %>%
  as.data.frame() %>%
  group_by(Class_name) %>%
  slice(chull(x=NMDS1,y=NMDS2))


# regular ggplot
p <- ggplot() +
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=Class_name),size=3) +
  scale_colour_brewer(palette = "Paired") +
  coord_equal()+
  theme_bw()
p

# plotly interaction
plotly::ggplotly(p,tooltip = "text")

# species scores
# selectedsp <- c("HIPPOGLOSSOIDES PLATESSOIDES",
#                 "GLYPTOCEPHALUS CYNOGLOSSUS",
#                 "SEBASTES",
#                 "AMBLYRAJA RADIATA",
#                 "GADUS MORHUA",
#                 "CHIONOECETES OPILIO")

selectedsp <- species.scores %>%
  mutate(vector=sqrt(NMDS1^2+NMDS2^2)) %>%
  arrange(desc(vector)) %>%
  head(10) %>%
  dplyr::select(species) %>%
  unlist()
sp_multfactor <- 0.7

p2 <- p +
  geom_segment(data=species.scores[species.scores$species %in% selectedsp,],
               aes(x=0, y=0, xend=NMDS1*sp_multfactor, yend=NMDS2*sp_multfactor),
               colour="red", arrow=arrow()) +
  geom_text(data=species.scores[species.scores$species %in% selectedsp,],
            aes(x=NMDS1*sp_multfactor, y=NMDS2*sp_multfactor, label=species),
            colour="red")
p2


# depth surface #TODO

# surf <- ordisurf(ord,)

# time as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=YEAR),size=3) +
  scale_colour_distiller(palette = "YlGnBu") +
  coord_equal()+
  theme_bw()
p

# MPA boundary as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=Boundary),size=3) +
  scale_colour_brewer(palette = "YlGnBu") +
  coord_equal()+
  theme_bw()
p

# z class boundary as colour
p <- ggplot() +
  # geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Class_name),alpha=0.30) + # add the hulls
  # geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the sp.labels
  geom_point(data=data.scores,aes(text=ID,x=NMDS1,y=NMDS2,shape=Class_name,colour=z_class),size=3) +
  scale_colour_brewer(palette = "YlGnBu") +
  coord_equal()+
  scale_shape_manual(values=c(15,7,18,16,10,8,17))+
  theme_bw()
p
