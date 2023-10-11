#######################################################
##Enhanced Snow Crab Survey Data in St. Anns Bank MPA##
#####################Fall 2023#########################
#######################################################

library(sf)
library(stringr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(Mar.datawrangling)
library(RODBC)
library(ggpubr)

sf_use_s2 = FALSE

#projections
latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#get the trawl data from our network
load("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/EnhancedStations2023/MPAFishMorph.RDATA")
load("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/EnhancedStations2023/MPATOWS.RDATA")

#Get the fish species code from ANDES
fishcodes <- read.csv("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv")
head(fishcodes) #fish codes, aphia-ID's, common and species names
head(MPAFishMorph) #fish SPECCD_ID, quantity, weight, length 
head(MPATOWS) #lots of variables about the set catch, including temperature, coordinates, # caught etc 

#First need to remove some redundant columns from MPAWTOWS
MPATOWS <- MPATOWS[,-c(12,13,17,36)]

#Merge species codes into MPAFishMorph
fishmorph <- merge(MPAFishMorph, fishcodes, by.x="SPECCD_ID", by.y="SPECCD_ID", all.x=T)


#load polygons
sab <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/StAnnsBank_MPA.shp")%>%
  st_transform(latlong)%>%
  mutate(name="St Anns Bank")

gully <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Gully/Gully_Boundary_Redone2014.shp")%>%
  st_transform(st_crs(sab))%>%
  mutate(name="Gully MPA")%>%
  dplyr::select(name,geometry)

mpas <- rbind(sab,gully)

#load basemap

admin<-rnaturalearth::ne_countries(scale="large", country = "Canada", returnclass="sf") #grabs all of Canada, we'll narrow this down to the northwest Atlantic when plotting


ggplot()+
  geom_sf(data=admin, fill=gray(.9),size=0)+
  geom_sf(data=sab,colour="red", fill=NA)+
  geom_sf(data=gully, colour="red", fill=NA)+
  coord_sf(xlim=c(-62, -58), ylim=c(43.5,48), expand=F)+
  geom_point(data=MPATOWS, aes(x=-LONGITUDE, y=LATITUDE), shape=21, fill="black", size=1.25)+
  theme_minimal()


#filter MPATOWS data by everything >45 degrees north to remove Gully stations 
SABTOWS <- MPATOWS %>% filter(LATITUDE>45.1)
SABTOWS$LONGITUDE <- SABTOWS$LONGITUDE*-1

#transform into an sf object and find stations inside vs outside the MPA
SABTOWS2 <- SABTOWS %>% st_as_sf(coords = c("LONGITUDE","LATITUDE"),crs=latlong) %>% 
  mutate(inside=as.numeric(st_intersects(.,sab, sparse=TRUE)), Year=format(as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y"))

SABTOWS2 <- merge(SABTOWS2, fishcodes,by.x="SPECCD_ID", by.y="SPECCD_ID", all.x=T)

SABTOWS2$COMM <- str_replace(SABTOWS2$COMM, pattern = "BASKET STARS; GORGONOCEPHALIDAE,ASTERONYCHIDAE", replacement = "BASKET STARS")

SABTOWS2 <- SABTOWS2  %>% filter(!is.na(COMM)) %>% filter(.,!grepl("SEAWEED, ALGAE ,KELP; THALLOPHYTA", "", "APHIA IS FOR GENUS; WAS TOSSIA"))

ggplot()+
  geom_sf(data=admin, fill=gray(.9),size=0)+
  geom_sf(data=sab,colour="blue", fill=NA, linewidth=1.25)+
  #geom_sf(data=gully2, colour="red", fill=NA)+
  coord_sf(xlim=c(-61, -58), ylim=c(45.25,47.1), expand=F)+
  geom_point(data=SABTOWS, aes(x=LONGITUDE, y=LATITUDE), shape=21, fill="black", size=1.25)+
  theme_minimal()

#Look at fish counts from the SAB tows
# Create boxplots for fish lengths by species
ggplot(SABTOWS2, aes(x = COMM, y = log(EST_NUM_CAUGHT))) +
  facet_wrap(vars(Year),nrow=7)+
  geom_boxplot() +
  labs(title = "Fish Lengths by Species",
       x = "Species",
       y = "# Caught") +
  theme_minimal()+
  theme(axis.text.x =  element_text(angle=90))

#####Nowlook at the morphology data (length and weight)
# Group the data by species and calculate summary statistics
summary_data <- fishmorph %>% na.omit(fishmorph) %>%
  group_by(COMM) %>%
  summarise(
    Mean_Length = mean(FISH_LENGTH, na.rm = TRUE),
    Median_Length = median(FISH_LENGTH, na.rm = TRUE),
    Mean_Weight = mean(MEASURED_WGT, na.rm = TRUE),
    Median_Weight = median(MEASURED_WGT, na.rm = TRUE)
  )

#Add a new column in fishmorph for year by extracting year from the BOARD DATE column
fishmorph <- na.omit(fishmorph) %>% mutate(Year = format
              (as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y"))
  
# Create boxplots for fish lengths by species
ggplot(fishmorph, aes(x = Year, y = FISH_LENGTH)) +
  facet_wrap(vars(COMM),nrow=8, scales="free_y")+
  geom_boxplot() +
  labs(title = "Fish Lengths by Species",
       x = "Species",
       y = "Length") +
  theme_minimal()+
  theme(axis.text.x =  element_text(angle=90))+
  stat_compare_means()

ggsave(filename = "CrabSurvey_SAB_FishLengths.png",plot = last_plot(),device = "png",width = 22, height=14, units = "in",dpi = 500, bg = "white")
# Create boxplots for fish weights by species
ggplot(fishmorph, aes(x = Year, y = MEASURED_WGT)) +
  facet_wrap(vars(COMM),nrow=6, scales="free_y")+
  geom_boxplot() +
  labs(title = "Fish Weights by Species",
       x = "Species",
       y = "Weight") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))

######################################################################################################################
#' For extractions only, you need all of this:

data.dir <- "r:/Science/CESD/HES_MPAGroup/Data/Crab Survey" 

#load the data unless you want it updated and if that is the case run/update code below
load("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/2022-10-31_SAB_crabsurvey_mpastations.RData") #or the most recent data
load("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/2022-10-31_SAB_crabsurvey_mpastations_lf.RData")

# get_data('isdb',usepkg = 'rodbc',
#          fn.oracle.username = "STANLEYR",
#          fn.oracle.password = "Homerun7",
#          fn.oracle.dsn = "PTRAN",
#          data.dir = data.dir)
#' #' Subsequent runs only need the line below, and you will never 
#' #' re-extract unless you add "force.extract = T" (and if you do 
#' #' want to re-extract, you need to include all of the oracle 
#' #' credential stuff from above)
#' 
#'     #' I see 2 potential triptype codes that I thought might be relevant
#'     #' 7061 = "SNOWCRAB SURVEY" 
#'     #' 7064 = "SNOWCRAB RESEARCH" - no data under this code
#'     ISTRIPTYPECODES <- ISTRIPTYPECODES[ISTRIPTYPECODES$TRIPCD_ID %in% c(7061), ]
#'     self_filter(keep_nullsets = T)
    
# get_data('isdb', data.dir = data.dir) 
# 
# #now select just the crab survey data
# ISTRIPTYPECODES <- ISTRIPTYPECODES[ISTRIPTYPECODES$TRIPCD_ID %in% c(7061), ]
# ISSETTYPECODES <- ISSETTYPECODES[ISSETTYPECODES$SETCD_ID == 11,] #from brent - note if you want all crab data you don't do this step
# self_filter(keep_nullsets = T)
# 
# crabdat <- summarize_catches()
#crab_lf <- ISFISH #this is trimped after the self_filter to the MPA stations

#save(crabdat,file=paste("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/",format(Sys.time(), "%Y-%m-%d"),"_SAB_crabsurvey_mpastations", ".RData", sep = ""))
#save(crab_lf,file=paste("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/",format(Sys.time(), "%Y-%m-%d"),"_SAB_crabsurvey_mpastations_lf", ".RData", sep = ""))





ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=mpas,fill=NA)+
  geom_sf(data=crab_MPA)
  theme_bw()+
  coord_sf(expand=0,ylim=c(43.5,47))
  