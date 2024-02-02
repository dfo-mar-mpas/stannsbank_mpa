#######################################################
##Enhanced Snow Crab Survey Data in St. Anns Bank MPA##
#####################Fall 2023#########################
#######################################################

#library(Mar.datawrangling)
library(sf)
library(stringr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(RODBC)
library(ggpubr)
library(tidyr)


#### 0. set projections and options
sf_use_s2 = FALSE
latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#### 1. get the trawl data from our network - this was pulled by Brent Wilson for us in summer 2023 and then again in Jan 2024 to add the 2023 fish morph data. See end of script for R code to pull data too. 
#Get the fish species code from ANDES - these get updated from time to time so may need to replace this csv eventually
fishcodes <- read.csv("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv")
head(fishcodes) #fish codes, aphia-ID's, common and species names

#These are OLD FILES only to 2022 - These files have ALL of the MPA stations (Gully and SAB) and the fish morphology data
#load("R:/Science/CESD/HES_MPAGroup/Data/Crab Survey/EnhancedStations2023/MPAFishMorph.RDATA")
load("data/CrabSurvey/OldData/MPATOWS.RDATA") #load this one so we can make a map of the Gully and SAB stations, this is only the enhanced stations

#Load the full fish morph data and then subset
load("data/CrabSurvey/FishMorph.RData")

#This is just restricted to SAB and includes all stations from 2015 to 2023
catchdat <- read.csv("data/CrabSurvey/SABMPA2023export.csv", header = T) 
catchdat<- catchdat %>% unite(TRIP_STATION, c("TRIP", "STATION"), remove = F)

#Load area swept data
arsw <- read.csv("C:/Users/JEFFERYN/Documents/GitHub/stannsbank_mpa/data/CrabSurvey/sabmpa_area_swept.csv", header = T)
arsw <- arsw %>% unite(TRIP_STATION, c("TRIP", "STATION"), remove = F)
arsw %>% summarise(mean=mean(AREA_SWEPT,na.rm = T), 
                   median=median(AREA_SWEPT, na.rm=T), 
                   stdev=sd(AREA_SWEPT, na.rm=T),
                   se=(sd(AREA_SWEPT)/sqrt(length(AREA_SWEPT))))
#mean        median        stdev          se 
#0.00435671 0.004279407 0.0008242672    6.730114-05

#merge catch data with species codes
catchdat <- merge(catchdat, fishcodes, by.x="SPECCD_ID", by.y="SPECCD_ID", all.x=T) #Amy Glass did a more recent pull in December 2023 so comparing this csv file to the MPATOWS RData object

#now join by area swept file - can use Trip and Station to join
catchdata <- left_join(arsw, catchdat, by = "TRIP_STATION", relationship = "many-to-many", keep = F)
#remove redundant columns
catchdata <- catchdata[,-c(8, 11)]
sab <- catchdata #make this compatible with Amy's annual report code just to make things easier
#head(MPAFishMorph) #fish SPECCD_ID, quantity, weight, length 
head(MPATOWS) #lots of variables about the set catch, including temperature, coordinates, # caught etc 

#First need to remove some redundant columns from MPAWTOWS
MPATOWS <- MPATOWS[,-c(12,13,17,36)]

#Merge species codes into fishmorphs object. This is a big file that needs to be subset by SAB TRIP and STATION
fishmorphs <- merge(fishmorphs, fishcodes, by.x="SPECCD_ID", by.y="SPECCD_ID", all.x=T)


##########################################################################################################
#### 2. Make some basic maps of the study areas
#load MPA polygons
sabshape <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/StAnnsBank_MPA.shp")%>%
  st_transform(latlong)%>%
  mutate(name="St Anns Bank")

gully <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Gully/Gully_Boundary_Redone2014.shp")%>%
  st_transform(st_crs(sabshape))%>%
  mutate(name="Gully MPA")%>%
  dplyr::select(name,geometry)


#basemap of Nova Scotia
novsco <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Coastline/NS_coastline_project_Erase1.shp")%>%st_transform(latlong)%>%
  mutate(name="Nova Scotia")%>%
  dplyr::select(name,geometry)

ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=0.75)+
  geom_sf(data=gully, colour="blue", fill=NA, linewidth=0.75)+
  coord_sf(xlim=c(-62, -58), ylim=c(43.5,47.5), expand=F)+
  geom_point(data=MPATOWS, aes(x=-LONGITUDE, y=LATITUDE), shape=21, fill="black", size=1.05)+
  #geom_point(data=set2023MPA, aes(x=LON, y=LAT), shape=21, fill="black", size=1.25)+
  labs(x="LONGITUDE")+
  theme_minimal()

#ggsave(filename = "output/CrabSurvey/EnhancedCrabStations.png", plot = last_plot(), device = "png", path = "output/", width = 6, height =10, units = "in", dpi=500)

###3. Do some data filtering and merging to just focus on the SAB inside and outside trawl stations
SABTOWS <- catchdata %>% filter(LATITUDE>45.1) #filter MPATOWS data by everything >45 degrees north to remove Gully stations
SABTOWS$LONGITUDE <- SABTOWS$LONGITUDE*-1

#transform into an sf object and find stations inside vs outside the MPA
SABTOWS2 <- SABTOWS %>% st_as_sf(coords = c("LONGITUDE","LATITUDE"),crs=latlong) %>% 
  mutate(Inside=as.logical(st_intersects(.,sabshape, sparse=TRUE)), Year=format(as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y"), LON=sf::st_coordinates(.)[,1], LAT=sf::st_coordinates(.)[,2])%>% 
  replace(is.na(.),FALSE)


#we'll make a better map later
ggplot()+
  geom_sf(data=novsco, fill=gray(.9),size=0)+
  geom_sf(data=sabshape,colour="blue", fill=NA, linewidth=1.25)+
  #geom_sf(data=gully2, colour="red", fill=NA)+
  coord_sf(xlim=c(-61, -58), ylim=c(45.25,47.1), expand=F)+
  geom_point(data=SABTOWS2, aes(x=LON, y=LAT,colour=Inside),size=3)+
  labs(y="Latitude",x="Longitude")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="lightblue1"), text=element_text(size=20))

#ggsave(filename = "SABCrabStations.png", plot = last_plot(), device = "png", path = "output/", width = 10, height =8, units = "in")


########################################################################################################
### 4. Plot up some of the fish catch and morpho data
SABTOWS2$COMM <- str_replace(SABTOWS2$COMM, pattern = "BASKET STARS; GORGONOCEPHALIDAE,ASTERONYCHIDAE", replacement = "BASKET STARS") #let's shorten this name so it plots nicer

SABTOWS3 <- SABTOWS2  %>% filter(!is.na(COMM) & !grepl("SEAWEED, ALGAE ,KELP; THALLOPHYTA|APHIA|BARNACLES|SARSI|BUCCINIDAE|SLUGS|ISOPOD1|PEACH|RUSSIAN|UNIDENTIFIED|SEA MOUSE|EGGS",COMM) & !grepl("ASTERIAS|LAMPRETAE|LYCODES|CARIS|HORMATHIA|BOLOCERA|CALATHURA|DECAPODA|SCLEROCRANGON", SPEC)) 

#Look at fish counts from the SAB tows
# Create boxplots for fish catch by species
ggplot(SABTOWS3, aes(x = Year, y = EST_NUM_CAUGHT, colour=Inside)) +
  facet_wrap(vars(SPEC),nrow=12, scales="free_y")+
  geom_boxplot() +
  labs(title = "# Caught by Species",
       x = "Species",
       y = "# Caught") +
  theme_minimal()+
  theme(axis.text.x =  element_text(angle=90),strip.text.x = element_text(size = 6))

#ggsave(filename = "FishCatchinSAB.png", bg="white", plot=last_plot(), device="png", path="output/CrabSurvey/", width=18, height = 12, units="in", dpi=500)


#####Now look at the morphology data (length and weight)
# Group the data by species and calculate summary statistics
#Also filter by TRIP_ID from the SABTOWS data 
fishmorph <- fishmorphs %>% filter(TRIP_ID %in% SABTOWS3$TRIP_ID)
fishmorph <- fishmorph[!(is.na(fishmorph$COMM) | fishmorph$COMM==""), ]

summary_data <- fishmorph %>% group_by(COMM) %>%
  summarise(
    Median_Length = median(FISH_LENGTH, na.rm = TRUE),
    Mean_Length = mean(FISH_LENGTH, na.rm = TRUE),
    Length_SD=sd(FISH_LENGTH, na.rm=TRUE),
    Median_Weight = median(FISH_WEIGHT, na.rm = TRUE),
    Mean_Weight = mean(FISH_WEIGHT, na.rm = TRUE),
    Weight_SD=sd(FISH_WEIGHT, na.rm=TRUE)
  )

#write.csv(x = summary_data, file = "output/CrabSurvey/Fish_weight_length_summary.csv", quote = F)
#Add a new column in fishmorph for year by extracting year from the BOARD DATE column
fishmorph2 <- fishmorph
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="YELLOWTAIL FLOUNDER; LIMANDA", replacement = "YELLOWTAIL")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="TURBOT,GREENLAND HALIBUT", replacement = "GREENLAND HALIBUT")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="MONKFISH,GOOSEFISH,ANGLER", replacement = "ANGLERFISH")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="HERRING ATLANTIC", replacement = "ATLANTIC HERRING")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="EELPOUT,NEWFOUNDLAND; INCL LYCODES ATLANTICUS; TERRAENOVA", replacement = "ATLANTIC EELPOUT")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="LONGFIN HAKE; UROPHYCIS", replacement = "LONGFIN HAKE")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="SNOW CRAB  QUEEN; INCL CHIONOECETES SP", replacement = "SNOW/QUEEN CRAB")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="SNAKE BLENNY; LUMPRETAEFORMIS", replacement = "SNAKE BLENNY")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="SNOWFLAKE HOOKEAR SCULPIN", replacement = "SNOWFLAKE SCULPIN")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="HOOKEAR SCULPIN,ATL.", replacement = "HOOKEAR SCULPIN")
fishmorph2$COMM <-str_replace(fishmorph2$COMM, pattern="SHORTTAILED EELPOUT VAHL", replacement = "LYCODES SP.")


#Have to replace 2020 with 2019 since there was no 2020 survey
#fishmorph2 <- fishmorph2 %>% mutate(Year=replace(Year, Year==2020, 2019)) %>% as.data.frame()
  
# Create boxplots for fish lengths by species
fishmorph3 <- fishmorph2[!(is.na(fishmorph2$FISH_LENGTH) | fishmorph2$FISH_LENGTH==""),]
fishmorph3 <- fishmorph3[!(is.na(fishmorph3$FISH_WEIGHT) | fishmorph3$FISH_WEIGHT==""),]

ggplot(fishmorph3 %>% filter(EST_NUM_CAUGHT > 2), aes(x = Year, y = log(FISH_LENGTH))) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(vars(COMM),nrow=5, scales="free_y")+
  labs(x = "Species",
       y = "Log(Length)") +
  theme_minimal()+
  theme(axis.text.x =  element_text(angle=90),strip.text.x = element_text(size = 8))#+
  #stat_compare_means(position = "jitter")

ggsave(filename = "CrabSurvey_SAB_FishLengths.png",plot = last_plot(),path="output/CrabSurvey/",device = "png",width = 12, height=10, units = "in",dpi = 500, bg = "white")

ggplot(fishmorph3 %>% filter(EST_NUM_CAUGHT > 40), aes(x = Year, y = FISH_LENGTH)) +
  geom_errorbar(width=0.1)+
  geom_point() +
  facet_wrap(vars(COMM),nrow=5, scales="free_y")+
  labs(x = "Species",
       y = "Length") +
  theme_minimal()+
  theme(axis.text.x =  element_text(angle=90),strip.text.x = element_text(size = 8))#+
#stat_compare_means(position = "jitter")

ggsave(filename = "CrabSurvey_SAB_FishLengths.png",plot = last_plot(),path="output/CrabSurvey/",device = "png",width = 12, height=10, units = "in",dpi = 500, bg = "white")
# Create boxplots for fish weights by species
ggplot(fishmorph3 %>% filter(EST_NUM_CAUGHT>2), aes(x = Year, y = FISH_WEIGHT)) +
  facet_wrap(vars(COMM),nrow=5, scales="free_y")+
  geom_boxplot() +
  labs(x = "Species",
       y = "Weight") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90), strip.text.x = element_text(size=7))

ggsave(filename = "CrabSurvey_SAB_FishWeights.png",plot = last_plot(), path = "output/CrabSurvey/",device = "png",width = 12, height=10, units = "in",dpi = 500, bg = "white")

save.image("data/CrabSurvey/2023CrabSurveyDat.RData")


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




  