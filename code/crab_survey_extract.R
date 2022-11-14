library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(Mar.datawrangling)
library(RODBC)

sf_use_s2 = FALSE

#projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load crab data
load("R:/Science/CESD/HES_MPAGroup/Projects/Gully/Data/SnowCrab.set.mpa.rdata")

crab_df <- set%>%
           dplyr::select(lon,lat,trip,set,station)%>%
           distinct(station,.keep_all = TRUE)%>%
           st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)

crab_dataframe <- crab_df%>%
                  data.frame()%>%
                  dplyr::select(-geometry)%>%
                  mutate(mpa = ifelse(lat>45,"St Anns Bank","Gully"),
                         mpa = factor(mpa,levels=c("St Anns Bank","Gully")))%>%
                  arrange(mpa,station)

crab_stations <- crab_dataframe$station

#write.csv(crab_dataframe,file = "R:/Science/CESD/HES_MPAGroup/Projects/St Anns Bank/data/snowcrab_stations.csv",row.names=F)

crab_df_long <- set%>%
              dplyr::select(lon,lat,trip,set,station,yr)%>%
              st_as_sf(coords=c("lon","lat"),crs=latlong)

#load polygons
sab <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/StAnnsBank_MPA.shp")%>%
  st_buffer(0.00005)%>%
  st_union()%>%
  st_as_sf()%>%
  st_transform(latlong)%>%
  mutate(name="St Anns Bank")%>%
  rename(geometry=x)%>%
  dplyr::select(name,geometry)

gully <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Gully/Gully_Boundary_Redone2014.shp")%>%
  st_transform(st_crs(sab))%>%
  mutate(name="Gully MPA")%>%
  dplyr::select(name,geometry)

mpas <- rbind(sab,gully)

#load basemap

basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="Canada")%>%
           st_intersection(.,mpas%>%st_bbox()%>%st_as_sfc()%>%st_as_sf()%>%st_buffer(2))

ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=mpas,fill=NA)+
  geom_sf(data=crab_df)+
  theme_bw()+
  coord_sf(expand=0,ylim=c(43.5,47))

ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=mpas,fill=NA)+
  geom_sf(data=crab_df_long)+
  facet_wrap(~yr,nrow=2)+
  theme_bw()+
  coord_sf(expand=0,ylim=c(43.5,47))



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
  