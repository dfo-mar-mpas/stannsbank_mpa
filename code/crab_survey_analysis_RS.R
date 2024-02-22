#load libraries
    library(sf)
    library(stringr)
    library(dplyr)
    library(ggplot2)
    library(rnaturalearth)
    library(RODBC)
    library(ggpubr)
    library(tidyr)

#### global options---------
    sf_use_s2 = FALSE
    
#### function -----
    
    capitalize_first <- function(text) { #from ChatGPT!
      words <- strsplit(text, " ")[[1]] 
      words[[1]] <- paste0(toupper(substr(words[[1]], 1, 1)), substr(words[[1]], 2, nchar(words[[1]])))
      capitalized_text <- paste(words, collapse = " ")  # Recombine the words into a single string
      return(capitalized_text)
    }

#projections -----------
    latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
    utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load polygons -----------
    sab <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
      st_transform(latlong)%>%
      mutate(name="St Anns Bank")
    
    sab_nozones  <- sab%>%
      st_transform(utm)%>%
      st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
      st_union()%>% #gets rid of the zones
      st_transform(latlong)%>%
      st_as_sf()%>%
      rename(geometry=x)%>%
      mutate(name="St Anns Bank")%>%
      dplyr::select(name,geometry)
    
    gully <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Gully/Gully_Boundary_Redone2014.shp")%>%
      st_transform(latlong)%>%
      mutate(name="Gully MPA")%>%
      dplyr::select(name,geometry)
    
    mpas <- rbind(sab_nozones,gully)

    #basemap of Nova Scotia
    novsco <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Coastline/NS_coastline_project_Erase1.shp")%>%st_transform(latlong)%>%
      mutate(name="Nova Scotia")%>%
      dplyr::select(name,geometry)

# fish taxonomy metadata -----------
fishcodes <- read.csv("data/CrabSurvey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv")

# Crab survey location metadata ---------------
stns<-read.csv("data/CrabSurvey/StationInfo.csv",header = T)
colnames(stns)<-c("STATION", "Enhanced", "Inside")

#Crab Survey data ----------

    #Load area swept data
    arsw <- read.csv("data/CrabSurvey/sabmpa_area_swept.csv", header = T)
    arsw <- arsw %>% unite(TRIP_STATION, c("TRIP", "STATION"), remove = F)

    #This is just restricted to SAB and includes all stations from 2015 to 2023
    catchdat <- read.csv("data/CrabSurvey/SABMPA2023export.csv", header = T)%>% 
                unite(TRIP_STATION, c("TRIP", "STATION"), remove = F)%>%
                left_join(.,fishcodes,by="SPECCD_ID")%>%
                left_join(arsw, catchdat, by = "TRIP_STATION", relationship = "many-to-many", keep = F)%>%
                dplyr::select(-grep(".y",names(.)))%>% #get rid of the duplicate names
                rename_with(~str_remove(., '.x'))


#Format the crab data
    load("data/CrabSurvey/OldData/MPATOWS.RDATA") #load this one so we can make a map of the Gully and SAB stations, this is only the enhanced stations
    
    MPATOWS <- MPATOWS%>%subset(., select = which(!duplicated(names(.))))
    
    MPATOWS_sf <- MPATOWS%>%
      mutate(trawlID = paste(TRIP,STATION,sep="_"))%>%
      group_by(trawlID)%>%
      summarise(lat_s = mean(START_LAT),
                lon_s = mean(START_LONG)*-1,
                lat_e = mean(END_LAT),
                lon_e = mean(END_LONG)*-1)%>%
      ungroup()%>%
      data.frame()
    
    #get the midpoints for plotting
    df_midpoint <- rbind(MPATOWS_sf%>%mutate(startend="start")%>%dplyr::select(trawlID,lat_s,lon_s)%>%rename(lat = lat_s, lon = lon_s),
                         MPATOWS_sf%>%mutate(startend="end")%>%dplyr::select(trawlID,lat_e,lon_e)%>%rename(lat = lat_e, lon = lon_e))%>%
      mutate(startend=rep(c("start","end"),each=n()/2))%>%
      st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
      group_by(trawlID)%>%
      summarise(geometry=st_union(geometry)%>%st_centroid())%>%
      ungroup()%>%
      mutate(inside=as.logical(st_intersects(.,sab, sparse=TRUE)))

#Load fish morphology data -----------

    #Load the full fish morph data and then subset
    load("data/CrabSurvey/FishMorph.RData")
    
    fishmorph_df <- fishmorphs%>%
                  left_join(.,fishcodes, by="SPECCD_ID")%>%
                  filter(COMM != "",!is.na(COMM),
                         !(is.na(FISH_LENGTH)), 
                         FISH_LENGTH!="")%>%
                  mutate(year=as.numeric(format(as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y")),
                         species = case_when(COMM == "YELLOWTAIL FLOUNDER; LIMANDA" ~ "Yellowtail",
                                          COMM=="YELLOWTAIL FLOUNDER; LIMANDA"~"Yellowtail",
                                          COMM=="TURBOT,GREENLAND HALIBUT"~"Greenland halibut",
                                          COMM=="MONKFISH,GOOSEFISH,ANGLER"~"Anglerfish",
                                          COMM=="HERRING ATLANTIC"~"Atlantic herring",
                                          COMM=="EELPOUT,NEWFOUNDLAND; INCL LYCODES ATLANTICUS; TERRAENOVA"~"Atlantic eelpout",
                                          COMM=="LONGFIN HAKE; UROPHYCIS"~"Longfin hake",
                                          COMM=="SNOW CRAB  QUEEN; INCL CHIONOECETES SP"~"Snow crab",
                                          COMM=="SNAKE BLENNY; LUMPRETAEFORMIS"~"Snake Blenny",
                                          COMM=="SNOWFLAKE HOOKEAR SCULPIN"~"Snowflake Sculpin",
                                          COMM=="HOOKEAR SCULPIN,ATL."~"Hookear sculpin",
                                          COMM=="SHORTTAILED EELPOUT VAHL"~"Lycodes sp.",
                                          COMM=="RIBBED SCULPIN; PINGELI"~"Ribbed sculpin",
                                          COMM=="STRIPED ATLANTIC WOLFFISH"~"Atlantic wolffish",
                                          COMM=="COD ATLANTIC"~"Atlantic cod",
                                          TRUE ~ NA), #these are the subset in the original plot
                        species2 = case_when(COMM!=species~species,
                                             COMM=="DAUBED SHANNY; LUMPENUS"~"Daubed shanny", #other ones to fix
                                             COMM=="WHITE BARRACUDINA; NOTOLEPIS RISSOI"~"White barracudina",
                                             COMM=="LITTLE SKATE; ERINACEA"~"Little skate",
                                             COMM=="SEASNAIL UNIDENTIFIED"~"Seasnail",
                                             COMM=="BRILL/WINDOWPANE"~"Windowpane flounder",
                                          TRUE ~ COMM))%>%
                  rowwise()%>%
                  mutate(species2 = capitalize_first(tolower(species2)))%>% #clean up the ugly all caps species
                  data.frame()
    

    
    bigfish_df <- fishmorph_df%>%
                  mutate(station=as.integer(STATION))%>%
                  filter(station %in% stns$STATION, 
                         EST_NUM_CAUGHT > 2)%>%
                  group_by(species2,station,year)%>%
                  summarise(med=median(FISH_LENGTH,na.rm=T),
                            mn=mean(FISH_LENGTH,na.rm=T),
                            bigfish = quantile(FISH_LENGTH,0.9))%>%
                  ungroup()%>%
                  left_join(.,stns%>%rename(station=STATION))%>%
                  group_by(species2,station)%>%
                  mutate(stand_length = bigfish/max(bigfish),
                         location = ifelse(Inside,"Inside","Outside"))%>%
                  ungroup()%>%
                  data.frame()
    
    yearrange <- range(bigfish_df$year)
    
    plotftr <- bigfish_df%>%
               group_by(species2)%>%
               mutate(inboth = ifelse(length(unique(location))>1,TRUE,FALSE))%>%
               ungroup()%>%
               group_by(species2,station)%>%
               mutate(enoughyears = ifelse(length(unique(year))>2,TRUE,FALSE),
                      prepost = ifelse(length(intersect(unique(year),2015:2016))>0 & length(intersect(unique(year),2017:2023))>0,TRUE,FALSE))%>%
               ungroup()%>%
               data.frame()
    
    
    plot_sp <- plotftr%>%
               distinct(species2,.keep_all=TRUE)%>%
               filter(inboth & enoughyears)%>%
               pull(species2)
    
    plot_df <- plotftr%>%filter(inboth & enoughyears & prepost)
               
                  
                  
    
    ggplot()+
      geom_vline(xintercept = 2017,lty=2)+
      geom_line(data=bigfish_df%>%filter(species2 %in% plot_sp),aes(x=year,y=stand_length,col=location,group=station))+
      geom_point(data=bigfish_df%>%filter(species2 %in% plot_sp),aes(x=year,y=stand_length,col=location,group=station),size=2)+
      theme_bw()+
      facet_wrap(~species,nrow=3)
    
    #Median size tracked over time with smooth
    p_med_a <- ggplot()+
      geom_vline(xintercept = 2017,lty=2)+
      geom_line(data=plot_df,aes(x=year,y=med,col=location,group=station),lwd=0.5,alpha=0.5)+
      geom_point(data=plot_df,aes(x=year,y=med,col=location,group=station),size=2,alpha=0.5)+
      theme_bw()+
      stat_smooth(data=plot_df,aes(x=year,y=med,col=location),lwd=2,se=FALSE)+
      facet_wrap(~species2,nrow=3,scales="free_y")+
      labs(col="",x="Year",y="Median recorded size");p_med_a
    
    p_med_b <- ggplot()+
      geom_vline(xintercept = 2017,lty=2)+
      geom_line(data=plot_df,aes(x=year,y=med,col=location,group=station),lwd=0.5)+
      geom_point(data=plot_df,aes(x=year,y=med,col=location,group=station),size=2)+
      theme_bw()+
      facet_wrap(~species2,nrow=3,scales="free_y")+
      labs(col="",x="Year",y="Median recorded size");p_med_b
    
    ggsave("output/CrabSurvey/MedianLength_inside-outside_wSmooth.png",p_med_a)
    ggsave("output/CrabSurvey/MedianLength_inside-outside.png",p_med_b)
    
    #90th percentile size tracked over time with smooth
    p_bigfish_a <- ggplot()+
      geom_vline(xintercept = 2017,lty=2)+
      geom_line(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),lwd=0.5,alpha=0.5)+
      geom_point(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),size=2,alpha=0.5)+
      theme_bw()+
      stat_smooth(data=plot_df,aes(x=year,y=bigfish,col=location),lwd=2,se=FALSE)+
      facet_wrap(~species2,nrow=3,scales="free_y")+
      labs(col="",x="Year",y="90th percentile size");p_bigfish_a
    
    p_bigfish_b <- ggplot()+
      geom_vline(xintercept = 2017,lty=2)+
      geom_line(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),lwd=0.5)+
      geom_point(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),size=2)+
      theme_bw()+
      facet_wrap(~species2,nrow=3,scales="free_y")+
      labs(col="",x="Year",y="90th percentile size");p_bigfish_b
      
    ggsave("output/CrabSurvey/LargeFish_inside-outside_wSmooth.png",p_bigfish_a)
    ggsave("output/CrabSurvey/LargeFish_inside-outside.png",p_bigfish_b)    
    
#### DIET ANALYSIS --------------
    
      
                  
