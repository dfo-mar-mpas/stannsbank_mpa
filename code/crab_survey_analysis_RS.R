#load libraries
    library(sf)
    library(stringr)
    library(dplyr)
    library(ggplot2)
    library(rnaturalearth)
    library(RODBC)
    library(ggpubr)
    library(tidyr)
    library(vegan)
    library(marmap)
    library(stars)
    library(tibble)
    library(ggrepel)
    library(patchwork)
    library(viridis)
    library(ggnewscale)
    library(ggridges)

#### Global options---------
    sf_use_s2 = FALSE
    
### Target species ---------
    target_sp <- data.frame(spec=c("ANARHICHAS LUPUS", "UROPHYCIS TENUIS", "HIPPOGLOSSOIDES PLATESSOIDES",
                                   "GLYPTOCEPHALUS CYNOGLOSSUS", "CHIONOECETES OPILIO", "GADUS MORHUA", 
                                   "SEBASTES SP.", "AMBLYRAJA RADIATA"),
                            common=c("Atlantic wolffish","White hake","American plaice",
                                     "Witch flounder","Snow crab","Atlantic cod",
                                     "Redfish sp","Thorny skate"))
#### Functions -----
    
    capitalize_first <- function(text) { #from ChatGPT!
      words <- strsplit(text, " ")[[1]] 
      words[[1]] <- paste0(toupper(substr(words[[1]], 1, 1)), substr(words[[1]], 2, nchar(words[[1]])))
      capitalized_text <- paste(words, collapse = " ")  # Recombine the words into a single string
      return(capitalized_text)
    }
    
    reverselog_trans <- function(base = exp(1)) { #https://stackoverflow.com/questions/20924705/plot-negative-values-in-logarithmic-scale-with-ggplot-2
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      trans_new(paste0("reverselog-", format(base)), trans, inv, 
                log_breaks(base = base), 
                domain = c(1e-100, Inf))
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
    
    #benthoscape
    benthoscape_classes <- data.frame(Assigned_c=c("A - Mud",
                                                   "Asp - Mud and seapens",
                                                   "B - Gravelly sand/mud <50% cobblles/gravel",
                                                   "C - Till >50% cobbles/gravel",
                                                   "D - Till with coraline algae",
                                                   "E - Gravel with crinoids",
                                                   "F - Sand with Sand dollars"),
                                      classn=c("Mud", #simplified names
                                               "Mud-Seapens",
                                               "Gravel-mud",
                                               "Till-gravel",
                                               "Till-coraline algae",
                                               "Gravel-crinoids",
                                               "Sand-SandDollars"))
    
    
    sab_benthoscape <- read_sf("data/Shapefiles/benthoscape.shp")%>%
                        st_transform(latlong)%>%
                        st_make_valid()%>%
                        left_join(.,benthoscape_classes)
    
    #load bathymetry
    bbsab <- st_bbox(sab_nozones)
    
    noaabathy <- getNOAA.bathy(bbsab[1]-1,bbsab[3]+1,bbsab[2]+1,bbsab[4]-1,resolution = 0.25,keep=T) %>%
                 fortify.bathy() %>%
                 st_as_stars() %>%
                 st_set_crs(4326)%>%
                 st_transform(latlong)
    

# fish taxonomy metadata -----------
    #fishcodes <- read.csv("data/CrabSurvey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv")
    load("data/CrabSurvey/GROUNDFISH.GSSPECIES.RData")#this is a pull by Mike McMahon that has the missing data items 
    fishcodes <- GSSPECIES%>%rename(SPECCD_ID = CODE); rm(GSSPECIES)# keep the code the same

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

    #standardize the data
    ggplot()+
      geom_boxplot(data=catchdat%>%distinct(TRIP_STATION,.keep_all=TRUE),aes(x=1,y=AREA_SWEPT))
    
    MedianAreaSwept <- catchdat%>%
                       filter(!is.na(AREA_SWEPT))%>%
                       distinct(TRIP_STATION,.keep_all=TRUE)%>%
                       pull(AREA_SWEPT)%>%
                       median()
    
    #standardize by the median area swept
    catchdat_stand <- catchdat%>%
                      filter(!is.na(AREA_SWEPT))%>%
                      mutate(number = EST_NUM_CAUGHT * MedianAreaSwept/AREA_SWEPT,
                             weight = EST_DISCARD_WT * MedianAreaSwept/AREA_SWEPT,
                             COMM = ifelse(COMM == "",NA,COMM),
                             SPEC = ifelse(SPEC =="",NA,SPEC),
                             species_filter = ifelse(is.na(SPEC),COMM,SPEC),
                             station = as.integer(STATION))%>% #get rid of unidentified species
                      filter(!is.na(species_filter))%>%
                      left_join(.,stns%>%rename(station=STATION))%>%
                      mutate(location = ifelse(Inside,"Inside","Outside"),
                             type = ifelse(Enhanced,"Enhanced","Standard"),
                             year=as.numeric(format(as.Date(BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),"%Y")))
 
    
    catchdat_sf <- catchdat_stand%>%
                   distinct(TRIP_STATION,.keep_all=TRUE)%>%
                   mutate(LONGITUDE = LONGITUDE *-1)%>%
                   st_as_sf(coords=c("LONGITUDE","LATITUDE"),crs=latlong)
                  
    ggplot()+geom_sf(data=mpas)+geom_sf(data=catchdat_sf) #all up by the MPA
    
#Save output for taxonomic cleaning -- run Crab_survey_taxonomy.R to get the full cleaned list
    # write.csv(catchdat_stand%>%
    #             select(COMM,SPEC,species_filter)%>%
    #             distinct(species_filter,.keep_all=TRUE),"output/taxonomic_raw.csv",row.names=FALSE)
    
#load the cleaned taxonomy
    PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    
    crab_tax_clean <- read.csv("output/crab_taxa_clean.csv")%>%
                      select(all_of(PhyloNames),aphiaID,common,species_filter)%>%
                      rowwise()%>%
                      mutate(taxcount = sum(c(is.na(Kingdom),is.na(Phylum),is.na(Class),
                                            is.na(Order),is.na(Family),is.na(Genus),is.na(Species))))%>%
                      data.frame()

    #get the minimum identification    
    crab_tax_clean$minID <- crab_tax_clean$Species
                      
    for(i in 1:nrow(crab_tax_clean)){
      
      if(crab_tax_clean[i,"taxcount"]!=0){crab_tax_clean[i,"minID"] <- crab_tax_clean[i,length(PhyloNames)-crab_tax_clean[i,"taxcount"]]}
      
    }
    
#link in the taxonomy data
    catchdat_stand <- catchdat_stand%>%left_join(.,crab_tax_clean)

#Format the crab data ---------
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
      mutate(inside=as.logical(st_intersects(.,sab, sparse=FALSE))) 

##### Fish morphology analysis -----------

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
               mutate(species2 = ifelse(species2=="Redfish unseparated","Redfish sp",species2))%>% #to match target sp
               data.frame()
    
    
    plot_df <- plotftr%>%
               filter(inboth & enoughyears & prepost,
                      species2 %in% target_sp$common)
               
     
    
    #Median size tracked over time with smooth
    p_med_a <- ggplot()+
               geom_vline(xintercept = 2017,lty=2)+
               geom_line(data=plot_df,aes(x=year,y=med,col=location,group=station),lty=2,lwd=0.25,alpha=0.5)+
               geom_point(data=plot_df,aes(x=year,y=med,col=location,group=station),size=2,alpha=0.5)+
               theme_bw()+
               theme(strip.background = element_rect(fill="white"),
                    legend.position = "bottom")+
               stat_smooth(data=plot_df,aes(x=year,y=med,col=location),lwd=2,se=FALSE)+
               facet_wrap(~species2,ncol=2,scales="free_y")+
               labs(col="",x="",y="Median recorded size")
    
    p_med_b <- ggplot()+
               geom_vline(xintercept = 2017,lty=2)+
               geom_line(data=plot_df,aes(x=year,y=med,col=location,group=station),lwd=0.5)+
               geom_point(data=plot_df,aes(x=year,y=med,col=location,group=station),size=2)+
               theme_bw()+
               theme(strip.background = element_rect(fill="white"),
                    legend.position = "bottom")+
               facet_wrap(~species2,nrow=3,scales="free_y")+
               labs(col="",x="",y="Median recorded size")
    
    ggsave("output/CrabSurvey/MedianLength_inside-outside_wSmooth.png",p_med_a,height=8,width=6,units="in",dpi=300)
    ggsave("output/CrabSurvey/MedianLength_inside-outside.png",p_med_b,height=8,width=6,units="in",dpi=300)
    
    #90th percentile size tracked over time with smooth
    p_bigfish_a <- ggplot()+
                   geom_vline(xintercept = 2017,lty=2)+
                   geom_line(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),lty=2,lwd=0.25,alpha=0.5)+
                   geom_point(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),size=2,alpha=0.55)+
                   theme_bw()+
                   theme(strip.background = element_rect(fill="white"),
                        legend.position = "bottom")+
                   stat_smooth(data=plot_df,aes(x=year,y=bigfish,col=location),lwd=2,se=FALSE)+
                   facet_wrap(~species2,ncol=2,scales="free_y")+
                   labs(col="",x="",y="90th percentile size")
    
    p_bigfish_b <- ggplot()+
                   geom_vline(xintercept = 2017,lty=2)+
                   geom_line(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),lwd=0.5)+
                   geom_point(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),size=2)+
                   theme_bw()+
                   theme(strip.background = element_rect(fill="white"),
                        legend.position = "bottom")+
                   facet_wrap(~species2,ncol=2,scales="free_y")+
                   labs(col="",x="",y="90th percentile size")
      
    ggsave("output/CrabSurvey/LargeFish_inside-outside_wSmooth.png",p_bigfish_a,height=8,width=6,units="in",dpi=300)
    ggsave("output/CrabSurvey/LargeFish_inside-outside.png",p_bigfish_b,height=8,width=6,units="in",dpi=300)    
    
    
##### Catch changes within stations --------------
    
    ##biodiversity per set
    biodiv_df <- catchdat_stand%>%
                 filter(type=="Enhanced")%>%
                 group_by(TRIP_STATION)%>%
                 summarise(richness=length(unique(species_filter)))%>%
                 ungroup()%>%
                 left_join(.,catchdat_stand%>%
                             distinct(TRIP_STATION,.keep_all=TRUE)%>%
                             select(TRIP_STATION,location,type,station,year))%>%
                 select(location,type,TRIP_STATION,station,richness,year)%>%
                 group_by(station)%>%
                 mutate(richness_stand = richness/max(richness))%>%
                 ungroup()%>%
                 data.frame()
                
    p_richness_a <- ggplot()+
                    geom_vline(xintercept = 2017,lty=2)+
                    geom_line(data=biodiv_df,aes(x=year,y=richness,group=station,col=location),lty=2,lwd=0.25)+
                    geom_point(data=biodiv_df,aes(x=year,y=richness,group=station,col=location),size=2,alpha=0.5)+
                    theme_bw()+
                    theme(strip.background = element_rect(fill="white"),
                          legend.position = "bottom")+
                    labs(col="",x="",y="Richness per station")+
                    theme(legend.position = "bottom")
    
    p_richness_b <- ggplot()+
                    geom_vline(xintercept = 2017,lty=2)+
                    geom_line(data=biodiv_df,aes(x=year,y=richness,group=station,col=location),lty=2,lwd=0.25,alpha=0.5)+
                    geom_point(data=biodiv_df,aes(x=year,y=richness,group=station,col=location),size=2,alpha=0.5)+
                    stat_smooth(data=biodiv_df,aes(x=year,y=richness,col=location),lwd=2)+
                    theme_bw()+
                    theme(strip.background = element_rect(fill="white"),
                          legend.position = "bottom")+
                    labs(col="",x="",y="Richness per station")+
                    theme(legend.position = "bottom")
    
    ggsave("output/CrabSurvey/Biodiversity_inside-outside.png",p_richness_a,height=8,width=6,units="in",dpi=300)
    ggsave("output/CrabSurvey/Biodiversity_inside-outside_wSmooth.png",p_richness_b,height=8,width=6,units="in",dpi=300)
    
    ##catch number per set
    
    num_df <- catchdat_stand%>%
              mutate(spec = ifelse(SPEC == "SEBASTES","SEBASTES SP.",SPEC))%>%
              filter(spec %in% target_sp$spec,
                     station %in% stns$STATION)%>%
              mutate(station=as.integer(STATION))%>%
              group_by(location,station,minID)%>%
              mutate(num_stand = number/max(number,na.rm=T))%>%
              ungroup()%>%
              data.frame()
      
    p_number_a <- ggplot()+
                  geom_vline(xintercept = 2017,lty=2)+
                  geom_line(data=num_df,aes(x=year,y=number,col=location,group=station),lty=2,lwd=0.25,alpha=0.5)+
                  geom_point(data=num_df,aes(x=year,y=number,col=location,group=station),size=2,alpha=0.55)+
                  theme_bw()+
                  theme(strip.background = element_rect(fill="white"),
                        legend.position = "bottom")+
                  stat_smooth(data=num_df,aes(x=year,y=number,col=location),lwd=2,se=FALSE)+
                  scale_y_log10()+
                  facet_wrap(~minID,ncol=2,scales="free_y")+
                  labs(col="",x="",y="Count per set")
                
    p_number_b <- ggplot()+
                  geom_vline(xintercept = 2017,lty=2)+
                  geom_line(data=num_df,aes(x=year,y=num_stand,col=location,group=station),lty=2,lwd=0.25,alpha=0.5)+
                  geom_point(data=num_df,aes(x=year,y=num_stand,col=location,group=station),size=2,alpha=0.55)+
                  theme_bw()+
                  theme(strip.background = element_rect(fill="white"),
                        legend.position = "bottom")+
                  stat_smooth(data=num_df,aes(x=year,y=num_stand,col=location),lwd=2)+
                  facet_wrap(~minID,ncol=2)+
                  labs(col="",x="",y=expression(paste("Count per set / Max count per set")))
                
    ggsave("output/CrabSurvey/CatchNumber_inside-outside.png",p_number_a,height=8,width=6,units="in",dpi=300)
    ggsave("output/CrabSurvey/CatchNumber_inside-outside_standardized.png",p_number_b,height=8,width=6,units="in",dpi=300)
    
    ##catch weight per set
    
    weight_sp <- catchdat_stand%>%
                 filter(station %in% stns$STATION,
                        !is.na(minID))%>%
                 group_by(minID)%>%
                 summarise(obs_cnt = n())%>%
                 ungroup()%>%
                 arrange(-obs_cnt)%>%
                 data.frame()%>%
                 slice(1:20)%>%
                 pull(minID)%>%
                 c(.,"Anarhichas lupus")
    
    weight_df <- catchdat_stand%>%
                 filter(minID %in% weight_sp)%>%
                 mutate(station=as.integer(STATION))%>%
                 group_by(location,station,minID)%>%
                 mutate(weight_stand = weight/max(weight,na.rm=T))%>%
                 ungroup()%>%
                 data.frame()
    
    p_weight_a <- ggplot()+
                  geom_vline(xintercept = 2017,lty=2)+
                  geom_line(data=weight_df,aes(x=year,y=weight,col=location,group=station),lty=2,lwd=0.25,alpha=0.5)+
                  geom_point(data=weight_df,aes(x=year,y=weight,col=location,group=station),size=2,alpha=0.55)+
                  theme_bw()+
                  stat_smooth(data=weight_df,aes(x=year,y=weight,col=location),lwd=2,se=FALSE)+
                  scale_y_log10()+
                  facet_wrap(~minID,nrow=3,scales="free_y")+
                  labs(col="",x="Year",y="Biomass per set (kg)");p_weight_a
    
    p_weight_b <- ggplot()+
                  geom_vline(xintercept = 2017,lty=2)+
                  geom_line(data=weight_df,aes(x=year,y=weight_stand,col=location,group=station),lty=2,lwd=0.25,alpha=0.5)+
                  geom_point(data=weight_df,aes(x=year,y=weight_stand,col=location,group=station),size=2,alpha=0.55)+
                  theme_bw()+
                  stat_smooth(data=weight_df,aes(x=year,y=weight_stand,col=location),lwd=2,se=FALSE)+
                  facet_wrap(~minID,ncol=6,scales="free_y")+
                  labs(col="",x="Year",y=expression(paste("Biomass per set / Max count per set")));p_weight_b
    
    ggsave("output/CrabSurvey/CatchNumber_inside-outside.png",p_weight_a)
    ggsave("output/CrabSurvey/CatchNumber_inside-outside_standardized.png",p_weight_b)
    
##### Diet analysis --------------
    
    diet_data <- read.csv("data/CrabSurvey/MPA.Diet.SnowCrabSurvey.Feb.2024.csv")%>%
                 filter(SLATDD>45) #filter to just SAB
    
    #There are some non-animal based things in the diet that can be flagged
    problem_things <- c(1099,2507,9500,9100,9000,7000,9001,3199,9700,4310,1180,4501,9620,9003,4223,2010,9002,2901,3031,9600)
    non_animal_prey <- c(9100,9000,9700,9600,9620)
    
    diet_df <- diet_data%>%
               rename(SPECCD_ID = PREYSPECCD,
                      PRED = SPEC)%>%
               left_join(.,fishcodes%>%select(SPEC,COMM, SPECCD_ID))%>%
               rename(prey = SPEC,
                      prey_comm = COMM,
                      prey_speccd_id = SPECCD_ID,
                      SPECCD_ID = PRED)%>%
               left_join(.,fishcodes%>%select(SPEC,COMM, SPECCD_ID))%>%
               rename(pred = SPEC,
                      pred_comm = COMM,
                      pred_speccd_id = SPECCD_ID)%>%
              rename(TRIP = MISSION)%>%
              mutate(TRIP_SET = paste(TRIP,SETNO,sep="_"))%>%
              left_join(.,arsw%>%mutate(TRIP_SET = paste(TRIP,SET,sep="_"))%>%select(-TRIP),by="TRIP_SET")%>%
              st_as_sf(coords=c("SLONGDD","SLATDD"),crs=latlong,remove=FALSE)%>%
              mutate(inside=as.logical(st_intersects(.,sab_nozones, sparse=FALSE)),
                     location=ifelse(inside,"Inside","Outside"),
                     prey=ifelse(prey %in% c("TEUTHOIDEA O."),"LOLIGINIDAE,OMMASTREPHIDAE F.",prey), #fix a double count for squid
                     prey_comm=ifelse(prey_speccd_id == 4501,"SQUID (NS)",prey_comm),
                     prey_speccd_id=ifelse(prey_speccd_id == 4501,4514,prey_speccd_id),
                     non_animal_prey=prey_speccd_id %in% non_animal_prey,
                     non_species=prey_speccd_id %in% problem_things,
                     non_species=ifelse(is.na(prey_speccd_id),TRUE,non_species))%>%
               data.frame()%>%
               select(-geometry)
    
        ## diet richness ~ location x year
    
    diet_richness <- diet_df%>%
                     group_by(pred,Year,location)%>%
                     summarise(rich = length(unique(prey_speccd_id)))%>%
                     ungroup()%>%
                     group_by(pred)%>%
                     mutate(rich_stand = rich/max(rich))%>% #this standardizes the richness so that average species trends can be shown
                     ungroup()%>%
                     data.frame()
    
    p_diet_richness_stand <- ggplot()+
                             geom_vline(xintercept = 2017,lty=2)+
                             geom_point(data=diet_richness,aes(x=Year,y=rich_stand,col=location,group=pred))+
                             stat_smooth(data=diet_richness,aes(x=Year,y=rich_stand,col=location))+
                             theme_bw()+
                             labs(x="",y="Standardized diet compositional richness",col="");p_diet_richness_stand
    
    diet_rich <- diet_df%>%
                 filter(pred %in% target_sp$spec,
                        !non_species)%>% #take out the non-species elements
                 group_by(pred,Year,location,TRIP_STATION)%>%
                 summarise(rich = length(unique(prey_speccd_id)))%>%
                 ungroup()%>%
                 group_by(pred,Year,location)%>%
                 summarise(mean_rich=mean(rich,na.rm=T),
                           sd_rich=sd(rich,na.rm=T))%>%
                 ungroup()%>%
                 data.frame()%>%
                 left_join(.,target_sp%>%rename(pred=spec))
              
    p_diet_richness <- ggplot(data=diet_rich,aes(x=Year,y=mean_rich,col=location))+
                        geom_vline(xintercept = 2017,lty=2)+
                        geom_line(lwd=0.25,alpha=0.25)+
                        geom_linerange(aes(ymin=mean_rich-sd_rich,ymax=mean_rich+sd_rich))+
                        geom_point(size=2, position = position_jitter(width = 0.1))+
                        facet_wrap(~common,ncol=2,scales="free_y")+
                        theme_bw()+
                        theme(strip.background = element_rect(fill="white"),
                              legend.position = "bottom")+
                        labs(x="",y=bquote(bar(x) ~ " diet richness ± sd"),col="");p_diet_richness
    
    ggsave("output/CrabSurvey/diet_richness.png",p_diet_richness,height=8,width=6,units="in",dpi=300)  
      
    #diet fullness
    
    diet_fullness <- diet_df%>%
                      filter(pred %in% target_sp$spec)%>%
                      group_by(pred,Year,location,TRIP_STATION)%>%
                      summarise(full = mean(FULLNESS,na.rm=T)/4)%>%
                      ungroup()%>%
                      group_by(pred,Year,location)%>%
                      summarise(mean_full=mean(full,na.rm=T),
                                sd_full=sd(full,na.rm=T))%>%
                      ungroup()%>%
                      data.frame()%>%
                      left_join(.,target_sp%>%rename(pred=spec))
    
    p_diet_fullness <- ggplot(data=diet_fullness,aes(x=Year,y=mean_full,col=location))+
                       geom_vline(xintercept = 2017,lty=2)+
                       geom_line(lwd=0.25,alpha=0.25)+
                       geom_linerange(aes(ymin=mean_full-sd_full,ymax=mean_full+sd_full))+
                       geom_point(size=2, position = position_jitter(width = 0.1))+
                       facet_wrap(~common,ncol=2)+
                       theme_bw()+
                       theme(strip.background = element_rect(fill="white"),
                             legend.position = "bottom")+
                       scale_y_continuous(labels=scales::percent_format(),breaks=c(0,0.5,1))+
                       labs(x="",y=bquote(bar(x) ~ " diet fullness ± sd"),col="");p_diet_fullness
    
    ggsave("output/CrabSurvey/diet_fullness.png",p_diet_fullness,height=8,width=6,units="in",dpi=300)    
                  
    
##### NMDS analysis (diet) ------------------
    
    diet_nmds_df <- diet_df%>%
                    filter(!non_species,
                          pred %in% target_sp$spec)%>%
                    group_by(pred,location,Year,prey)%>%
                    summarise(prey_weight=mean(PWT,na.rm=T))%>%
                    ungroup()
    
    diet_nmds_list <- list()
    
    #dummy dataframes to hold the outputs -- can be loaded after first run with load("output/CrabSurvey/diet_nmds_analysis.RData")
    nmds_list <- list()
    nmds_pa_list <- list()
    nmds_species_list <- list()
    nmds_data <- NULL
    
    #run the nmds analyses per species. 
    for(i in 1:length(unique(diet_nmds_df$pred))){
      
      sp <- unique(diet_nmds_df$pred)[i] #species of the nmds
      
      clusters <- ifelse(sp == "ANARHICHAS LUPUS",2,ifelse(sp=="SEBASTES SP.",3,5)) # note that there are few data points for wolffish and Sebastes
      
      message(paste0("Working on ",sp," ",i, " of ",length(unique(diet_nmds_df$pred)))) #progress message
      
      #create the nmds sample by species data.frame - standardize
      
      nmds_df <- diet_nmds_df%>%
                 filter(pred==sp)%>%
                 mutate(id = paste(location,Year,sep="_"))
      
      nmds_temp <- diet_nmds_df%>% #if wolffish it returns an issue
                   filter(pred==sp)%>%
                   mutate(id = paste(location,Year,sep="_"))%>%
                   select(id,prey,prey_weight)%>%
                   spread(prey,prey_weight,fill=0)%>%
                   column_to_rownames('id')%>%
                   mutate(tot = rowSums(.,na.rm=T))%>%
                   filter(tot>0)%>% #filter out any zero catches
                   select(-tot)%>%
                   decostand(method = "total")%>%
                   metaMDS(.,k=clusters)
      
      nmds_temp_pa <- diet_nmds_df%>% #presence-absence - diet composition 
                      filter(pred==sp)%>%
                      mutate(id = paste(location,Year,sep="_"))%>% #this converts to binary - technically the condition isn't needed since there is a value for each, but keep it in there in case
                      select(id,prey,prey_weight)%>%
                      spread(prey,prey_weight,fill=0)%>%
                      column_to_rownames('id')%>%
                      mutate(tot = rowSums(.,na.rm=T))%>%
                      filter(tot>0)%>% #filter out any zero catches
                      select(-tot)%>%
                      decostand(method = "total")%>%
                      dist(., method="binary")%>%
                      metaMDS(.,k=clusters)
        
      #data processing
      data.scores <- as.data.frame(scores(nmds_temp,"sites"))%>%
                     mutate(id=rownames(.),
                            method="bray",
                            stress=round(nmds_temp$stress,3))%>%
                     separate(id,c("location","year"),sep="_")%>%
                     rbind(.,
                           as.data.frame(scores(nmds_temp_pa,"sites"))%>%
                             mutate(id=rownames(.),
                                    method="pa",
                                    stress=round(nmds_temp_pa$stress,3))%>%
                             separate(id,c("location","year"),sep="_"))%>%
                    select(NMDS1,NMDS2,location,year,method,stress)%>%
                    mutate(year=as.numeric(year),pred=sp)
      
      species.scores <- as.data.frame(scores(nmds_temp,"species"))%>%
                        mutate(method="bray")%>%
                        rbind(.,
                              as.data.frame(scores(nmds_temp_pa,"species"))%>%
                                mutate(method="pa"))%>%
                        mutate(pred=sp)
      
      #looped outputs
      nmds_data <- rbind(nmds_data,data.scores)
      nmds_species_list[[sp]] <- species.scores
      nmds_list[[sp]] <- nmds_temp
      nmds_temp_pa[[sp]] <- nmds_temp_pa
      
      rm(data.scores,species.scores,nmds_temp,nmds_temp_pa,nmds_df) #clean up 
      
    }
    
    save(nmds_data,nmds_species_list,nmds_list,file="output/CrabSurvey/diet_nmds_analysis.RData")# interim save. 
  
  #Assemble the convex hull data for plotting    
    hull_data <- nmds_data%>%
                 group_by(pred,location,method)%>%
                 slice(chull(x=NMDS1,y=NMDS2))%>%
                 rename(spec=pred)%>%
                 left_join(.,target_sp)%>%
                 filter(spec !="ANARHICHAS LUPUS") #wolffish do not have enough 'outside' recoreds to make a hull

    #Convex hull plots of diet
    nmds_plot_bray <- ggplot(data=hull_data%>%filter(method=="bray"))+
                      geom_polygon(aes(x=NMDS1,y=NMDS2,fill=location),alpha=0.30)+
                      labs(fill="")+
                      new_scale_fill()+
                      geom_point(aes(x=NMDS1,y=NMDS2,shape=location,fill=year),size=3)+
                      scale_shape_manual(values=c(21,23),guide="none")+
                      theme_bw()+
                      theme(axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            legend.position = "bottom",
                            strip.background = element_rect(fill="white"))+
                      facet_wrap(~common,ncol=2,scales="free")+
                      labs(shape="",fill="",col="Sample year")+
                      scale_fill_viridis()
    
    nmds_plot_pa <- ggplot(data=hull_data%>%filter(method=="pa"))+
                    geom_polygon(aes(x=NMDS1,y=NMDS2,fill=location),alpha=0.30)+
                    labs(fill="")+
                    new_scale_fill()+
                    geom_point(aes(x=NMDS1,y=NMDS2,shape=location,fill=year),size=3)+
                    scale_shape_manual(values=c(21,23),guide="none")+
                    theme_bw()+
                    theme(axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          legend.position = "bottom",
                          strip.background = element_rect(fill="white"))+
                    facet_wrap(~common,ncol=2,scales="free")+
                    labs(shape="",fill="",col="Sample year")+
                    scale_fill_viridis()
    
    #save plots
    ggsave("output/CrabSurvey/diet_hull_bray.png",nmds_plot_bray,height=6,width=6,units="in",dpi=300)
    ggsave("output/CrabSurvey/diet_hull_pa.png",nmds_plot_pa,height=6,width=6,units="in",dpi=300)
    
     ##Assemble yearly progression plots --- 
     plot_sp <- unique(hull.data$common)%>%sort()
     
     for(i in c("bray","pa")){
       
       for(j in target_sp$common){
         
         temp_plotdf <- hull_data%>%
                        filter(common==j,
                               method==i)%>%
                        arrange(location,year)%>%
                        data.frame()
         
         #set up common range limits for the axes
         
         range_offset <- 0.05 # this is the % buffer to add to the plot
         
         nmds1_range <- range(temp_plotdf$NMDS1)
         nmds1_offset <- (nmds1_range[2]-abs(nmds1_range[1]))*range_offset
         nmds1_range[1] <- nmds1_range[1] - nmds1_offset
         nmds1_range[2] <- nmds1_range[2] + nmds1_offset
         
         nmds2_range <- range(temp_plotdf$NMDS2)
         nmds2_offset <- (nmds2_range[2]-abs(nmds2_range[1]))*range_offset
         nmds2_range[1] <- nmds2_range[1] - nmds2_offset
         nmds2_range[2] <- nmds2_range[2] + nmds2_offset
     
    
    p1 <- ggplot(data=temp_plotdf%>%filter(location=="Inside"),aes(x=NMDS1,y=NMDS2,group=location))+
          geom_segment(aes(xend=c(tail(NMDS1, n=-1), NA),yend=c(tail(NMDS2, n=-1), NA),group=location),
          arrow=arrow(length=unit(0.25,"cm"),type="closed"),lty=3)+
          geom_point(aes(fill=year),size=3,shape=21)+
          geom_point(aes(fill=year),data=temp_plotdf%>%filter(location=="Inside",year==max(year)),size=5,shape=21)+
          geom_text_repel(aes(label=year))+
          scale_x_continuous(limits=nmds1_range)+
          scale_y_continuous(limits=nmds2_range)+
          theme_bw()+
          facet_grid(common~location)+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                strip.background.y = element_blank(),
                strip.text.y = element_blank(),
                axis.title = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text=element_text(size=18),
                strip.background.x = element_rect(fill="white"),
                legend.position = "none")+
          labs(fill="")
    
    p2 <- ggplot(data=temp_plotdf%>%filter(location=="Outside"),aes(x=NMDS1,y=NMDS2,group=location))+
          geom_segment(aes(xend=c(tail(NMDS1, n=-1), NA),yend=c(tail(NMDS2, n=-1), NA),group=location),
                       arrow=arrow(length=unit(0.25,"cm"),type="closed"),lty=3)+
          geom_point(aes(fill=year),size=3,shape=21)+
          geom_point(aes(fill=year),data=temp_plotdf%>%filter(location=="Outside",year==max(year)),size=5,shape=21)+
          geom_text_repel(aes(label=year))+
          scale_x_continuous(limits=nmds1_range)+
          scale_y_continuous(limits=nmds2_range)+
          theme_bw()+
          facet_grid(common~location)+
          theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.title = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text=element_text(size=18),
                strip.background = element_rect(fill="white"),
                legend.position = "none")+
          labs(fill="")
      
      
      #some hard coding based on on a set diet species list based on the focal species for the power analysis. 
      if(j != "American plaice") {p3 <- p1+theme(strip.background.x = element_blank(),strip.text.x = element_blank()) + 
                                        p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()) + plot_layout(ncol=2)}
    
      if(j == "American plaice") {p3 <- p1 + p2 + plot_layout(ncol=2)}
      
      assign(paste0(gsub(" ","_",j),"_plot"),p3)
    
       } #end j loop (species)
       
       group_plot <- American_plaice_plot / Atlantic_cod_plot / Redfish_sp_plot / Thorny_skate_plot / White_hake_plot / Witch_flounder_plot 
       
       ggsave(paste0("output/CrabSurvey/diet_year_",i,".png"),group_plot,height=12,width=7,units="in",dpi=300)
       
       rm(group_plot,American_plaice_plot,Atlantic_cod_plot,Redfish_sp_plot,Thorny_skate_plot,White_hake_plot,Witch_flounder_plot) #clean up ws
       
     } #end i loop (method)
    
    

### diet ~ benthoscape ----
     
     
     #get the environmental covariates for each station. 
      diet_nmds_env <- diet_df%>%
                       filter(!non_species)%>%
                       st_as_sf(coords=c("SLONGDD","SLATDD"),crs=latlong,remove=FALSE)%>%
                       group_by(STATION)%>%
                       summarise(geometry=st_union(geometry)%>%st_centroid())%>% #get a consensus station location
                       ungroup()%>%
                       left_join(.,diet_df%>%distinct(STATION,.keep_all = TRUE)%>%select(location,TRIP_SET,TRIP_STATION,STATION))%>%
                       st_intersection(.,sab_benthoscape%>%dplyr::select(classn,geometry))%>%
                       rbind(.,diet_df%>% #coordinate in the top right is just adjacent the mud-seapens shelf break. TO get replication we will assign this one manually. 
                               filter(!non_species)%>%
                               st_as_sf(coords=c("SLONGDD","SLATDD"),crs=latlong,remove=FALSE)%>%
                               group_by(STATION)%>%
                               summarise(geometry=st_union(geometry)%>%st_centroid())%>% #get a consensus station location
                               ungroup()%>%
                               left_join(.,diet_df%>%distinct(STATION,.keep_all = TRUE)%>%
                                           select(location,TRIP_SET,TRIP_STATION,STATION))%>%
                               mutate(lat=st_coordinates(.)[,2])%>%
                               arrange(-lat)%>%
                               slice(2)%>%
                               mutate(classn="Mud-Seapens")%>%
                               select(STATION,location,TRIP_SET,TRIP_STATION,classn,geometry))%>%
                       select(STATION,location,classn,geometry)%>%
                       st_join(noaabathy %>% st_as_sf())%>% # add in the bathymetry 
                       data.frame()%>%
                       select(-geometry)
                
     diet_benthoscape <- diet_df%>%
                         filter(!non_species,
                                pred %in% target_sp$spec,
                                STATION %in% diet_nmds_env$STATION)%>%
                         left_join(.,diet_nmds_env%>%select(-location))%>%
                         left_join(.,target_sp%>%rename(pred=spec))%>%
                         group_by(common,classn)%>%
                         summarise(diet_rich = length(unique(prey)))%>%
                         ungroup()%>%
                         group_by(common)%>%
                         mutate(diet_rich_stand = diet_rich/max(diet_rich))%>%
                         ungroup()%>%
                         data.frame()
                      
     #focal species analysis
     
         #Simple boxplot
         # ggplot(data=diet_benthoscape,aes(x=classn,y=diet_rich))+
         #    geom_jitter()+
         #    geom_boxplot(alpha=0.5)+
         #    theme_bw()+
         #    labs(x="",y="Diet richness")
         # 
     
     bentho_richness <- ggplot(data=diet_benthoscape,aes(x=common,y=diet_rich,fill=classn))+
                       geom_bar(position = "dodge", stat = "identity",col="black")+
                       theme_bw()+
                       facet_wrap(~common,nrow=1,scales="free_x")+
                       scale_y_continuous(expand= expansion(mult = c(0, .05)))+
                       scale_fill_viridis(discrete=TRUE)+
                       theme(legend.position = "bottom",
                             axis.text.x=element_blank(),
                             strip.background = element_rect(fill="white"))+
                       labs(x="",y="Diet richness",fill="")
     
     ggsave("output/CrabSurvey/benthoscape_diet_richness.png",bentho_richness,width=8,height=5,units="in",dpi=300)
     
     #comparison of all species by benthoscape class
     diet_benthoscape_all <- diet_df%>%
                           filter(!non_species,
                                  #pred %in% target_sp$spec,
                                  STATION %in% diet_nmds_env$STATION)%>%
                           left_join(.,diet_nmds_env%>%select(-location))%>%
       
                           group_by(pred,classn)%>%
                           summarise(diet_rich = length(unique(prey)),
                                     pred_size = mean(FLEN,na.rm=T),
                                     pred_mass = mean(FWT,na.rm=T),
                                     pred_full = mean(FULLNESS,na.rm=T))%>%
                           ungroup()%>%
                           
                           group_by(pred)%>%
                           mutate(diet_rich_stand = diet_rich/max(diet_rich))%>%
                           ungroup()%>%
                           
                           group_by(classn)%>%
                           mutate(num_pred = length(unique(pred)))%>%
                           ungroup()%>%
                           data.frame()
     
     bentho_richness2 <- ggplot(data=diet_benthoscape_all,aes(x=classn,y=diet_rich,fill=num_pred))+
                         geom_point(position = position_jitter(),aes(size=pred_mass),alpha=0.5,shape=21,fill="grey30",col="black",show.legend = FALSE)+
                         geom_boxplot(alpha=0.5)+
                         theme_bw()+
                         labs(x="",y="Diet richness",fill="Predator richness")+
                         scale_fill_viridis()+
                         theme(legend.position = "bottom")
     
     ggsave("output/CrabSurvey/benthoscape_diet_richness_all.png",bentho_richness2,width=5,height=5,units="in",dpi=300)
     
     
##### benthoscape nmds -----
     
     benthodiet_nmds_bray <- diet_df%>%
                           filter(!non_species,
                                  STATION %in% diet_nmds_env$STATION)%>%
                           left_join(.,diet_nmds_env%>%select(-location))%>%
                           group_by(pred,classn,prey)%>%
                           summarise(prey_weight=mean(PWT,na.rm=T))%>%
                           ungroup()%>%
                           mutate(id = paste(pred,classn,sep="_"))%>%
                           select(id,prey,prey_weight)%>%
                           spread(prey,prey_weight,fill=0)%>%
                           column_to_rownames('id')%>%
                           mutate(tot = rowSums(.,na.rm=T))%>%
                           filter(tot>0)%>% #filter out any zero catches
                           select(-tot)%>%
                           decostand(method = "total")%>%
                           metaMDS(.,k=5)
     
     benthodiet_nmds_pa <- diet_df%>%
                           filter(!non_species,
                                  STATION %in% diet_nmds_env$STATION)%>%
                           left_join(.,diet_nmds_env%>%select(-location))%>%
                           group_by(pred,classn,prey)%>%
                           summarise(prey_weight=mean(PWT,na.rm=T))%>%
                           ungroup()%>%
                           mutate(id = paste(pred,classn,sep="_"))%>%
                           select(id,prey,prey_weight)%>%
                           spread(prey,prey_weight,fill=0)%>%
                           column_to_rownames('id')%>%
                           mutate(tot = rowSums(.,na.rm=T))%>%
                           filter(tot>0)%>% #filter out any zero catches
                           select(-tot)%>%
                           decostand(method = "total")%>%
                           dist(., method="binary")%>%
                           metaMDS(.,k=clusters)
     
     bentho_data_scores <- as.data.frame(scores(benthodiet_nmds_bray,"sites"))%>%
                           mutate(id=rownames(.),
                                  method="bray",
                                  stress=round(benthodiet_nmds_bray$stress,3))%>%
                           separate(id,c("pred","classn"),sep="_")%>%
                           rbind(.,
                                 as.data.frame(scores(benthodiet_nmds_pa,"sites"))%>%
                                   mutate(id=rownames(.),
                                          method="pa",
                                          stress=round(benthodiet_nmds_pa$stress,3))%>%
                                   separate(id,c("pred","classn"),sep="_"))%>%
                           select(NMDS1,NMDS2,pred,classn,method,stress)
     
     hull_data_bentho <- bentho_data_scores%>%
                         group_by(classn,method)%>%
                         slice(chull(x=NMDS1,y=NMDS2))
     
     
     
     #Convex hull plots of diet
      p1_bray <- ggplot(data=hull_data_bentho%>%filter(method=="bray"))+
                 geom_polygon(aes(x=NMDS1,y=NMDS2,fill=classn),alpha=0.30,col="black",lwd=0.1)+
                 geom_point(aes(x=NMDS1,y=NMDS2,fill=classn),shape=21,size=2)+
                 theme_bw()+
                 theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       legend.position = "none",
                       strip.background = element_rect(fill="white"))+
                 labs(shape="",fill="",col="Sample year",title="Bray-Curtis")+
                 geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(benthodiet_nmds_bray$stress,3),"k =",benthodiet_nmds_bray$ndim)))+
                 scale_fill_viridis(discrete=TRUE)
      
      p1_pa <- ggplot(data=hull_data_bentho%>%filter(method=="pa"))+
                geom_polygon(aes(x=NMDS1,y=NMDS2,fill=classn),alpha=0.30,col="black",lwd=0.1)+
                geom_point(aes(x=NMDS1,y=NMDS2,fill=classn),shape=21,size=2)+
                theme_bw()+
                theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.position = "right",
                      strip.background = element_rect(fill="white"))+
                labs(shape="",fill="",col="Sample year",title="Euclidean")+
                geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(benthodiet_nmds_pa$stress,3),"k =",benthodiet_nmds_pa$ndim)))+
                scale_fill_viridis(discrete=TRUE)
      
      #combine using patchwork
      p1_combo <- p1_bray + p1_pa + plot_layout(ncol=2)
                     
      #save plot
      ggsave("output/CrabSurvey/Benthoscape_nmds_comp.png",p1_combo,width=7,height=5,units="in",dpi=300)
      
      
      #now pool the species and repliation by  year
      benthodiet_nmds_bray_all <- diet_df%>%
                                  filter(!non_species,
                                         STATION %in% diet_nmds_env$STATION)%>%
                                  left_join(.,diet_nmds_env%>%select(-location))%>%
            
                                  group_by(Year,classn,pred,prey)%>%
                                  summarise(prey_weight=mean(PWT,na.rm=T))%>% #mean prey weight
                                  ungroup()%>%
                                  
                                  group_by(pred,classn,prey)%>%
                                  mutate(prey_weight = prey_weight/max(prey_weight,na.rm=T))%>% #standardized prey weight amongst years by predator 
                                  ungroup()%>%
            
                                  group_by(Year,classn,prey)%>%
                                  summarise(prey_weight=mean(prey_weight,na.rm=T))%>% #mean prey weight, standarized amongst years, averaged amongst predators
                                  ungroup()%>%
            
                                  mutate(id = paste(Year,classn,sep="_"))%>%
                                  select(id,prey,prey_weight)%>%
                                  spread(prey,prey_weight,fill=0)%>%
                                  column_to_rownames('id')%>%
                                  mutate(tot = rowSums(.,na.rm=T))%>%
                                  filter(tot>0)%>% #filter out any zero catches
                                  select(-tot)%>%
                                  decostand(method = "total")%>%
                                  metaMDS(.,k=5)
      
      benthodiet_nmds_pa_all <- diet_df%>%
                                filter(!non_species,
                                       STATION %in% diet_nmds_env$STATION)%>%
                                left_join(.,diet_nmds_env%>%select(-location))%>%
                                group_by(Year,classn,prey)%>%
                                summarise(prey_weight=mean(PWT,na.rm=T))%>%
                                ungroup()%>%
                                mutate(id = paste(Year,classn,sep="_"))%>%
                                select(id,prey,prey_weight)%>%
                                spread(prey,prey_weight,fill=0)%>%
                                column_to_rownames('id')%>%
                                mutate(tot = rowSums(.,na.rm=T))%>%
                                filter(tot>0)%>% #filter out any zero catches
                                select(-tot)%>%
                                decostand(method = "total")%>%
                                dist(., method="binary")%>%
                                metaMDS(.,k=clusters)
      
      
      bentho_data_scores_all <- as.data.frame(scores(benthodiet_nmds_bray_all,"sites"))%>%
                                mutate(id=rownames(.),
                                       method="bray",
                                       stress=round(benthodiet_nmds_bray_all$stress,3))%>%
                                separate(id,c("year","classn"),sep="_")%>%
                                rbind(.,
                                      as.data.frame(scores(benthodiet_nmds_pa_all,"sites"))%>%
                                        mutate(id=rownames(.),
                                               method="pa",
                                               stress=round(benthodiet_nmds_pa_all$stress,3))%>%
                                        separate(id,c("year","classn"),sep="_"))%>%
                                select(NMDS1,NMDS2,year,classn,method,stress)
                              
      hull_data_bentho_all <- bentho_data_scores_all%>%
                              group_by(classn,method)%>%
                              slice(chull(x=NMDS1,y=NMDS2))
      
      
      #Convex hull plots of diet
      p1_bray_all <- ggplot(data=hull_data_bentho_all%>%filter(method=="bray"))+
                    geom_polygon(aes(x=NMDS1,y=NMDS2,fill=classn),alpha=0.30,col="black",lwd=0.1)+
                    geom_point(aes(x=NMDS1,y=NMDS2,fill=classn),shape=21,size=2)+
                    theme_bw()+
                    theme(axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          legend.position = "none",
                          strip.background = element_rect(fill="white"))+
                    labs(shape="",fill="",col="Sample year",title="Bray-Curtis")+
                    geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,
                                  label=paste("Stress =",round(benthodiet_nmds_bray_all$stress,3),"k =",benthodiet_nmds_bray_all$ndim)))+
                    scale_fill_viridis(discrete=TRUE)
      
      p1_pa_all <- ggplot(data=hull_data_bentho_all%>%filter(method=="pa"))+
                    geom_polygon(aes(x=NMDS1,y=NMDS2,fill=classn),alpha=0.30,col="black",lwd=0.1)+
                    geom_point(aes(x=NMDS1,y=NMDS2,fill=classn),shape=21,size=2)+
                    theme_bw()+
                    theme(axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          legend.position = "right",
                          strip.background = element_rect(fill="white"))+
                    labs(shape="",fill="",col="Sample year",title="Euclidean")+
                    geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(benthodiet_nmds_pa_all$stress,3),"k =",benthodiet_nmds_pa_all$ndim)))+
                    scale_fill_viridis(discrete=TRUE)
      
      #combine using patchwork
      p1_combo_all <- p1_bray_all + p1_pa_all + plot_layout(ncol=2)
      
      #save plot
      ggsave("output/CrabSurvey/Benthoscape_nmds_comp_all.png",p1_combo_all,width=7,height=5,units="in",dpi=300)
      
      
      ##By individual sdm (no pooling) - have to fill in the sparse matrix to make this work 
      
      prey_sp_list <- diet_df%>%filter(!non_species,STATION %in% diet_nmds_env$STATION)%>%pull(prey)%>%unique() # list of species to populate the sparse matrix (all prey items amongst stations)
      
      #assemble the nmds dataset - where each row is a species and each column is all prey items consumed. 
     diet_nmds_indv1 <- diet_df%>%
                        filter(!non_species,
                               STATION %in% diet_nmds_env$STATION)%>%
                        left_join(.,diet_nmds_env%>%select(-location))%>%
                        select(PRED_SEQ,prey,PWT,classn,pred)
     
     diet_nmds_indv2 <- diet_nmds_indv1%>%
                        group_by(PRED_SEQ)%>%
                        complete(prey=prey_sp_list)%>%
                        ungroup()%>%
                        select(-c(classn,pred))%>%#because the 'complete' populates with NAs along with PWT
                        left_join(.,diet_nmds_indv1%>%distinct(PRED_SEQ,.keep_all=TRUE)%>%select(PRED_SEQ,pred,classn)) 
     
     diet_nmds_indv <- diet_nmds_indv2%>%
                       group_by(PRED_SEQ,prey)%>%
                       summarise(PWT = sum(PWT,na.rm=T))%>% # some stomachs have two weights for two of a prey item
                       ungroup()%>%
                       left_join(.,diet_nmds_indv2%>%select(PRED_SEQ,pred,classn)%>%distinct(PRED_SEQ,.keep_all=TRUE))%>%
                       mutate(id = paste(pred,PRED_SEQ,classn,sep="_"))%>%
                       select(-c(classn,PRED_SEQ,pred))%>%
                       spread(prey,PWT,fill=0)%>%
                       column_to_rownames('id')%>%
                       mutate(tot = rowSums(.,na.rm=T))%>%
                       filter(tot>0)%>% #filter out any zero catches
                       select(-tot)%>%
                       decostand(method = "total")
     
     
     #run the metaMDS using bray (prey weight, weighted) and euclidean (PA) - takes time for each 
     diet_nmds_indv_bray <- diet_nmds_indv%>%
                            metaMDS(.,k=5)
     
     diet_nmds_indv_pa <- diet_nmds_indv%>%
                         dist(., method="binary")%>%
                         metaMDS(.,k=5)
     
     bentho_data_scores_indv <- as.data.frame(scores(diet_nmds_indv_bray,"sites"))%>%
                               mutate(id=rownames(.),
                                      method="bray",
                                      stress=round(diet_nmds_indv_bray$stress,3))%>%
                               separate(id,c("pred","PRED_SEQ","classn"),sep="_")%>%
                               rbind(.,
                                     as.data.frame(scores(diet_nmds_indv_pa,"sites"))%>%
                                       mutate(id=rownames(.),
                                              method="pa",
                                              stress=round(diet_nmds_indv_pa$stress,3))%>%
                                       separate(id,c("pred","PRED_SEQ","classn"),sep="_"))%>%
                               select(NMDS1,NMDS2,pred,PRED_SEQ,classn,method,stress)
     
     #get convex hulls for each 
     hull_data_bentho_indv <- bentho_data_scores_indv %>%
                              group_by(classn,method)%>%
                              slice(chull(x=NMDS1,y=NMDS2))
     
     hull_data_bentho_indv_pred <- bentho_data_scores_indv %>%
                                   group_by(pred,method)%>%
                                   slice(chull(x=NMDS1,y=NMDS2))
     
     #plot it up
     pa_all_indv <- ggplot(data=hull_data_bentho_indv%>%filter(method=="pa"))+
                     geom_polygon(aes(x=NMDS1,y=NMDS2,fill=classn),alpha=0.30,col="black",lwd=0.1)+
                     geom_point(aes(x=NMDS1,y=NMDS2,fill=classn),shape=21,size=2)+
                     theme_bw()+
                     theme(axis.text.x=element_blank(),
                           axis.ticks.x=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.position = "right",
                           strip.background = element_rect(fill="white"))+
                     labs(shape="",fill="",col="Sample year",title="Euclidean")+
                     geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(diet_nmds_indv_pa$stress,3),"k =",diet_nmds_indv_pa$ndim)))+
                     scale_fill_viridis(discrete=TRUE)
     
     bray_all_indv <- ggplot(data=hull_data_bentho_indv%>%filter(method=="bray"))+
                     geom_polygon(aes(x=NMDS1,y=NMDS2,fill=classn),alpha=0.30,col="black",lwd=0.1)+
                     geom_point(aes(x=NMDS1,y=NMDS2,fill=classn),shape=21,size=2)+
                     theme_bw()+
                     theme(axis.text.x=element_blank(),
                           axis.ticks.x=element_blank(),
                           axis.text.y=element_blank(),
                           axis.ticks.y=element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           legend.position = "right",
                           strip.background = element_rect(fill="white"))+
                     labs(shape="",fill="",col="Sample year",title="Bray-Curtis")+
                     geom_text(aes(x=Inf,y=-Inf,hjust=1.05,vjust=-0.5,label=paste("Stress =",round(diet_nmds_indv_bray$stress,3),"k =",diet_nmds_indv_bray$ndim)))+
                     scale_fill_viridis(discrete=TRUE)
     
     #combine using patchwork
     indv_combo_all <- bray_all_indv + pa_all_indv + plot_layout(ncol=2)
     
     #save plot
     ggsave("output/CrabSurvey/Benthoscape_nmds_comp_indv.png",indv_combo_all,width=7,height=5,units="in",dpi=300)
     
     
     ##investigate year of year changes by target species 
     
     predyr_mean <- bentho_data_scores_indv %>%
                    filter(pred %in% target_sp$spec)%>%
                    left_join(.,diet_df%>%select(PRED_SEQ,Year)%>%distinct(PRED_SEQ,.keep_all=TRUE)%>%mutate(PRED_SEQ = as.character(PRED_SEQ)))%>%
                    group_by(pred,Year,method)%>%
                    summarise(NMDS1_sd=sd(NMDS1),
                              NMDS2_sd=sd(NMDS2),
                              NMDS1=mean(NMDS1),
                              NMDS2=mean(NMDS2))%>%
                    ungroup()%>%
                    data.frame()
     
     predyr_df <- bentho_data_scores_indv %>%
                  filter(pred %in% target_sp$spec)%>%
                  left_join(.,diet_df%>%select(PRED_SEQ,Year)%>%distinct(PRED_SEQ,.keep_all=TRUE)%>%mutate(PRED_SEQ = as.character(PRED_SEQ)))%>%
                  data.frame()
     
     
     ggplot(data=predyr_df%>%filter(method=="bray",pred=="GADUS MORHUA"),aes(x=NMDS1,y=NMDS2,fill=Year))+
       # geom_segment(aes(xend=c(tail(NMDS1, n=-1), NA),yend=c(tail(NMDS2, n=-1), NA)),
       #              arrow=arrow(length=unit(0.25,"cm"),type="closed"),lty=3)+
       geom_point(shape=21,size=2,col="black")+
       theme_bw()+
       facet_wrap(~pred,ncol=2,scales="free")
       
