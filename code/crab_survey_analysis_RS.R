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
    
### Target species ---------
    target_sp <- data.frame(spec=c("ANARHICHAS LUPUS", "UROPHYCIS TENUIS", "HIPPOGLOSSOIDES PLATESSOIDES",
                                   "GLYPTOCEPHALUS CYNOGLOSSUS", "CHIONOECETES OPILIO", "GADUS MORHUA", 
                                   "SEBASTES SP.", "AMBLYRAJA RADIATA"),
                            common=c("Atlantic wolffish","White hake","American plaice",
                                     "Witch flounder","Snow crab","Atlantic cod",
                                     "Redfish sp","Thorny skate"))
    
#### function -----
    
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
    write.csv(catchdat_stand%>%
                select(COMM,SPEC,species_filter)%>%
                distinct(species_filter,.keep_all=TRUE),"output/taxonomic_raw.csv",row.names=FALSE)
    
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
      mutate(inside=as.logical(st_intersects(.,sab, sparse=TRUE))) #this doesn't work

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
               labs(col="",x="Year",y="Median recorded size");p_med_a
    
    p_med_b <- ggplot()+
               geom_vline(xintercept = 2017,lty=2)+
               geom_line(data=plot_df,aes(x=year,y=med,col=location,group=station),lwd=0.5)+
               geom_point(data=plot_df,aes(x=year,y=med,col=location,group=station),size=2)+
               theme_bw()+
               theme(strip.background = element_rect(fill="white"),
                    legend.position = "bottom")+
               facet_wrap(~species2,nrow=3,scales="free_y")+
               labs(col="",x="Year",y="Median recorded size");p_med_b
    
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
                   labs(col="",x="Year",y="90th percentile size");p_bigfish_a
    
    p_bigfish_b <- ggplot()+
                   geom_vline(xintercept = 2017,lty=2)+
                   geom_line(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),lwd=0.5)+
                   geom_point(data=plot_df,aes(x=year,y=bigfish,col=location,group=station),size=2)+
                   theme_bw()+
                   theme(strip.background = element_rect(fill="white"),
                        legend.position = "bottom")+
                   facet_wrap(~species2,ncol=2,scales="free_y")+
                   labs(col="",x="Year",y="90th percentile size");p_bigfish_b
      
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
      labs(col="",x="Year",y="Richness per station")+
      theme(legend.position = "bottom");p_richness_a
    
    p_richness_b <- ggplot()+
      geom_vline(xintercept = 2017,lty=2)+
      geom_line(data=biodiv_df,aes(x=year,y=richness,group=station,col=location),lty=2,lwd=0.25,alpha=0.5)+
      geom_point(data=biodiv_df,aes(x=year,y=richness,group=station,col=location),size=2,alpha=0.5)+
      stat_smooth(data=biodiv_df,aes(x=year,y=richness,col=location),lwd=2)+
      theme_bw()+
      labs(col="",x="Year",y="Richness per station")+
      theme(legend.position = "bottom");p_richness_b
    
    ggsave("output/CrabSurvey/Biodiversity_inside-outside.png",p_richness_a)
    ggsave("output/CrabSurvey/Biodiversity_inside-outside_wSmooth.png",p_richness_b)
    
    ##catch number per set
    
    freq_sp <- catchdat_stand%>%
               filter(minID != "Actiniaria")%>% #sea anemones are captures but not really a species ae can tell a story around 
               group_by(minID)%>%
               summarise(count=n())%>%
               ungroup()%>%
               arrange(-count)%>%
               slice(1:10)%>%
               pull(minID)%>%
               c(.,"Anarhichas lupus","Merluccius bilinearis")
    
    num_df <- catchdat_stand%>%
              filter(minID %in% freq_sp,
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
                  stat_smooth(data=num_df,aes(x=year,y=number,col=location),lwd=2,se=FALSE)+
                  scale_y_log10()+
                  facet_wrap(~minID,nrow=3,scales="free_y")+
                  labs(col="",x="Year",y="Count per set");p_number_a
                
    p_number_b <- ggplot()+
                  geom_vline(xintercept = 2017,lty=2)+
                  geom_line(data=num_df,aes(x=year,y=num_stand,col=location,group=station),lty=2,lwd=0.25,alpha=0.5)+
                  geom_point(data=num_df,aes(x=year,y=num_stand,col=location,group=station),size=2,alpha=0.55)+
                  theme_bw()+
                  stat_smooth(data=num_df,aes(x=year,y=num_stand,col=location),lwd=2)+
                  facet_wrap(~minID,nrow=3,scales="free_y")+
                  labs(col="",x="Year",y=expression(paste("Count per set / Max count per set")));p_number_b
                
    ggsave("output/CrabSurvey/CatchNumber_inside-outside.png",p_number_a)
    ggsave("output/CrabSurvey/CatchNumber_inside-outside_standardized.png",p_number_b)
    
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
    
#### DIET ANALYSIS --------------
    
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
              st_as_sf(coords=c("SLONGDD","SLATDD"),crs=latlong)%>%
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
                  
