#Figure 1 for the Acoustic Report

#load libraries ----
library(dplyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(lubridate)
library(arcpullr)
library(marmap)
library(stars)
library(ggrepel)
library(worrms)
library(patchwork)
library(raster)
library(rphylopic)
library(ggspatial)
library(viridis)
library(terra)
library(tidyterra)

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load St Anns Bank MPA shapefile 
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

sab_banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%st_transform(latlong)

#load bioregion and get plotting limits----
bioregion <- read_sf("data/Shapefiles/MaritimesPlanningArea.shp")%>%
             st_transform(latlong)

mar_net <- read_sf("data/Shapefiles/networksites_proposed_OEM_MPA_20220617.shp")%>%
            st_transform(latlong)

plot_lim <- bioregion%>% #bounding area for the inset (zoom out)
            st_transform(utm)%>%
            st_buffer(50)%>%
            st_transform(latlong)%>%
            st_bbox()

#get the banks
bank_df <- data.frame(bank=c("Curdo Bank","Scatarie Bank"),
                      lon=c(-59.598049,-59.213127),
                      lat=c(45.893003,45.988401))%>%
            st_as_sf(coords=c("lon","lat"),crs=latlong)

bank_range <- bank_df%>%
              filter(bank=="Scatarie Bank")%>%
              st_transform(utm)%>%
              st_buffer(10)%>%
              st_transform(latlong)%>%
              st_bbox()

#for inset plot 
plot_boundaries <- c(c(-60,-58,45.6,46.6)) #for inset plot
inset_box <- sab%>%st_bbox()
inset_box[c(1,3,2,4)] <- plot_boundaries     
inset_box <- inset_box%>%st_as_sfc()%>%st_as_sf()

basemap_inset <- rbind(
                  ne_states(country = "Canada",returnclass = "sf")%>%
                  dplyr::select(name_en,geometry)%>%
                  st_as_sf()%>%
                  st_union()%>%
                  st_transform(latlong)%>%
                  st_as_sf()%>%
                  mutate(country="Canada"),
                  ne_states(country = "United States of America",returnclass = "sf")%>%
                    dplyr::select(name_en,geometry)%>%
                    st_as_sf()%>%
                    st_union()%>%
                    st_transform(latlong)%>%
                    st_as_sf()%>%
                    mutate(country="US")
                  )

global_basemap <- ne_states() #this is a very large scale map


#bathymetry 
    bbsab <- st_bbox(sab)
    
    #extract bathymetric data
    noaabathy_sab <- getNOAA.bathy(bbsab[1]-1,bbsab[3]+1,bbsab[2]+1,bbsab[4]-1,resolution = 0.25,keep=T)
    
    #create xyz dataframes that can be used by geom_contour
    isobath_sab_df <- as.xyz(noaabathy_sab)%>%
                      rename(lon=1,lat=2,depth=3)
    
    
## station depth extract ------
    
    #SAB dem
    dem_sab <- raster("data/Bathymetry/sab_dem.tif")
    
    reciever_depth <- read.csv("data/Acoustic/SAB_deployments_recovered_allFeb2024.csv")%>%
                      rename(lon=DEPLOY_LONG,lat=DEPLOY_LAT)%>%
                      dplyr::select(-geometry)%>%
                      st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                      st_transform(proj4string(dem_sab))%>%
                      mutate(depth = round(raster::extract(dem_sab,as_Spatial(.)),1))%>%#extract depth
                      st_transform(latlong)%>%
                      st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
                      data.frame()%>%
                      dplyr::select(-geometry)
                      
    #write.csv(reciever_depth,file="data/Acoustic/Acoustic_Stations.csv",row.names=FALSE)

#load coastline and make basemap ------
coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")

#download OTN data --  

  #code from https://rpubs.com/MaritimesMSP/OTN-DFO-Summary
  
  proj_start <- ymd("20190101")
  proj_end <- ymd("20230909")
  proj_long_upp <- -40.00
  proj_long_low <- -70.00
  proj_lat_upp <- 60.00
  proj_lat_low <- 40.00
  
  geoserver_receivers <- readr::read_csv('https://members.oceantrack.org/geoserver/otn/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=otn:stations_receivers&outputFormat=csv', guess_max = 13579)
  
  otn_stations <- geoserver_receivers %>%
                  filter(!is.na(deploy_date)) %>% # Remove rows with missing deploy_date values
                  filter((deploy_date > proj_start & deploy_date < proj_end) | # Select stations that fall within the project timeframe defined above
                           (recovery_date < proj_end & recovery_date > proj_start) |
                           (deploy_date < proj_end & is.na(recovery_date) & deploy_date > proj_start - duration(18, 'months')) |
                           (grepl('VR3', model) & deploy_date < proj_end & is.na(recovery_date) & deploy_date > proj_start - duration(4, 'years')) | # Select specific models within certain date ranges
                           (grepl('VR4', model) & deploy_date < proj_end & is.na(recovery_date) & deploy_date > proj_start - duration(6, 'years'))) %>%
                  filter(stn_lat >= proj_lat_low & stn_lat <= proj_lat_upp &# Filter stations based on latitude and longitude bounds
                           stn_long >= proj_long_low & stn_long <= proj_long_upp)%>%
                  st_as_sf(coords=c("stn_long","stn_lat"),crs=latlong)%>%
                  mutate(year=year(deploy_date),
                         array = ifelse(year<2021,"2015-2020","2020-2024"))%>%
                  filter(year>2014)#just within the SAB window
                  
  
  url <- "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Oceans_Act_Marine_Protected_Areas/MapServer/0"
  
  MPAs <- get_spatial_layer(url) %>% 
    st_make_valid() %>% 
    st_crop(xmin=proj_long_low,
            ymin=proj_lat_low,
            xmax=proj_long_upp,
            ymax=proj_lat_upp)%>%
    st_transform(latlong)
    
## OTN tag deployment locations ---------
  tag_url <- 'https://members.oceantrack.org/geoserver/otn/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=otn:animals&outputFormat=csv'
  
  geoserver_tag_releases <- readr::read_csv(tag_url, guess_max = 13579) %>%
                            filter(collectioncode == "SABMPA")%>%
                            mutate(id = paste(abs(longitude),latitude,sep="_"))%>%
                            distinct(id,.keep_all=TRUE)%>%
                            st_as_sf(coords=c("longitude","latitude"),crs=latlong)
  
  tag_tax <- data.frame(aphiaid = unique(geoserver_tag_releases$aphiaid))%>%
             rowwise()%>%
             mutate(sp=wm_id2name(aphiaid))%>%
             data.frame()%>%
             mutate(common=case_when(sp == "Amblyraja radiata" ~ "Thorny skate",
                                     sp == "Anarhichas lupus" ~ "Atlantic wolffish",
                                     sp == "Chionoecetes opilio" ~ "Snow crab",
                                     sp == "Gadus morhua" ~ "Atlantic cod",
                                     sp == "Hippoglossus hippoglossus" ~ "Atlantic halibut",
                                     sp == "Myoxocephalus scorpius" ~ "Shorthorn sculpin"))
  
  geoserver_tag_releases <- geoserver_tag_releases%>%left_join(.,tag_tax)
  
  #some of the geoserver extracts are missing information. Use this csv instead
  tag_df <- read.csv("data/Acoustic/SABMPA_tagrelease_Dec2023.csv")%>%
            rename(common=COMMON_NAME_E)%>%
            mutate(id = paste(abs(RELEASE_LONGITUDE),RELEASE_LATITUDE,common,sep="_"),
                   year=year(as.POSIXct(UTC_RELEASE_DATE_TIME)),
                   common=case_when(common == "ATLANTIC COD" ~ "Atlantic cod",
                                    common == "ATLANTIC HALIBUT" ~ "Atlantic halibut",
                                    common == "Atlantic Striped wolffish" ~ "Atlantic wolffish",
                                    TRUE ~ common),
                   array = ifelse(year<2021,"2015-2020","2020-2024"))%>% #only need unique locations
            distinct(id,.keep_all=TRUE)%>%
            st_as_sf(coords=c("RELEASE_LONGITUDE","RELEASE_LATITUDE"),crs=latlong)%>%
            dplyr::select(common,year,array,geometry)%>%
            rbind(.,geoserver_tag_releases%>% #add the snow crab data
                    filter(common == "Snow crab")%>%
                    rename(year=yearcollected)%>%
                    mutate(array = ifelse(year<2021,"2015-2020","2020-2024"))%>%
                    dplyr::select(common,year,array,geometry))
  
 #bounding area for the primary plot----------
      bounding_area <- sab%>%
        st_transform(utm)%>%
        st_buffer(dist = 100)%>%
        st_bbox()%>%
        st_as_sfc()%>%
        st_as_sf()%>%
        st_transform(latlong)%>%
        st_bbox()%>%
        st_as_sfc()%>%
        st_as_sf()%>%
        suppressWarnings()%>%
        suppressMessages()
      
      basemap <- coast_hr%>%
        st_intersection(bounding_area%>%st_transform(st_crs(coast_hr)))%>%
        st_transform(latlong)%>%
        suppressWarnings()
  
      
#make primary map
      recievers_sf <- otn_stations%>%
                      filter(collectioncode=="SABMPA")%>%
                      mutate(array = ifelse(year<2021,"2015-2020","2020-2024"),
                             lon = st_coordinates(.)[,1],
                             lat = st_coordinates(.)[,2])%>%
                      st_transform(proj4string(dem_sab))%>%
                      mutate(depth = round(raster::extract(dem_sab,as_Spatial(.)),1))%>%#extract depth
                      st_transform(latlong)
      
      recievers_sf_simple <- recievers_sf%>%filter(year %in% c(2019,2021))%>%
                             mutate(group=ifelse(lat>46.13,"north array",array))
     
      recievers_lines <- recievers_sf_simple%>%
                                  group_by(group)%>%
                                  arrange(lon)%>%
                                  summarise(do_union=FALSE)%>%
                                  st_cast("LINESTRING")%>%
                                  mutate(array=ifelse(group=="2020-2024","2020-2024","2015-2020"))
      banks_df <- sab_banks%>%
                  st_centroid()
      
      primary_plot <- ggplot()+
                      geom_contour(data=isobath_sab_df,aes(x=lon,y=lat,z=depth),breaks=-250,color = "grey80", size = 0.5)+
                      geom_sf(data=sab_banks,fill="forestgreen")+
                      # ggrepel::geom_label_repel(data=sab_banks,aes(label = name, geometry = geometry),
                      #                           stat = "sf_coordinates",
                      #                           min.segment.length = 0)+
                      geom_sf(data=otn_stations,size=0.1)+
                      geom_sf(data=coast_hr,fill="grey70")+
                      geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
                      geom_sf(data=mar_net_df%>%filter(TYPE == "TBD"),fill="grey70",alpha=0.25,lty=3)+
                      geom_sf(data=recievers_sf,aes(fill=array),shape=21)+
                      theme_bw()+
                      theme(legend.position = "bottom")+
                      guides(shape = guide_legend(override.aes = list(size=6)))+
                      labs(fill="",x="",y="")+
                      coord_sf(expand=0,xlim=plot_boundaries[1:2],ylim=plot_boundaries[3:4])
      
      ggsave("output/Acoustic/primary_plot.png",primary_plot,height=6,width=6,units="in",dpi=300)
      
      
#inset map
      
      inset_labels <- otn_stations%>%
                      filter(collectioncode == "HFX",year==2020)%>%
                      group_by(collectioncode)%>%
                      st_combine()%>%
                      st_cast("LINESTRING")%>%
                      st_centroid()%>%
                      st_as_sf()%>%
                      mutate(lab="HFX")%>%
                      rbind(.,otn_stations%>%
                              filter(collectioncode == "CBS",year==2020)%>%
                              group_by(collectioncode)%>%
                              st_combine()%>%
                              st_cast("LINESTRING")%>%
                              st_centroid()%>%
                              st_as_sf()%>%
                              mutate(lab="CBS"))%>%
                      mutate(lon=st_coordinates(.)[,1],
                             lat=st_coordinates(.)[,2])
      
      #modified MPA network shapefile
      mar_net_df <- mar_net%>%
                    mutate(NAME = case_when(NAME == "St Anns Bank Marine Protected Area" ~ "St. Anns Bank Marine Protected Area", #to match the FGP pull
                                            NAME == "The Gully Marine Protected Area" ~ "Gully Marine Protected Area",
                                            TRUE ~ NAME),
                           TYPE = case_when(NAME %in% c("Corsair/Georges Canyons Conservation Area","Western Emerald Bank Conservation Area",
                                                        "Emerald Basin Sponge Conservation Area","Sambro Bank Sponge Conservation Area",
                                                        "Jordan Basin Conservation Area","Eastern Canyons Marine Refuge",
                                                        "Northeast Channel Coral Conservation Area")~ "MR",
                                            TRUE ~ TYPE))%>%
                    filter(!NAME %in% MPAs$NAME_E) #these are plotted as part of the GIS pull, which has zones. 
      
      inset_plot <- ggplot()+
                    geom_sf(data=bioregion,fill=NA)+
                    geom_contour(data=isobath_df,aes(x=lon,y=lat,z=depth),breaks=-250,color = "grey80", size = 0.5)+
                    geom_sf(data=mar_net_df%>%filter(TYPE == "TBD"),fill="grey70",alpha=0.25,lty=3)+
                    geom_sf(data=mar_net_df%>%filter(TYPE == "MR"),fill="grey10",alpha=0.3,lty=3)+
                    geom_sf(data=mar_net_df%>%filter(TYPE == "AOI"),fill="salmon2",alpha=0.3,lty=3)+
                    geom_sf(data=otn_stations,size=0.25,pch=20)+
                    geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
                    geom_sf(data=otn_stations%>%filter(collectioncode %in% c("HFX","CBS")),col="red",size=0.7)+
                    geom_sf(data=basemap_inset,fill="grey70")+
                    geom_sf(data=inset_box,fill=NA)+ #this is the zoomed out
                    coord_sf(expand=0,xlim=plot_lim[c(1,3)],ylim=plot_lim[c(2,4)])+
                    theme_bw()+
                    labs(x="",y="")+
                    theme(axis.text=element_blank())
      
      ggsave("output/Acoustic/inset_plot.png",inset_plot,height=6,width=6,units="in",dpi=300)
      
      #tag map
      plot_list <- NULL
      
      for(i in unique(tag_df$common)){
        
        array_df <- tag_df%>%filter(common == i)%>%pull(array)%>%unique()
        
        for(j in array_df){
        
        temp_tag <- tag_df%>%filter(common == i,array==j)
        
                        p1 <- ggplot()+
                        geom_sf(data=shelfbreak,col="grey20",lty=2,fill=NA)+
                        geom_sf(data=otn_stations%>%filter(array==j),size=0.1)+
                        geom_sf(data=basemap)+
                        geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
                        geom_sf(data=recievers_sf%>%filter(array==j),fill="white",shape=21)+
                        geom_sf(data=temp_tag,fill="coral2",shape=21,size=2)+
                        geom_text(data=temp_tag%>%slice(1),aes(label=array),x=(-59.9),y=46.5,fontface = "bold",hjust=0)+
                        theme_bw()+
                        theme(legend.position = "bottom",
                              strip.background = element_rect(fill="white"),
                              axis.text=element_blank(),
                              axis.title = element_blank())+
                        labs(fill="")+
                        coord_sf(expand=0,xlim=plot_boundaries[1:2],ylim=plot_boundaries[3:4])+
                        facet_wrap(~common)
                        
                        assign(paste(gsub(" ","_",i),gsub("-","_",j),"plot",sep="_"),p1)
                        
                        plot_list <- c(plot_list,paste(gsub(" ","_",i),gsub("-","_",j),"plot",sep="_"))
                    
        } #end array loop
      } #end common loop
      

      combo_tag_plot <- (Atlantic_cod_2015_2020_plot + Atlantic_cod_2020_2024_plot)/
                        (Atlantic_tomcod_2015_2020_plot + Atlantic_tomcod_2020_2024_plot)/
                        (Atlantic_halibut_2020_2024_plot + Atlantic_wolffish_2015_2020_plot)/
                        (Shorthorn_Sculpin_2015_2020_plot+Shorthorn_sculpin_2020_2024_plot)/
                        (Thorny_Skate_2015_2020_plot+Snow_crab_2015_2020_plot)
      
      ggsave("output/Acoustic/tagging_plot.png",combo_tag_plot,height=9,width=5,units="in",dpi=300)
      
      
### Tagging plot ----- 
      
    tag_df <- read.csv("data/Acoustic/qdet_releaselocs.csv")%>%
              rename(lon=RELEASE_LONGITUDE,lat=RELEASE_LATITUDE)%>%
              filter(!is.na(lon))%>%
              st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
              mutate(group=case_when(Common.Name %in% c("White shark","Shortfin mako","Blue shark","Blue shark/Mako/Bluefin tuna","Porbeagle shark") ~ "Sharks",
                                     Common.Name %in% c("Atlantic sturgeon","Bluefin tuna")~"Pelagics",
                                     Common.Name %in% c("Atlantic cod","Snow crab","Atlantic halibut") ~ "Atlantic fisheries",
                                     Common.Name %in% c("Atlantic salmon","American eel") ~ "Diadramous",
                                     Common.Name %in% c("Grey seal","Leatherback turtle") ~ "Large pelagics",
                                     TRUE ~ NA),
                     Common.Name = factor(Common.Name,levels=c("American eel","Atlantic cod","Atlantic halibut","Atlantic salmon","Atlantic sturgeon",
                                                               "Blue shark","Blue shark/Mako/Bluefin tuna","Bluefin tuna","Grey seal","Leatherback turtle",           
                                                               "Porbeagle shark","Shortfin mako","Snow crab","White shark")))
              
      
    
    tag_bound <- tag_df%>%
                 st_bbox()%>%
                 st_as_sfc()%>%
                 st_buffer(0.5)%>% #0.5 degree buffer
                 st_bbox()
    
    tag_df_sp <- data.frame(sp = unique(tag_df$Common.Name))%>%
                 mutate(search=case_when(sp == "Atlantic salmon"~"Salmo trutta", #closest match in phylopic
                                         sp == "Atlantic halibut" ~"Hippoglossoides platessoides", #closest match in phylopic
                                         sp == "Leatherback turtle"~"Dermochelys coriacea",
                                         sp == "Porbeagle shark"~"Lamna nasus",
                                         sp == "Snow crab"~"Pugettia quadridens", #closest I could find on phylopic
                                         sp == "Blue shark/Mako/Bluefin tuna" ~ "Isurus oxyrinchus", #general shark pic for this species group
                                         TRUE ~sp))%>%
                 rowwise()%>%
                 mutate(uid = get_uuid(search),
                        sp = ifelse(sp == "Blue shark/Mako/Bluefin tuna",gsub("/","-",sp),sp))%>%
                 data.frame()
    
    
    
    # five scales
    
    scale0 <- sab%>%
              st_transform(utm)%>%
              st_buffer(7)%>%
              st_transform(latlong)%>%
              st_bbox()
    
    scale0_box <- scale0%>%st_as_sfc()%>%st_as_sf()
    
    scale0_pts <- tag_df%>%st_intersection(scale0_box)
    
    scale0_plot <- ggplot()+
                    geom_sf(data=shelfbreak,col="grey20",lty=2,fill=NA)+
                    geom_sf(data=basemap_inset)+
                    geom_sf(data=sab_zones,fill="cornflowerblue",alpha=0.3)+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=3,pch=21)+
                    coord_sf(expand=0,xlim=scale0[c(1,3)],ylim=scale0[c(2,4)])+
                    theme_bw()+
                    annotation_scale(location="br")+
                    theme(legend.position = "none",
                          axis.text=element_blank())
    
    
    #scale 1 (fine scale)
    scale1 <- sab%>%
              st_transform(utm)%>%
              st_buffer(300)%>%
              st_transform(latlong)%>%
              st_bbox()
    
    scale1[2] <- 44
    
    scale1_box <- scale1%>%st_as_sfc()%>%st_as_sf()
    
    scale1_pts <- tag_df%>%
                  st_intersection(scale1_box)%>%
                  filter(!animal_id %in% scale0_pts$animal_id)
    
    scale1_plot <- ggplot()+
                    geom_sf(data=basemap_inset)+
                    geom_sf(data=sab_zones,fill="cornflowerblue",alpha=0.3)+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=1,pch=21)+
                    geom_sf(data=tag_df%>%filter(animal_id %in% scale1_pts$animal_id),aes(fill=Common.Name),size=4,pch=21)+
                    coord_sf(expand=0,xlim=scale1[c(1,3)],ylim=scale1[c(2,4)])+
                    theme_bw()+
                    annotation_scale(location="br")+
                    theme(legend.position = "none",
                          axis.text=element_blank())
    
    #scale 2 fine-medium scale 
    scale2 <- sab%>%
              st_transform(utm)%>%
              st_buffer(650)%>%
              st_transform(latlong)%>%
              st_bbox()
    
    scale2_box <- scale2%>%st_as_sfc()%>%st_as_sf()
    
    scale2_pts <- tag_df%>%
                  st_intersection(.,scale2_box)%>%
                  filter(!animal_id %in% c(scale0_pts$animal_id,scale1_pts$animal_id))
    
    scale2_plot <- ggplot()+
                    geom_sf(data=basemap_inset)+
                    geom_sf(data=sab_zones,fill="cornflowerblue",alpha=0.3)+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=1,pch=21)+
                    geom_sf(data=tag_df%>%filter(animal_id %in% scale2_pts$animal_id),aes(fill=Common.Name),size=3,pch=21)+
                    coord_sf(expand=0,xlim=scale2[c(1,3)],ylim=scale2[c(2,4)])+
                    theme_bw()+
                    annotation_scale(location="br")+
                    theme(legend.position = "none",
                          axis.text=element_blank())
    
    #scale 3 medium-large scale
    scale3 <- sab%>%
              st_bbox()
    
    scale3[c(3,4)] <- scale2[c(3,4)]
    scale3[3] <- -51
    scale3[2] <- 38.5
    scale3[1] <- -78
    
    scale3_box <- scale3%>%st_as_sfc()%>%st_as_sf()
    
    scale3_pts <- tag_df%>%
                  st_intersection(.,scale3_box)
                
    scale3_plot <- ggplot()+
                   geom_sf(data=global_basemap)+
                   geom_sf(data=sab,fill="cornflowerblue",alpha=0.3)+
                   geom_sf(data=tag_df,aes(fill=Common.Name),size=2,pch=21)+
                   geom_sf(data=tag_df%>%filter(animal_id %in% scale3_pts$animal_id),aes(fill=Common.Name),size=4,pch=21)+
                   coord_sf(expand=0,xlim=scale3[c(1,3)],ylim=scale3[c(2,4)])+
                   theme_bw()+
                   annotation_scale(location="br")+
                   theme(legend.position = "none",
                         axis.text=element_blank())
    
    #scale 4 largest
    scale4 <- sab%>%
              st_bbox()
    
    scale4[c(3,4)] <- scale3[c(3,4)]
    scale4[3] <- -47
    scale4[2] <- 4
    scale4[1] <- -85
    
    #Sab distance circles
    sab_centre <- sab%>%
      st_centroid()
    
    sab_centre_terra <- sab_centre%>%vect()
    
    rings <- NULL
    
    for(i in c(100,500,1000,2500,5000)){rings <- rbind(rings,sab_centre_terra%>%buffer(i*1000)%>%st_as_sf())}
    
    
    scale4_box <- scale4%>%st_as_sfc()%>%st_as_sf()
    
    scale4_pts <- tag_df%>%
                  st_intersection(.,scale4_box)%>%
                  filter(!animal_id %in% c(scale0_pts$animal_id,scale1_pts$animal_id,scale2_pts$animal_id,scale3_pts$animal_id))
    
    scale4_plot <- ggplot()+
                  geom_sf(data=global_basemap)+
                  geom_sf(data=sab,fill="cornflowerblue",alpha=0.3)+
                  geom_sf(data=tag_df,aes(fill=Common.Name),size=3,pch=21)+
                  geom_sf(data=tag_df%>%filter(animal_id %in% scale4_pts$animal_id),aes(fill=Common.Name),size=3,pch=21)+
                  geom_sf(data=rings,lty=2,lwd=0.5,fill=NA)+
                  coord_sf(expand=0,xlim=scale4[c(1,3)],ylim=scale4[c(2,4)],label_axes = "-NE")+
                  theme_bw()+
                  annotation_scale(location = "br")+
                  annotation_north_arrow(location = "tr")+
                  theme(legend.position = "none")
    
    #make a legend plot
    legend_dummy <- ggplot()+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=5,pch=21)+
                    theme_bw()+
                    labs(fill="")+
                    coord_sf(expand=0,label_axes = "-NE")
    
    #save plot - note used Scale 3 and 4 in the end plot to make it simple
    ggsave("output/Acoustic/tagmap_scale0.png",scale0_plot+scale_fill_viridis(discrete=T),height=5,width=5,units="in",dpi=300)
    ggsave("output/Acoustic/tagmap_scale1.png",scale1_plot+scale_fill_viridis(discrete=T),height=5,width=6,units="in",dpi=300)
    ggsave("output/Acoustic/tagmap_scale2.png",scale2_plot+scale_fill_viridis(discrete=T),height=6,width=8,units="in",dpi=300)
    ggsave("output/Acoustic/tagmap_scale3.png",scale3_plot+scale_fill_viridis(discrete=T),height=6,width=10,units="in",dpi=300)
    ggsave("output/Acoustic/tagmap_scale4.png",scale4_plot+scale_fill_viridis(discrete=T),height=10,width=6,units="in",dpi=300) 
    ggsave("output/Acoustic/tagmap_legend.png",legend_dummy+scale_fill_viridis(discrete=T),height=5,width=5,units="in",dpi=300)

#### Interpolated species track plot -----------
    
    qdets <- read.csv("data/Acoustic/qdets_pos.csv")
    pos_cod <- read.csv("data/Acoustic/pos_acod.csv")%>%mutate(sp="Atlantic cod")
    pos_halibut <- read.csv("data/Acoustic/pos_halibut.csv")%>%mutate(sp="Atlantic halibut")
    pos_wolffish <- read.csv("data/Acoustic/pos_wolffish.csv")%>%mutate(sp="Atlantic wolffish")
    
    pos_df <- rbind(pos_cod,pos_halibut,pos_wolffish)%>%
              mutate(timestamp=as.POSIXct(bin_timestamp),
                     julian=yday(timestamp),
                     year=year(timestamp),
                     id=gsub("SABMPA-","",animal_id), #individual fish id
                     id=gsub("MPA-SAB-","",id),
                     id=sub("-.*", "", id),
                     id=as.numeric(id),
                     array=ifelse(year>2020,"2020-2024","2015-2020"))%>%
              st_as_sf(coords=c("longitude","latitude"),crs=latlong,remove=FALSE)
    
    pos_filter <- pos_df%>%
                  mutate(position_tag = paste(longitude,latitude,sep="-"))%>% #some points have multiple days but one coordinate
                  group_by(id,year)%>%
                  summarise(cnt=length(unique(julian)),
                            cnt2=length(unique(position_tag)))%>%
                  ungroup()%>%
                  arrange(cnt)%>%
                  filter(cnt>1 & cnt2>1)%>%
                  mutate(id2 = paste(id,year,sep="_"))%>%
                  pull(id2)
                  
    
    pos_df_line <-  pos_df%>%
                    mutate(id2 = paste(id,year,sep="_"))%>%
                    filter(id2 %in% pos_filter)%>%
                    arrange(sp,year,id,timestamp)%>%
                    group_by(id,year)%>%
                    mutate(geometry=st_union(geometry))%>%
                    st_cast("LINESTRING")
    
    plot_line_lim <- pos_df%>%
                     st_bbox()%>%
                     st_as_sfc()%>%
                     st_transform(utm)%>%
                     st_buffer(10)%>%
                     st_transform(latlong)%>%
                     st_bbox()
       
    
    p1 <- ggplot()+
      geom_sf(data=sab)+
      geom_sf(data=basemap_inset,fill="grey60")+
      geom_sf(data=sab_banks,fill="forestgreen",alpha=0.5)+
      geom_sf(data=recievers_sf%>%filter(year %in% c(2019,2021)))+
      geom_sf(data=pos_df_line)+
      coord_sf(expand=0,xlim=plot_line_lim[c(1,3)],ylim=plot_line_lim[c(2,4)])+
      theme_bw()+
      facet_grid(sp~array)
    
    ggsave("test.png",height=6,width=6,dpi=)
    
    #count raster
    points_extent <- pos_df%>%
                     st_intersection(.,sab)%>%
                     st_transform(utm)%>%
                     st_buffer(10)%>%
                     st_bbox()%>%
                     st_as_sfc()%>%
                     st_as_sf()
    
    raster_empty <-  raster(extent(points_extent),res=2,crs=utm) # 1 km resolution 
    
    raster_list <- list()
    
    for(i in unique(pos_df$sp)){
      for(j in unique(pos_df$array)){
        
        temp_df <- pos_df%>%
                    st_transform(utm)%>%
                    filter(sp==i,array==j)
        
        if(nrow(temp_df)>0){
    
            point_counts <- rasterize(temp_df%>%dplyr::select(geometry), raster_empty, fun = "count")%>%
                            projectRaster(.,crs=latlong)%>%
                            crop(.,extent(sab))%>%
                            mask(.,as_Spatial(sab))
            
            raster_list[[paste(i,j,sep="_")]] <- point_counts
    
        } #nend if check
      } # end j loop
    } # end i loop
    
    plot_df <- data.frame(ras_name=names(raster_list))%>%separate(ras_name,c("sp","array"),sep="_",remove = FALSE)
    
    
   
    #other tagged
    qdets <- read.csv("data/Acoustic/qdets_pos.csv")%>%
             mutate(sp = case_when(grepl("BLUEFIN",animal_id) ~ "Bluefin tuna",
                                   grepl("BLUESH",animal_id) ~ "Blue shark",
                                   grepl("COD",animal_id) ~ "Atlantic cod",
                                   grepl("CRAB",animal_id) ~ "Snow crab",
                                   grepl("EEL",animal_id) ~ "Atlantic eel",
                                   grepl("GSEAL",animal_id) ~ "Grey seal",
                                   grepl("HAL",animal_id) ~ "Atlantic halibut",
                                   grepl("LHBACK",animal_id) ~ "Leatherback sea turtle",
                                   grepl("MAKO",animal_id) ~ "Mako shark",
                                   grepl("PORB",animal_id) ~ "Porbeagle shark",
                                   grepl("SALM",animal_id) ~ "Atlantic salmon",
                                   grepl("STURG",animal_id) ~ "Atlantic sturgeon",
                                   grepl("WHITE",animal_id) ~ "White shark",
                                   TRUE ~ NA),
                    timestamp = as.POSIXct(bin_timestamp),
                    year=year(timestamp),
                    month=month(timestamp),
                    season=case_when(month %in% 1:3 ~ "Winter",
                                     month %in% 4:6 ~ "Spring",
                                     month %in% 7:9 ~ "Summer",
                                     month %in% 10:12 ~ "Winter"),
                    array=ifelse(year>2020,"2020-2024","2015-2020"))%>%
              st_as_sf(coords=c("longitude","latitude"),crs=latlong,remove=FALSE)
    
    qdet_rasts <- NULL
    
    for(i in unique(qdets$sp)){
      for(j in unique(qdets$array)){
        
        temp_df <- rbind(qdets%>%dplyr::select(record_type,sp,array,geometry),
                         pos_df%>%dplyr::select(record_type,sp,array,geometry))%>%
          filter(record_type=="detection")%>%
          st_transform(utm)%>%
          filter(sp==i,array==j)
        
        if(nrow(temp_df)>0){
          
          point_counts <- rasterize(temp_df%>%dplyr::select(geometry), raster_empty, fun = "count")%>%
            projectRaster(.,crs=latlong)%>%
            crop(.,extent(sab))%>%
            mask(.,as_Spatial(sab))
          
          qdet_rasts[[paste(i,j,sep="_")]] <- point_counts
          
        } #nend if check
      } # end j loop
    } # end i loop
    
   
    rast_df_list <- lapply(seq_along(qdet_rasts), function(i) {
                          rast <- qdet_rasts[[i]]
                          rast_df <- as.data.frame(rast,xy=TRUE)
                          rast_df$Raster <- names(qdet_rasts)[i]
                          return(rast_df)
                            })%>%
                      do.call(rbind, .)%>%
                      separate(Raster,c("sp","array"),sep="_",remove = FALSE)
    
    ggplot() +
      geom_sf(data=sab)+
      geom_raster(aes(x = x, y = y, fill = layer),data=rast_df_list%>%filter(sp=="Atlantic cod")) +
      geom_sf(data=recievers_sf_simple%>%filter(array=="2015-2020"),size=0.5,pch=19)+
      facet_wrap(~ Raster) +
      coord_sf(crs = st_crs(sab))+
      scale_fill_viridis_c(na.value = "transparent") +  # Example color scale
      theme_bw()
    
    
    plot_rast_df <- data.frame(sp=c("Atlantic cod","Snow crab", "Atlantic eel","Atlantic salmon","Atlantic sturgeon","Atlantic halibut",
                                    "Bluefin tuna","Blue shark","Grey seal","Leatherback sea turtle","Mako shark","Porbeagle shark","White shark"),
                               group=c(rep("short",6),rep("long",7)),
                               first=c("first",rep("other",5),"first",rep("other",6)),
                               bottom_axis=c(rep("no",5),"yes",rep("no",6),"yes"))%>%# to flag the one that has the top facet label
                    mutate(plot_name = paste0(tolower(gsub(" ","_",sp)),"_plot"))
    
    for(i in c("short","long")){
      
      target_sp <- plot_rast_df%>%filter(group==i)%>%pull(sp)
      
      temp_df <- rast_df_list%>%
                 filter(sp %in% target_sp)%>%
                 left_join(.,plot_rast_df%>%dplyr::select(sp,first))
      
      for(j in target_sp){
        
        temp_plot_df <- temp_df%>%
                        filter(sp == j)
        
        p1 <-  ggplot() +
                geom_sf(data=sab_zones,fill=NA,lty=2)+
                geom_raster(aes(x = x, y = y, fill = layer),data=temp_plot_df) +
                geom_sf(data=recievers_lines,lwd=0.25,alpha=0.8)+
                geom_sf(data=sab_banks,alpha=0.1)+
                scale_x_continuous(breaks=seq(-59.6,-58.4,0.4))+
                scale_y_continuous(breaks=seq(45.8,46.4,0.2))+
                facet_grid(sp ~ array) +
                coord_sf(crs = st_crs(sab))+
                scale_fill_viridis_c(na.value = "transparent") +  # Example color scale
                theme_bw()+
                theme(strip.background = element_rect(fill="white"),
                      axis.title = element_blank(),
                      legend.position = "none")
        
        # if(!(plot_rast_df%>%filter(sp==j)%>%pull(first) == "first")){p1 <- p1 + theme(strip.background.x = element_blank(),
        #                                                                               strip.text.x = element_blank())}
        # 
        # if(plot_rast_df%>%filter(sp==j)%>%pull(bottom_axis) == "no"){p1 <- p1 + theme(axis.text.x=element_blank())}
        
        assign(plot_rast_df%>%filter(sp==j)%>%pull(plot_name),p1)
        
        
    } #end sp j loop
      
    } #end of i group loop
    
    
    #assemble plots based on plot_rast_df
    
    short_plot1 <- (atlantic_cod_plot+theme(axis.text.x=element_blank()))/
                  (snow_crab_plot+theme(strip.background.x = element_blank(),
                                        strip.text.x = element_blank(),
                                        axis.text.x=element_blank()))/
                  (atlantic_halibut_plot+theme(strip.background.x = element_blank(),
                                           strip.text.x = element_blank()))
    
    short_plot2 <- (atlantic_salmon_plot+theme(axis.text.x=element_blank()))/
                  (atlantic_sturgeon_plot+theme(strip.background.x = element_blank(),
                                        strip.text.x = element_blank(),
                                        axis.text.x=element_blank()))/
                  (atlantic_eel_plot+theme(strip.background.x = element_blank(),
                                           strip.text.x = element_blank()))
    
   long_plot1 <- (blue_shark_plot+theme(axis.text.x=element_blank()))/
                 (mako_shark_plot+theme(strip.background.x = element_blank(),
                                        strip.text.x = element_blank(),
                                        axis.text.x=element_blank()))/
                 (porbeagle_shark_plot+theme(strip.background.x = element_blank(),
                                        strip.text.x = element_blank(),
                                        axis.text.x=element_blank()))/
                 (white_shark_plot+theme(strip.background.x = element_blank(),
                                               strip.text.x = element_blank()))
    
    long_plot2 <- (bluefin_tuna_plot+theme(axis.text.x=element_blank()))/
                  (grey_seal_plot+theme(strip.background.x = element_blank(),
                                                strip.text.x = element_blank(),
                                                axis.text.x=element_blank()))/
                  (leatherback_sea_turtle_plot+theme(strip.background.x = element_blank(),
                                           strip.text.x = element_blank()))
                    
     
    #save plots
    ggsave("output/Acoustic/detectionraster_fisheries_sp.png",short_plot1,height=8,width=6,units="in",dpi=300)
    ggsave("output/Acoustic/dtectionraster_anadramous_sp.png",short_plot2,height=8,width=6,units="in",dpi=300)
    
    ggsave("output/Acoustic/detectionraster_sharks.png",long_plot1,height=10,width=6,units="in",dpi=300)
    ggsave("output/Acoustic/detectionraster_lg_plegics.png",long_plot2,height=8,width=6,units="in",dpi=300)
   
    