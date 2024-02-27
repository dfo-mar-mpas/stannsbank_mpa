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

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load St Anns Bank MPA shapefile 
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

#load bioregion and get plotting limits----
bioregion <- read_sf("data/Shapefiles/MaritimesPlanningArea.shp")%>%
             st_transform(latlong)

plot_lim <- bioregion%>% #bounding area for the inset (zoom out)
            st_transform(utm)%>%
            st_buffer(50)%>%
            st_transform(latlong)%>%
            st_bbox()

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

#bathymetry 
    bbsab <- st_bbox(sab)
    
    noaabathy <- getNOAA.bathy(bbsab[1]-1,bbsab[3]+1,bbsab[2]+1,bbsab[4]-1,resolution = 0.25,keep=T)
    
    sab_bathy <- noaabathy%>%
      fortify.bathy() %>%
      st_as_stars() %>%
      st_set_crs(4326)%>%
      st_transform(latlong)
    
    load('data/Bathymetry/250m_stars_object.RData')
    dem_sab <- raster("data/Bathymetry/sab_dem.tif")

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
                  filter(year>2014)%>%#just within the SAB window
                  
  
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
      
      plot_boundaries <- c(c(-60,-58,45.6,46.6))
      
      inset_box <- sab%>%st_bbox()
      inset_box[c(1,3,2,4)] <- plot_boundaries     
      inset_box <- inset_box%>%st_as_sfc()%>%st_as_sf()
      
#make primary map
      
      recievers_sf <- otn_stations%>%
                      filter(collectioncode=="SABMPA")%>%
                      mutate(array = ifelse(year<2021,"2015-2020","2020-2024"),
                             lon = st_coordinates(.)[,1],
                             lat = st_coordinates(.)[,2])%>%
                      st_transform(proj4string(dem_sab))%>%
                      mutate(depth = round(raster::extract(dem_sab,as_Spatial(.)),1))#extract depth
      
      write.csv(recievers_sf%>%data.frame()%>%dplyr::select(-geometry),file="data/Acoustic/Acoustic_Stations.csv",row.names=FALSE)
      
      
      primary_plot <- ggplot()+
                      geom_sf(data=shelfbreak,col="grey20",lty=2,fill=NA)+
                      geom_sf(data=otn_stations,size=0.1)+
                      geom_sf(data=basemap)+
                      geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
                      geom_sf(data=recievers_sf,aes(fill=array),shape=21)+
                      geom_sf(data=geoserver_tag_releases,col="red")+
                      theme_bw()+
                      theme(legend.position = "bottom")+
                      guides(shape = guide_legend(override.aes = list(size=6)))+
                      labs(fill="")+
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
      
      inset_plot <- ggplot()+
                    geom_sf(data=bioregion,fill=NA)+
                    geom_sf(data=otn_stations,size=0.05,pch=3)+
                    geom_sf(data=basemap_inset)+
                    geom_sf(data=coast_hr)+
                    geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
                    geom_sf(data=coast_hr)+
                    #geom_text_repel(data=inset_labels,aes(label=lab,x=lon,y=lat),nudge_x=0.1,nudge_y=0.1)+
                    geom_sf(data=otn_stations%>%filter(collectioncode %in% c("HFX","CBS")),col="red")+
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
      