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

noaabathy <- getNOAA.bathy(bbsab[1]-1,bbsab[3]+1,bbsab[2]+1,bbsab[4]-1,resolution = 0.25,keep=T) %>%
              fortify.bathy() %>%
              st_as_stars() %>%
              st_set_crs(4326)%>%
              st_transform(latlong)

#load coastline and make basemap ------
coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")

#load the reciever coordinates
recievers <- read.csv("data/Acoustic/OTN_redesign_coords.csv")%>%
             st_as_sf(coords=c("long","lat"),crs=latlong)

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
                  mutate(year=year(deploy_date))%>%
                  filter(year>2014) #just within the SAB window
  
  url <- "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Oceans_Act_Marine_Protected_Areas/MapServer/0"
  
  MPAs <- get_spatial_layer(url) %>% 
    st_make_valid() %>% 
    st_crop(xmin=proj_long_low,
            ymin=proj_lat_low,
            xmax=proj_long_upp,
            ymax=proj_lat_upp)%>%
    st_transform(latlong)
    

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
      
      load('data/Bathymetry/250m_stars_object.RData')
      
      recievers_sf <- otn_stations%>%
                     filter(collectioncode=="SABMPA")%>%
                     mutate(array = ifelse(year<2021,"2015-2020","2020-2024"))
      
      primary_plot <- ggplot()+
                      geom_sf(data=shelfbreak,col="grey20",lty=2,fill=NA)+
                      geom_sf(data=otn_stations,size=0.1)+
                      geom_sf(data=basemap)+
                      geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
                      geom_sf(data=recievers_sf,aes(fill=array),shape=21)+
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
        