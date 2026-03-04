#Figure 1 for the Acoustic Report

#load libraries ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(lubridate)
library(arcpullr)
library(marmap)
library(stars)
library(ggrepel)
library(ggspatial)
library(ggnewscale)
library(worrms)
library(patchwork)
library(raster)
library(rphylopic)
library(ggspatial)
library(viridis)
library(terra)
library(tidyterra)
library(MarConsNetData)

s2_as_sf = FALSE

#source functions
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/trim_img_ws.R")
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/marmap_to_isobath.R")
source("code/Acoustic/create_gate_raster.R")

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load shapefiles -----

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

#Canadian EEZ 
can_eez <- read_sf("data/Shapefiles/can_eez.shp")%>%st_transform(latlong)
us_eez <- read_sf("data/Shapefiles/US EEZ/eez.shp")%>%st_transform(latlong)

#load bioregion and get plotting limits 
bioregion <- data_planning_areas()%>%
             st_transform(latlong)

mar_net <- read_sf("data/Shapefiles/networksites_proposed_OEM_MPA_20220617.shp")%>%
            st_transform(latlong)

plot_lim <- bioregion%>% #bounding area for the inset (zoom out)
            st_transform(utm)%>%
            st_buffer(10)%>%
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

      #for sab plot 
      # plot_boundaries <- c(c(-60,-58,45.6,46.6)) #for inset plot
      # inset_box <- sab%>%st_bbox()
      # inset_box[c(1,3,2,4)] <- plot_boundaries     
      # inset_box <- inset_box%>%st_as_sfc()%>%st_as_sf()

plot_boundaries <- sab%>%
                   st_transform(utm)%>%
                   st_buffer(5)%>%
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


bbox_vals <- c( #these are the values of scale4
                xmin = -85- 2,
                ymin = 4- 2,
                xmax = -47+ 2,
                ymax = 52.3+ 2)%>%
              st_bbox(crs=latlong)%>%
              st_as_sfc()

global_basemap <- ne_states(returnclass = "sf")%>%
  st_transform(latlong)

focal_countries <- global_basemap%>%
                   st_make_valid()%>%
                   st_intersection(bbox_vals)%>%
                   group_by(admin) %>% 
                   st_make_valid()%>%
                   summarise(geometry = st_union(geometry), .groups = "drop")


#bathymetry ----
    bbsab <- st_bbox(sab)
    
    #extract bathymetric data
    
    # noaabathy_sab <- getNOAA.bathy(bbsab[1]-1,bbsab[3]+1,bbsab[2]-1,bbsab[4]+1,resolution = 0.25,keep=T)
    # noaabathy_plotregion <- getNOAA.bathy(plot_lim[1]-1,plot_lim[3]+1,plot_lim[2]-1,plot_lim[4]+1,resolution = 0.25,keep=T)
    #  
    # cont_250_sab <- marmap_to_isobath(noaabathy_sab,-250,latlong)
    # cont_250_plotregion <- marmap_to_isobath(noaabathy_plotregion,-250,latlong)
    # 
    # write_sf(cont_250_sab,dsn="data/Bathymetry/contour_250_sab.shp")
    # write_sf(cont_250_plotregion,dsn="data/Bathymetry/contour_250_plotregion.shp")
    
    cont_250_sab <- read_sf("data/Bathymetry/contour_250_sab.shp")
    cont_250_plotregion <- read_sf("data/Bathymetry/contour_250_plotregion.shp")
    
    #create xyz dataframes that can be used by geom_contour
    # isobath_sab_df <- as.xyz(noaabathy_sab)%>%
    #                   rename(lon=1,lat=2,depth=3)
    # 
    # isobath_df <- as.xyz(noaabathy_plotregion )%>%
    #   rename(lon=1,lat=2,depth=3)
    
    load("data/Bathymetry/gebco_plotregion.RData") #this is a trimmed GEBCO for the plot region (+ 100 km buffer)
    
    
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
  proj_end <- ymd("20260101")
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
                         array = ifelse(year<2021,"2015-2020","2020-2025"))%>%
                  filter(year>2014)#just within the SAB window
                  
  
  can_mcn <- read_sf("data/Shapefiles/All_GoC_MCAs.shp")%>%
             st_transform(latlong)
  
  MPAs <- can_mcn%>%filter(grepl("Marine Protected Area",TYPE_E))
    
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
  # tag_df <- read.csv("data/Acoustic/SABMPA_tagrelease_Dec2023.csv")%>%
  #           rename(common=COMMON_NAME_E)%>%
  #           mutate(id = paste(abs(RELEASE_LONGITUDE),RELEASE_LATITUDE,common,sep="_"),
  #                  year=year(as.POSIXct(UTC_RELEASE_DATE_TIME)),
  #                  common=case_when(common == "ATLANTIC COD" ~ "Atlantic cod",
  #                                   common == "ATLANTIC HALIBUT" ~ "Atlantic halibut",
  #                                   common == "Atlantic Striped wolffish" ~ "Atlantic wolffish",
  #                                   TRUE ~ common),
  #                  array = ifelse(year<2021,"2015-2020","2020-2024"))%>% #only need unique locations
  #           distinct(id,.keep_all=TRUE)%>%
  #           st_as_sf(coords=c("RELEASE_LONGITUDE","RELEASE_LATITUDE"),crs=latlong)%>%
  #           dplyr::select(common,year,array,geometry)%>%
  #           rbind(.,geoserver_tag_releases%>% #add the snow crab data
  #                   filter(common == "Snow crab")%>%
  #                   rename(year=yearcollected)%>%
  #                   mutate(array = ifelse(year<2021,"2015-2020","2020-2024"))%>%
  #                   dplyr::select(common,year,array,geometry))
  
  
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
  
      
#make primary map
      recievers_sf <- otn_stations%>%
                      filter(collectioncode=="SABMPA")%>%
                      mutate(array = ifelse(year<2021,"2015-2020","2020-2025"),
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
      
      #high res extraction of bathymetry (combo of multibeam and GEBCO)
      
      sab_dem_hr <- terra::rast("data/Bathymetry/TotalSABbathy/wholeBathy.tif")
      
      rast_proj <- st_crs(sab_dem_hr)
      
      sab_dem_hr_trimmed <- sab_dem_hr%>%mask(.,sab%>%st_transform(rast_proj))

## Figure 1 -----------    
      primary_plot <- ggplot()+
                      geom_spatraster(data = sab_dem_hr_trimmed)+
                      geom_sf(data=cont_250_sab,col="grey35",linewidth=0.5)+
                      geom_sf(data=MPAs%>%filter(NAME_E == "St. Anns Bank Marine Protected Area"),fill=NA,linewidth=1.1,col="black")+
                      geom_sf(data=sab_banks%>%st_transform(rast_proj),fill=NA,col="black")+
                      geom_sf(data=otn_stations%>%filter(!grepl("_SB",station_name)),size=1.5,fill = "grey70", colour="black",shape=21,
                              stroke = 0.1)+
                      scale_fill_viridis(na.value = "transparent")+
                      geom_sf(data=coast_hr,fill="grey70")+
                      labs(fill="Depth (m)")+
                      new_scale_fill()+
                      geom_sf(data=recievers_sf%>%filter(!grepl("_SB",station_name)),aes(fill=array),col="black",shape=21,size=2,show.legend = FALSE)+
                      scale_fill_manual(values=c("2015-2020" = "#D73027","2020-2025"="white"))+
                      annotation_scale(location="bl")+
                      theme_bw()+
                      theme(legend.position = "inside",
                            legend.position.inside = c(0.12,0.85),
                            legend.background = element_blank(),
                            axis.text=element_blank())+
                      guides(shape = guide_legend(override.aes = list(size=6)))+
                      labs(x="",y="")+
                      coord_sf(expand=0,xlim=plot_boundaries[c(1,3)],ylim=plot_boundaries[c(2,4)])
      
      ggsave("output/Acoustic/primary_plot.png",primary_plot,height=7.5,width=7.5,units="in",dpi=300)
      trim_img_ws("output/Acoustic/primary_plot.png")
      
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
      
      inset_otn_stations <- otn_stations%>%
                            #filter(station_type %in% c("Acoustic")%>%
                            st_difference(sab)
     
      
      inset_plot <- ggplot()+
                    geom_sf(data=can_eez,fill=NA)+
                    geom_sf(data=cont_250_plotregion,color = "grey80",linewidth = 0.6)+
                    #geom_spatraster_contour(data = gebco_region,breaks = -250,color = "grey80",linewidth = 0.6)+
                    #geom_contour(data=isobath_df,aes(x=lon,y=lat,z=depth),breaks=-250,color = "grey80", size = 0.5)+
                    geom_sf(data=mar_net_df%>%filter(TYPE == "TBD"),fill="grey70",alpha=0.25,lty=3)+
                    geom_sf(data=mar_net_df%>%filter(TYPE == "MR"),fill="grey10",alpha=0.3,lty=3)+
                    geom_sf(data=mar_net_df%>%filter(TYPE == "AOI"),fill="salmon2",alpha=0.3,lty=3)+
                    #geom_sf(data=otn_stations,size=0.25,pch=20)+
                    geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
                    geom_sf(data=recievers_sf%>%filter(!grepl("_SB",station_name)),aes(fill=array),shape=21,size=1.1,stroke = 0.1)+
                    scale_fill_manual(values=c("2015-2020" = "#D73027","2020-2025"="white"))+
                    #geom_sf(data=otn_stations%>%filter(collectioncode %in% c("HFX","CBS")),col="red",size=0.7)+
                    geom_sf(data=inset_otn_stations%>%filter(!grepl("_SB",station_name)),size=1.1,fill = "grey70", colour="black",shape=21,
                            stroke = 0.1)+
                    geom_sf(data=basemap_inset,fill="grey70")+
                    geom_sf(data=plot_boundaries%>%st_as_sfc(),fill=NA)+ 
                    coord_sf(expand=0,xlim=plot_lim[c(1,3)],ylim=plot_lim[c(2,4)])+
                    theme_bw()+
                    labs(x="",y="")+
                    theme(
                          legend.position = "none")
      
      ggsave("output/Acoustic/inset_plot.png",inset_plot,height=6,width=6,units="in",dpi=300)
      trim_img_ws("output/Acoustic/inset_plot.png")
      
      #global inset ... inset --
      
      global_scale <- st_bbox(
        c(
          xmin = -85.00000,
          ymin =  4.00000,
          xmax = -47.00000,
          ymax = 52.25748
        ),
        crs = latlong)
      
      global_inset <- ggplot()+
        geom_sf(data=can_eez,fill=NA,colour="black")+
        geom_sf(data=focal_countries)+
        geom_sf(data=basemap_inset%>%filter(country=="Canada"),fill="grey70")+
        geom_sf(data=cont_250_plotregion,color = "grey80",linewidth = 0.6)+
        geom_sf(data=plot_lim%>%st_as_sfc(),fill=NA)+
        geom_sf(data=mar_net_df%>%filter(TYPE == "TBD"),fill="grey70",alpha=0.25)+
        geom_sf(data=mar_net_df%>%filter(TYPE == "MR"),fill="grey10",alpha=0.3)+
        geom_sf(data=mar_net_df%>%filter(TYPE == "AOI"),fill="salmon2",alpha=0.3)+
        geom_sf(data=MPAs,fill="cornflowerblue",alpha=0.5)+
        theme_bw()+
        coord_sf(xlim=global_scale[c(1,3)],ylim=global_scale[c(2,4)],expand=0)+
        annotation_north_arrow(location = "tl",
                               height = unit(0.7, "cm"),
                               width  = unit(0.7, "cm"))
        
        
        ggsave("output/Acoustic/global_inset.png",global_inset,height=6,width=5,units="in",dpi=300)
        trim_img_ws("output/Acoustic/global_inset.png")
      
      #globe version of the map
        
        #download the world globe basemap
        world_globe <- ne_countries(
          scale = "medium",
          returnclass = "sf"
        ) %>%
          st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
          st_transform(globe_crs)
        
        #Define a global region that is centered on the study area
        
        center_pt <- plot_lim%>%
                     st_as_sfc()%>%
                     st_centroid()
        
        lon0 <- st_coordinates(center_pt)[1]
        lat0 <- st_coordinates(center_pt)[2]
        
        globe_crs <- sprintf("+proj=ortho +lat_0=%s +lon_0=%s",lat0, lon0)
        
        #define the box denoting the study region you want to highlight
        global_box <- plot_lim%>%
                      st_as_sfc() %>%
                      st_transform(globe_crs)
        
        #make a circle to wrap the globe plot
        
        globe_circle <- st_sfc(
          st_buffer(
            st_point(c(0, 0)),   # center of orthographic projection
            dist =  6378137  # meters
          ),
          crs = globe_crs
        )
        
        #crudgy way to make it so that the oceans are white in the plot
        globe_disc <- st_sfc(
          st_point(c(0, 0)),  # center in projected coords
          crs = globe_crs
        ) %>%
          st_buffer(dist = 6378137) %>%   # Earth radius in meters
          st_as_sf()
       
        
        global_inset2 <- ggplot() +
          
          #global background
          geom_sf(data = globe_disc, fill = "white", colour = "black", linewidth = 0.4)+
          
          # Land
          geom_sf(
            data = world_globe,
            colour = "grey20",
            linewidth = 0.2
          ) +
          
          geom_sf(
            data = world_globe%>%filter(formal_en == "Canada"),
            fill = "grey70",
            colour = "grey20",
            linewidth = 0.2
          ) +
          
          # Study region box
          geom_sf(
            data = global_box,
            fill = NA,
            colour = "black",
            linewidth = 0.9
          ) +
          
          geom_sf(data = globe_circle,
                  fill = NA,
                  colour = "grey30",
                  linewidth = 0.4)+
          
          coord_sf(crs = globe_crs) +
          
          theme_void() +
          
          theme(
            panel.background = element_rect(fill = NA, colour = NA),
            plot.background  = element_rect(fill = NA, colour = NA)
          )
        
        ggsave(
          "output/Acoustic/global_inset2.png",
          plot = global_inset2,
          width = 4,
          height = 4,
          dpi = 600,
          bg = "transparent"
        )
      
#tag maps ---------------
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
      
      tag_df <- read.csv("data/Acoustic/CJFAS paper/SAB_filt_Q25_releaselocs1.csv")%>%
                rename(Common.Name=common_name)%>%
                filter(!is.na(mean_lon))%>%
                st_as_sf(coords=c("mean_lon","mean_lat"),crs=latlong)
      
      #Colour scale to match other plots
      tag_colours <- read.csv("data/Acoustic/CJFAS paper/tag_palette.csv")%>%
                     mutate(Common.Name = tag_df%>%arrange(Common.Name)%>%distinct(Common.Name)%>%pull(Common.Name))%>%
                     rename(colour=2)%>%
                     dplyr::select(Common.Name,colour)
      
      fill_vals <- setNames(tag_colours$colour,
                            tag_colours$Common.Name)
      
    # tag_df <- read.csv("data/Acoustic/qdet_releaselocs.csv")%>%
    #           rename(lon=RELEASE_LONGITUDE,lat=RELEASE_LATITUDE)%>%
    #           filter(!is.na(lon))%>%
    #           st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
    #           mutate(group=case_when(Common.Name %in% c("White shark","Shortfin mako","Blue shark","Blue shark/Mako/Bluefin tuna","Porbeagle shark") ~ "Sharks",
    #                                  Common.Name %in% c("Atlantic sturgeon","Bluefin tuna")~"Pelagics",
    #                                  Common.Name %in% c("Atlantic cod","Snow crab","Atlantic halibut") ~ "Atlantic fisheries",
    #                                  Common.Name %in% c("Atlantic salmon","American eel") ~ "Diadramous",
    #                                  Common.Name %in% c("Grey seal","Leatherback turtle") ~ "Large pelagics",
    #                                  TRUE ~ NA),
    #                  Common.Name = factor(Common.Name,levels=c("American eel","Atlantic cod","Atlantic halibut","Atlantic salmon","Atlantic sturgeon",
    #                                                            "Blue shark","Blue shark/Mako/Bluefin tuna","Bluefin tuna","Grey seal","Leatherback turtle",           
    #                                                            "Porbeagle shark","Shortfin mako","Snow crab","White shark")))
    #           
      
    tag_bound <- tag_df%>%
                 st_bbox()%>%
                 st_as_sfc()%>%
                 st_buffer(0.5)%>% #0.5 degree buffer
                 st_bbox()
    
    #Set up the scales -----
    
    #scale 0
    scale0 <- sab%>%
              st_transform(utm)%>%
              st_buffer(7)%>%
              st_transform(latlong)%>%
              st_bbox()
    
    scale0_box <- scale0%>%st_as_sfc()%>%st_as_sf()
    #scale 1
    scale1 <- sab%>%
      st_transform(utm)%>%
      st_buffer(300)%>%
      st_transform(latlong)%>%
      st_bbox()
    
    scale1[2] <- 44
    
    scale1_box <- scale1%>%st_as_sfc()%>%st_as_sf()
    
    #scale 2
    scale2 <- sab%>%
      st_transform(utm)%>%
      st_buffer(650)%>%
      st_transform(latlong)%>%
      st_bbox()
    
    scale2_box <- scale2%>%st_as_sfc()%>%st_as_sf()
    
    #scale 3
    scale3 <- sab%>%
      st_bbox()
    
    scale3[c(3,4)] <- scale2[c(3,4)]
    scale3[3] <- -51
    scale3[2] <- 38.5
    scale3[1] <- -78
    
    scale3_box <- scale3%>%st_as_sfc()%>%st_as_sf()
    
    #scale 4 
    scale4 <- sab%>%
      st_bbox()
    
    scale4[c(3,4)] <- scale3[c(3,4)]
    scale4[3] <- -47
    scale4[2] <- 4
    scale4[1] <- -85
    
    #download the bathymetry data
    # noaabathy_fine <- getNOAA.bathy(scale1[1]-1,scale1[3]+1,scale1[2]-1,scale1[4]+1,resolution = 0.25,keep=T)
    # noaabathy_med <- getNOAA.bathy(scale2[1]-1,scale2[3]+1,scale2[2]-1,scale2[4]+1,resolution = 0.5,keep=T)
    # noaabathy_large <- getNOAA.bathy(scale4[1]-1,scale4[3]+1,scale4[2]-1,scale4[4]+1,resolution = 1,keep=T)
    # 
    # cont_250_fine <- marmap_to_isobath(noaabathy_fine,-250,latlong)
    # cont_250_med <- marmap_to_isobath(noaabathy_med,-250,latlong)
    # cont_250_large <- marmap_to_isobath(noaabathy_large,-250,latlong)
    # 
    # write_sf(cont_250_fine,dsn="data/Bathymetry/contour_250_fine.shp")
    # write_sf(cont_250_med,dsn="data/Bathymetry/contour_250_medium.shp")
    # write_sf(cont_250_large,dsn="data/Bathymetry/contour_250_large.shp")
    
    cont_250_fine <- read_sf("data/Bathymetry/contour_250_fine.shp")
    cont_250_med <- read_sf("data/Bathymetry/contour_250_medium.shp")
    cont_250_large <- read_sf("data/Bathymetry/contour_250_large.shp")
      
    
    #make plots  
    
    #full scale - 
    scale0_pts <- tag_df%>%st_intersection(scale0_box)
    
    scale0_plot <- ggplot()+
                    geom_sf(data=cont_250_sab,col="grey20",lty=2,fill=NA)+
                    geom_sf(data=basemap_inset)+
                    geom_sf(data=sab_zones,fill="cornflowerblue",alpha=0.3)+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=3,pch=21)+
                    coord_sf(expand=0,xlim=scale0[c(1,3)],ylim=scale0[c(2,4)])+
                    theme_bw()+
                    scale_fill_manual(values=tag_colours$colour)
                    annotation_scale(location="br")+
                    theme(legend.position = "none",
                          axis.text=element_blank())
    
    
    #scale 1 (fine scale)
    
    scale1_pts <- tag_df%>%
                  st_intersection(scale1_box)%>%
                  filter(!animal_id %in% scale0_pts$animal_id)
    
    scale1_plot <- ggplot()+
                    geom_sf(data=basemap_inset)+
                    geom_sf(data=cont_250_fine,col="grey60",linewidth=0.6)+
                    geom_sf(data=sab_zones,fill="cornflowerblue",alpha=0.3)+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=1,pch=21)+
                    geom_sf(data=tag_df%>%filter(animal_id %in% scale1_pts$animal_id),aes(fill=Common.Name),size=4,pch=21)+
                    coord_sf(expand=0,xlim=scale1[c(1,3)],ylim=scale1[c(2,4)])+
                    theme_bw()+
                    annotation_scale(location="br")+
                    theme(legend.position = "none",
                          axis.text=element_blank())
    
    #scale 2 fine-medium scale 
    
    scale2_pts <- tag_df%>%
                  st_intersection(.,scale2_box)%>%
                  filter(!animal_id %in% c(scale0_pts$animal_id,scale1_pts$animal_id))
    
    scale2_plot <- ggplot()+
                    geom_sf(data=basemap_inset)+
                    geom_sf(data=cont_250_med,col="grey60",linewidth=0.6)+
                    geom_sf(data=sab_zones,fill="cornflowerblue",alpha=0.3)+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=1,pch=21)+
                    geom_sf(data=tag_df%>%filter(animal_id %in% scale2_pts$animal_id),aes(fill=Common.Name),size=3,pch=21)+
                    coord_sf(expand=0,xlim=scale2[c(1,3)],ylim=scale2[c(2,4)])+
                    theme_bw()+
                    annotation_scale(location="br")+
                    theme(legend.position = "none",
                          axis.text=element_blank())
    
    #scale 3 medium-large scale
   
    
    scale3_pts <- tag_df%>%
                  st_intersection(.,scale3_box)
                
    scale3_plot <- ggplot()+
                   geom_sf(data=can_eez,fill=NA,colour="black")+
                   geom_sf(data=basemap_inset%>%filter(country!="Canada"))+
                   geom_sf(data=basemap_inset%>%filter(country=="Canada"),fill="grey70")+
                   geom_sf(data=cont_250_large,col="grey65",linewidth=0.4)+
                   geom_sf(data=sab,fill="cornflowerblue",alpha=0.3)+
                   geom_sf(data=tag_df,aes(fill=Common.Name),size=4,pch=21)+
                   #geom_sf(data=tag_df%>%filter(animal_id %in% scale3_pts$animal_id),aes(fill=Common.Name),size=4,pch=21)+
                   coord_sf(expand=0,xlim=scale3[c(1,3)],ylim=scale3[c(2,4)])+
                   scale_fill_manual(values=fill_vals)+
                   theme_bw()+
                   annotation_scale(location="br")+
                   theme(legend.position = "none",
                         axis.text=element_blank())
    
    #scale 3 with labels ----
    
    tag_proj <- st_transform(tag_df, 3347)
    
    near_list <- st_is_within_distance(tag_proj, dist = 20000)
    
    # Convert neighbour list to cluster IDs
    cluster_id <- rep(NA_integer_, length(near_list))
    current_id <- 1
    
    for (i in seq_along(near_list)) {
      if (is.na(cluster_id[i])) {
        members <- unlist(near_list[i])
        cluster_id[members] <- current_id
        current_id <- current_id + 1
      }
    }
    
    tag_clustered <- tag_proj %>%
      mutate(cluster = cluster_id) %>%
      group_by(cluster) %>%
      slice(1) %>%
      ungroup()%>%
      st_transform(latlong)%>%
      st_intersection(scale3_box)%>%
      mutate(Project = gsub("ACT.","",Project_full),
             Project = gsub("FACT.","",Project_full),
             Project = sub(" -.*", "",Project))
    
    
    scale3_plot_labs <- ggplot()+
                        geom_sf(data=basemap_inset %>% filter(country!="Canada"))+
                        geom_sf(data=basemap_inset %>% filter(country=="Canada"),fill="grey70")+
                        geom_sf(data=sab,fill="cornflowerblue",alpha=0.3)+
                        geom_sf(data=tag_df,aes(fill=Common.Name),
                                                  size=4,pch=21)+
                        geom_label_repel(
                          data = tag_clustered,
                          aes(label = Project, geometry = geometry),
                          stat = "sf_coordinates",
                          size = 3,
                          max.overlaps = Inf,
                          box.padding = 0.4,
                          point.padding = 0.3,
                          fill = "white",        # box fill
                          color = "black",       # text colour
                          label.size = 0.2       # border thickness
                        )+
      
                          coord_sf(expand=0,
                                 xlim=scale3[c(1,3)],
                                 ylim=scale3[c(2,4)])+
                        scale_fill_manual(values=fill_vals)+
                        theme_bw()+
                        annotation_scale(location="br")+
                        theme(legend.position = "bottom")+
                        labs(x="",y="")
    
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
    
    # scale4_pts <- tag_df%>%
    #               st_intersection(.,scale4_box)%>%
    #               filter(!animal_id %in% c(scale0_pts$animal_id,scale1_pts$animal_id,scale2_pts$animal_id,scale3_pts$animal_id))
    # 
    scale4_plot <- ggplot()+
                  geom_sf(data=can_eez,fill=NA,colour="black")+
                  geom_sf(data=focal_countries)+
                  geom_sf(data=basemap_inset%>%filter(country=="Canada"),fill="grey70")+
                  geom_sf(data=cont_250_large,col="grey65",linewidth=0.4)+
                  geom_sf(data=sab,fill="cornflowerblue",alpha=0.3,colour="black")+
                  geom_sf(data=tag_df,aes(fill=Common.Name),size=3,pch=21)+
                  #geom_sf(data=tag_df%>%filter(animal_id %in% scale4_pts$animal_id),aes(fill=Common.Name),size=3,pch=21)+
                  geom_sf(data=rings,lty=2,lwd=0.5,fill=NA)+
                  scale_fill_manual(values=fill_vals)+
                  coord_sf(expand=0,xlim=scale4[c(1,3)],ylim=scale4[c(2,4)],label_axes = "-NE")+
                  theme_bw()+
                  annotation_scale(location = "br")+
                  annotation_north_arrow(location = "tl",
                                         height = unit(0.7, "cm"),
                                         width  = unit(0.7, "cm"))+
                  theme(legend.position = "none")
    
    #make a legend plot
    legend_dummy <- ggplot()+
                    geom_sf(data=tag_df,aes(fill=Common.Name),size=5,pch=21)+
                    theme_bw()+
                    labs(fill="")+
                    coord_sf(expand=0,label_axes = "-NE")
    
    #save plot - note used Scale 3 and 4 in the end plot to make it simple
    # ggsave("output/Acoustic/tagmap_scale0.png",scale0_plot+scale_fill_viridis(discrete=T),height=5,width=5,units="in",dpi=300)
    # ggsave("output/Acoustic/tagmap_scale1.png",scale1_plot+scale_fill_viridis(discrete=T),height=5,width=6,units="in",dpi=300)
    #ggsave("output/Acoustic/tagmap_scale2.png",scale2_plot+scale_fill_viridis(discrete=T),height=6,width=8,units="in",dpi=300)
    ggsave("output/Acoustic/tagmap_scale3.png",scale3_plot,height=6,width=10,units="in",dpi=300)
    ggsave("output/Acoustic/tagmap_scale3_labs.png",scale3_plot_labs,height=6,width=10,units="in",dpi=300)
    ggsave("output/Acoustic/tagmap_scale4.png",scale4_plot,height=10,width=6,units="in",dpi=300) 
     
    trim_img_ws("output/Acoustic/tagmap_scale3.png")
    trim_img_ws("output/Acoustic/tagmap_scale3_labs.png")
    trim_img_ws("output/Acoustic/tagmap_scale4.png")

#### Interpolated species track plot -----------
    
    #build a raster for the cells -----------
    
    #identify the between receiver spacing
        
    # receiver_spacing <- recievers_sf%>%
    #                     mutate(datetime=as.POSIXct(deploy_date),
    #                            year=year(datetime))%>%
    #                     filter(year %in% c(2018,2023))%>%
    #                     mutate(gate=ifelse(year==2018,"Double","Single"),
    #                            gate2=case_when(lat>46.15~"Double_top",
    #                                            lat<46.15 & gate == "Double" ~ "Double_bottom",
    #                                            TRUE ~ "Single"))%>%
    #                     arrange(gate,gate2,-lon)%>%
    #                     group_by(gate2)%>%
    #                     mutate(station_id = paste(gate2,1:n(),sep="-"),
    #                            .groups = "drop")%>%
    #                     group_by(gate2) %>%
    #                     mutate(
    #                       station_id = paste(gate2, row_number(), sep = "-"),
    #                       spacing_m = as.numeric(st_distance(geometry, lag(geometry), by_element = TRUE))
    #                     ) %>%
    #                     ungroup()%>%
    #                     st_transform(utm)
    # 
    # round(median(receiver_spacing$spacing_m,na.rm=T)/1000,2) #roughly 1 km spacing


    #now create rasters for the cells and plotting. this will be done by creating a linestring converted from the points,
    #sampling along that lines string at even intervals (~ 1 km nominal spacing) and then converting to cells and a raster

    # gate_cells <- NULL
    # gate_rasters <- list()
    # 
    # for(i in unique(receiver_spacing$gate2)){
    # 
    #   input <- receiver_spacing%>%filter(gate2==i)
    # 
    #   gate_single_r <- create_gate_raster(input, spacing_m = 1000, cell_width = 1000)
    # 
    #   gate_cells <- rbind(gate_cells,gate_single_r[[2]]%>%mutate(gate2=i))
    #   gate_rasters[[i]] <- gate_single_r[[1]]
    # 
    # }
    # 
    # gate_cells2 <- gate_cells%>%
    #   left_join(receiver_spacing%>%
    #               st_drop_geometry() %>%
    #               distinct(gate2,.keep_all=TRUE) %>%
    #               dplyr::select(gate,gate2))%>%
    #   st_transform(latlong)%>%
    #   mutate(cell_id = 1:n(),
    #          array=ifelse(gate=="Double","2015-2021","2021-2025"))
    # 
    # write_sf(gate_cells2,"data/Acoustic/CJFAS paper/reciever_cells.shp")
    
    gate_cells <- read_sf("data/Acoustic/CJFAS paper/reciever_cells.shp")
    
    gate_cells_cents <- gate_cells%>%
                        st_centroid()%>%
                        mutate(lon=st_coordinates(.)[,1],
                               lat=st_coordinates(.)[,2])%>%
                        data.frame()%>%
                        dplyr::select(cell_id,lon,lat) #for coordinates
    
    
    #benchmark transition date
    transition_date <- read.csv("data/Acoustic/CJFAS paper/SAB_Q25dets_array_stats.csv")%>%
                       filter(grepl("Double",array_design))%>%
                       mutate(timestamp=as.POSIXct(Latest_Detection),
                              transition=ceiling_date(timestamp, unit = "day") %>% as.Date()) #round up the the next day
    
    #process fish we tagged 
    pos_cod <- read.csv("data/Acoustic/CJFAS paper/cod_final_M25.csv")%>%mutate(sp="Atlantic cod")
    pos_halibut <- read.csv("data/Acoustic/CJFAS paper/hal_final_M25.csv")%>%mutate(sp="Atlantic halibut")
    pos_wolffish <- read.csv("data/Acoustic/CJFAS paper/wolf_final_M25.csv")%>%mutate(sp="Atlantic wolffish")
    
    pos_df_step1 <- rbind(pos_cod,pos_halibut,pos_wolffish)%>%
                filter(glatos_array == "OTN.SABMPA")%>%
               mutate(timestamp=as.POSIXct(detection_timestamp_utc),
                     julian=yday(timestamp),
                     year=year(timestamp),
                     id=gsub("SABMPA-","",animal_id), #individual fish id
                     id=gsub("MPA-SAB-","",id),
                     id=sub("-.*", "", id),
                     id=as.numeric(id),
                     #array=ifelse(year>2020,"2021-2025","2015-2021"),
                     array = ifelse(
                       as.Date(timestamp) >= transition_date$transition,
                       "2021-2025",
                       "2015-2021"
                     ),
                     stn_id=paste(receiver_sn,year,sep="-"))%>%
              st_as_sf(coords=c("deploy_long","deploy_lat"),crs=latlong)
    
    pos_df1 <- pos_df_step1%>%
               filter(array=="2021-2025")%>%
               st_join(gate_cells%>%
                       filter(array=="2021-2025")%>%
                       dplyr::select(cell_id),
                       join=st_within)%>% #now you have a unique cell_id
               data.frame()%>%
               rbind(., #have to rbind so you don't miss assign detections to the wrong yearly grid cell
                     pos_df_step1%>%
                       filter(array=="2015-2021")%>%
                       st_join(gate_cells%>%
                                 filter(array=="2015-2021")%>%
                                 dplyr::select(cell_id),
                               join=st_within)%>% 
                       data.frame())%>%
               left_join(.,gate_cells_cents)%>%
               st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
               mutate(tag_type="SABMPA")
    
    #process qualified detections
    qdets_step1 <- read.csv("data/Acoustic/CJFAS paper/SAB_12hr_events_Q25.csv")%>%
             mutate(timestamp=as.POSIXct(first_detection),
                    TEMPID=1:n(),
                    julian=yday(timestamp),
                    year=year(timestamp),
                    array = ifelse(
                      as.Date(timestamp) >= transition_date$transition,
                      "2021-2025",
                      "2015-2021"
                    ),
                    stn_id=location, #not quite the same as the fish we tagged
                    sp = case_when(grepl("BLUEFIN",animal_id) ~ "Bluefin tuna",
                                   grepl("BLUESH",animal_id) ~ "Blue shark",
                                   grepl("COD",animal_id) ~ "Atlantic cod",
                                   grepl("CRAB",animal_id) ~ "Snow crab",
                                   grepl("EEL",animal_id) ~ "Atlantic eel",
                                   grepl("GSEAL",animal_id) ~ "Grey seal",
                                   grepl("HAL",animal_id) ~ "Atlantic halibut",
                                   grepl("LHBACK",animal_id) ~ "Leatherback sea turtle",
                                   grepl("MAKO",animal_id) ~ "Shortfin mako",
                                   grepl("PORB",animal_id) ~ "Porbeagle shark",
                                   grepl("SALM",animal_id) ~ "Atlantic salmon",
                                   grepl("STURG",animal_id) ~ "Atlantic sturgeon",
                                   grepl("WHITE",animal_id) ~ "White shark",
                                   grepl("MAC",animal_id) ~ "Atlantic mackerel",
                                   TRUE ~ NA))%>%
             st_as_sf(coords=c("mean_longitude","mean_latitude"),crs=latlong,remove=FALSE)
    
    qdets <- qdets_step1%>%
             filter(array=="2021-2025")%>%
             st_join(gate_cells%>%
                       filter(array=="2021-2025")%>%
                       dplyr::select(cell_id),join=st_within)%>%
             data.frame()%>%
             rbind(.,
                   qdets_step1%>%
                     filter(array=="2015-2021")%>%
                     st_join(gate_cells%>%
                               filter(array=="2015-2021")%>%
                               dplyr::select(cell_id),join=st_within)%>%
                     data.frame())%>%
             filter(!is.na(cell_id))%>%
             left_join(.,gate_cells_cents)%>%
             st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
             mutate(tag_type="qdet")
    

    #sum of total daily unique detection per species, per station, over the deployment period
    pos_df_daily <- pos_df1%>%
                  data.frame()%>%
                  dplyr::select(sp,year,cell_id,julian,animal_id)%>%
                  rbind(.,qdets%>%
                          data.frame()%>%
                          dplyr::select(sp,year,cell_id,julian,animal_id))%>%
                  group_by(sp,cell_id,year,julian)%>%
                  summarise(total_count = n_distinct(animal_id),
                            .groups = "drop")%>%
                  left_join(pos_df1%>%
                              data.frame()%>%
                              dplyr::select(cell_id,array)%>%
                              rbind(.,qdets%>%data.frame()%>%dplyr::select(cell_id,array))%>%
                              distinct(cell_id,.keep_all=TRUE))%>%
                  mutate(gate=ifelse(array=="2015-2020","Double","Single"))%>%
                  left_join(gate_cells_cents)%>%
                  st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)
    
    pos_df_unique <- pos_df1%>%
                    data.frame()%>%
                    dplyr::select(sp,year,cell_id,julian,animal_id)%>%
                    rbind(.,qdets%>%
                            data.frame()%>%
                            dplyr::select(sp,year,cell_id,julian,animal_id))%>%
                    group_by(sp,cell_id)%>%
                    summarise(total_count = n_distinct(animal_id),
                              .groups = "drop")%>%
                    left_join(pos_df1%>%
                                data.frame()%>%
                                dplyr::select(cell_id,array)%>%
                                rbind(.,qdets%>%data.frame()%>%dplyr::select(cell_id,array))%>%
                                distinct(cell_id,.keep_all=TRUE))%>%
                    mutate(gate=ifelse(array=="2015-2020","Double","Single"))%>%
                    left_join(gate_cells_cents)%>%
                    st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)
    
    #count raster
    points_extent <- gate_cells%>%
                     st_transform(utm)%>%
                     st_buffer(2)%>% # 2 km buffer
                     st_bbox()%>%
                     st_as_sfc()%>%
                     st_as_sf()
    
     
    #Group and sum by species, cell, and array (if needed by array)
    detection_summary <- pos_df_daily %>%
                        rename(species=sp) %>%
                        st_drop_geometry() %>%
                        group_by(species, array, cell_id ) %>%  
                        summarise(total_detections = sum(total_count, na.rm = TRUE),
                                  .groups = 'drop')
    
    detection_summary_unique <- pos_df_unique%>%
                                rename(species=sp) %>%
                                st_drop_geometry()
    
    
    #make species specific cell data.frames -- total daily unique sum
    cell_df <- list()
    
    for(i in unique(detection_summary$species)){
      for(j in unique(detection_summary$array)){
        
        temp <- detection_summary%>%
                filter(species == i,array==j)%>%
                right_join(.,gate_cells%>%
                             filter(array==j)%>%
                             dplyr::select(cell_id))%>%
                st_as_sf()%>%
                st_transform(latlong)%>%
                mutate(species=i,
                       array=j)
        
        cell_df[[paste(i,j,sep="_")]] <- temp
        
      }# end j loop
    }# end i loop
    
    cell_df <- do.call("rbind",cell_df)
    
    #for unique individuals 
    
    cell_df_unique <- list()
    
    for(i in unique(detection_summary_unique$species)){
      for(j in unique(detection_summary_unique$array)){
        
        temp <- detection_summary_unique%>%
          filter(species == i,array==j)%>%
          right_join(.,gate_cells%>%
                       filter(array==j)%>%
                       dplyr::select(cell_id))%>%
          st_as_sf()%>%
          st_transform(latlong)%>%
          mutate(species=i,
                 array=j)
        
        cell_df_unique[[paste(i,j,sep="_")]] <- temp
        
      }# end j loop
    }# end i loop
    
    cell_df_unique <- do.call("rbind",cell_df_unique)
    
    #set up plotting limits
    grid_lims <- cell_df%>%
                 st_bbox()%>%
                 st_as_sfc()%>%
                 st_transform(utm)%>%
                 st_buffer(1.5)%>%
                 st_transform(latlong)%>%
                 st_bbox()
    
    #Focal species - these are the most common species
    focal_sp <- c("Atlantic cod","Atlantic salmon","Atlantic halibut",
                  "Bluefin tuna","Blue shark","White shark")
    
    label_df <- cell_df %>%
      filter(species %in% focal_sp) %>%
      distinct(species) %>%
      mutate(
        x = grid_lims["xmax"] - 0.01,
        y = grid_lims["ymax"] - 0.01
      )
    
    #make a general map for the inset
    p_raster_inset <- ggplot()+
                      geom_sf(data=sab_zones,fill="white")+
                      geom_sf(data = sab_banks %>% filter(name == "Scatarie Bank"),
                              fill = "grey90")+
                      geom_sf(data=gate_cells,fill=NA,linewidth=0.005)+
                      geom_sf(data=gate_cells%>%st_centroid(),size=0.0005,pch=20)+
                      geom_sf(data=grid_lims%>%st_as_sfc(),linetype=2,fill=NA)+
                      theme_void()+
                      theme(panel.background = element_rect(fill = NA, colour = NA),
                            plot.background  = element_rect(fill = NA, colour = NA))+
                      coord_sf(expand=0)
    
    
    #count range per species
    range_df <- cell_df %>%
      filter(species %in% focal_sp) %>%
      group_by(species) %>%
      summarize(
        min = min(total_detections, na.rm = TRUE),
        max = max(total_detections, na.rm = TRUE),
        range = paste(min, max, sep = "-"),
        .groups = "drop"
      )
    
    #plotting function (this makes it all work together better with patchwork)
    plot_species <- function(sp) {
      
      df <- cell_df %>%
        filter(species == sp,
               !(species == "Atlantic halibut" & array == "2015-2021"))
      
      # Get the range string for this species
      range_label <- range_df %>% filter(species == sp) %>% pull(range)
      
      p <- ggplot(df) +
        geom_sf(data = sab_zones, fill = NA) +
        geom_sf(data = sab_banks %>% filter(name == "Scatarie Bank"),
                fill = "grey90") +
        geom_sf(aes(fill = total_detections)) +
        scale_fill_viridis(
          trans = "log10",
          limits = c(1, 200),
          oob = scales::squish,
          na.value = NA
        ) +
        coord_sf(expand = 0,
                 xlim = grid_lims[c(1,3)],
                 ylim = grid_lims[c(2,4)]) +
        theme_bw() +
        theme(
          legend.position = "none",
          axis.text = element_blank(),
          panel.grid = element_blank()
        ) +
        # Add the species name and range to the title
        labs(title = paste0(sp, " (", range_label, ")"),
             fill = "Σ unique daily individuals")
      
      # Optional scale bar for Bluefin tuna
      if(sp=="Bluefin tuna"){
        p <- p +
          annotation_scale(
            location = "tl",      # top-left
            width_hint = 0.25
          ) 
      }
      
      # If halibut, insert inset
      if (sp == "Atlantic halibut") {
        p <- p +
          inset_element(
            p_raster_inset,
            left   = 0.02,   
            bottom = 0.55,
            right  = 0.55,   
            top    = 0.98    
          )
      }
      
      return(p)
    }
  
    #create the plots for each of the focal spacies
    plots <- lapply(focal_sp, plot_species)
    
    #combine the plot with patchwork
    final_plot <- wrap_plots(plots, ncol = 3) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    #save the outputs
    ggsave("output/Acoustic/Figure4.png",final_plot,width=8,height=6,units="in",dpi=300)
    trim_img_ws("output/Acoustic/Figure4.png")
    
    #now repeat for the unique detections
    range_df_unique <- cell_df_unique %>%
      filter(species %in% focal_sp) %>%
      group_by(species) %>%
      summarize(
        min = min(total_count, na.rm = TRUE),
        max = max(total_count, na.rm = TRUE),
        range = paste0("n=",max), #just the total number
        .groups = "drop"
      )
    
    global_min <- min(range_df_unique$min)
    global_max <- max(range_df_unique$max)
    
    plot_species_unique <- function(sp) {
      
      df <- cell_df_unique %>%
        filter(species == sp,
               !(species == "Atlantic halibut" & array == "2015-2021"))
      
      # Get the range string for this species
      range_label <- range_df_unique %>% filter(species == sp) %>% pull(range)
      
      p <- ggplot(df) +
        geom_sf(data = sab_zones, fill = NA) +
        geom_sf(data = sab_banks %>% filter(name == "Scatarie Bank"),
                fill = "grey90") +
        geom_sf(aes(fill = total_count)) +
        scale_fill_viridis(
          limits = c(global_min, global_max),   # <-- unified across panels
          breaks = c(1, 12, 24, 36, 48),       # adjust to taste
          na.value = NA
        ) +
        coord_sf(expand = 0,
                 xlim = grid_lims[c(1,3)],
                 ylim = grid_lims[c(2,4)]) +
        theme_bw() +
        theme(
          legend.position = "none",
          axis.text = element_blank(),
          panel.grid = element_blank()
        ) +
        labs(title = paste0(sp, " (", range_label, ")"),
             fill = "Total unique individuals")
      
      # Optional scale bar for Bluefin tuna
      if(sp=="Bluefin tuna"){
        p <- p +
          annotation_scale(
            location = "tl",      # top-left
            width_hint = 0.25
          ) 
      }
      
      # If halibut, insert inset
      if (sp == "Atlantic halibut") {
        p <- p +
          inset_element(
            p_raster_inset,
            left   = 0.02,   
            bottom = 0.55,
            right  = 0.55,   
            top    = 0.98    
          )
      }
      
      return(p)
    }
    
    #create the plots for each of the focal spacies
    plots_unique <- lapply(focal_sp, plot_species_unique)
    
    #combine the plot with patchwork
    final_plot_unique <- wrap_plots(plots_unique, ncol = 3) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    
    #save the outputs
    ggsave("output/Acoustic/Figure4_unique_ids.png",final_plot_unique,width=8,height=6,units="in",dpi=300)
    trim_img_ws("output/Acoustic/Figure4_unique_ids.png")
    
    

    