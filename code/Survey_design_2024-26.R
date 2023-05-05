## map of recievers in ESI 2022

#load libraries ----
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)
library(ggplot2)
library(raster)
library(stars)
library(geosphere)
library(ggspatial)
library(ggnewscale)
library(viridis)
library(patchwork)
library(tmaptools)
library(patchwork)

s2_as_sf = FALSE

source("r:/Science/CESD/HES_MPAGroup/R/Functions/box_grid.R")

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load bathymetry ------
dem_sab <- raster("data/Bathymetry/sab_dem.tif") # SAB based on the 'predSpace' data
dem_esi <- raster("data/Bathymetry/esi_dem.tif") # ESI based on the Greenlaw 35 dem
load("data/Bathymetry/250m_stars_object.RData") # STARS object based on the 250m depth contour extracted (converted) from GEBCO

#simplify the names for COIP application
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

#load benthoscape classification
sab_benthoscape <- read_sf("data/Shapefiles/benthoscape.shp")%>%
                   st_transform(latlong)%>%
                   st_make_valid()%>%
                   left_join(.,benthoscape_classes)


#load St Anns Bank MPA shapefile 
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
              st_transform(latlong)

sab  <- sab_zones%>%
        st_transform(utm)%>%
        st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
        st_union()%>% #gets rid of the zones
        st_transform(latlong)%>%
        st_as_sf()

#load coastline and make basemap ------
    coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")
    
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
    
    plot_boundaries <- c(c(-61,-58,45.6,46.6))

#waypoints for leaving Cape Breton -------------
    
    km_nm <- 1.852 #coefficent to go from nautical miles to km
    nm_km <- 1/km_nm #coefficient to go from km to nautical miles
    perley_speed <- 9.3 #knots - roughly the max speed. 
    perley_transit_time <- 5 #time permitted for traveling to station. 

      waypoints_syd <- data.frame(name=c("sydney","WP1","WP2","WP3","WP4"),
                              lon=c(-60.200447,-60.211195,-60.220671,-60.181529,-60.042532),
                              lat=c(46.141232,46.156192,46.191570,46.276581, 46.315330))%>%
                      mutate(location="sydney")%>%
                      st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)
        
      waypoints_lb <- data.frame(name=c("louisbourg","WP1","WP2","WP3"),
                                 lon=c(-59.970314,-59.974171,-59.959054,-59.941645),
                                 lat=c(45.917173,45.905090,45.902111,45.900265))%>%
                       mutate(location="louisbourg")%>%
                       st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)
      
      #distance to outer waypoint
      dist_outer_syd <- as.numeric(waypoints_syd%>%summarise(do_union=FALSE)%>%st_cast("LINESTRING")%>%st_length()/1000*nm_km)
      dist_outer_lb <- as.numeric(waypoints_lb%>%summarise(do_union=FALSE)%>%st_cast("LINESTRING")%>%st_length()/1000*nm_km)
      
      #time to outer waypoint
      time_outer_syd <- dist_outer_syd/perley_speed
      time_outer_lb <- dist_outer_lb/perley_speed
      
      
      #max distance (km) from outer way point based on perley_speed, perley_transit_time and the amount of time to get from the dock to the outer waypoint
      maxdist_syd <- (perley_transit_time-time_outer_syd)*perley_speed*km_nm
      maxdist_lb <- (perley_transit_time-time_outer_lb)*perley_speed*km_nm
      
      #Sydney operational window
      perley_buffer_syd <- waypoints_syd[nrow(waypoints_syd),]%>%
                            st_transform(utm)%>% #km based planar projection
                            st_buffer(dist=maxdist_syd)%>%
                            st_as_sf()%>%
                            st_transform(latlong)%>%
                            mutate(name="Sydney")
      
      #Louisbourg operational window
      perley_buffer_lb <- waypoints_lb[nrow(waypoints_lb),]%>%
                          st_transform(utm)%>% #km based planar projection
                          st_buffer(dist=maxdist_lb)%>%
                          st_as_sf()%>%
                          st_transform(latlong)%>%
                          mutate(name="Louisbourg")
      
      #combined and intersected with the MPA
      perley_buffers <- rbind(perley_buffer_syd,perley_buffer_lb)%>%
                        group_by(name)%>%
                        summarise(geometry=st_intersection(geometry,sab))%>%
                        st_make_valid()%>%
                        ungroup()%>%
                        st_as_sf()
      
      #the total operational window from both ports
      perley_net_buffer <- rbind(perley_buffer_syd,perley_buffer_lb)%>%
                           st_union()%>%
                           st_as_sf()%>%
                           mutate(buffer="St Anns Bank MPA")
                    
        
      #plot it out - shows that for the most part Louisbourg is the closer of the two stations in terms of steam time, but that all of the Sydny                 
      buffer_plot <- ggplot()+
                    geom_sf(data=basemap)+
                    geom_sf(data=waypoints_lb[1,])+
                    geom_sf(data=waypoints_syd[1,])+
                    geom_sf(data=sab,fill=NA)+
                    geom_sf(data=perley_buffers,aes(fill=name),alpha=0.6)+
                    scale_fill_manual(values=c("cornflowerblue","coral"))+
                    theme_bw()+
                    coord_sf(xlim=plot_boundaries[1:2],ylim=plot_boundaries[3:4],expand=0)+
                    labs(fill="")+
                    theme(legend.position = "bottom",
                          plot.margin=grid::unit(c(0,0,0,0), "mm"));buffer_plot
      
      ggsave("output/2023_mission/2023_buffer_plot.png",buffer_plot,width=(get_asp_ratio(sab)*6)+0.5,height=6,dpi=300)
      
      #load stations ---------
      
      #load the CAMPOD locations 
      photo_coords <- read_sf("data/shapefiles/Kenchington2023_SAB_MPA_Photo_Archive_Metadata.shp")%>%
                      st_transform(latlong)%>%
                      group_by(Station)%>%
                      summarise(geometry = st_union(geometry)) %>%
                      st_centroid()%>% # get the centroid of point clusters
                      st_as_sf()%>%
                      mutate(lon=st_coordinates(.)[,1],
                             lat=st_coordinates(.)[,2],
                             type="Campod photos")%>%
                      rename(station=Station)%>%
                      suppressWarnings() #taking centriod of non-planar coordinates is ok since these are just approximate
      
      video_coords <- read_sf("data/shapefiles/Kenchington2023_SAB_MPA_Video_Transect_Metadata.shp")%>%
                      st_transform(latlong)%>%
                      group_by(Station)%>%
                      summarise(geometry = st_union(geometry)) %>%
                      st_centroid()%>% # get the centroid of point clusters
                      st_as_sf()%>%
                      mutate(lon=st_coordinates(.)[,1],
                             lat=st_coordinates(.)[,2],
                             type="Campod videos")%>%
                      rename(station=Station)%>%
                      dplyr::select(station,type,lon,lat,geometry)%>%
                      suppressWarnings() #taking centriod of non-planar coordinates is ok since these are just approximate
      
      camera_locations <- rbind(photo_coords%>%
                                  filter(!station %in% unique(video_coords$station)), #take the centriod from the videos and not the photos 
                                video_coords)%>%
                          st_intersection(.,sab_benthoscape%>%dplyr::select(classn,geometry))%>%
                          st_intersection(.,perley_net_buffer)%>% #cut out the ones not in the operational window of the Perley
                          st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))
      
      camera_locations_df <- data.frame(camera_locations)%>%dplyr::select(-geometry)%>%
                              rename(latitude=lat,longitude=lon,name=station)%>%
                              mutate(type="Camera drop locaiton")
            
#load the 2022 stations of the fixed design from 2022 ----
    edna_stations <- read.csv("data/sampling_2022_coords.csv")%>%
                     filter(!grepl("SR_",station),!grepl("MR_",station),!grepl("ESI",station))%>% #these are the receiver coordinates
                     st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                     mutate(type=case_when(grepl("trans",station) ~ "eDNA depth transect",
                                           TRUE ~ "eDNA zonal comparision"))%>%
                     dplyr::select(station,lon,lat,type)%>%
                     st_intersection(.,sab_benthoscape%>%dplyr::select(classn,geometry))%>%
                     st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
                     filter(type == "eDNA depth transect") # the 2022 design with the structured zonal comparisons are being covered via the random stratified samples. 
    
    edna_stations_df <- edna_stations%>%
                        data.frame()%>%
                        rename(latitude=lat,longitude=lon,name=station)%>% # to fit the COIP input template
                         dplyr::select(-geometry)

    #cut out the benthoscape areas that are too far to get too
    buffered_benthoscape <- sab_benthoscape%>%
                            st_intersection(.,perley_buffer_lb)%>%
                            st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
                            group_by(Zone,classn)%>%
                            summarise(geometry=st_union(geometry))%>%
                            ungroup()%>%
                            st_as_sf()%>%
                            st_make_valid()

    #benthoscape area constrained by operation distance of the Perley (from Louisbourg) and intersected with the MPA Zones
    sab_benthoscape_zone <- buffered_benthoscape%>%
                            dplyr::select(classn,geometry)%>%
                            st_intersection(.,perley_buffer_lb)%>%
                            st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
                            st_make_valid()
    
    #Identify the area of each benthoscape classication within the MPA within each zone. This will be used to conduct the stratified sampling. 
    sab_zone_bentho_areas <- sab_benthoscape_zone%>%
                              group_by(Zone,classn)%>%
                              summarise(geometry=st_union(geometry),
                                        area=as.numeric(st_area(geometry)/1000/1000))%>%
                              ungroup()%>%
                              data.frame()%>%
                              dplyr::select(-geometry)%>%
                              arrange(Zone);sab_zone_bentho_areas
    
  #Identify zones within each area with zone 2b being the smallest with the fewest zones. 
    
    #which classes are in each zone
    three_zone_classes <- Reduce(intersect,split(sab_zone_bentho_areas$classn,sab_zone_bentho_areas$Zone))
    
    #which zones are shared amongst Zones 2a and 1 that aren't in zone 2b (smallest zone)
    two_zone_classes <- Reduce(intersect,split(sab_zone_bentho_areas[sab_zone_bentho_areas$Zone!="2b","classn"],
                                               sab_zone_bentho_areas[sab_zone_bentho_areas$Zone!="2b","Zone"]))%>%
      setdiff(.,three_zone_classes)
    
    #class just in zone 1
    zone1_classes <-  setdiff(sab_zone_bentho_areas[sab_zone_bentho_areas$Zone=="1","classn"],c(two_zone_classes,three_zone_classes))
    
    #check to see if any are missing
    if(length(unique(sab_zone_bentho_areas$classn)) == length(unique(c(zone1_classes,two_zone_classes,three_zone_classes)))){message("Good no missing classes!")} #no missing zones
    
    #compile a sampling scheme
    
    sample_size <- 3 #target sample size
    
    #this is our target sampling frame so that each benthoscape class within each zone is sampled 'sample_size' times. 
    strat_sampling <-   data.frame(classn=c(zone1_classes,rep(two_zone_classes,each=2),rep(three_zone_classes,each=3)),
                                   Zone = c("1",rep(c("1","2a"),2),rep(c("1","2a","2b"),3)))%>%
                          mutate(sample_size=sample_size)
    
    
    #Exiting stations - locations we will go too or want to get too (Campod) for another question (depth transect) 
    existing_targets <- edna_stations%>%
      filter(Zone == "1")%>%
      st_intersection(.,perley_buffers%>%dplyr::select(geometry))%>%
      rbind(.,camera_locations%>%dplyr::select(names(edna_stations)))%>%
      suppressWarnings()
    
    
    #Do the random stratfied sampling for each benthoscape classification wihtin the MPA zones that aren't captured by presumptive camera locations
    #or locations required to repeat the depth transect radiating from Scatterie Bank. 
    strat_samples <- NULL
    
    for(i in unique(strat_sampling$Zone)){
      
      message(paste0("Working on Zone ",i))
      
      if(i %in% existing_targets$Zone){
        
        temp_existing <- existing_targets%>%
          filter(Zone == i)%>%
          pull(classn)%>%
          table()%>%
          data.frame()%>%
          rename(classn=1)
        
        temp <- strat_sampling%>%
          filter(Zone == i)%>%
          left_join(temp_existing)%>%
          mutate(Freq = ifelse(is.na(Freq),0,Freq),
                 sample_size = sample_size - Freq,
                 sample_size = ifelse(sample_size<0,0,sample_size),
                 classn=factor(classn,levels=unique(buffered_benthoscape$classn)))%>%
          arrange(classn)
      }else{temp <- strat_sampling%>%
        filter(Zone == i)%>%
        mutate(classn=factor(classn,levels=unique(buffered_benthoscape$classn)))%>%
        arrange(classn)}
      
      set.seed(23) #this will ensure that the analysis is repeatable. COmment out this code to resample the MPA
      
      strat_out <- buffered_benthoscape%>%
        filter(Zone==i)%>%
        st_sample(size=temp$sample_size,by_polygon=TRUE,exact=T,type="random")%>% #take 5 of each minus the locations of campod
        st_as_sf()%>%
        st_intersection(.,sab_benthoscape)%>%
        mutate(lon=st_coordinates(.)[,1],
               lat=st_coordinates(.)[,2],
               type="St Anns Bank Stratified eDNA")%>%
        dplyr::select(lon,lat,type,classn)%>%
        st_as_sf()%>%
        rename(geometry=x)%>%
        rbind(.,existing_targets%>%
                filter(Zone==i)%>%
                #mutate(type="St Anns Bank Stratified eDNA - existing")%>%
                dplyr::select(lon,lat,type,classn,geometry))%>%
        group_by(classn)%>%
        mutate(name=paste(unique(classn),1:n(),sep="_"))%>%
        ungroup()%>%
        arrange(name)%>%
        st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
        suppressWarnings()
      
      strat_samples <- rbind(strat_samples,strat_out)
      
    }
    
    
    strat_samples_format <- strat_samples%>%
                            rename(latitude=lat,longitude=lon)%>%
                            mutate(name=case_when(grepl("depth",type) ~ paste(name,"DT",sep=" - "), #DT - 'depth transect"
                                                  grepl("Campod",type) ~ paste(name,"CP",sep=" - "), #CP - Campod
                                                  TRUE ~ name),
                                   type=ifelse(grepl("Campod",type),"Campod",type))%>%
                            dplyr::select(name,longitude,latitude,type)
    
    strat_samples_df <- strat_samples_format%>%data.frame()
    
    #make the plot
    
    #boundaries for zoomed in - inset
    edna_boundaries <- strat_samples_format%>%
                       st_bbox()%>%
                       st_as_sfc()%>%
                       st_as_sf()%>%
                       st_transform(utm)%>%
                       st_buffer(10)%>%
                       st_transform(latlong)%>%
                       st_bbox()%>% #make it square
                       st_as_sfc()%>%
                       st_as_sf()
    
    sample_plot <- ggplot()+
                    geom_sf(data=basemap)+
                    geom_sf(data=buffered_benthoscape,aes(fill=classn))+
                    geom_sf(data=sab_zones,fill=NA)+
                    geom_sf(data=perley_buffers,fill=NA,lty=2)+
                    geom_sf(data=strat_samples_format,aes(shape=type),fill="white",col="grey20",size=1.5)+
                    geom_sf(data=strat_samples_format%>%filter(type=="Campod"),aes(shape=type),fill="white",col="grey20",size=3)+
                    coord_sf(xlim=st_bbox(edna_boundaries)[c(1,3)],ylim=st_bbox(edna_boundaries)[c(2,4)],expand=0)+
                    scale_shape_manual(values=c(21,25,19))+
                    theme_bw()+
                    #scale_fill_viridis(discrete = T,option="C")+
                    labs(fill="",shape="",title="Proposed sample locations")
    
    outplot <- ggplot()+
                geom_sf(data=basemap)+
                geom_sf(data=buffered_benthoscape,aes(fill=classn),show.legend = FALSE)+
                geom_sf(data=sab_zones,fill=NA)+
                geom_sf(data=edna_boundaries,fill=NA)+
                geom_sf_text(data=waypoints_lb[1,],label="Louisbourg",nudge_x=-0.2,nudge_y=0.05)+
                geom_sf(data=waypoints_lb[1,],pch=19,size=3)+
                geom_sf_text(data=waypoints_syd[1,],label="Sydney",nudge_x=0.1,nudge_y=0.05)+
                geom_sf(data=waypoints_syd[1,],pch=19,size=3)+
                theme_bw()+
                coord_sf(xlim=plot_boundaries[c(1,2)],ylim=plot_boundaries[c(3,4)],expand=0)+
                labs(title="Study Area",x="",y="")
      
    
    study_plot <- outplot + sample_plot + plot_layout(guides = "collect",ncol = 2) & theme(legend.position = "bottom")       
    
    
    ggsave("output/2023_mission/sab_survey_design_2024-26.png",study_plot,width=10,height=4,units="in",dpi=300)
    
    
## now format the final station list ---------
    
    station_output <- rbind(strat_samples_format%>%
                              mutate(type=ifelse(grepl("Campod",type),"St Anns Bank Stratified eDNA",type))%>%
                              dplyr::select(latitude,longitude,name,type,geometry),
                            camera_locations%>%
                              rename(name=station,latitude=lat,longitude=lon)%>%
                              mutate(type="Camera drop")%>%
                              dplyr::select(latitude,longitude,name,type,geometry))%>%
                       st_transform(proj4string(dem_sab))%>%
                       mutate(depth = round(raster::extract(dem_sab,as_Spatial(.)),1))%>%
                       arrange(type,name)%>%
                       st_transform(latlong)
    
    #for COIP
    write.csv(x=station_output%>%data.frame()%>%dplyr::select(-c(geometry,depth)),file="output/2023_mission/2024-26_StationList.csv",row.names = FALSE)
    
    #for activity approval
    write.csv(x=station_output%>%
                st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
                data.frame()%>%
                rename(station=name)%>%
                mutate(latitude = round(latitude,2),
                       longitude = round(longitude,2))%>%
                arrange(type,Zone,station)%>%
                dplyr::select(Zone,type,station,depth,longitude,latitude),file="output/2023_mission/2023_StationList_activityapproval.csv",row.names = FALSE)
    
    
    edna_stations_ESI <- read.csv("data/sampling_2022_coords.csv")%>%
                        filter(grepl("ESI",station))%>% 
                        mutate(type="ESI eDNA",
                               name=paste(type,1:n(),sep="-"))%>%
                        rename(latitude=lat,longitude=lon)%>%
                        dplyr::select(latitude,longitude,name,type)
    
    write.csv(x=rbind(station_output,edna_stations_ESI),file="output/2023_mission/COIP_2024-26_StationList.csv")
    
    
#make a map of the proposed survey locations for the 2023 Activity approval  ----------
    
    sample_plot_AA <- ggplot()+
                      geom_sf(data=shelfbreak,col="grey20",lty=2,fill=NA)+
                      geom_sf(data=basemap)+
                      geom_sf(data=buffered_benthoscape,aes(fill=classn))+
                      geom_sf(data=sab_zones,fill=NA)+
                      geom_sf(data=strat_samples_format,aes(shape=type),fill="white",col="grey20",size=1.5)+
                      geom_sf(data=strat_samples_format%>%filter(type=="Campod"),aes(shape=type),fill="white",col="grey20",size=3)+
                      coord_sf(xlim=st_bbox(edna_boundaries)[c(1,3)],ylim=st_bbox(edna_boundaries)[c(2,4)],expand=0)+
                      scale_shape_manual(values=c(21,25,19))+
                      theme_bw()+
                      #scale_fill_viridis(discrete = T,option="C")+
                      labs(fill="",shape="",title="Proposed sample locations")
    
    ggsave("output/2023_mission/proposed_sample_locations.png",sample_plot_AA,width=7,height=5,units="in",dpi=300)

##extra code  -------------

#create a random stratified sampling of the MPA based on benthoscape to structure the COIP request
# 
# #set the seed for the random sampling so that this will generate the same results each time
# set.seed(23)
# 
# #make a data-frame that counts the benthoscape classifications associated with the general camera locations so that each class as a count including 0s
# camera_strata <- camera_locations%>%
#   pull(classn)%>%
#   table(.)%>%
#   data.frame()%>%
#   rename(classn=1)%>%
#   rbind(.,data.frame(classn=setdiff(unique(buffered_benthoscape$classn),unique(camera_locations$classn)),Freq=0))%>%
#   mutate(classn=factor(classn,levels=unique(buffered_benthoscape$classn)))%>% #arrange them the same way that the stratification by polygon conducts the sampling
#   arrange(classn)
# 
# strat_samples <- buffered_benthoscape%>%
#   #st_sample(size=c(5,5),by_polygon=TRUE,exact=T,type="random")%>%
#   st_sample(size=5-camera_strata$Freq,by_polygon=TRUE,exact=T,type="random")%>% #take 5 of each minus the locations of campod
#   st_as_sf()%>%
#   st_intersection(.,sab_benthoscape)%>%
#   mutate(lon=st_coordinates(.)[,1],
#          lat=st_coordinates(.)[,2],
#          type="St Anns Bank Stratified eDNA")%>%
#   dplyr::select(lon,lat,type,classn)%>%
#   st_as_sf()%>%
#   rename(geometry=x)%>%
#   rbind(.,camera_locations%>%
#           mutate(type="St Anns Bank Stratified eDNA - camera")%>%
#           dplyr::select(lon,lat,type,classn,geometry))%>%
#   group_by(classn)%>%
#   mutate(name=paste(unique(classn),1:5,sep="_"))%>%
#   ungroup()%>%
#   arrange(name)%>%
#   st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))



