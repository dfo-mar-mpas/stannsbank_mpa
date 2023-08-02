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
                                             "Mud and Seapens",
                                             "Gravelly sand/mud <50% cobbles/gravel",
                                             "Till >50% cobbles/gravel",
                                             "Till with coraline algae",
                                             "Gravel with crinoids",
                                             "Sand with Sand dollars"))

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
  
#load the Eastern Shore Islands Area of Interest
  esi <- read_sf("data/Shapefiles/EasternShoreIslands_networksite.shp")%>%
         st_transform(latlong)

#load coastline and make basemap ------
    coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")
    
  #sab map coast
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
    
  #esi map coast
    bounding_area_esi <- esi%>%
      st_transform(utm)%>%
      st_buffer(dist = 25)%>%
      st_bbox()%>%
      st_as_sfc()%>%
      st_as_sf()%>%
      st_transform(latlong)%>%
      st_bbox()%>%
      st_as_sfc()%>%
      st_as_sf()%>%
      suppressWarnings()%>%
      suppressMessages()
    
    basemap_esi <- coast_hr%>%
      st_intersection(bounding_area_esi%>%st_transform(st_crs(coast_hr)))%>%
      st_transform(latlong)%>%
      suppressWarnings()
    
#load stations ---------
    cam_locations <- read.csv("data/Camera/SAB_Camera_targets.csv")%>%
                     rename(lon=Longitude,lat=Latitude)%>%
                     st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                     mutate(type="Proposed camera transect location",
                            station=Camera_WP_name)%>%
                     rename(classn=Benthoscape)%>%
                     st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))#new coordinates from Laura                     
    
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
                        rbind(.,cam_locations%>%dplyr::select(names(edna_stations)))%>%
                        suppressWarnings()
    
    
    #Do the random stratfied sampling for each benthoscape classification within the MPA zones that aren't captured by presumptive camera locations
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
      
      set.seed(23) #this will ensure that the analysis is repeatable. Comment out this code to re-sample the MPA
      
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
                            dplyr::select(longitude,latitude,type,name,Zone)
    
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
    
    plot_df <- strat_samples_format%>%
               mutate(type=case_when(grepl("camera",type) ~ "Camera/eDNA",
                                     TRUE ~ type))
    
    sample_plot <- ggplot()+
                    geom_sf(data=basemap)+
                    geom_sf(data=buffered_benthoscape,aes(fill=classn),show.legend = FALSE)+
                    geom_sf(data=sab_zones,fill=NA)+
                    geom_sf(data=perley_buffers,fill=NA,lty=2)+
                    geom_sf(data=strat_samples_format,aes(shape=type),fill="grey90",col="black",size=1.5)+
                    #geom_sf(data=plot_df%>%filter(type=="Camera/eDNA"),aes(shape=type),fill="grey90",col="black",size=3)+
                    coord_sf(xlim=st_bbox(edna_boundaries)[c(1,3)],ylim=st_bbox(edna_boundaries)[c(2,4)],expand=0)+
                    scale_shape_manual(values=c(21,22,25,19))+
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
    
    station_output <- strat_samples_format%>%
                              dplyr::select(Zone,name,type,latitude,longitude,geometry)%>%
                       st_transform(proj4string(dem_sab))%>%
                       mutate(depth = round(raster::extract(dem_sab,as_Spatial(.)),1))%>%
                       st_transform(latlong)%>%
                        mutate(lon = st_coordinates(.)[,1]*-1,
                              lat = st_coordinates(.)[,2],
                              deg_lon = floor(lon),
                              deg_lat = floor(lat),
                              dm_lon = (lon-deg_lon)*60,
                              dm_lat = (lat-deg_lat)*60,
                              deg_lon = deg_lon*-1,
                              activity = case_when(grepl("camera",type) ~ "Camera/eDNA",
                                                   TRUE ~ "eDNA"))%>%
                       arrange(Zone,name,type)
    
    station_df <- station_output%>%
                  mutate(longitude=st_coordinates(.)[,1],
                         latitude=st_coordinates(.)[,2])%>%
                  data.frame()%>%
                  mutate(id = paste(longitude,latitude,sep="-"))%>%
                  left_join(.,cam_locations%>% #merge with the camera locations
                              data.frame()%>%
                              mutate(id = paste(lon,lat,sep="-"))%>%
                              dplyr::select(station,id))%>%
                  mutate(NA_id = cumsum(is.na(station)), #create station names for the eDNA stations
                         station = ifelse(is.na(station),paste("eDNA",NA_id,sep="_"),station),
                         location = "St Anns Bank")%>%
                  dplyr::select(-geometry)
   
    ##station grouping --- 
    
    ##south west group (upper)
    sw_group_a <- station_df%>%
                filter(latitude<45.98,latitude>45.85,longitude<(-59.4))%>%
                mutate(group="Western cluster 2")%>%
                st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    sw_group_b <- station_df%>%
                  filter(latitude<45.85,longitude<(-59.4))%>%
                  mutate(group="Western cluster 1")%>%
                  st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    ##north west group
    nw_group <- station_df%>%
                filter(latitude<46.09,Zone == "2a",!station %in% sw_group$station)%>%
                mutate(group="Western cluster 3")%>%
                st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    nw_group <- nw_group%>%
                rbind(., #there is one point in Zone 1 that can be grouped
                      station_df%>%
                      filter(!station%in%c(sw_group$station,nw_group$station))%>%
                      mutate(group="Western cluster 3")%>%
                      st_as_sf(coords=c("longitude","latitude"),crs=latlong)%>%
                      slice(nw_group%>%st_union()%>%st_nearest_feature(station_df%>%
                                                                         filter(!station%in%c(sw_group$station,nw_group$station))%>%
                                                                         st_as_sf(coords=c("longitude","latitude"),crs=latlong))))
    
    ##far field outside of benthoscape cluster (best appraoched from Sydney)
    out_bentho <- station_df%>%
                  filter(latitude>46.09,Zone == "2a")%>%
                  mutate(group="Western cluster 4")%>%
                  st_as_sf(coords=c("longitude","latitude"),crs=latlong)
  
    #scattarie Bank cluster
    scatterie_bank <- station_df%>%
                      filter(Zone == "1",
                             !station %in% c(sw_group_a$station,sw_group_b$station,out_bentho$station,nw_group$station),
                             latitude<46.07,longitude<(-59.1))%>%
                      mutate(group="Scattarie Bank")%>%
                      st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    #northeast cluster
    ne_group_a <- station_df%>%
                  filter(Zone == "1",
                         latitude<46.2, 
                         longitude<(-59),
                         !station %in% c(sw_group_a$station,sw_group_b$station,out_bentho$station,nw_group$station,scatterie_bank$station),
                         latitude>46.07)%>%
                  mutate(group="Eastern cluster 3")%>%
                  st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    ne_group_b <- station_df%>%
                  filter(Zone == "1",
                         latitude>46.2 | latitude < 46.2 & longitude>(-59),
                         !station %in% c(sw_group_a$station,sw_group_b$station,out_bentho$station,nw_group$station,scatterie_bank$station),
                         latitude>46.07)%>%
                  mutate(group="Eastern cluster 4")%>%
                  st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    #southeast cluster (top group) - split due to distance from port
    se_group_a <- station_df%>%
                  filter(!station %in% c(sw_group_a$station,sw_group_b$station,out_bentho$station,nw_group$station,scatterie_bank$station,ne_group_a$station,ne_group_b$station),
                         latitude>45.95)%>%
                  mutate(group="Eastern cluster 2")%>%
                  st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    #southeast cluster (bottom group)
    se_group_b <- station_df%>%
                    filter(!station %in% c(sw_group_a$station,sw_group_b$station,out_bentho$station,nw_group$station,scatterie_bank$station,ne_group_a$station,ne_group_b$station),
                           latitude<45.95)%>%
                    mutate(group="Eastern cluster 1")%>%
                    st_as_sf(coords=c("longitude","latitude"),crs=latlong)
    
    grouped_stations <- rbind(sw_group_a,sw_group_b,nw_group,out_bentho,scatterie_bank,ne_group_a,ne_group_b,se_group_a,se_group_b)
    
    grouped_hulls <- grouped_stations%>%
                     group_by(group)%>%
                     summarise()%>%
                     st_convex_hull()
    
    grouped_hull_centroids <- grouped_hulls%>%
                              group_by(group)%>%
                              summarise()%>%
                              st_centroid()
    
    grouped_plot <- ggplot()+
                    geom_sf(data=basemap)+
                    geom_sf(data=sab_zones,fill=NA)+
                    geom_sf(data=perley_buffers,fill=NA,lty=2)+
                    geom_sf(data=strat_samples_format)+
                    geom_sf(data=grouped_hulls,aes(fill=group),alpha=0.5,show.legend = FALSE)+
                    geom_sf(data=grouped_stations,aes(fill=group),pch=21,size=1.5)+
                    coord_sf(xlim=st_bbox(edna_boundaries)[c(1,3)],ylim=st_bbox(edna_boundaries)[c(2,4)],expand=0)+
                    scale_shape_manual(values=c(21,22,25,19))+
                    theme_bw()+
                    guides(fill = guide_legend(override.aes = list(size = 3)))+
                    labs(fill="",shape="")+
                    annotation_scale();grouped_plot
    
    ggsave("output/2023_mission/SAB_grouped_station_Plot.png",grouped_plot,width=10,height=6,units="in",dpi=300)
    
    #for form b - Cruise plan. ----------
    edna_esi <- read.csv("data/SamplingDesign_65by20km.csv")%>%
                mutate(lon = lon *-1,
                       deg_lon = floor(lon),
                       deg_lat = floor(lat),
                       dm_lon = (lon-deg_lon)*60,
                       dm_lat = (lat-deg_lat)*60,
                       deg_lon = deg_lon*-1,
                       activity="eDNA",
                       location="ESI",
                       id=c(16,12,8,4,15,11,7,3,14,10,6,2,13,9,5,1))%>%#reorder because they were done along shore not perpendicular
                arrange(id)%>%
                mutate(station=paste0("ESI_",1:16),
              transect = rep(paste0("Transect_",1:4),each=4),
              Zone=" ")
    
    esi_sf <- edna_esi%>%
              mutate(lon=lon*-1)%>%
              st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
              rename(name=station)%>%
              dplyr::select(name,geometry)
    
    #write output for formb combining eSI and SAB stations. -------------
    
    station_output <- rbind(edna_esi%>%dplyr::select(station,activity,depth,deg_lon,dm_lon,deg_lat,dm_lat,lon,lat),
                            station_df%>%dplyr::select(station,activity,depth,deg_lon,dm_lon,deg_lat,dm_lat,lon,lat))%>%
                      mutate(depth = round(depth,1),
                               lon=lon*-1)%>%
                      left_join(.,grouped_stations%>%data.frame()%>%dplyr::select(station,group))%>%
                      mutate(group = case_when(station %in% paste0("ESI_",1:8) ~ "Day 2",
                                               station %in% paste0("ESI_",9:16) ~ "Day 3",
                                               TRUE ~ group))%>%
                      arrange(group,station)
                      
    write.csv(station_output,
              file="data/eDNA/2023_station_list_formb.csv",row.names=FALSE) ## need to create station names that match to the camera file
    
##distances for form b 
    
    waypoints <- data.frame(name=c("WP1","WP2","WP3"),
                            lon=c(-62.483136,-62.497729,-62.475491),
                            lat=c(44.891778,44.835977,44.803352))%>%
      st_as_sf(coords=c("lon","lat"),crs=latlong)
    
    sheet_harbour_port <- data.frame(lon=-62.504091,lat=44.903582,name="Port of Sheet Harbour")%>%
      st_as_sf(coords=c("lon","lat"),crs=latlong)
    
    waypoints_esi <- rbind(sheet_harbour_port,waypoints)
    
    louisbourg <- data.frame(lon=-59.970314,lat=45.917173,name="Louisbourg")%>%
      st_as_sf(coords=c("lon","lat"),crs=latlong)
    
    waypoints_cb <- rbind(louisbourg,data.frame(name=c("WP1","WP2","WP3"),
                               lon=c(-59.974171,-59.959054,-59.941645),
                               lat=c(45.905090,45.902111,45.900265))%>%
      st_as_sf(coords=c("lon","lat"),crs=latlong))
    
    #Day 2 
    
    day2_stations <- c("ESI_1","ESI_2","ESI_3","ESI_4","ESI_8","ESI_7","ESI_6","ESI_5")
    
    day2 <- rbind(waypoints_esi,
                  esi_sf%>%
                    filter(name %in% day2_stations)%>%
                    mutate(order=match(name,day2_stations))%>%
                    arrange(order)%>%
                    dplyr::select(-order),
                  waypoints_esi[rev(1:nrow(waypoints_esi)),])%>%
      mutate(day="Day 2")
    
    day2_line <- day2%>%
                summarise(do_union = FALSE)%>%
                st_cast("LINESTRING")%>%
                mutate(day="Day 2",
                       group="eDNA")
    
    day2_dist <- day2_line%>%
                  st_length()%>%
                  as.numeric()/1000/1.852 
    
    day_2_plot <- ggplot()+
                  geom_sf(data=basemap_esi)+
                  geom_sf(data=esi,fill=NA)+
                  geom_sf(data=day2)+
                  geom_sf(data=day2_line)
    
    #Day 3
    day3_stations <- c("ESI_13","ESI_14","ESI_15","ESI_16","ESI_12","ESI_11","ESI_10","ESI_9")
    
    day3 <- rbind(waypoints_esi,
                  esi_sf%>%
                    filter(name %in% day3_stations)%>%
                    mutate(order=match(name,day3_stations))%>%
                    arrange(order)%>%
                    dplyr::select(-order),
                  waypoints_esi[rev(1:nrow(waypoints_esi)),])%>%
             mutate(day="Day 3")
    
    day3_line <- day3%>%
                  summarise(do_union = FALSE)%>%
                  st_cast("LINESTRING")%>%
                  mutate(day="Day 3",
                         group="eDNA")
    
    day3_dist <- day3_line%>%
      st_length()%>%
      as.numeric()/1000/1.852 
    
    #eastern shore islands map combo
    esi_points <- rbind(day2,day3)
    esi_lines <- rbind(day2_line,day3_line)
    
    esi_plot <- ggplot()+
                geom_sf(data=basemap_esi,col="black")+
                geom_sf(data=esi,fill=NA)+
                geom_sf(data=esi_lines,aes(col=day),show.legend = FALSE)+
                geom_sf(data=esi_points,aes(fill=day),pch=21,col="black")+
                geom_sf(data=waypoints_esi,pch=19,col="black")+
                coord_sf(expand=0)+
                geom_sf_label(data=sheet_harbour_port,aes(label=name),nudge_x = 0.02,nudge_y = 0.06)+
                theme_bw()+
                labs(fill="",x="",y="")+
                theme(legend.position = c(0.9,0.1),
                      legend.title = element_blank())
    
    ggsave("output/2023_mission/per_2023_767_esi.png",esi_plot,width=7,height=5,units="in",dpi=300)