## Survey Design for the 2025 Field Season in St. Anns Bank

#load libraries ----
library(tidyr)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(terra)
library(ggspatial)
library(ggnewscale)
library(viridis)
library(patchwork)

s2_as_sf = FALSE

#Functions ------

ddm_to_dd <- function(ddm_string) {
  # Clean input string
  ddm_string <- trimws(ddm_string)
  
  # Extract components using regex
  # Pattern matches: [degrees]° [decimal minutes]' [direction]
  pattern <- "^\\s*(\\d+)°\\s*(\\d+\\.\\d+)'\\s*([NSEW])\\s*$"
  
  # Check if string matches the expected pattern
  if (!grepl(pattern, ddm_string)) {
    stop("Input must be in format: DD° MM.MMM' D (e.g., '59° 35.186' W')")
  }
  
  # Extract components
  degrees <- as.numeric(gsub(pattern, "\\1", ddm_string))
  decimal_minutes <- as.numeric(gsub(pattern, "\\2", ddm_string))
  direction <- gsub(pattern, "\\3", ddm_string)
  
  # Convert to decimal degrees
  decimal_degrees <- degrees + (decimal_minutes / 60)
  
  # Apply sign based on direction
  if (direction %in% c("W", "S")) {
    decimal_degrees <- -decimal_degrees
  }
  
  return(decimal_degrees)
}


convert_to_ddm <- function(decimal_degrees, is_longitude = TRUE) {
  sign <- ifelse(decimal_degrees < 0, -1, 1)
  abs_dd <- abs(decimal_degrees)
  degrees <- floor(abs_dd)
  minutes <- round((abs_dd - degrees) * 60, 3)
  if(minutes >= 60) {
    degrees <- degrees + 1
    minutes <- 0
  }
  # Determine direction based on coordinate type
  if(is_longitude) {
    direction <- ifelse(decimal_degrees < 0, "W", "E")
  } else {
    direction <- ifelse(decimal_degrees < 0, "S", "N")
  }
  return(sprintf("%d° %06.3f' %s", degrees, minutes, direction))
}

determine_transect_group <- function(lon, threshold) {
  diff_lon <- c(0, diff(sort(lon)))
  cumsum(diff_lon > threshold) + 1
}

generate_points <- function(center, radius, num_points) {
  angles <- seq(0, 2 * pi, length.out = num_points + 1)[-1]  # Exclude the last point as it duplicates the first
  x <- st_coordinates(center)[1] + radius * cos(angles)
  y <- st_coordinates(center)[2] + radius * sin(angles)
  points <- st_sfc(st_multipoint(cbind(x, y)), crs = st_crs(center))
  points_sf <- st_as_sf(data.frame(radius = radius), geom = st_cast(points, "POINT"))
  return(points_sf)
}

create_circle <- function(center, radius) {
  circle <- st_buffer(center, dist = radius)
  return(circle)
}

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load benthoscape 

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


#load St Anns Bank MPA shapefile  -------------
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
              st_transform(latlong)

sab  <- sab_zones%>%
        st_transform(utm)%>%
        st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
        st_union()%>% #gets rid of the zones
        st_transform(latlong)%>%
        st_as_sf()

#load bathymetric data
sab_dem <- rast("data/Bathymetry/sab_dem.tif")%>%project(latlong)

contour_sab <- as.contour(sab_dem, levels = 150)%>% # could also use 237 from 'deep' but just round out to 250
  st_as_sf()%>%
  st_transform(latlong)


#load the banks polygons -------
banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%st_transform(latlong)

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
    

#Load the 2024 SAB field season coordinates. 
    
    edna_coords <- read.csv("output/2024_mission/station_coords.csv")%>%
                   dplyr::select(Station,Depth,Longitude,Latitude)%>%
                   rename(name=Station,depth=Depth,long=Longitude,lat=Latitude)%>%
                   st_as_sf(coords=c("long","lat"),crs=latlong)%>%
                   mutate(sample="eDNA")
    
    #acoustic receivers
    rec_coords_new <- read.csv("data/Acoustic/SABMPA_NEW_2025_deploys.csv")%>%
                      rename(name=Station)%>%
                      st_as_sf(coords=c("Long","Lat"),crs=latlong)%>%
                      mutate(type="New deployments")%>%
                      dplyr::select(name,type)
    
    rec_coords_recov <- read.csv("data/Acoustic/2024_Additional_Station_Design.csv")%>%
                        st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                        mutate(type="Recovery")%>%
                        dplyr::select(name,type)
    
    rec_coords_recov_depth <- cbind(rec_coords_recov%>%pull(name),
                                    terra::extract(sab_dem,rec_coords_recov%>%st_transform(st_crs(sab_dem)))*-1)
    
    array_coords <- read.csv("data/Acoustic/OTN_redesign_coords.csv")%>%
                    st_as_sf(coords=c("long","lat"),crs=latlong)%>%
                    mutate(name=paste("SAB_",Id),
                           type="Existing array")%>%
                    dplyr::select(name,type)
    
    rov_recov_coords <- read.csv("data/Acoustic/2025_rov_stations.csv")%>%
                        st_as_sf(coords=c("long","lat"),crs=latlong)%>%
                        dplyr::select(name,type)
    
    acoustic_coords <- rbind(rec_coords_new,rec_coords_recov,rov_recov_coords,array_coords)%>%
                       mutate(sample="Acoustic recievers")
    
    
    
#evaluate operational buffer for 2024 season  -------------
    
    km_nm <- 1.852 #coefficent to go from nautical miles to km
    nm_km <- 1/km_nm #coefficient to go from km to nautical miles
    perley_speed <- 8.9 #knots - roughly the max speed with some buffer
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
    
    
#Create a geographic stratification of the MPA  ----------------
    

    
#make main map figure for formb
    
    curdo_plot_bound <- edna_coords%>%
                        filter(grepl("Curdo",name))%>%
                        st_transform(utm)%>%
                        st_buffer(0.5)%>%
                        st_transform(latlong)%>%
                        st_bbox()
    
    curdo_box <- curdo_plot_bound%>%st_as_sfc()
    
    point_scaler <- 0.5
                      
    p1_curdo <- ggplot()+
                geom_sf(data=banks%>%filter(name=="Curdo Bank"),fill="#F8766D")+
                geom_sf(data=edna_coords%>%filter(grepl("Curdo",name)),aes(fill=sample),shape=21,size=2.5*point_scaler)+
                geom_sf(data=acoustic_coords%>%mutate(sample="AR: New deployments"),aes(fill=sample),shape=24,size=2.5*point_scaler)+
                theme_bw()+
                theme(axis.text=element_blank(),
                      plot.margin = margin(0, 0, 0, 0),
                      plot.title = element_text(size=10,vjust=-2),
                      # legend.position = "inside",
                      # legend.position.inside = c(0.80,0.9),
                      # legend.title = element_blank(),
                      # legend.background = element_blank(),
                      legend.position = "none")+
                coord_sf(xlim=curdo_plot_bound[c(1,3)],ylim=curdo_plot_bound[c(2,4)],expand=0)+
                labs(title="Curdo Bank",fill="")+
                annotation_scale(location="br")+
                scale_fill_manual(values=c("AR: Recovery" = "aquamarine3",
                                           "AR: Existing array" = "grey80",
                                           "AR: New deployments"="blue2",
                                           "eDNA" = "white"))
  
    
    
    #Sactarie Bank
    scatarie_coords <- rbind(
                      edna_coords%>%
                        filter(!grepl("Cam",name))%>%
                        dplyr::select(sample),
                      acoustic_coords%>%
                        mutate(sample=paste0("AR: ",type))%>%
                        dplyr::select(sample)
                    )
    
    
    scatarie_plot_bound <- banks%>%
                          filter(name=="Scatarie Bank")%>%
                           st_transform(utm)%>%
                           st_buffer(3)%>%
                           st_transform(latlong)%>%
                           st_bbox()
    
    scatarie_box <- scatarie_plot_bound%>%st_as_sfc()
    
    p1_scatarie <- ggplot()+
                    geom_sf(data=banks%>%filter(name=="Scatarie Bank"),fill="cornflowerblue")+
                    geom_sf(data=scatarie_coords%>%filter(sample=="AR: Existing array"),aes(fill=sample),shape=23,size=2.5*point_scaler)+
                    geom_sf(data=scatarie_coords%>%filter(sample=="AR: New deployments"),aes(fill=sample),shape=24,size=2.5*point_scaler)+
                    geom_sf(data=scatarie_coords%>%filter(sample=="AR: Recovery"),aes(fill=sample),shape=25,size=2.5*point_scaler)+
                    geom_sf(data=scatarie_coords%>%filter(sample=="eDNA"),aes(fill=sample),shape=21,size=2.5*point_scaler)+
                    geom_sf(data=scatarie_coords%>%filter(sample=="AR: ROV recovery"),aes(fill=sample),shape=21,size=3.25*point_scaler)+              
                    theme_bw()+
                    theme(axis.text=element_blank(),
                          plot.margin = margin(0, 0, 0, 0),
                          plot.title = element_text(size=10,vjust=-2),
                          # legend.position = "inside",
                          # legend.position.inside = c(0.80,0.9),
                          # legend.title = element_blank(),
                          # legend.background = element_blank(),
                          legend.position = "none")+
                    coord_sf(xlim=scatarie_plot_bound [c(1,3)],ylim=scatarie_plot_bound [c(2,4)],expand=0)+
                    labs(title="Scatarie Bank",fill="")+
                    annotation_scale(location="br")+
                    scale_fill_manual(values=c("AR: Recovery" = "aquamarine3",
                                               "AR: Existing array" = "grey80",
                                               "AR: New deployments"="blue2",
                                               "eDNA" = "white",
                                               "AR: ROV recovery" = "red"))
    
    
    sab_points <- scatarie_coords%>%
                  mutate(sample=ifelse(grepl("AR",sample),"Acoustic receiver",sample))
    
      p1 <- ggplot()+
              geom_sf(data=basemap)+
              geom_sf(data=contour_sab,linewidth=0.5,linetype=2,col="grey")+
              geom_sf(data=curdo_box,fill=NA)+
              geom_sf(data=scatarie_box,fill=NA)+
              geom_sf(data=sab_zones,fill=NA)+
              geom_sf(data=banks,fill=c("#F8766D","#00BFC4"))+
              geom_sf(data=scatarie_coords%>%filter(sample=="AR: Existing array"),aes(fill=sample),shape=23,size=1)+
              geom_sf(data=scatarie_coords%>%filter(sample=="AR: New deployments"),aes(fill=sample),shape=24,size=1)+
              geom_sf(data=scatarie_coords%>%filter(sample=="AR: Recovery"),aes(fill=sample),shape=25,size=1)+
              geom_sf(data=scatarie_coords%>%filter(sample=="AR: ROV recovery"),aes(fill=sample),shape=21,size=1.75)+
              geom_sf(data=sab_points%>%filter(sample=="eDNA"),aes(fill=sample),shape=21,size=1.5)+
              theme_bw()+
              coord_sf(xlim=c(-60.3,-58.3),ylim=c(45.75,46.5),expand=0)+
              theme(plot.margin = margin(0, 0, 0, 0),
                    legend.position="inside",
                    legend.position.inside = c(0.25,0.8),
                    legend.background = element_blank(),
                    legend.title = element_blank(),
                    legend.key = element_blank(),
                    legend.text = element_text(size = 6))+
                annotation_scale()+
                scale_fill_manual(values=c("AR: Recovery" = "aquamarine3",
                                           "AR: Existing array" = "grey80",
                                           "AR: New deployments"="blue2",
                                           "eDNA" = "white",
                                           "AR: ROV recovery" = "red"))+
                guides(fill = guide_legend(override.aes = list(size = 2.75)))
    
    mosaic_map <- p1 + (p1_scatarie/p1_curdo) + plot_layout(ncol=2,widths=c(4,1))
    
    ggsave("output/2025_mission/mission_map.png",mosaic_map,width=8,height=4,units="in",dpi=300)

    
### now compile the coordinates for the map
    
    form_coords <- rbind(
                  edna_coords%>%
                    filter(!grepl("Cam",name))%>%
                    dplyr::select(name,sample),
                  acoustic_coords%>%
                    filter(type!="Existing array")%>%
                    mutate(sample=paste0("AR: ",type))%>%
                    dplyr::select(name,sample)
                   )%>%
                   mutate(Longitude=st_coordinates(.)[,1],
                          Latitude=st_coordinates(.)[,2])%>%
                  rowwise()%>%
                  mutate(lon_ddm = convert_to_ddm(Longitude, is_longitude = TRUE),
                          lat_ddm = convert_to_ddm(Latitude, is_longitude = FALSE))%>%
                  ungroup()%>%
                  mutate(depth = round(terra::extract(sab_dem,.),1))%>%
                   st_transform(latlong)%>%
                   data.frame()%>%
                   dplyr::select(sample,name,depth,lon_ddm,lat_ddm)
    
    write.csv(form_coords,file="output/2025_mission/station_coords.csv",row.names=FALSE)
    