## Survey Design for the 2026 Field Season in St. Anns Bank

#load libraries ----
library(tidyr)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(terra)
library(tidyterra)
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
#utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- 32620 

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
banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%
  st_transform(latlong)

#load the Offshore wind polygons -----
osw <- read_sf("C:/Users/stanleyr/Documents/GitHub/offfshore_wind/data/shapefiles/Designated_WEAs_25_07_29.shp")%>%
        st_transform(latlong)

#load coastline and make basemap ------
    coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")%>%
                st_transform(latlong)
    
    bounding_area <- sab%>%
      st_transform(utm)%>%
      st_buffer(dist = 100*1000)%>%
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
    

#Assemble coordinates ----
    
    #eDNA coords
    edna_coords <- read.csv("output/2024_mission/station_coords.csv")%>%
                   filter(grepl("eDNA",Activity))%>%
                   dplyr::select(Station,Depth,Longitude,Latitude)%>%
                   rename(station=Station,depth=Depth,lon=Longitude,lat=Latitude)%>%
                   st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                   mutate(type="eDNA: SAB")%>%
                   dplyr::select(type,station,depth,lon,lat)
    
    osw_edna_coords <- read.csv("data/2026 Mission/osw_stations.csv")%>%
                       st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                       mutate(type="eDNA: OSW")%>%
                       dplyr::select(type,station,depth,lon,lat)
    
    #acoustic coordiantes
    sab_stations <- read.csv("data/2026 Mission/otn-instrument-deployment-SABMPAnew_June2025 3.csv")%>%
                    rename(lon=DEPLOY_LONG_DD,lat=DEPLOY_LAT_DD,depth=BOTTOM_DEPTH)%>%
                    st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                    mutate(station=gsub("SB25","",STATION_NO),
                           station=gsub("SAB_CB25","CB_",station),
                           station=gsub("SAB_","SCAT_",station),
                           type="AR: rollover")%>%
                    dplyr::select(type,station,depth,lon,lat)
    
    kipca_stations <- read.csv("data/2026 Mission/otn-instrument-deployments-KIPCA_June2025 1.csv")%>%
                      rename(lon=DEPLOY_LONG,lat=DEPLOY_LAT,depth=BOTTOM_DEPTH,station=STATION_NO)%>%
                      st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                      mutate(type="AR: rollover")%>%
                      dplyr::select(type,station,depth,lon,lat)
    
    otn_sab_array <- read.csv("data/Acoustic/OTN_redesign_coords.csv")%>%
                     rename(lon=long,station=Id)%>%
                     mutate(paste0("SAB_",station))%>%
                     st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                     mutate(station=paste("SAB_",station),
                            type="AR: array")%>%
                     dplyr::select(type,station,depth,lon,lat)
    
    #ROV stations
    rov_stations <- read.csv("data/2026 Mission/rov_coords_decimal.csv")%>%
                    st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                    mutate(type="AR: recovery",
                           station=paste0("ROV_",1:n()))%>%
                    dplyr::select(type,station,depth,lon,lat)
    
    
  #combined all together
    all_coords <- rbind(edna_coords,osw_edna_coords,sab_stations,kipca_stations,rov_stations,otn_sab_array)
      
            
#evaluate operational buffer for 2025 season  -------------
    
    km_nm <- 1.852 #coefficent to go from nautical miles to km
    nm_km <- 1/km_nm #coefficient to go from km to nautical miles
    perley_speed <- 8.9 #knots - roughly the max speed with some buffer
    perley_transit_time <- 5 #time permitted for traveling to station. 
    
    
    waypoints_syd <- data.frame(name=c("North Sydney","WP1","WP2"),
                                lon=c(-60.241412,-60.212059,-60.179058),
                                lat=c(46.210767,46.223997,46.254106))%>%
      mutate(location="North Sydney")%>%
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
      st_buffer(dist=maxdist_syd*1000)%>%
      st_as_sf()%>%
      st_transform(latlong)%>%
      mutate(name="Sydney")
    
    #Louisbourg operational window
    perley_buffer_lb <- waypoints_lb[nrow(waypoints_lb),]%>%
      st_transform(utm)%>% #km based planar projection
      st_buffer(dist=maxdist_lb*1000)%>%
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
    

#make main map figure for formb ------------
    
    curdo_plot_bound <- edna_coords%>%
                        filter(grepl("Curdo",station))%>%
                        st_transform(utm)%>%
                        st_buffer(0.5*1000)%>%
                        st_transform(latlong)%>%
                        st_bbox()
    
    curdo_box <- curdo_plot_bound%>%st_as_sfc()
    
    point_scaler <- 0.5
                      
    p1_curdo <- ggplot()+
                geom_sf(data=banks%>%filter(name=="Curdo Bank"),fill="#F8766D")+
                geom_sf(data=all_coords%>%filter(grepl("Curdo",station),type!="eDNA: OSW"),aes(fill=type),shape=21,size=2.5*point_scaler)+
                geom_sf(data=all_coords%>%filter(!grepl("eDNA",type),type!="eDNA: OSW"),aes(fill=type),shape=24,size=2.5*point_scaler)+
                theme_bw()+
                theme(
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  #plot.margin = margin(0, 0, 0, 0),  # top, right, bottom, left
                  plot.margin = margin(0,-10,0,-10, unit="pt"),
                  plot.title = element_text(size = 10, vjust = -2),
                  legend.position = "none",
                  panel.spacing = margin(0, 0, 0, 0)
                )+
                coord_sf(xlim=curdo_plot_bound[c(1,3)],ylim=curdo_plot_bound[c(2,4)],expand=FALSE)+
                labs(title="Curdo Bank",fill="")+
                annotation_scale( location = "br",
                                  pad_x = unit(0.02, "in"),
                                  pad_y = unit(0.02, "in"))+
                scale_fill_manual(values=c("AR: array" = "aquamarine3",
                                           "AR: recovery" = "red",
                                           "AR: rollover"="blue2",
                                           "eDNA: SAB" = "white"))
    
    
    #Sactarie Bank
    scatarie_plot_bound <- banks%>%
                          filter(name=="Scatarie Bank")%>%
                           st_transform(utm)%>%
                           st_buffer(3*1000)%>%
                           st_transform(latlong)%>%
                           st_bbox()
    
    scatarie_box <- scatarie_plot_bound%>%st_as_sfc()
    
    p1_scatarie <- ggplot()+
                    geom_sf(data=banks%>%filter(name=="Scatarie Bank"),fill="cornflowerblue")+
                    geom_sf(data=all_coords%>%filter(type=="AR: array"),aes(fill=type),shape=23,size=2.5*point_scaler)+
                    geom_sf(data=all_coords%>%filter(type=="AR: rollover"),aes(fill=type),shape=24,size=2.5*point_scaler)+
                    geom_sf(data=all_coords%>%filter(type=="eDNA: SAB"),aes(fill=type),shape=21,size=2.5*point_scaler)+
                    geom_sf(data=all_coords%>%filter(type=="AR: recovery"),aes(fill=type),shape=21,size=3.25*point_scaler)+              
                    theme_bw()+
                    theme(
                      axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      axis.title = element_blank(),
                      #plot.margin = margin(0, 0, 0, 0),  # top, right, bottom, left
                      plot.margin = margin(0,-10,0,0, unit="pt"),
                      plot.title = element_text(size = 10, vjust = -2),
                      legend.position = "none",
                      panel.spacing = margin(0, 0, 0, 0)
                    )+
                    coord_sf(xlim=scatarie_plot_bound [c(1,3)],ylim=scatarie_plot_bound [c(2,4)],expand=FALSE)+
                    labs(title="Scatarie Bank",fill="")+
                    annotation_scale(location="br",
                                     pad_x = unit(0.02, "in"),
                                     pad_y = unit(0.02, "in"))+
                    scale_fill_manual(values=c("AR: array" = "aquamarine3",
                                               "AR: recovery" = "red",
                                               "AR: rollover"="blue2",
                                               "eDNA: SAB" = "white"))
    
    #offshore wind plot
    
    osw_lims <- osw%>%
                filter(WEA == "Sydney Bight")%>%
                st_transform(utm)%>%
                st_buffer(10*1000)%>%
                st_transform(latlong)%>%
                st_bbox()
    
    kipca_lims <- kipca_stations%>%
                  st_transform(utm)%>%
                  st_buffer(10*1000)%>%
                  st_transform(latlong)%>%
                  st_bbox()
    
    outside_lims <- osw_lims%>%
                    st_as_sfc()%>%
                    st_union(kipca_lims%>%st_as_sfc())%>%
                    st_bbox()
      
    
    p1_osw <- ggplot()+
              geom_sf(data=coast_hr)+
              geom_sf(data=osw%>%
                        filter(WEA == "Sydney Bight"),fill="green3")+
              geom_sf(data=all_coords%>%filter(type=="eDNA: OSW"),aes(fill=type),shape=21,size=2.5*point_scaler)+
              geom_sf(data=all_coords%>%filter(type=="AR: rollover"),aes(fill=type),shape=24,size=2.5*point_scaler)+             
              theme_bw()+
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                #plot.margin = margin(0, 0, 0, 0),  # top, right, bottom, left
                plot.margin = margin(0,0,0,-10, unit="pt"),
                plot.title = element_text(size = 10, vjust = -2),
                legend.position = "none",
                panel.spacing = margin(0, 0, 0, 0)
              )+
              coord_sf(xlim=outside_lims[c(1,3)],ylim=outside_lims[c(2,4)],expand=FALSE)+
              labs(title="Sydney Bite",fill="")+
              annotation_scale( location = "br",
                                pad_x = unit(0.02, "in"),
                                pad_y = unit(0.02, "in"))+
              scale_fill_manual(values=c("AR: array" = "aquamarine3",
                                         "AR: recovery" = "red",
                                         "AR: rollover"="blue2",
                                         "eDNA: OSW" = "white"))
            
    
    #zoomed out plot
    sab_points <- all_coords%>%
                  mutate(type=ifelse(grepl("eDNA",type),"eDNA",type))
    
    #set regional plotting limits
    total_plot <- sab%>%
                  st_transform(utm)%>%
                  st_buffer(10*1000)%>%
                  st_transform(latlong)%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_union(outside_lims%>%st_as_sfc())%>%
                  st_bbox()
    
    #download bathymetry 
    bathy <- getNOAA.bathy(
      lon1 = total_plot["xmin"],
      lon2 = total_plot["xmax"],
      lat1 = total_plot["ymin"],
      lat2 = total_plot["ymax"],
      resolution = 1
    )
    
    bathy_r <- marmap::as.raster(bathy)
    
    bathy_terra <- terra::rast(bathy_r)
    
    contours <- terra::as.contour(
      bathy_terra,
      levels = -200   # marmap uses negative depths
    )
  
    p1 <- ggplot() +
      geom_sf(data = basemap) +
      geom_sf(data = sab_zones, fill = NA) +
      geom_spatvector(data = contours, colour = "grey70", linewidth = 0.4,linetype=2) +
      geom_sf(data = osw %>% filter(WEA == "Sydney Bight"), fill = "green3") +
      geom_sf(data = banks, fill = c("#F8766D", "cornflowerblue")) +
      geom_sf(data = sab_points %>% filter(type == "AR: array"),
              aes(fill = type), shape = 23, size = 1) +
      geom_sf(data = sab_points %>% filter(type == "AR: rollover"),
              aes(fill = type), shape = 25, size = 1) +
      geom_sf(data = sab_points %>% filter(type == "AR: recovery"),
              aes(fill = type), shape = 21, size = 1.75) +
      geom_sf(data = sab_points %>% filter(grepl("eDNA", type)),
              aes(fill = type), shape = 21, size = 1.5) +
      coord_sf(xlim = total_plot[c(1,3)],
               ylim = total_plot[c(2,4)],
               expand = FALSE) +
      theme_bw() +
      theme(
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 6)
      ) +
      annotation_scale() +
      scale_fill_manual(values = c(
        "AR: array" = "aquamarine3",
        "AR: recovery" = "red",
        "AR: rollover" = "blue2",
        "eDNA" = "white"
      )) +
      guides(fill = guide_legend(override.aes = list(size = 2.75)))
    
    
    #mosaic patchwork plot
    
    bottom_row <- p1_scatarie + p1_curdo + p1_osw +
      plot_layout(
        ncol = 3,
        widths = c(1, 1, 1.25)
      ) &
      theme(
        plot.margin = margin(0, 0, 0, 0)
      )
    
    mosaic_map <- p1 / bottom_row +
      plot_layout(
        heights = c(3.8, 1)
      )
    
    ggsave(
      "output/2026_mission/mission_map.png",
      mosaic_map,
      width = 8,
      height = 6.2,
      units = "in",
      dpi = 300,
      bg = "white"
    )

    
### now compile the coordinates for the map
    
    form_coords <- all_coords%>%
                   filter(type!="AR: array")%>%
                   st_drop_geometry()%>%
                   mutate(POC = case_when(grepl("KIPCA",station) ~ "North Sydney",
                                          grepl("OSW",type)~ "North Sydney",
                                          TRUE ~ "Louisbourg"),
                          POC = factor(POC,levels=c("North Sydney","Louisbourg")),
                          type=factor(type,levels=c("eDNA: SAB","eDNA: OSW","AR: rollover","AR: recovery")))%>%
                  arrange(POC,type)%>%
                   rename(Longitude=lon,
                          Latitude=lat)%>%
                  rowwise()%>%
                  mutate(lon_ddm = convert_to_ddm(Longitude, is_longitude = TRUE),
                          lat_ddm = convert_to_ddm(Latitude, is_longitude = FALSE))%>%
                  ungroup()%>%
                  data.frame()%>%
                  dplyr::select(type,station,depth,lon_ddm,lat_ddm)
    
    write.csv(form_coords,file="output/2026_mission/station_coords_2026.csv",row.names=FALSE)
    