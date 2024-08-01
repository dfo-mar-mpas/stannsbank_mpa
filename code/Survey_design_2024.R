## map of recievers in ESI 2024

#load libraries ----
library(tidyr)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
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

#Functions ------

convert_to_decimal_degrees <- function(degree_minute) {
  parts <- strsplit(degree_minute, " ")[[1]]
  degrees <- as.numeric(parts[1])
  decimal_minutes <- as.numeric(parts[2])
  decimal_degrees <- degrees + (decimal_minutes / 60)
  return(decimal_degrees)
}

convert_to_ddm <- function(decimal_degrees) {
  degrees <- floor(decimal_degrees)  # Extract degrees
  minutes <- round((decimal_degrees - degrees) * 60,3)  # Calculate decimal minutes
  return(sprintf("%d° %f'", degrees, minutes))  # Format as "D° M.M'"
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

#load bathymetry ------
dem_sab <- raster("data/Bathymetry/sab_dem.tif") # SAB based on the 'predSpace' data
dem_esi <- raster("data/Bathymetry/esi_dem.tif") # ESI based on the Greenlaw 35 dem
load("data/Bathymetry/250m_stars_object.RData") # STARS object based on the 250m depth contour extracted (converted) from GEBCO

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
    grid_size <- 12.5 #12.5 grid which generates ~ 19 stations
    
    sab_grid <- sab%>%
                st_transform(utm)%>%
                st_make_grid(cellsize=grid_size)%>%
                st_intersection(sab%>%st_transform(utm))%>%
                st_intersection(perley_buffer_lb%>%st_transform(utm))
    
    sab_grid_sf <- grid%>%
                   st_as_sf()%>%
                   mutate(area=as.numeric(st_area(.)))%>%
                   filter(area>15)%>%
                   st_centroid()%>%
                   mutate(type="eDNA Sample")%>%
                   rename(geometry = 1)%>%
                   mutate(lon = st_coordinates(.)[, "X"]) %>%
                   arrange(lon) %>%
                   mutate(transect_group = determine_transect_group(lon, 2.5))%>%
                   group_by(transect_group) %>%
                   arrange(st_coordinates(.)[, "Y"], .by_group = TRUE) %>%
                   ungroup() %>%
                   mutate(name = paste0("SAB_", row_number()))%>%
                   dplyr::select(name,area,type,geometry)
    
    ggplot()+
      geom_sf(data=sab_grid,fill=NA)+
      geom_sf(data=sab_zones,fill=NA)+
      geom_sf(data=banks,aes(fill=name))+
      geom_sf(data=sab_grid_sf)+
      theme_bw()
    
    
#make a survey design for the banks ----
    curdo_eDNA <- banks%>%
                  filter(name=="Curdo Bank")%>%
                  st_transform(utm)%>%
                  st_centroid()
    
    scatarie_eDNA <- banks%>%
                      filter(name=="Scatarie Bank")%>%
                      st_transform(utm)%>%
                      st_centroid()
    #make a survey around the bank
    
    distances <- c(0.5,1,1.5)
    num_points <- 3
    
    curdo_points <- lapply(distances, function(radius) generate_points(curdo_eDNA, radius, num_points)) %>%
                    bind_rows()%>%
                    rename(geometry=geom)%>%
                    rbind(.,curdo_eDNA%>%
                            mutate(radius = 0)%>%
                            st_as_sf()%>%
                            dplyr::select(radius,geometry))%>%
                    mutate(name="Curdo Bank")
    
    curdo_circles <- NULL
   for(i in distances){
     
     curdo_circles <- rbind(curdo_circles,curdo_eDNA%>%
                                          st_buffer(i)%>%
                                          mutate(dist=i,name="Curdo Bank"))
     
   }
    
    #scatarie Bank
    
    distances_sb <- c(1,3,5)
    
    scatarie_points <- lapply(distances_sb, function(radius) generate_points(scatarie_eDNA, radius, 3)) %>%
                        bind_rows()%>%
                        rename(geometry=geom)%>%
                        rbind(.,scatarie_eDNA%>%
                                mutate(radius = 0)%>%
                                st_as_sf()%>%
                                dplyr::select(radius,geometry))%>%
                        mutate(name="Scatarie Bank")
    
    #scatarie circles
    scatarie_circles <- NULL
    for(i in distances_sb){
      
      scatarie_circles <- rbind(scatarie_circles,scatarie_eDNA%>%
                               st_buffer(i)%>%
                               mutate(dist=i,name="Scatarie Bank"))
      
    }
    
    
    sub_samp_pts <- rbind(curdo_points,scatarie_points)
    sub_samp_cir <- rbind(curdo_circles,scatarie_circles)
      
    ggplot()+
      geom_sf(data=banks,aes(fill=name))+
      geom_sf(data=sub_samp_cir,fill=NA)+
      geom_sf(data=sub_samp_pts)+
      theme_bw()
     
## load priority locations for camera drops
    cam_coords <- read_sf("data/2024 Mission/SAB_2024_cameraTargets.csv")%>%
                  st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                  mutate(type="Camera Survey")
    # %>%
    #               mutate(group = case_when(lat<45.88 ~ "Southwest_1",
    #                                        lat<45.9 & lon >(-59.45) ~ "Curdo_sw",
    #                                        )

#load coordinates sampled last year 
    
    stations_2023 <- read.csv("data/eDNA/PER-2023-767/SampleDataSheet.csv")%>%
                     filter(grepl("SAB",Station),
                            !is.na(Lat),
                            !is.na(Long),
                            Lat != "BLANK",
                            Long != "BLANK")%>%
                     mutate(Lat = trimws(Lat),
                            Long = trimws(Long))%>%
                     distinct(Station,.keep_all=TRUE)%>%
                     separate(Lat,c("dlat","ddlat"),sep=" ")%>%
                     separate(Long,c("dlong","ddlong"),sep=" ")%>%
                     mutate(dlat = as.numeric(dlat),
                            ddlat = as.numeric(ddlat),
                            dlong = as.numeric(dlong),
                            ddlong = as.numeric(ddlong),
                            lat=dlat + ddlat/60,
                            lon=(dlong + ddlong/60)*-1)%>%
                     st_as_sf(coords=c("lon","lat"),crs=latlong)
      

    
#make main map figure for formb
    
    curdo_plot_bound <- sub_samp_cir%>%
                        filter(name=="Curdo Bank",
                               dist==distances)%>%
                        st_buffer(0.25)%>%
                        st_transform(latlong)%>%
                        st_bbox()
    
    curdo_box <- curdo_plot_bound%>%st_as_sfc()
                      
    p1_curdo <- ggplot()+
                geom_sf(data=banks%>%filter(name=="Curdo Bank"),fill="#F8766D")+
                geom_sf(data=sub_samp_cir%>%filter(name=="Curdo Bank"),fill=NA)+
                geom_sf(data=sub_samp_pts%>%filter(name=="Curdo Bank"),fill="grey80",shape=21,size=2)+
                geom_sf(data=cam_coords,size=0.9,shape=21,fill="cornflowerblue")+
                theme_bw()+
                theme(axis.text=element_blank(),
                      plot.margin = margin(0, 0, 0, 0),
                      plot.title = element_text(size=10,vjust=-2))+
                coord_sf(xlim=curdo_plot_bound[c(1,3)],ylim=curdo_plot_bound[c(2,4)],expand=0)+
                labs(title="Curdo Bank")
  
    
    scatarie_plot_bound <- sub_samp_cir%>%
                            filter(name=="Scatarie Bank",
                                   dist==distances_sb)%>%
                            st_buffer(0.25)%>%
                            st_transform(latlong)%>%
                            st_bbox()
    
    scatarie_box <- scatarie_plot_bound%>%st_as_sfc()
    
    p1_scatarie <- ggplot()+
      geom_sf(data=banks%>%filter(name=="Scatarie Bank"),fill="#00BFC4")+
      geom_sf(data=sub_samp_cir%>%filter(name=="Scatarie Bank"),fill=NA)+
      geom_sf(data=sub_samp_pts%>%filter(name=="Scatarie Bank"),fill="grey80",shape=21,size=2)+
      geom_sf(data=cam_coords,size=0.9,shape=21,fill="cornflowerblue")+
      theme_bw()+
      theme(axis.text=element_blank(),
            plot.margin = margin(0, 0, 0, 0),
            plot.title = element_text(size=10,vjust=-2))+
      coord_sf(xlim=scatarie_plot_bound[c(1,3)],ylim=scatarie_plot_bound[c(2,4)],expand=0)+
      labs(title="Scatarie Bank")
    
    sample_points <- rbind(cam_coords%>%dplyr::select(type,geometry),
                           sab_grid_sf%>%st_transform(latlong)%>%dplyr::select(type,geometry))
    
    p1 <- ggplot()+
      geom_sf(data=basemap)+
      geom_sf(data=waypoints_lb%>%filter(name=="louisbourg"),size=3)+
      geom_sf(data=curdo_box,fill=NA)+
      geom_sf(data=scatarie_box,fill=NA)+
      geom_sf(data=sab_zones,fill=NA)+
      geom_sf(data=banks,fill=c("#F8766D","#00BFC4"))+
      geom_sf(data=sample_points,aes(fill=type),shape=21,size=2)+
      theme_bw()+
      coord_sf(xlim=c(-60.3,-58.3),ylim=c(45.75,46.5),expand=0)+
      theme(plot.margin = margin(0, 0, 0, 0),
            legend.position="inside",
            legend.position.inside = c(0.1,0.92),
            legend.background = element_blank(),
            legend.title = element_blank())+
      scale_fill_manual(values=c("cornflowerblue","aquamarine2"))
    
    mosaic_map <- p1 + (p1_scatarie/p1_curdo) + plot_layout(ncol=2,widths=c(4,1))
    
    ggsave("output/2024_mission/mission_map.png",mosaic_map,width=8,height=4,units="in",dpi=300)

    
### now compile the coordinates for the map
    
  form_b_coords <- cam_coords%>%
                   mutate(Activity="Camera Survey")%>%
                   rename(Attribute = rationale)%>%
                   dplyr::select(name,Activity,Attribute,geometry)%>%
                   rbind(.,
                         sab_grid_sf%>%
                         st_transform(latlong)%>%
                         mutate(Activity="eDNA Sample",
                                Attribute="Representative grid sample")%>%
                         dplyr::select(name,Activity,Attribute,geometry),
                         sub_samp_pts%>%
                         group_by(name)%>%
                         mutate(name=paste(gsub(" Bank","",name),row_number(),sep="_"),
                                Activity="eDNA sample",
                                Attribute="Species transition")%>%
                         ungroup()%>%
                         st_transform(latlong)%>%
                           dplyr::select(name,Activity,Attribute,geometry))%>%
                    mutate(Longitude=st_coordinates(.)[,1],
                           Latitude=st_coordinates(.)[,2],
                           lon_ddm=paste0("'",convert_to_ddm(Longitude)),
                           lat_ddm=convert_to_ddm(Latitude))%>%
                    st_transform(proj4string(dem_sab))%>%
                    mutate(Depth = round(raster::extract(dem_sab,as_Spatial(.)),1))%>%
                    st_transform(latlong)%>%
                    data.frame()%>%
                    rename(Station=name)%>%
                    dplyr::select(Station,Activity,Depth,lon_ddm,lat_ddm,
                                  Longitude,Latitude,Attribute)
            
  write.csv(form_b_coords,file = "output/2024_mission/station_coords.csv",row.names=FALSE)
    