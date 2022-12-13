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
library(viridis)

s2_as_sf = FALSE

source("r:/Science/CESD/HES_MPAGroup/R/Functions/box_grid.R")

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load bathymetry ------
bathy <- raster("r:/Science/CESD/HES_MPAGroup/Data/Bathymetry/dem35c_5c.tif")
bathy_zoom <- raster("r:/Science/CESD/HES_MPAGroup/Data/Bathymetry/NONNA 10/Moser Shallow geotiff/NONNA100_4400N06300W.tiff") #nonna 100 ** the nonna 10 didn't have data

# load("r:/Science/CESD/HES_MPAGroup/Data/Bathymetry/PredSpaceData.RData")
# bathydat <- predSpace%>%dplyr::select("plon","plat","z")
# colnames(bathydat) <- c("long","lat","depth")
# coordinates(bathydat) <- ~long+lat
# proj4string(bathydat) <- utm_mar
# gridded(bathydat) <- TRUE
# dem <- crop(raster(bathydat),bounding_area%>%st_transform(utm_mar)%>%as_Spatial())
# writeRaster(dem,"r:/Science/CESD/HES_MPAGroup/Projects/eDNA Genome Quebec/data/sab_dem.tif")

dem <- raster("r:/Science/CESD/HES_MPAGroup/Projects/eDNA Genome Quebec/data/sab_dem.tif")


#load stations ---------

#load the 2022 stations -- note that these coordinates are derived from below. 
esi_stations <- read.csv("r:/Science/CESD/HES_MPAGroup/Projects/Acoustic/ESI/2022_stations_depth.csv")%>%
                rename(lon=lon_dec_deg,lat=lat_dec_deg)%>%
                st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                mutate(line=ifelse(grepl("SR",station),"Shag Roost","Moser River"),
                       retrieve = ifelse(station %in% paste0("SR_",c("01","02", "04", "05", "08", "09", 10, 11, 14,15)),FALSE,TRUE),#these were retrieved prior to the 2022 mission on the PackCat
                       retrieve = ifelse(station %in% paste0("MR_",c("01","02","03")),FALSE,retrieve)) #too shallow for the Perley will be retrieved in fall 2022  
              
edna_stations <- read.csv("r:/Science/CESD/HES_MPAGroup/Projects/Eastern Shore Islands/Water Chemistry/SamplingDesign_65by20km.csv")%>%
                 st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                 mutate(id=c(16,12,8,4,15,11,7,3,14,10,6,2,13,9,5,1))%>%#reorder because they were done along shore not perpendicular
                 arrange(id)%>%
                 mutate(station=paste0("ESI_",1:16),
                        transect = rep(paste0("Transect_",1:4),each=4))


# esi_stations <- read.csv("r:/Science/CESD/HES_MPAGroup/Projects/Acoustic/ESI/AcousticStations.csv")%>%
#                 filter(grepl("SR",Name)|grepl("MR",Name))%>%
#                 st_as_sf(coords=c("Lon","Lat"),crs=latlong,remove=FALSE)%>%
#                 mutate(line = ifelse(grepl("MR",Name),"Moser River","Shag Roost"),
#                        status = "Core stations")%>%
#                 rename(lon=Lon,lat=Lat,name=Name)%>%
#                 dplyr::select(line,name,status,lon,lat,geometry)

#way points to the mouth of sheet harbour

waypoints <- data.frame(name=c("WP1","WP2","WP3"),
                        lon=c(-62.483136,-62.497729,-62.475491),
                        lat=c(44.891778,44.835977,44.803352))%>%
             st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)


sheet_harbour_port <- data.frame(lon=-62.504091,lat=44.903582,name="Port of Sheet Harbour")%>%
                      st_as_sf(coords=c("lon","lat"),crs=latlong)

louisbourg <- data.frame(lon=-59.970314,lat=45.917173,name="Louisbourg")%>%
        st_as_sf(coords=c("lon","lat"),crs=latlong)

waypoints_cb <- data.frame(name=c("WP1","WP2","WP3"),
                           lon=c(-59.974171,-59.959054,-59.941645),
                           lat=c(45.905090,45.902111,45.900265))%>%
        st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)

#Shape files --------

#ESI
esi_poly <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/EasternShoreIslands_networksite.shp")%>%
            st_transform(latlong)

##SAB
sab_zones <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
        st_transform(latlong)

sab  <- sab_zones%>%
        st_transform(utm)%>%
        st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
        st_union()%>% #gets rid of the zones
        st_transform(latlong)%>%
        st_as_sf()

sab_bounds <- st_bbox(sab)

#SAB Benthoscape
sab_benthoscape <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/StAnns_Benthoscape/benthoscape/benthoscape.shp")%>%
        st_transform(latlong)%>%
        rename(class=Assigned_c)%>%
        dplyr::select(class,geometry)%>%
        st_make_valid()%>%
        group_by(class)%>%
        summarise()
                             

#load coastline ------
coast_hr <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/NS_coastline_project_Erase1.shp")

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

# basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
#                    dplyr::select(name_en,geometry)%>%
#                    st_as_sf()%>%
#                    st_union()%>%
#                    st_transform(latlong)%>%
#                    st_as_sf()%>%
#                    mutate(country="Canada")%>%
#                    st_intersection(bounding_area)%>%
#                   suppressWarnings()%>%
#                   suppressMessages()


basemap <- coast_hr%>%
        st_intersection(bounding_area%>%st_transform(st_crs(coast_hr)))%>%
        st_transform(latlong)%>%
        suppressWarnings()

#load in 


#find the distance between sequential stations ---------

# mr_stations <- esi_stations%>%filter(line=="Moser River")%>%pull(name)
# 
# dist_mr <- NULL
# for(i in 2:length(mr_stations)){
#   
#   dist_mr <- c(dist_mr,esi_stations%>%
#                        filter(name %in% mr_stations[c(i-1,i)])%>%
#                        st_distance()%>%
#                        max()%>%
#                        as.numeric())
# }
# 
# mean_spacing <- mean(dist_mr[dist_mr<1000]) #spacing (m) between mr_03 and mr_04 is higher because it snakes around an island. 
# 
# ##calculate the bearing of the coordinates ------
# 
# x <- esi_stations%>%
#      filter(name %in% c("MR_04","MR_21"))%>%
#      as_Spatial()
# 
# bear <- geosphere::bearing(x[1,],x[2,])
# 
# #estimate locations of next 5 stations ----------
# 
# newstations <- geosphere::destPoint(x[2,],bear,mean_spacing*1:5)%>%
#                data.frame()%>%
#                st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
#                mutate(name = paste0("MR_",22:26),
#                       line = "Moser River",
#                       status = "2022 station additions")%>%
#                dplyr::select(names(esi_stations))
#             
# 
# #bind it back together ----------
# esi_stations_revised <- rbind(esi_stations,newstations)
# 
# #check the distance again
# mr_stations_new <- esi_stations%>%filter(line=="Moser River")%>%pull(name)
# 
# dist_mr <- NULL
# for(i in 2:length(mr_stations_new)){
#   
#   dist_mr <- c(dist_mr,esi_stations%>%
#                  filter(name %in% mr_stations[c(i-1,i)])%>%
#                  st_distance()%>%
#                  max()%>%
#                  as.numeric())
# }
# 
# #check the distances are about the same (with some rounding errors) as the mean_spacing
# round(mean(dist_mr[dist_mr<1000]),1) == round(mean_spacing,1)      

# #now we extract the depths for the stations -----------
# esi_extent <- esi_stations_revised%>%
#               st_bbox()%>%
#               st_as_sfc()%>%
#               st_as_sf()%>%
#               st_transform(utm)%>%
#               st_buffer(2)%>% #2 km buffer
#               st_transform(proj4string(bathy))%>%
#               extent()
# 
# 
# bathy_crop <- crop(bathy,esi_extent)
# plot(bathy_crop)
# plot(esi_stations_revised%>%st_transform(proj4string(bathy))%>%as_Spatial(),add=T)
# 
# esi_stations_revised <- esi_stations_revised%>%
#                         st_transform(proj4string(bathy_crop))%>%
#                         mutate(depth=extract(bathy_crop,as_Spatial(.)))%>%
#                         st_transform(latlong)%>%
#                         dplyr::select(line,name,depth,status,lon,lat,geometry)
# 
# ## get the depth extracts for the shallowest MR line stations
# moser_shallow <- c("MR_01","MR_02","MR_03")
# 
# moser_shallow_extract <- esi_stations%>%
#                           filter(name %in% moser_shallow)%>%
#                           st_transform(proj4string(bathy_zoom))%>%
#                           mutate(depth=extract(bathy_zoom,as_Spatial(.)))%>%
#                           st_transform(latlong)%>%
#                           dplyr::select(line,name,depth,status,lon,lat,geometry)
# 
# comparison_df <- cbind(moster_shallow_extract%>%
#                          rename(depth_nonna=depth)%>%
#                          data.frame()%>%
#                          dplyr::select(name,depth_nonna),
#                        esi_stations_revised%>%
#                          rename(depth_greenlaw = depth)%>%
#                          data.frame()%>%
#                          mutate(depth_greenlaw = depth_greenlaw*-1)%>%
#                          filter(name %in% moser_shallow)%>%
#                          dplyr::select(depth_greenlaw))


# #save the revised stations
# write.csv(esi_stations_revised%>%
#             data.frame()%>%
#             dplyr::select(line,name,depth,lon,lat),file="ESI/2022_stations_depth.csv",row.names=F)


#Make the plot ---------

esi_stations_revised <- esi_stations #to keep code consistent

basemap_extent <- esi_stations_revised%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_as_sf()%>%
                  st_transform(utm)%>%
                  st_buffer(5)%>%
                  st_transform(st_crs(coast_hr))%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_as_sf()

aoi_extent <- esi_poly%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()%>%
                st_transform(utm)%>%
                st_buffer(5)%>%
                st_transform(st_crs(coast_hr))%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()

basemap_lims <- basemap_extent%>%st_transform(latlong)%>%st_bbox()
aoi_lims <- aoi_extent%>%st_transform(latlong)%>%st_bbox()

coast_hr_crop <- coast_hr%>%
                 st_intersection(basemap_extent)%>%
                 st_transform(latlong)

coast_hr_crop_esi <- coast_hr%>%
                     st_intersection(aoi_extent)%>%
                     st_transform(latlong)

p1 <- ggplot()+
        geom_sf(data=esi_poly,fill=NA)+
        geom_sf(data=coast_hr_crop,fill="darkolivegreen3")+
        geom_sf(data=sheet_harbour_port,size=3)+
        geom_sf_label(data=sheet_harbour_port,aes(label=name),nudge_x = 0.02,nudge_y = 0.02)+
        geom_sf(data=esi_stations_revised%>%mutate(Retrieval = ifelse(retrieve,'Yes','No')),aes(fill=Retrieval),shape=21,size=2)+
        scale_fill_manual(values=c("white","red"))+
        theme_bw()+
        coord_sf(expand=0,xlim = basemap_lims[c(1,3)],ylim=basemap_lims[c(2,4)])+
        theme(legend.position = "none")+
        labs(x="",y="",col="")+
        annotation_scale();p1

ggsave("2022 Water Sampling/2022_station_design_formb.png",p1,width=7,height=5,units="in",dpi=300)

#Calculate distances for the travel

#Leg 1 = Sheet Harbour - SR_03 - SR_13 - MR_04-MR_26
leg1 <- rbind(sheet_harbour_port,
              waypoints%>%
                 dplyr::select(name)%>%
                 filter(name !="WP3"), # this is more direct to not go to mouth and go west from WP2
              esi_stations%>%
                 filter(line=="Shag Roost",retrieve)%>%
                 dplyr::select(station)%>%
                 rename(name=station),
              esi_stations%>%
                 filter(line=="Moser River",retrieve)%>%
                 dplyr::select(station)%>%
                 rename(name=station),
              waypoints%>%
                 mutate(order=1:n())%>%
                 arrange(-order)%>%
                 dplyr::select(name),
              sheet_harbour_port)%>%
        mutate(leg="Leg 1",
               group="Acoustic")

leg1_line <- leg1%>%
        summarise(do_union = FALSE)%>%
        st_cast("LINESTRING")%>%
        mutate(leg="Leg 1",
               group="Acoustic")

leg1_dist <- leg1_line%>%
             st_length()%>%
             as.numeric()/1000/1.852 #convert to nautical miles from 'km' -- base calculation is in 'm'

leg1_dist/9.5 

leg1_map <- ggplot()+
        geom_sf(data=esi_poly,fill=NA)+
        geom_sf(data=coast_hr_crop_esi)+
        geom_sf(data=leg1)+
        geom_sf(data=leg1_line)+
        geom_sf(data=sheet_harbour_port,size=3,pch=21)+
        coord_sf(expand=0,xlim = aoi_lims[c(1,3)],ylim=aoi_lims[c(2,4)])+
        theme_bw()+
        annotation_scale();leg1_map



#Leg 2 Sheet Harbour ESI_01 - 04
        
leg2 <- rbind(sheet_harbour_port,
              waypoints%>%
                      dplyr::select(name), # this is more direct to not go to mouth and go west from WP2
              edna_stations%>%
                      filter(station %in% paste0("ESI_",1:4))%>%
                      rename(name=station)%>%
                      dplyr::select(name),
              waypoints%>%
                      mutate(order=1:n())%>%
                      arrange(-order)%>%
                      dplyr::select(name),
              sheet_harbour_port)%>%
        mutate(leg="Leg 2",
               group="eDNA Day 1")

leg2_line <- leg2%>%
        summarise(do_union = FALSE)%>%
        st_cast("LINESTRING")%>%
        mutate(leg="Leg 2",
               group="eDNA Day 1")

leg2_dist <- leg2_line%>%
        st_length()%>%
        as.numeric()/1000/1.852 #convert to nautical miles from 'km' -- base calculation is in 'm'

outer_dist_leg2 <- rbind(sheet_harbour_port,
                         waypoints%>%
                                 dplyr::select(name)%>%
                                 filter(name !="WP3"), # this is more direct to not go to mouth and go west from WP2
                         edna_stations%>%
                                 filter(station == "ESI_4")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name))%>%
        summarise(do_union=FALSE)%>%
        st_cast("LINESTRING")%>%
        st_length()%>%
        as.numeric()/1000/1.852

inner_dist_leg2 <- rbind(edna_stations%>%
                                 filter(station == "ESI_1")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name),
                         waypoints%>%
                                 mutate(order=1:n())%>%
                                 arrange(-order)%>%
                                 dplyr::select(name),
                         sheet_harbour_port)%>%
        summarise(do_union=FALSE)%>%
        st_cast("LINESTRING")%>%
        st_length()%>%
        as.numeric()/1000/1.852

message(paste0("Distance to inner station ESI 05 is ",round(inner_dist_leg3,2)," nautical miles and distance to outer station ESI 08 is ",round(outer_dist_leg3,2),". Total distance is ",round(leg3_dist,2)," nautical miles."))


leg2_map <- ggplot()+
        geom_sf(data=esi_poly,fill=NA)+
        geom_sf(data=coast_hr_crop_esi)+
        geom_sf(data=leg2)+
        geom_sf(data=leg2_line)+
        geom_sf(data=sheet_harbour_port,size=3,pch=21)+
        coord_sf(expand=0,xlim = aoi_lims[c(1,3)],ylim=aoi_lims[c(2,4)])+
        theme_bw()+
        annotation_scale();leg2_map

#Leg 3 Sheet Harbour ESI_05 - 08

leg3 <- rbind(sheet_harbour_port,
              waypoints%>%
                      dplyr::select(name)%>%
                      filter(name !="WP3"), # this is more direct to not go to mouth and go west from WP2
              edna_stations%>%
                      filter(station %in% paste0("ESI_",5:8))%>%
                      rename(name=station)%>%
                      dplyr::select(name)%>%
                      mutate(id=1:n())%>%
                      arrange(-id)%>%
                      dplyr::select(-id),
              edna_stations%>%
                      filter(station=="ESI_5")%>%
                      rename(name=station)%>%
                      dplyr::select(name),
              waypoints%>%
                      filter(name !="WP3")%>%
                      mutate(order=1:n())%>%
                      arrange(-order)%>%
                      dplyr::select(name),
              sheet_harbour_port)%>%
        mutate(leg="Leg 3",
               group="eDNA Day 2")

leg3_line <- leg3%>%
        summarise(do_union = FALSE)%>%
        st_cast("LINESTRING")%>%
        mutate(leg="Leg 3",
               group="eDNA Day 2")

leg3_dist <- leg3_line%>%
        st_length()%>%
        as.numeric()/1000/1.852 #convert to nautical miles from 'km' -- base calculation is in 'm'

outer_dist_leg3 <- rbind(sheet_harbour_port,
                         waypoints%>%
                                 dplyr::select(name)%>%
                                 filter(name !="WP3"), # this is more direct to not go to mouth and go west from WP2
                         edna_stations%>%
                                 filter(station == "ESI_8")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name))%>%
        summarise(do_union=FALSE)%>%
        st_cast("LINESTRING")%>%
        st_length()%>%
        as.numeric()/1000/1.852

inner_dist_leg3 <- rbind(edna_stations%>%
                                 filter(station == "ESI_5")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name),
                         waypoints%>%
                                 mutate(order=1:n())%>%
                                 arrange(-order)%>%
                                 dplyr::select(name),
                         sheet_harbour_port)%>%
        summarise(do_union=FALSE)%>%
        st_cast("LINESTRING")%>%
        st_length()%>%
        as.numeric()/1000/1.852

message(paste0("Distance to inner station ESI 05 is ",round(inner_dist_leg3,2)," nautical miles and distance to outer station ESI 08 is ",round(outer_dist_leg3,2),". Total distance is ",round(leg3_dist,2)," nautical miles."))


leg3_map <- ggplot()+
        geom_sf(data=esi_poly,fill=NA)+
        geom_sf(data=coast_hr_crop_esi)+
        geom_sf(data=leg3)+
        geom_sf(data=leg3_line)+
        geom_sf(data=sheet_harbour_port,size=3,pch=21)+
        coord_sf(expand=0,xlim = aoi_lims[c(1,3)],ylim=aoi_lims[c(2,4)])+
        theme_bw()+
        annotation_scale();leg3_map


#Leg 4 Sheet Harbour ESI_13 - 16 (East)

leg4 <- rbind(sheet_harbour_port,
              waypoints%>%
                      dplyr::select(name), # this is more direct to not go to mouth and go west from WP2
              edna_stations%>%
                      filter(station %in% paste0("ESI_",13:16))%>%
                      rename(name=station)%>%
                      dplyr::select(name),
              waypoints%>%
                      mutate(order=1:n())%>%
                      arrange(-order)%>%
                      dplyr::select(name),
              sheet_harbour_port)%>%
        mutate(leg="Leg 4",
               group="eDNA Day 3")

leg4_line <- leg4%>%
        summarise(do_union = FALSE)%>%
        st_cast("LINESTRING")%>%
        mutate(leg="Leg 4",
               group="eDNA Day 3")

leg4_dist <- leg4_line%>%
        st_length()%>%
        as.numeric()/1000/1.852 #convert to nautical miles from 'km' -- base calculation is in 'm'

outer_dist_leg4 <- rbind(sheet_harbour_port,
                         waypoints%>%
                                 dplyr::select(name), # this is more direct to not go to mouth and go west from WP2
                         edna_stations%>%
                                 filter(station == "ESI_16")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name))%>%
                         summarise(do_union=FALSE)%>%
                         st_cast("LINESTRING")%>%
                         st_length()%>%
                         as.numeric()/1000/1.852

inner_dist_leg4 <- rbind(edna_stations%>%
                                 filter(station == "ESI_13")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name),
                         waypoints%>%
                                 mutate(order=1:n())%>%
                                 arrange(-order)%>%
                                 dplyr::select(name),
                         sheet_harbour_port)%>%
                        summarise(do_union=FALSE)%>%
                        st_cast("LINESTRING")%>%
                        st_length()%>%
                        as.numeric()/1000/1.852

message(paste0("Distance to inner station ESI 13 is ",round(inner_dist_leg4,2)," nautical miles and distance to outer station ESI 16 is ",round(outer_dist_leg4,2),". Total distance is ",round(leg4_dist,2)," nautical miles."))

leg4_map <- ggplot()+
        geom_sf(data=esi_poly,fill=NA)+
        geom_sf(data=coast_hr_crop_esi)+
        geom_sf(data=leg4)+
        geom_sf(data=leg4_line)+
        geom_sf(data=sheet_harbour_port,size=3,pch=21)+
        coord_sf(expand=0,xlim = aoi_lims[c(1,3)],ylim=aoi_lims[c(2,4)])+
        theme_bw()+
        annotation_scale();leg4_map

#Leg 5 Sheet Harbour ESI_09 - 12

leg5 <- rbind(sheet_harbour_port,
              waypoints%>%
                      dplyr::select(name), # this is more direct to not go to mouth and go west from WP2
              edna_stations%>%
                      filter(station %in% paste0("ESI_",9:12))%>%
                      rename(name=station)%>%
                      dplyr::select(name),
              waypoints%>%
                      mutate(order=1:n())%>%
                      arrange(-order)%>%
                      dplyr::select(name),
              sheet_harbour_port)%>%
        mutate(leg="Leg 5",
               group="eDNA Day 4")

leg5_line <- leg5%>%
        summarise(do_union = FALSE)%>%
        st_cast("LINESTRING")%>%
        mutate(leg="Leg 5",
               group="eDNA Day 4")

leg5_dist <- leg5_line%>%
        st_length()%>%
        as.numeric()/1000/1.852 #convert to nautical miles from 'km' -- base calculation is in 'm'

outer_dist_leg5 <- rbind(sheet_harbour_port,
                         waypoints%>%
                                 dplyr::select(name), # this is more direct to not go to mouth and go west from WP2
                         edna_stations%>%
                                 filter(station == "ESI_12")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name))%>%
        summarise(do_union=FALSE)%>%
        st_cast("LINESTRING")%>%
        st_length()%>%
        as.numeric()/1000/1.852

inner_dist_leg5 <- rbind(edna_stations%>%
                                 filter(station == "ESI_9")%>%
                                 rename(name=station)%>%
                                 dplyr::select(name),
                         waypoints%>%
                                 mutate(order=1:n())%>%
                                 arrange(-order)%>%
                                 dplyr::select(name),
                         sheet_harbour_port)%>%
        summarise(do_union=FALSE)%>%
        st_cast("LINESTRING")%>%
        st_length()%>%
        as.numeric()/1000/1.852

message(paste0("Distance to inner station ESI 09 is ",round(inner_dist_leg5,2)," nautical miles and distance to outer station ESI 12 is ",round(outer_dist_leg5,2),". Total distance is ",round(leg5_dist,2)," nautical miles."))


leg5_map <- ggplot()+
        geom_sf(data=esi_poly,fill=NA)+
        geom_sf(data=coast_hr_crop_esi)+
        geom_sf(data=leg5)+
        geom_sf(data=leg5_line)+
        geom_sf(data=sheet_harbour_port,size=3,pch=21)+
        coord_sf(expand=0,xlim = aoi_lims[c(1,3)],ylim=aoi_lims[c(2,4)])+
        theme_bw()+
        annotation_scale();leg5_map




#Map all of the legs ---------


#assemble plotting data

                #labels
                Facet_labs <- data.frame(leg=paste0("Leg ",1:5),
                                         group=c("Acoustic",paste0("eDNA Day ",1:4)))%>%
                                  mutate(lab=paste0(group," (~",c(round(leg1_dist,2),
                                                         round(leg2_dist,2),
                                                         round(leg3_dist,2),
                                                         round(leg4_dist,2),
                                                         round(leg5_dist,2)),
                                                    " nm)"))%>%
                              rbind(.,data.frame(leg = "Eastern Shore Islands 2022",group="Eastern Shore Islands 2022", lab="Eastern Shore Islands 2022"))
                
                #lines
                esi_leg_lines <- rbind(leg1_line,leg2_line,leg3_line,leg4_line,leg5_line)%>%
                                 left_join(.,Facet_labs)
                
                # esi_leg_lines <- rbind(esi_leg_lines,
                #                        esi_leg_lines%>%mutate(lab="Eastern Shore Islands 2022"))
                 
                 esi_leg_lines$lab <- factor(esi_leg_lines$lab,levels=Facet_labs$lab)
                #       
                #points  
                esi_leg_pts <- rbind(leg1,leg2,leg3,leg4,leg5)%>%
                                left_join(.,Facet_labs)
                
                esi_leg_pts <- rbind(esi_leg_pts,
                                       esi_leg_pts%>%mutate(group="Eastern Shore Islands 2022",lab="Eastern Shore Islands 2022"))
                
                esi_leg_pts$lab <- factor(esi_leg_pts$lab,levels=Facet_labs$lab)
   #plot it
        ESI_plan <- ggplot()+
         geom_sf(data=esi_poly,fill=NA,col="black")+
         geom_sf(data=coast_hr_crop_esi,fill="darkolivegreen3",col="black")+
         geom_sf(data=esi_leg_lines)+
         geom_sf(data=esi_leg_pts,aes(col=leg))+
         facet_wrap(~lab,nrow=3)+
         theme_bw()+
         coord_sf(expand=0,xlim = aoi_lims[c(1,3)],ylim=aoi_lims[c(2,4)])+
         theme(legend.position="none",
               strip.background = element_rect(fill="white"))
        
        ggsave("2022 Water Sampling//Cruise_Plan_ESI_2022.png",ESI_plan, width=8, height=8, units="in",dpi=300)



## St. Anns Bank -------------
        
        #Create a 6.5 hour radius around the entry to Louisburg to define the maximum sample distance assuming travel at 9.3 knots
        maxdist <- 5.5*9.3*1.852
        
        perley_buffer <- waypoints_cb%>%
                          filter(name=="WP3")%>%
                          st_transform(utm)%>%
                          st_buffer(dist=maxdist)%>%
                          st_as_sf()%>%
                          st_transform(latlong)
                          
        
        ggplot()+
                geom_sf(data=sab_benthoscape,aes(fill=class))+
                geom_sf(data=sab,fill=NA)+
                geom_sf(data=perley_buffer,fill=NA)+
                theme_bw()+
                theme(legend.position = "bottom")


        sab_dem <- dem%>%
                   crop(.,extent(sab%>%st_transform(proj4string(dem))))%>%        
                   mask(.,sab%>%st_transform(proj4string(dem))%>%as_Spatial())%>%
                   projectRaster(.,crs=latlong)
        
        #find Scattarie 
        #coordinates from charts - https://www.marineregions.org/gazetteer.php?p=details&id=34185
        scattarie_bank <- data.frame(lon=-59.25,lat=45.966667,name="scattarie_bank",remove=FALSE)%>%
                                st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                                st_transform(proj4string(dem))%>%
                                mutate(depth=extract(dem,as_Spatial(.)))%>%
                                st_transform(latlong)
        
        #create a buffer to find the shallow point
        scattarie_shallow_buff <- scattarie_bank%>%
                                 st_transform(utm)%>%
                                 st_buffer(10)%>%
                                 st_transform(proj4string(dem))%>%
                                 st_as_sf()
        
        scattarie_shallow_dem <- dem%>%
                                 crop(.,extent(scattarie_shallow_buff))%>%
                                 mask(.,scattarie_shallow_buff%>%as_Spatial())
        
        scattarie_shallow <- xyFromCell(scattarie_shallow_dem,which(scattarie_shallow_dem[] == min(values(scattarie_shallow_dem),na.rm=T)))%>%
                             data.frame()%>%
                             rename(lon=x,lat=y)%>%
                             st_as_sf(coords=c("lon","lat"),crs=proj4string(dem))%>%
                             mutate(depth=extract(dem,as_Spatial(.)))%>%
                             st_transform(latlong)%>%
                             mutate(lon=st_coordinates(.)[,1],
                                    lat=st_coordinates(.)[,2],
                                    name="Scattarie Bank")
        
        #Now construct transects from shallow to deep starting with anchor points in the LC to the north and northeast 
        sab_deep <- data.frame(lon=c(-58.96635,-58.81155),
                               lat=c(46.3443,46.21858),
                               name=c("northeast anchor","north anchor"))%>%
                    st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)
        
        #Depth transition transect 1 
        trans1 <- rbind(scattarie_shallow%>%dplyr::select(name,lon,lat,geometry),
                             sab_deep[1,]%>%dplyr::select(name,lon,lat,geometry))
        
        trans1_dist <- st_distance(trans1)%>%as.numeric()%>%max()/1000
        
        trans1_bear <- bearing(scattarie_shallow%>%as_Spatial,sab_deep[1,]%>%as_Spatial)
        
        trans2 <- rbind(scattarie_shallow%>%dplyr::select(name,lon,lat,geometry),
                        sab_deep[2,]%>%dplyr::select(name,lon,lat,geometry))
        
        trans2_dist <- st_distance(trans2)%>%as.numeric()%>%max()/1000
        
        trans2_bear <- bearing(scattarie_shallow%>%as_Spatial,sab_deep[2,]%>%as_Spatial)
        
        #both close to 45 km so we will extend to that distance with equally spaced with closer within the bank
        n <- 5
        transect_spaceing <- c(1,2.5,45/n*1:n)*1000
        
        coord_names <- c("name","id","depth","lon","lat")
        
        transect_coords <- rbind(geosphere::destPoint(p=scattarie_shallow%>%as_Spatial(),
                                                b=trans1_bear,
                                                d=transect_spaceing)%>%
                                                data.frame()%>%
                                                st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                                                st_transform(proj4string(dem))%>%
                                                mutate(name="Transect 1",
                                                       depth=extract(dem,as_Spatial(.)),
                                                       id=1:n())%>%
                                                st_transform(latlong)%>%
                                         dplyr::select(all_of(coord_names)),
                                   geosphere::destPoint(p=scattarie_shallow%>%as_Spatial(),
                                                 b=trans2_bear,
                                                 d=transect_spaceing)%>%
                                                data.frame()%>%
                                                st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                                                st_transform(proj4string(dem))%>%
                                                mutate(name="Transect 2",
                                                       depth=extract(dem,as_Spatial(.)),
                                                       id=1:n())%>%
                                                st_transform(latlong)%>%
                                                dplyr::select(all_of(coord_names)),
                                 scattarie_shallow%>%
                                         mutate(id=NA)%>%
                                         dplyr::select(all_of(coord_names)))%>%
                            mutate(class=sab_benthoscape$class[as.numeric(st_within(.,sab_benthoscape))],
                                   dist_scat = as.numeric(round(st_distance(.,scattarie_shallow)/1000)),
                                   day_group = case_when(dist_scat<10 ~ "Day 1",
                                                         dist_scat>30 & name == "Transect 1" ~ "Day 2",
                                                         dist_scat>30 & name == "Transect 2" ~ "Day 3",
                                                         dist_scat>10 & dist_scat<40 & name == "Transect 1" ~ "Day 4",
                                                         dist_scat>10 & dist_scat<40 & name == "Transect 2" ~ "Day 5"))
        
        
        
        ggplot()+
        geom_sf(data=sab_benthoscape,aes(fill=class))+
        geom_sf(data=sab_zones,fill=NA)+
        geom_sf(data=transect_coords,aes(col=name))+
        geom_sf(data=perley_buffer,fill=NA)+
        coord_sf(expand=0,xlim=sab_bounds[c(1,3)],ylim=sab_bounds[c(2,4)])+
        theme_bw()+
        scale_color_manual(values=c("white","darkblue","black"))
        
        
        scatterie_bounds <- scattarie_shallow%>%
                            st_transform(utm)%>%
                            st_buffer(10)%>%
                            st_transform(latlong)%>%
                            st_bbox()
        
        
        ggplot()+
                geom_sf(data=sab_zones)+
                geom_sf(data=sab_benthoscape,aes(fill=class))+
                geom_sf(data=transect_coords,aes(fill=name))+
                coord_sf(expand=0,xlim=scatterie_bounds[c(1,3)],ylim=scatterie_bounds[c(2,4)])+
                theme_bw()
        
        #how much area for each class
        sab_benthoscape%>%
                mutate(area=round(as.numeric(st_area(.)/1000/1000),2))%>%
                data.frame()%>%
                dplyr::select(class,area)
                
         

#St. Anns Bank Zonal comparison -------
        
zone_spacing <- 5   
        
zone2b <- sab_zones%>%filter(Zone == "2b")

zone2b_centre <- st_centroid(zone2b)%>% #shift up a bit
                 as_Spatial()%>%
                 geosphere::destPoint(p=.,b=0,d=1500)%>%
                 data.frame()%>%
                 st_as_sf(coords=c("lon","lat"),crs=latlong)


zone2b_centre_reference <- geosphere::destPoint(p=as_Spatial(zone2b_centre),b=270,d=c(23*1000,46*1000))%>% #this puts coordinates within zones 1 and 2a
        data.frame()%>%
        st_as_sf(coords=c("lon","lat"),crs=latlong)
        

zone_grid <- rbind(
        box_grid(zone2b_centre,dist=zone_spacing ,full=FALSE)%>%
                mutate(zone="2b",
                       type="fished"),
        box_grid(zone2b_centre_reference[1,],dist=zone_spacing ,full=FALSE)%>%
                mutate(zone="1",
                       type="core"),
        box_grid(zone2b_centre_reference[2,],dist=zone_spacing ,full=FALSE)%>%
                mutate(zone="2a",
                       type="fished")
        )%>%
        mutate(class=sab_benthoscape$class[as.numeric(st_within(.,sab_benthoscape))])%>%
        st_transform(proj4string(dem))%>%
        mutate(depth=extract(dem,as_Spatial(.)),
               day_group = case_when(zone=="2b" ~ "Day 6",
                                     zone=="2a" ~ "Day 7",
                                     zone=="1"~ "Day 8"))%>%
        st_transform(latlong)

#evaluate the coverage of each benthoscape classification
table(zone_grid$class,zone_grid$zone)


#plot it out in detail overlaid on the benthoscape
zone_bounds <- zone_grid%>%
        st_transform(utm)%>%
        st_bbox()%>%
        st_as_sfc()%>%
        st_as_sf()%>%
        st_buffer(10)%>%
        st_as_sf()%>%
        st_transform(latlong)%>%
        st_bbox()

inside_outside_comp <- ggplot()+
        geom_sf(data=sab_benthoscape,aes(fill=class))+
        geom_sf(data=sab_zones,fill=NA,lwd=1.2)+
        geom_sf(data=zone_grid)+
        theme_bw()+
        theme(legend.position="bottom")+
        coord_sf(xlim=zone_bounds[c(1,3)],ylim=zone_bounds[c(2,4)],expand=0)

ggsave("2022 Water Sampling/inside_outside_comp.png",inside_outside_comp,height=6,width=8,units="in",dpi=300)


#Distances for each leg ----------------

#Inside outside comparison 


sab_map_bound <- rbind(louisbourg%>%dplyr::select(geometry),
                       zone_grid%>%dplyr::select(geometry),
                       transect_coords%>%dplyr::select(geometry))%>%
              st_transform(utm)%>%
              st_bbox()%>%
              st_as_sfc()%>%
              st_as_sf()%>%
              st_buffer(3)%>%
              st_as_sf()%>%
              st_transform(latlong)%>%
              rbind(.,sab)%>%
              st_bbox()

#inside_outside comparisons

inside_outside_lines <- NULL

for(i in unique(zone_grid$day_group)){
  
  temp <- zone_grid%>%
          filter(day_group==i)%>%
          mutate(lon=st_coordinates(.)[,1],
                 lat=st_coordinates(.)[,2])%>%
          arrange(-lat)%>%
          mutate(name = day_group)
  
  #gets them in the right order

  if(temp$lon[1]>temp$lon[2]){temp <- temp[c(2,1,3,4,5),]}
  if(temp$lon[4]<temp$lon[5]){temp <- temp[c(1,2,3,5,4),]}
  
  temp <- temp%>%dplyr::select(c("name","geometry"))
  
  
  #fix the ordering here **** day 8 isn't working correctly ****************
  
  temp_line <- rbind(louisbourg%>%dplyr::select(c("name","geometry")),
                 waypoints_cb%>%dplyr::select(c("name","geometry")),
                 temp%>%dplyr::select(c("name","geometry")),
                 waypoints_cb%>%
                   mutate(order=1:n())%>%
                   arrange(-order)%>%
                   dplyr::select(c("name","geometry")),
                 louisbourg%>%dplyr::select(c("name","geometry")))%>%
            summarise(do_union=FALSE)%>%
            st_cast("LINESTRING")%>%
            mutate(length=as.numeric(st_length(.))/1000/1.852,
                   day_group=i)
  
  inside_outside_lines <- rbind(inside_outside_lines,temp_line)
  
}

ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=louisbourg)+
  geom_sf(data=waypoints_cb)+
  geom_sf(data=sab_zones)+
  geom_sf(data=inside_outside_lines)+
  coord_sf(expand=0,xlim=sab_map_bound[c(1,3)],ylim=sab_map_bound[c(2,4)])


#St. Anns Bank Depth transects -

#Day 1
trans_line_day1 <- transect_coords%>%
                   filter(day_group=="Day 1")

start_day1 <- st_nearest_feature(louisbourg,trans_line_day1)

trans_line1_day1 <- trans_line_day1%>%
                    filter(name=="Transect 1")%>%
                    arrange(lat)

trans_line2_day1 <- trans_line_day1%>%
                    filter(name=="Transect 2")%>%
                    arrange(-lat)

transect_line_day1 <- rbind(louisbourg%>%dplyr::select(c("name","geometry")),
                         waypoints_cb%>%dplyr::select(c("name","geometry")),
                         trans_line_day1[start_day1,]%>%dplyr::select(c("name","geometry")),
                         trans_line1_day1%>%dplyr::select(c("name","geometry")),
                         trans_line2_day1%>%dplyr::select(c("name","geometry")),
                         waypoints_cb%>%
                           mutate(order=1:n())%>%
                           arrange(-order)%>%
                           dplyr::select(c("name","geometry")),
                         louisbourg%>%dplyr::select(c("name","geometry")))%>%
                    summarise(do_union=FALSE)%>%
                      st_cast("LINESTRING")%>%
                      mutate(length=as.numeric(st_length(.))/1000/1.852,
                             day_group="Day 1")


#Day 2
trans_line_day2 <- transect_coords%>%
  filter(day_group %in% c("Day 4","Day 5"))

start_day2 <- st_nearest_feature(louisbourg,trans_line_day2)

trans_line1_day2 <- trans_line_day2%>%
  filter(name=="Transect 1")%>%
  slice(-start_day2)

trans_line2_day2 <- trans_line_day2%>%
  filter(name=="Transect 2")%>%
  arrange(-lat)

transect_line_day2 <- rbind(louisbourg%>%dplyr::select(c("name","geometry")),
                            waypoints_cb%>%dplyr::select(c("name","geometry")),
                            trans_line_day2[start_day1,]%>%dplyr::select(c("name","geometry")),
                            trans_line1_day2%>%dplyr::select(c("name","geometry")),
                            trans_line2_day2%>%dplyr::select(c("name","geometry")),
                            waypoints_cb%>%
                              mutate(order=1:n())%>%
                              arrange(-order)%>%
                              dplyr::select(c("name","geometry")),
                            louisbourg%>%dplyr::select(c("name","geometry")))%>%
  summarise(do_union=FALSE)%>%
  st_cast("MULTIPOINT")%>%
  st_cast("LINESTRING")%>%
  mutate(length=as.numeric(st_length(.))/1000/1.852,
         day_group="Day 2")

#Day 3
transect_line_day3 <- rbind(louisbourg%>%dplyr::select(c("name","geometry")),
                            waypoints_cb%>%dplyr::select(c("name","geometry")),
                            transect_coords%>%
                              filter(day_group=="Day 2")%>% #old naming
                              arrange(-lat)%>%
                              dplyr::select(c("name","geometry")),
                            waypoints_cb%>%
                              mutate(order=1:n())%>%
                              arrange(-order)%>%
                              dplyr::select(c("name","geometry")),
                            louisbourg%>%dplyr::select(c("name","geometry")))%>%
  summarise(do_union=FALSE)%>%
  st_cast("LINESTRING")%>%
  mutate(length=as.numeric(st_length(.))/1000/1.852,
         day_group="Day 3")
                   

#Day 4
transect_line_day4 <- rbind(louisbourg%>%dplyr::select(c("name","geometry")),
                            waypoints_cb%>%dplyr::select(c("name","geometry")),
                            transect_coords%>%
                              filter(day_group=="Day 3")%>% #old naming
                              arrange(-lat)%>%
                              dplyr::select(c("name","geometry")),
                            waypoints_cb%>%
                              mutate(order=1:n())%>%
                              arrange(-order)%>%
                              dplyr::select(c("name","geometry")),
                            louisbourg%>%dplyr::select(c("name","geometry")))%>%
  summarise(do_union=FALSE)%>%
  st_cast("LINESTRING")%>%
  mutate(length=as.numeric(st_length(.))/1000/1.852,
         day_group="Day 4")

#bind them together
transect_lines <- rbind(transect_line_day1,
                        transect_line_day2,
                        transect_line_day3,
                        transect_line_day4)


#combine SAB transect and inside_outside comparison 

sab_sample_plan_lines <- rbind(transect_lines%>%
                                 mutate(type="Depth transect",
                                        lab=paste0("Depth transect ",day_group," (~",round(length,2),"nm)"))%>%
                                 dplyr::select(type,lab,length,geometry),
                               inside_outside_lines%>%
                           mutate(day_group <- c("Day 5","Day 6","Day 7"),
                                  type="Zonal comparison",
                                  lab=paste0("Zonal comparison ",day_group," (~",round(length,2),"nm)"))%>%
                           dplyr::select(type,lab,length,geometry))

#construct the points
sab_sample_plan_points=NULL
for(i in unique(sab_sample_plan_lines$lab)){
  
  temp <- sab_sample_plan_lines%>%
           filter(lab==i)%>%
           st_cast("MULTIPOINT")%>%
           st_intersection(sab)%>%suppressWarnings()
  
  sab_sample_plan_points <- rbind(sab_sample_plan_points,temp)
  
  
}

all_points <- sab_sample_plan_points%>%
              st_union()%>%
              st_as_sf()%>%
              rename(geometry=x)%>%
              mutate(lab="St Anns Bank 2022")


sab_sample_plan_points_map <- rbind(sab_sample_plan_points%>%dplyr::select(lab,geometry),
                                    all_points%>%dplyr::select(lab,geometry))%>%
                              mutate(lab=factor(lab,levels=c(sab_sample_plan_points$lab,all_points$lab)))

sab_sample_plan_lines$lab <- factor(sab_sample_plan_lines$lab,levels=c(sab_sample_plan_points$lab,all_points$lab))

sab_edna_map <- ggplot()+
                geom_sf(data=basemap,fill="darkolivegreen3")+
                #geom_sf(data=sab_benthoscape,aes(fill=class))+
                geom_sf(data=louisbourg)+
                geom_sf(data=waypoints_cb)+
                geom_sf(data=sab_zones,fill=NA)+
                geom_sf(data=sab_sample_plan_points_map)+
                geom_sf(data=sab_sample_plan_lines)+
                facet_wrap(~lab,nrow=4)+
                theme_bw()+
                coord_sf(expand=0,xlim=sab_map_bound[c(1,3)],ylim=sab_map_bound[c(2,4)])+
                theme(legend.position="none",
                      strip.background = element_rect(fill="white"),
                      axis.text=element_text(size=6))

ggsave("2022 Water Sampling/sab_edna_map.png",sab_edna_map,height=8,width=8,units="in",dpi=300)


                
##plot them all ----------- 

sab_samples <- ggplot()+
        geom_sf(data=basemap)+
        geom_sf(data=louisbourg,pch=21,fill="white")+
        geom_sf(data=sab_zones,fill=NA,lwd=1.05)+
        geom_sf(data=zone_grid,aes(fill=day_group),pch=21,size=2)+
        geom_sf(data=transect_coords,aes(fill=day_group),pch=21,size=2)+
        theme_bw()+
        theme()+
        coord_sf(expand=0,xlim=sab_map_bound[c(1,3)],ylim=sab_map_bound[c(2,4)])+
        scale_fill_viridis(discrete = TRUE)+
        labs(fill="")

ggsave("2022 Water Sampling/sab_samples_collective_2022.png",sab_samples,width=8,height=6,units="in",dpi=300)


#Coordinate sheet ----------------

acoustic <- esi_stations%>%
            filter(retrieve)%>%
            data.frame()%>%
            mutate(deg_lon = trunc(lon),
                   dm_lon = abs(lon-deg_lon)*60,
                   deg_lat = trunc(lat),
                   dm_lat = (lat-deg_lat)*60)%>%
            dplyr::select(station,lon,deg_lon,dm_lon,lat,deg_lat,dm_lat)

edna <- edna_stations%>%
        mutate(lon=st_coordinates(.)[,1],
               lat=st_coordinates(.)[,2],
               deg_lon = trunc(lon),
               dm_lon = abs(lon-deg_lon)*60,
               deg_lat = trunc(lat),
               dm_lat = (lat-deg_lat)*60)%>%
        data.frame()%>%
        dplyr::select(station,lon,deg_lon,dm_lon,lat,deg_lat,dm_lat)
        

sab_trans <- transect_coords%>%
             mutate(id2 = ifelse(is.na(id),1,id+1),
                    id2 = ifelse(name=="Transect 2",id+8,id2),
                    station = paste0("SAB_trans_",id2),
                    deg_lon = trunc(lon),
                    dm_lon = abs(lon-deg_lon)*60,
                    deg_lat = trunc(lat),
                    dm_lat = (lat-deg_lat)*60)%>%
             arrange(id2)%>%
             data.frame()%>%
             dplyr::select(station,lon,deg_lon,dm_lon,lat,deg_lat,dm_lat)

sab_zone <- zone_grid%>%
            mutate(station = paste0("SAB_zone_",1:15),
                   lon = st_coordinates(.)[,1],
                   lat = st_coordinates(.)[,2],
                   deg_lon = trunc(lon),
                   dm_lon = abs(lon-deg_lon)*60,
                   deg_lat = trunc(lat),
                   dm_lat = (lat-deg_lat)*60)%>%
            data.frame()%>%
            dplyr::select(station,lon,deg_lon,dm_lon,lat,deg_lat,dm_lat)
             

sampling_2022_coords <- rbind(acoustic,edna,sab_trans,sab_zone)%>%
                        mutate(lon=round(lon,4),
                               lat=round(lat,4),
                               dm_lon=round(dm_lon,4),
                               dm_lat=round(dm_lat,4))

write.csv(x=sampling_2022_coords,file="2022 Water Sampling/sampling_2022_coords.csv",row.names=F)




#making a map for the approval email to SAB

sab_dem <- dem%>%
          crop(.,extent(sab%>%st_transform(proj4string(dem))))%>%
          disaggregate(.,fact=8)%>% #increase the resolution so it doesn't look so choppy with the mask
          mask(.,sab%>%st_transform(proj4string(dem))%>%as_Spatial())%>%
          projectRaster(.,crs=latlong)%>%
          st_as_stars()

sab_map_bound2 <- sab%>%
                 st_bbox()%>%
                 st_as_sfc()%>%
                 st_as_sf()%>%
                 st_buffer(10)%>%
                 st_bbox()%>%
                 suppressWarnings()

sab_map_bound2[1] <- -59.9

sab_dem_form <- sab_dem
sab_dem_form$sab_dem <- sab_dem$sab_dem*-1
names(sab_dem_form) <- "Depth (m)"


contour_250 <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Bathymetry/Countour_250.shp")%>%
               st_transform(latlong)%>%
               st_make_valid()%>%
               st_intersection(sab_bound)

coords <- read.csv("r:/Science/CESD/HES_MPAGroup/Projects/eDNA/2022 Water Sampling/sampling_2022_coords.csv")%>%
          filter(grepl("SAB",station))%>%
          dplyr::select(station,lon,lat)%>%
          st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
          st_transform(proj4string(dem))%>%
          mutate(depth=round(extract(dem,as_Spatial(.)),2))%>%
          st_transform(latlong)%>%
          mutate(type=ifelse(grepl("trans",station),"Depth transect","Zonal comparison"))

sab_edna_plan

sab_approval_plot <- ggplot()+
                      geom_sf(data=contour_250,fill=NA,col="grey50")+
                      geom_sf(data=basemap,fill="darkolivegreen3")+
                      geom_sf(data=louisbourg)+
                      geom_sf(data=sab_zones,fill=NA)+
                      geom_sf(data=coords,aes(fill=type),pch=21,col="black",size=1.5)+
                      theme_bw()+
                      coord_sf(xlim=sab_map_bound[c(1,3)],ylim=sab_map_bound[c(2,4)])+
                      labs(fill="")+
                      theme(legend.position = c(0.15,0.85),
                            legend.background = element_rect(fill="white"),
                            legend.box.background = element_rect(colour = "black"),
                            legend.title = element_blank())

ggsave("R:/Science/CESD/HES_MPAGroup/Projects/eDNA/sab_approval_plot.png",sab_approval_plot,height=5,width=7,dpi=300,units="in")


sab_tweet_plot <- ggplot()+
                  geom_stars(data=sab_dem_form)+
                  scale_fill_viridis(option="C",na.value="white")+
                  geom_sf(data=basemap,fill="darkolivegreen3")+
                  geom_sf(data=sab_zones,fill=NA)+
                  new_scale_fill() + 
                  geom_sf(data=coords,aes(fill=type),pch=21,col="black",size=2.25)+
                  scale_fill_manual(values=c("brown2","chartreuse1"))+
                  theme_bw()+
                  theme(#axis.text = element_blank(),
                        axis.title = element_blank())+
                  coord_sf(expand=0,xlim=sab_map_bound2[c(1,3)],ylim=sab_map_bound2[c(2,4)])+
                  labs(fill="")+
                  annotation_scale(location="br")+
                  annotation_north_arrow(location="tl")

ggsave("r:/Science/CESD/HES_MPAGroup/Projects/eDNA/sab_map_bathy.png",sab_tweet_plot,width=8,height=8,units="in",dpi=300)

#output coordinates for sab
sab_coords_output <- coords%>%
                     data.frame()%>%
                     dplyr::select(station,lon,lat,depth)

write.csv(sab_coords_output,"r:/Science/CESD/HES_MPAGroup/Projects/eDNA/2022 Water Sampling/sab_coordinates_approval.csv",row.names=FALSE)

#ESI plot

#code to have two colour scales https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/
devtools::source_gist("https://gist.github.com/eliocamp/eabafab2825779b88905954d84c82b32")

coords_ESI <- rbind(read.csv("r:/Science/CESD/HES_MPAGroup/Projects/eDNA/2022 Water Sampling/sampling_2022_coords.csv")%>%
                      filter(grepl("ESI",station))%>%
                      mutate(type="eDNA")%>%
                      dplyr::select(station,type,lon,lat),
                    read.csv("r:/Science/CESD/HES_MPAGroup/Projects/Acoustic/ESI/2022_stations_depth.csv")%>%
                      mutate(type="Acoustic")%>%
                      rename(lon=lon_dec_deg,lat=lat_dec_deg)%>%
                      dplyr::select(station,type,lon,lat))%>%
              st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)
  
bathy_esi <- bathy%>%
             crop(.,extent(esi_poly%>%st_transform(proj4string(bathy))))%>%
             mask(.,esi_poly%>%st_transform(proj4string(bathy))%>%as_Spatial())%>%
             projectRaster(.,crs=latlong)%>%
             st_as_stars()

bathy_esi_form <- bahthy_esi
bathy_esi_form$dem35c_5c <- bathy_esi$dem35c_5c*-1 # only do this once!
names(bathy_esi_form) <- "Depth (m)"

coast_hr_crop2 <- coast_hr%>%
  st_intersection(st_bbox(esi_poly)%>%st_as_sfc()%>%st_as_sf()%>%st_transform(st_crs(coast_hr)))%>%
  st_transform(latlong)

p1 <- ggplot()+
  geom_sf(data=esi_poly,fill=NA)+
  geom_stars(data=bathy_esi_form)+
  scale_fill_viridis(option="C",na.value="white")+
  new_scale_fill() + 
  geom_sf(data=coast_hr_crop2,fill="darkolivegreen3")+
  geom_sf(data=coords_ESI,aes(fill=type),shape=21,size=2)+
  scale_fill_manual(values=c("brown2","chartreuse1"))+
  theme_bw()+
  coord_sf(expand=0)+
  theme(axis.text=element_blank())+
  labs(x="",y="",col="")+
  annotation_scale(location="br")+
  annotation_north_arrow(location="tl")

ggsave("r:/Science/CESD/HES_MPAGroup/Projects/eDNA/2022 Water Sampling/sample_map.png",p1,width=8,height=8,units="in",dpi=300)
