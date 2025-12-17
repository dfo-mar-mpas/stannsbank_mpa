#load libraries ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(stars)
library(ggspatial)
library(ggnewscale)
library(viridis)
library(patchwork)
library(dbscan)
library(cmocean)

s2_as_sf = FALSE


#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load bathymetry ------
dem_sab <- rast("data/Bathymetry/sab_dem.tif") # SAB based on the 'predSpace' data
load("data/Bathymetry/250m_stars_object.RData") # STARS object based on the 250m depth contour extracted (converted) from GEBCO

#load St Anns Bank MPA shapefile 
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()

#Load the bigh resolution bathymetry 
sab_highres <- rast("data/Bathymetry/SAB_highres/fullBathy_projected_latlong.tif")%>%
               crop(.,st_bbox(sab))

bathy_proj <- st_crs(sab_highres)

#load the banks shapfile
sab_banks <- read_sf("data/Shapefiles/sab_banks_40m.shp")%>%
             st_transform(latlong)%>%
             st_polygonize()

curdo <- sab_banks%>%
         filter(name=="Curdo Bank")%>%
         st_simplify()%>%
         st_make_valid()

scatarie <- sab_banks%>%
  filter(name=="Scatarie Bank")%>%
  st_simplify()%>%
  st_make_valid()

#bank-based bathymetry
curdo_bathy <- sab_highres%>%
               crop(.,st_bbox(curdo))%>%
               terra::mask(.,curdo)

scatarie_bathy <- sab_highres%>%
                  crop(.,st_bbox(scatarie))%>%
                  terra::mask(.,scatarie)

#trim to 35m for the plot
curdo_35m_contour <- as.contour(curdo_bathy, levels = -35)%>%
                      st_as_sf()%>%
                      st_set_crs(latlong)%>%
                      st_cast("POLYGON")%>%
                      mutate(name="Curdo Bank")

scatarie_35m_contour <- as.contour(scatarie_bathy, levels = -35)%>%
                        st_as_sf()%>%
                        st_set_crs(latlong)%>%
                        st_cast("POLYGON")%>%
                        mutate(name="Scatarie Bank")

curdo_bathy_35 <- curdo_bathy%>%
                  terra::mask(.,curdo_35m_contour)

scatarie_bathy_35 <- scatarie_bathy%>%
                    terra::mask(.,scatarie_35m_contour)




#load the reciever locations
reciever_locs_old <- read.csv("data/Acoustic/OTN_redesign_coords.csv")%>%
                 st_as_sf(coords=c("long","lat"),crs=latlong)%>%
                  st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
                  mutate(lon=round(st_coordinates(.)[,1],3),
                         lat=round(st_coordinates(.)[,2],3),
                         type="Acoustic Reciever Locations",)%>%
                  dplyr::select(type,Zone,lon,lat,geometry)

reciever_locs_new <- read.csv("data/Acoustic/otn-instrument-deployment-SABMPAnew_June2025.csv")%>%
                     st_as_sf(coords=c("DEPLOY_LONG_DD","DEPLOY_LAT_DD"),crs=latlong)%>%
                      st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
                      mutate(lon=round(st_coordinates(.)[,1],3),
                             lat=round(st_coordinates(.)[,2],3),
                             type="New Acoustic Reciever Locations")%>%
                      dplyr::select(type,Zone,lon,lat,geometry)

reciever_locs <- rbind(reciever_locs_old,reciever_locs_new)%>%mutate(type="Acoustic Receiver")

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

#plot_boundaries <- c(c(-61,-58,45.6,46.6))



#tagging locations -----

tag_2022 <- read.csv("data/Acoustic/Halibut tagging sheet Nov 2022 OTN.csv")%>%
                filter(!is.na(TAG_SERIAL_NUMBER))%>% #tailing blank spaces in the Excel sheet - this will be inclusive of cod and halibut where sets were deployed
                mutate(lon.deg=as.numeric(substr(RELEASE_LONGITUDE, 1, 2)),
                       lon.min=RELEASE_LONGITUDE - (lon.deg*100),
                       lon = (lon.deg + lon.min/60)*-1,
                       lat.deg=as.numeric(substr(RELEASE_LATITUDE, 1, 2)),
                       lat.min=RELEASE_LATITUDE - (lat.deg*100),
                       lat = lat.deg + lat.min/60,
                       year=2022)%>%
                st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                rename(species=COMMON_NAME_E)%>%
                dplyr::select(lon,lat,year,species,geometry)

tag_2023_25 <- read.csv("data/Acoustic/SABMPA_tagging_basic_metadata_all_QCApril2025.csv")%>%
                  st_as_sf(coords=c("long","lat"),crs=latlong,remove=FALSE)%>%
                  rename(lon=long)%>%
                  dplyr::select(lon,lat,year,species,geometry)

tag_df <- rbind(tag_2022,tag_2023_25)


#create a tag table (this will summarize the points)
coords <- tag_df%>%
          st_transform(utm)%>%
          st_coordinates()
          
clusters <- dbscan(coords, eps = 1, minPts = 1)

tag_df_summary <- tag_df %>%
  mutate(cluster_id = clusters$cluster) %>%
  group_by(cluster_id, year) %>%
  summarise(
    n_points = n(),
    geometry = st_centroid(st_union(geometry)),
    .groups = "drop"
  ) %>%
  st_as_sf(keep=T)%>%
  st_intersection(sab_zones)%>%
  arrange(year,Zone)%>%
  mutate(Longitude = st_coordinates(.)[,1],
         Latitude = st_coordinates(.)[,2])%>%
  rename(Year=year)

write.csv(tag_df_summary%>%
            data.frame()%>%
            dplyr::select(Year,Zone,Longitude,Latitude),
          "output/Acoustic/2025_tag_summary_table.csv")



#save coordinates for the activity approval applcation - to simplify the included table, the 2022 tag locations will be presented as centriods of the three clusters
mod_sets <- tag_2022%>%
            mutate(cluster = case_when(lon> -59.26 ~ 1, #rougly 4 clusters
                                       lon> -59.3 & lon < (-59.26) ~ 2,
                                       lon< (-59.5) ~ 3,
                                       TRUE ~ 4))%>%
            group_by(cluster)%>%
            summarise(geometry = st_union(geometry))%>%
            st_centroid()%>%
            ungroup()%>%
            st_as_sf()%>%
            mutate(lon=round(st_coordinates(.)[,1],3),
                   lat=round(st_coordinates(.)[,2],3),
                   year=2022)%>%
            st_intersection(.,sab_zones%>%dplyr::select(Zone,geometry))%>%
            suppressWarnings()

tag_df <- rbind(tag_2022%>%dplyr::select(year,species,geometry),tag_2023_25%>%dplyr::select(year,species,geometry))%>%
          st_intersection(.,sab)%>%
          mutate(type="Tagging location")

#make a curdo and a scatarie plot 
curdo_bound <- sab_banks%>%
  filter(name=="Curdo Bank")%>%
  st_transform(utm)%>%
  st_buffer(0.5)%>%
  st_transform(latlong)%>%
  st_bbox()

scatarie_bound <- sab_banks%>%
  filter(name=="Scatarie Bank")%>%
  st_transform(utm)%>%
  st_buffer(0.75)%>%
  st_transform(latlong)%>%
  st_bbox()


curdo_plot <- ggplot() +
  geom_sf(data=shelfbreak,fill=NA,lwd=0.2)+
  geom_sf(data=basemap)+
  geom_spatraster(data=curdo_bathy_35)+
  geom_sf(data=reciever_locs,fill="white",pch=21,size=2)+
  geom_sf(data=tag_df_summary,shape=21,col="coral")+
  coord_sf(expand=0,xlim=curdo_bound[c(1,3)],ylim=curdo_bound[c(2,4)])+
  scale_fill_cmocean(name = "deep", direction = -1, na.value = NA)+
  theme_bw()+
  theme(axis.text=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),
        legend.background = element_rect(color="black"),
        legend.key.height = unit(0.3/2, "cm"),
        legend.key.width = unit(0.25/2, "cm"),
        legend.text = element_text(size = 5*0.8),
        legend.title = element_text(size = 5*0.7))+
  annotation_scale(location="bl")+
  labs(title="Curdo Bank",
       fill="Depth (m)")

scatarie_plot <- ggplot() +
                geom_sf(data=shelfbreak,fill=NA,lwd=0.2)+
                geom_sf(data=basemap)+
                #geom_sf(data=sab_banks,fill=NA)+
                geom_spatraster(data=scatarie_bathy_35)+
                geom_sf(data=reciever_locs,fill="white",pch=21,size=2)+
                #geom_sf(data=tag_df,size=0.5,pch=3)+
                geom_sf(data=tag_df_summary,shape=21,col="coral")+
                coord_sf(expand=0,xlim=scatarie_bound[c(1,3)],ylim=scatarie_bound[c(2,4)])+
                scale_fill_cmocean(name = "deep", direction = -1, na.value = NA)+
                theme_bw()+
                theme(axis.text=element_blank(),
                      legend.position = "none")+
                annotation_scale(location="bl")+
                labs(title="Scatarie Bank")

plot_boundaries <- c(c(-59.7,-58.3,45.7,46.5))

#make plot

# Calculate hillshade from your DEM
# hill <- terrain(sab_highres, v = "slope", neighbors = 8)
# shade <- terrain(sab_highres, v = "aspect", neighbors = 8)
# 
# # Create hillshade (0-1 scale)
# hillshade <- shade(hill, shade)
# 
# # Convert to data frame for plotting
# hillshade_df <- as.data.frame(hillshade, xy = TRUE)


sab_plot <- ggplot()+
  geom_sf(data=basemap)+
  geom_spatraster(data=sab_highres)+
  geom_sf(data=sab_zones,fill=NA,linewidth=1.1,col="grey60")+
  geom_sf(data=sab_banks,fill="grey",alpha=0.5)+
  geom_sf(data=curdo_bound%>%st_as_sfc(),fill=NA)+
  geom_sf(data=scatarie_bound%>%st_as_sfc(),fill=NA)+
  geom_sf(data=reciever_locs,fill="white",pch=21)+
  geom_sf(data=tag_df,size=0.5,pch=3)+
  geom_sf(data=tag_df_summary,size=2,col="grey60",pch=21)+
  scale_fill_cmocean(name = "deep", direction = -1, na.value = NA)+
  coord_sf(xlim=plot_boundaries[c(1,2)],ylim=plot_boundaries[c(3,4)],expand=0)+
  labs(fill="Depth (m)")+
  theme_bw()+
  theme(
        legend.position = "inside",
        legend.position.inside = c(0.9,0.22),
        legend.background = element_rect(color="black")
        )+
  annotation_scale(location="bl")

mosaic_map <- sab_plot + (scatarie_plot/curdo_plot) + plot_layout(ncol=2,widths=c(3.5,1.8))

ggsave("output/Acoustic/2025_ActivityApproval_taglocations.png",mosaic_map,height=5,width=8,units="in",dpi=300)


