#code to estimate distance from tag location to the St. Anns Bank MPA

#load libraries ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(lubridate)
library(stars)
library(ggspatial)
library(patchwork)
library(ggspatial)
library(viridis)
library(terra)
library(tidyterra)
library(worrms)
library(osmdata)
library(nngeo)
library(sfnetworks)
library(fasterize)
library(gdistance)
library(ggridges)
library(readxl)

s2_as_sf = FALSE

#source helper functions
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/trim_img_ws.R")
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/coord_bump.R")
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/lcp_site.R")

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


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

#basemap of the coast - high enough resolution for this type of analysis -----
basemap <- rbind(
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


#load the reciever array ---- 
reciever_df <- read.csv("data/Acoustic/SAB_deployments_recovered_allFeb2024.csv")%>%
                  rename(lon=DEPLOY_LONG,lat=DEPLOY_LAT)%>%
                  dplyr::select(-geometry)%>%
                  st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                  st_transform(latlong)%>%
                  filter(OTN_ARRAY_2 == "SAB_NEW")%>%
                  distinct(STATION_NO,.keep_all=TRUE)%>%
                  mutate(id = as.numeric(gsub(".*_", "", STATION_NO)))%>%
                  arrange(id)

sab_array_centre <- reciever_df%>%
                    filter(id == 23)


#Qualified Detection Locations -----
tag_df <- read.csv("data/Acoustic/qdet_releaselocs.csv")%>%
  rename(lon=RELEASE_LONGITUDE,lat=RELEASE_LATITUDE)%>%
  filter(!is.na(lon))%>%
  mutate(lon = ifelse(RELEASE_LOCATION == "Bay of Quinte, Lake Ontario",-69.787025,lon), #bay of Quinte shows up in 'water' but the great lakes are landlocked with the basemap we use. 
         lat = ifelse(RELEASE_LOCATION == "Bay of Quinte, Lake Ontario",47.773176,lat))%>%
  st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
  mutate(group=case_when(Common.Name %in% c("White shark","Shortfin mako","Blue shark","Blue shark/Mako/Bluefin tuna","Porbeagle shark") ~ "Sharks",
                         Common.Name %in% c("Atlantic sturgeon","Bluefin tuna")~"Pelagics",
                         Common.Name %in% c("Atlantic cod","Snow crab","Atlantic halibut") ~ "Atlantic fisheries",
                         Common.Name %in% c("Atlantic salmon","American eel") ~ "Diadramous",
                         Common.Name %in% c("Grey seal","Leatherback turtle") ~ "Large pelagics",
                         TRUE ~ NA),
         Common.Name = factor(Common.Name,levels=c("American eel","Atlantic cod","Atlantic halibut","Atlantic salmon","Atlantic sturgeon",
                                                   "Blue shark","Blue shark/Mako/Bluefin tuna","Bluefin tuna","Grey seal","Leatherback turtle",           
                                                   "Porbeagle shark","Shortfin mako","Snow crab","White shark")),
         id = 1:n())





#simple coastaline for the distance analysis ----

tag_bound <- tag_df%>%
  st_transform(3857) %>% #  web mercator (m)
  st_bbox()%>%
  st_as_sfc()%>%
  st_buffer(500*1000)%>%
  st_make_valid()%>%
  st_transform(latlong)%>%
  st_bbox()

coastline_tag <- basemap%>%
                 st_intersection(tag_bound%>%st_as_sfc())%>%
                 st_make_valid()%>%
                 st_boundary()

#get a global basemap
global_basemap <- ne_countries(returnclass = "sf",scale="medium")%>%
                  st_make_valid() %>%                 # fix invalid geometries first
                  st_transform(latlong) %>%           # match CRS
                  st_intersection(tag_bound%>%st_as_sfc()) %>%      # trim to extent
                  st_union()
                  


# identify which are on land and how far they are to land 
tag_land <- tag_df%>%
            st_intersection(basemap)

#for the land-based tags. Identify the common marine starting points. Distance to SAB will be based on the marine travel. Gettting accurate distances up river is difficult to get for smaller rivers. 
marine_offsets <- data.frame(location = c("Bras d'Or Lake","St. Mary's River","Morell River","West River","Magaguadavic Basin","St. John Harbour","Cascapedia River","Penobscot River","Kennebec River","LaHave River","Bay of Quinte"),
                             lat = c(46.212952,45.097969,46.427347,44.903232,45.121248, 45.260275,48.175941,44.454093,43.746538,44.277128,47.773176),
                             lon = c(-60.214464,-61.835046,-62.685043,-62.497673,-66.906926,-66.057328,-65.925574,-68.813083,-69.773837,-64.348486, -69.787025)) #this will need to be further 'bumped' but they are now at the outfall

tag_land_marine <- tag_land%>%
                   data.frame()%>%
                   dplyr::select(-geometry,-lon,-lat)%>%
                   mutate(location = case_when(grepl("Bras",RELEASE_LOCATION) ~ "Bras d'Or Lake",
                                               grepl("Mary",RELEASE_LOCATION) ~ "St. Mary's River",
                                               grepl("Morell",RELEASE_LOCATION) ~ "Morell River",
                                               grepl("West River",RELEASE_LOCATION) ~ "West River",
                                               grepl("magag",tolower(RELEASE_LOCATION))~ "Magaguadavic Basin",
                                               grepl("Tobique",RELEASE_LOCATION)~ "St. John Harbour",
                                               grepl("Bucksport",RELEASE_LOCATION) ~ "Penobscot River",
                                               grepl("Penobscot",RELEASE_LOCATION) ~ "Penobscot River",
                                               grepl("Kennebec",RELEASE_LOCATION)~ "Kennebec River",
                                               grepl("Cascap",RELEASE_LOCATION)~ "Cascapedia River",
                                               grepl("LeHave",RELEASE_LOCATION) ~ "LaHave River", #spelling error in the data
                                               grepl("LaHave",RELEASE_LOCATION) ~ "LaHave River",
                                               TRUE ~ "Didn't work"
                          ))%>%
                    left_join(marine_offsets)%>% #add back the marine offset points that can be used in the analysis now. 
                    st_as_sf(coords=c("lon","lat"),crs=latlong,remove=FALSE)%>%
                    dplyr::select(names(tag_df))

lcp_dataframe <- rbind(tag_df%>%filter(!id %in% tag_land_marine$id),
                       tag_land_marine)%>%
                 mutate(site_id = paste("site",id,sep="-"),
                        dist_sab = apply(st_distance(geometry, sab_array_centre), 1, min) / 1000) 



#Now do the least cost path analysis 

tag_df_bbox <- lcp_dataframe%>%
               st_bbox()%>%
               st_as_sfc()%>%
               st_buffer(0.25)%>% # 0.25 degree buffer
               st_bbox()%>%
               st_as_sfc()%>%
               st_as_sf()

#will break this up by distance

 med_scale_lcp <- lcp_site(x=lcp_dataframe%>%filter(dist_sab<=750),#high resolution distance finding for those within 750 km
                           target = sab_array_centre,
                                basemap=global_basemap,
                                bound_box = tag_df_bbox,
                                resolution = 10,
                                transition_name = "tag_df",
                                rad=20,
                                lines=TRUE,
                                dirs=16,
                                recalculate = FALSE,
                                midpoint=FALSE)
 
 large_scale_lcp <- lcp_site(x=lcp_dataframe%>%filter(dist_sab > 750),#high resolution distance finding for those within 750 km
                           target = sab_array_centre,
                           basemap=global_basemap,
                           bound_box = tag_df_bbox,
                           resolution = 10,
                           transition_name = "tag_df",
                           rad=20,
                           lines=TRUE,
                           dirs=16, # just a lower resolution for the distance seeking at the large scale
                           recalculate=FALSE,
                           midpoint=FALSE)

#save interim outputs
save(med_scale_lcp,large_scale_lcp,file="output/Acoustic/lcp_outputs.RData")

#bring the data together

lcp_distances <- lcp_dataframe%>%
                 left_join(.,rbind(med_scale_lcp[[1]],large_scale_lcp[[1]])%>%
                 rename(lcp_dist = dist))

lcp_paths <- rbind(med_scale_lcp[[2]],large_scale_lcp[[2]])%>%st_transform(latlong)


p1 <- ggplot()+
  geom_sf(data=global_basemap)+
  geom_sf(data=sab,fill="cornflowerblue")+
  geom_sf(data=reciever_df)+
  geom_sf(data=sab_array_centre,size=3)+
  geom_sf(data=lcp_paths,col="grey20",linetype=2)+
  geom_sf(data=lcp_dataframe,aes(fill=Common.Name),size=1.5,pch=21)+
  theme_bw()+
  coord_sf(expand=0,xlim=tag_bound[c(1,3)],ylim=tag_bound[c(2,4)])+
  labs(fill="")+
  annotation_scale()+
  theme(
        # legend.position="inside",
        # legend.position.inside = c(0.85,0.5),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.spacing.x = unit(0.05, "cm"));p1

ggsave("output/Acoustic/distance_analysis_qdet_lines.png",p1,height=7,width=5,units="in")
trim_img_ws("output/Acoustic/distance_analysis_qdet_lines.png")  
                 
#plot the distances

breaks = c(10, 30, 100, 300, 1000, 3000)

p2 <- ggplot(lcp_distances,
       aes(x = lcp_dist,
           y = Common.Name)) +
  
  geom_density_ridges_gradient(
    aes(fill = after_stat(x)),   # ← raw distance
    scale = 1.2,
    rel_min_height = 0.01,
    alpha = 0.8
  ) +
  
  geom_point(
    position = position_jitter(height = 0.05),
    size = 1.1,
    alpha = 0.5,
    shape = 21,
    fill = "grey30",
    colour = "black"
  ) +
  
  scale_fill_viridis_c(
    trans = "log10",   # ← log transform here instead
    breaks = breaks,
    labels = scales::label_number()
  ) +
  
  scale_x_log10(
    breaks = breaks,
    labels = scales::label_number()
  ) +
  
  annotation_logticks(sides = "b",
                      short = unit(0.1, "cm"),
                      mid   = unit(0.15, "cm"),
                      long  = unit(0.2, "cm")) +
  theme_bw()+
  labs(x="Marine distance (km)",
       y="")+
  theme(legend.position="none")

ggsave("output/Acoustic/distance_ridgelines_qdet.png",p2,height=6,width=6,dpi=300,units="in")

#save the csv output

write.csv(lcp_distances%>%data.frame()%>%dplyr::select(-c(dist_sab,id,site_id,geometry)),"output/Acoustic/lcp_distances.csv",row.names=FALSE)


## with our tags --- 

#read in the animal tag metadata from OTN
otn_meta <- read.csv("data/Acoustic/SABMPAtagging_release_locs_5March2024.csv")%>% #note this gives the same outputs as the 'animals.csv' from OTN
  rename(lon_dep=deploy_long,lat_dep=deploy_lat,date_dep=datelastmodified)%>%
  dplyr::select(animal_id,date_dep,lon_dep,lat_dep)%>%
  distinct(animal_id,.keep_all=TRUE)

#read in where the tags have been detected
det_sheets <- excel_sheets("data/Acoustic/det_evts_match_bind_old.xlsm")

det_df <- NULL
for(i in det_sheets){det_df <- rbind(det_df,read_excel("data/Acoustic/det_evts_match_bind_old.xlsm",sheet=i))}

#merge in the deployment information
det_df <- det_df%>%
          left_join(otn_meta)%>%
          rename(lon=mean_longitude,
                 lat=mean_latitude)%>%
          mutate(event_id = 1:n())%>%
          filter(res_time_sec>0,
                 num_detections>2)

det_boundbox <- det_df%>%
                dplyr::select(lon,lat)%>%
                st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                rbind(.,det_df%>%
                        dplyr::select(lon_dep,lat_dep)%>%
                        st_as_sf(coords=c("lon_dep","lat_dep"),crs=latlong))%>%
                st_transform(utm)%>%
                st_buffer(50)%>% # buffer around the points
                st_transform(latlong)%>%
                st_bbox()

#generate a fine(r) scale transition object for the tag detections around St. Anns Bank
tag_trans <- trans_gen(basemap = basemap,
                       bound_box = det_boundbox,
                       resolution = 2, #2 km resolution. For this smaller map, transition object takes a while to be created but it will give more accurate assessments
                       transition_name = "tag_df_sab",
                       dirs=16,
                       recalculate = TRUE)



#check for errors - this was an issue pre-screening for res_time_sec and num_detections
# det_issues <- det_df%>%
#               mutate(diff_lat = lat-lat_dep, #this means the tag deployment location matches also where it was detected
#                      diff_lon=lon-lon_dep)%>%
#               filter(diff_lat == 0 | diff_lon == 0)
# 
# det_df <- det_df%>%
#           filter(!event_id %in% det_issues$event_id)

#load the transition object. Note that this is made during another analysis for the qualified detections using lcp_site.R function. These can be
#manually made if needed.

load("data/transition objects/tag_df_16_10km.RData") #large scale 10km resolution
load('data/transition objects/')

basemap_d <- basemap%>%st_intersection(det_boundbox%>%st_as_sfc()) #used for the straight line check in lcp_otn_pair

lcp_tags <- lcp_otn_pair(x=det_df,trans_d=tag_trans,lines=TRUE,basemap_d = basemap_d,sld_check=TRUE) #takes a while to get through the tags

#unpack the outputs

#lines
spec_ord <- c("Atlantic Halibut","Atlantic Cod","Atlantic Striped Wolffish")

lcp_tags_lines_df <- lcp_tags[[2]]%>%
                     st_set_crs(latlong)%>%
                     left_join(det_df%>%distinct(animal_id,.keep_all=TRUE)%>%dplyr::select(animal_id,species))%>%
                     mutate(species = factor(species,levels = spec_ord))

tag_lims <- lcp_tags_lines_df%>%
            st_transform(utm)%>%
            st_bbox()%>%
            st_as_sfc()%>%
            st_buffer(20)%>%
            st_transform(latlong)%>%
            st_bbox()

spec_ord <- c("Atlantic Halibut","Atlantic Cod","Atlantic Striped Wolffish")

release_locs <- det_df%>%
                st_as_sf(coords=c("lon_dep","lat_dep"),crs=latlong)%>%
                mutate(species = factor(species,levels = spec_ord))

detection_locs <- det_df%>%
                  st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                  mutate(species = factor(species,levels = spec_ord))


lcp_plot <- ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=sab,fill="cornflowerblue")+
  geom_sf(data=lcp_tags_lines_df,linewidth=1,alpha=0.3,col="grey")+
  geom_sf(data=release_locs,aes(fill=species),size=2,shape=21)+
  geom_sf(data=detection_locs,aes(fill=species),size=2,shape=21)+
  facet_wrap(~species)+
  theme_bw()+
  coord_sf(xlim=tag_lims[c(1,3)],ylim=tag_lims[c(2,4)],expand=0)+
  theme(legend.position="none",
        strip.background = element_rect(fill="white"))+
  scale_fill_viridis(discrete=T)+
  annotation_scale(
    data=data.frame(species="Atlantic Halibut"),
    location="br"
  )

ggsave("output/Acoustic/lcp_analysis_species.png",lcp_plot,width=10,height=7,units="in",dpi=300)
trim_img_ws("output/Acoustic/lcp_analysis_species.png")


#save the distances
sab_tags_dists <- lcp_tags[[1]]

write.csv(lcp_tags[[1]],"output/Acoustic/lcp_sab_tags_distances.csv",row.names=FALSE)



