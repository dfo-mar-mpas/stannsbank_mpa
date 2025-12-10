#Looking at the high resolution raster

#load libraries
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggspatial)
library(patchwork)

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

sab_highres <- rast("data/Bathymetry/SAB_highres/fullBathy.tif")

bathy_proj <- st_crs(sab_highres)

#load the St. Anns Bank MPA shapefiles
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(bathy_proj)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(bathy_proj)%>%
  st_as_sf()

sab_banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%
  st_transform(bathy_proj)

#crop the extent of the raster
sab_bound <- st_bbox(sab)%>%
             st_as_sfc()



#large scale bathy cropped to SAB -----
bathy_cropped <- sab_highres %>%
  crop(., sab_bound) %>%
  mask(., sab) %>%
  rename(depth=1)


#Banks zoomed in -----
curdo_bound <- sab_banks%>%
               filter(name=="Curdo Bank")%>%
               st_transform(utm)%>%
               st_buffer(0.5)%>%
               st_transform(bathy_proj)%>%
               st_bbox()

curdo_bathy <- bathy_cropped%>%
               crop(.,curdo_bound%>%st_as_sfc())

curdo_40m_contour <- as.contour(curdo_bathy, levels = -40)%>%
                     st_as_sf()%>%
                     st_set_crs(bathy_proj)%>%
                     mutate(name="Curdo Bank")

  
#scatarie bank
scatarie_bound <- sab_banks%>%
                  filter(name=="Scatarie Bank")%>%
                  st_transform(utm)%>%
                  st_buffer(0.75)%>%
                  st_transform(bathy_proj)%>%
                  st_bbox()

scatarie_bathy <- bathy_cropped%>%
  crop(.,scatarie_bound%>%st_as_sfc())

scatarie_40m_contour <- as.contour(scatarie_bathy , levels = -40)%>%
                        st_as_sf()%>%
                        st_set_crs(bathy_proj)%>%
                        mutate(name="Scatarie Bank")

#combine the banks
sab_banks_40m <- rbind(curdo_40m_contour,scatarie_40m_contour)

#write_sf(sab_banks_40m,"data/Shapefiles/sab_banks_40m.shp")

#make the plots 

depths <- values(bathy_cropped)%>%
          as.data.frame()%>%
          filter(!is.nan(depth))


depth_limits <- c(min(depths), 
                  max(depths))

curdo_plot <- ggplot() +
              geom_spatraster(data = curdo_bathy, aes(fill = depth)) +
              geom_sf(data=curdo_40m_contour,fill=NA)+
              scale_fill_viridis_c(name = "Depth (m)",na.value = NA) +
              annotation_scale(location="tl")+
              theme_bw()+
              coord_sf(expand=0)+
              theme(axis.text=element_blank())+
              labs(title="Curdo Bank")

scatarie_plot <- ggplot() +
                  geom_spatraster(data = scatarie_bathy, aes(fill = depth)) +
                  geom_sf(data=scatarie_40m_contour,fill=NA)+
                  scale_fill_viridis_c(name = "Depth (m)",na.value = NA) +
                  # scale_fill_viridis_c(name = "Depth (m)", na.value = NA,
                  #                      limits = depth_limits)+
                  annotation_scale(location="tl")+
                  theme_bw()+
                  coord_sf(expand=0)+
                  theme(axis.text=element_blank())+
                  labs(title="Scatarie Bank")
    #buffer for the mosaic plotting             
sab_plot_bound <- sab%>%
  st_transform(utm)%>%
  st_buffer(5)%>%
  st_transform(bathy_proj)%>%
  st_bbox()


sab_plot <- ggplot() +
            geom_spatraster(data = bathy_cropped, aes(fill = depth)) +
            geom_sf(data=sab_banks_40m,fill=NA,linewidth=0.25)+
            geom_sf(data=curdo_bound%>%st_as_sfc(),fill=NA,linewidth=0.2)+
            geom_sf(data=scatarie_bound%>%st_as_sfc(),fill=NA,linewidth=0.2)+
            geom_sf(data = sab_zones, fill = NA, color = "black", linewidth = 1) +
            scale_fill_viridis_c(name = "Depth (m)",na.value = NA) +
            annotation_scale(location="tl")+
            labs(title="St Anns Bank MPA - multibeam bathymetry coverage")+
            theme_bw()+
            theme(legend.position = "inside",
                  legend.position.inside = c(0.93,0.25),
                  legend.background = element_blank(),
                  legend.text = element_text(size=8),
                  legend.title = element_text(size=9))+
            coord_sf(expand=0,xlim=sab_plot_bound[c(1,3)],ylim=sab_plot_bound[c(2,4)])

mosaic_map <- sab_plot + (scatarie_plot/curdo_plot) + plot_layout(ncol=2,widths=c(3.5,1.4))

ggsave("output/multibeam_bathymap_stannsbank.png",mosaic_map,width=8,height=4,units="in",dpi=600)
