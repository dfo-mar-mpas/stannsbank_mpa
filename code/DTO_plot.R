# Code to support Lazin et al. DTO Tech rep and paper 

#load libraries
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(MarConsNetData)
library(patchwork)

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#source functions
source("https://raw.githubusercontent.com/dfo-mar-mpas/MCRG_functions/refs/heads/main/code/trim_img_ws.R")

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

#Maritimes draft network 
mar_net <- data_draft_areas()%>% # note zones are missing
  st_transform(latlong)

#canada MPAs
can_mcn <- read_sf("data/Shapefiles/All_GoC_MCAs.shp")%>%
           mutate(name = case_when(grepl("Strait of Georgia and Howe Sound Glass Sponge Reef Closure",NAME_E) ~ "Strait of Georgia and Howe Sound Glass Sponge Reef Closure",
                          grepl("Eastport",NAME_E) ~ "Easport MPA",
                          TRUE ~ NAME_E))%>%
           st_transform(latlong)

MPAs <- can_mcn%>%filter(grepl("Marine Protected Area",TYPE_E))

#basemap 

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

coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")%>%st_transform(latlong)

#250 m contour for the Maritimes region 
cont_250_plotregion <- read_sf("data/Bathymetry/contour_250_plotregion.shp")
cont_250_sab <- read_sf("data/Bathymetry/contour_250_sab.shp")

#Canadian EEZ
can_eez <- read_sf("data/Shapefiles/can_eez.shp")%>%st_transform(latlong)

#bathymetry for the St. Anns bank MPA
sab_dem_hr_trimmed <- terra::rast("data/Bathymetry/TotalSABbathy/sab_dem_hr_trimmed.tif")%>%terra::project(latlong)


#make a map of the regional network
mar_net_df <- mar_net%>%
  mutate(NAME = case_when(grepl("anns",tolower(SiteName_E)) ~ "St. Anns Bank Marine Protected Area", #to match the FGP pull
                          SiteName_E == "The Gully Marine Protected Area" ~ "Gully Marine Protected Area",
                          TRUE ~ SiteName_E),
         TYPE = case_when(grepl("refuge",tolower(SiteName_E))~ "MR",
                          grepl("Tier",Classification_E)~ "Draft",
                          grepl("Interest",Classification_E) ~ "AOI",
                          Classification_E == "Existing site" & LeadAgency_E != "Fisheries and Oceans Canada" ~ "OECM",
                          TRUE ~ "MPA"))%>%
  filter(!NAME %in% MPAs$NAME_E) #these are plotted as part of the GIS pull, which has zones. 

plot_polys <- MPAs%>%
  mutate(TYPE="MPA")%>%
  dplyr::select(TYPE)%>%
  rbind(.,mar_net_df%>%
          dplyr::select(TYPE)%>%
          rename(geometry=2))%>%
  mutate(TYPE = factor(TYPE,levels=c("MPA","MR","AOI","OECM","Draft")))

#make the plots
primary_plot <- ggplot()+
                geom_spatraster(data = sab_dem_hr_trimmed)+
                geom_sf(data=cont_250_sab,col="grey35",linewidth=0.3)+
                geom_sf(data=MPAs%>%filter(NAME_E == "St. Anns Bank Marine Protected Area"),fill=NA,linewidth=1.1,col="black")+
                scale_fill_viridis(na.value = "transparent")+
                geom_sf(data=coast_hr,fill="#8AB58A")+
                labs(fill="Depth (m)")+
                # annotation_scale(location="br",
                #                  width_hint = 0.7 * 0.25,  
                #                  text_cex = 0.7)+
                theme_bw()+
                theme(legend.position = "inside",
                      legend.position.inside = c(0.12,0.85),
                      legend.background = element_blank(),
                      legend.text = element_text(size=6),
                      legend.title = element_text(size = 7),       
                      legend.key.size = unit(0.3, "cm"),
                      axis.text=element_blank())+
                guides(shape = guide_legend(override.aes = list(size=6)))+
                labs(x="",y="")+
                coord_sf(expand=0,xlim=plot_boundaries[c(1,3)],ylim=plot_boundaries[c(2,4)])


p1 <- ggplot() +
  geom_sf(data = can_eez, fill = NA) +
  geom_sf(data = cont_250_plotregion, color = "grey80", linewidth = 0.3) +
  
  geom_sf(
    data = plot_polys%>%filter(TYPE!="Draft"),
    aes(
      fill = TYPE,
    ),
    color = "black",
    alpha = 0.3
  ) +
  geom_sf(
    data = plot_polys%>%filter(TYPE=="Draft"),
    aes(
      fill = TYPE,
    ),
    linetype=2, #just emphasize 'draftiness'
    color = "black",
    alpha = 0.3
  ) +
  
  # Basemap for countries
  geom_sf(data = basemap_inset %>% filter(country == "Canada"), fill = "grey70") +
  geom_sf(data = basemap_inset %>% filter(country != "Canada"), fill = "grey90") +
  # Plot boundaries
  geom_sf(data = plot_boundaries %>% st_as_sfc(), fill = NA) + #For St. Anns Bank emphasis
  # Coordinate limits
  coord_sf(expand = 0,
           xlim = plot_lim[c(1, 3)],
           ylim = plot_lim[c(2, 4)]) +
  
  # Theme
  theme_bw() +
  labs(x = "", y = "") +
  
  # Optional: manually set scales for TYPE
  scale_fill_manual(
    values = c(
      "Draft" = "grey70",
      "MR"    = "grey10",
      "AOI"   = "salmon2",
      "MPA"   = "cornflowerblue",
      "OECM" = "springgreen2"
    )
  ) +
  annotation_scale(location="br")+
                   #width_hint = 0.7 * 0.25,  
                   #text_cex = 0.5)+
  labs(fill="")+
  theme(legend.title=element_blank())


combo_plot <- primary_plot + p1 + plot_layout(ncol=2)

ggsave("output/DTO_network_plot.png",combo_plot,width=10,height=10,units="in",dpi=300)
trim_img_ws("output/DTO_network_plot.png")


#Summaries of the Canadian MCN ----

can_mcn_union <- can_mcn %>% #get rid of the zones
  group_by(name) %>%
  st_make_valid()%>%
  summarise(.groups = "drop")%>%
  mutate(area_km=st_area(.)/1e6)%>%
  left_join(.,can_mcn%>%
              data.frame()%>%
              distinct(name,.keep_all=TRUE)%>%
              dplyr::select(name,MGMT_E))

median(can_mcn_union$area_km)
mean(can_mcn_union$area_km)

can_mcn_union%>%
  filter(#MGMT_E == "Fisheries And Oceans Canada",
    name != "Tuvaijuittuq Marine Protected Area")%>% #this one is so skewing
  arrange(-area_km)%>%
  pull(area_km)%>%
  mean()
