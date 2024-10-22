### cluster analysis 

#load libraries -----
library(vegan)
library(tidyverse)
library(sf)
library(arcpullr)
library(stars)
library(marmap)
library(patchwork)
library(gridExtra)
library(cowplot)
library(grid)

nmdstheme <- theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(size=18))

# get SAB shapefile
url <- "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Oceans_Act_Marine_Protected_Areas/MapServer/0"
sab <- get_spatial_layer(url,where="NAME_E='St. Anns Bank Marine Protected Area'")

#get the RV survey strata
rv_strata <- get_spatial_layer("https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Maritimes_Summer_Research_Vessel_Survey_En/MapServer/1")
rv_strata <- rv_strata%>%st_transform(4326)

rv_sab <- rv_strata%>%st_intersection(sab%>%st_union())

ggplot()+
  geom_sf(data=rv_sab,aes(fill=STRAT))+
  geom_sf(data=sab,col="black",lwd=1.01,fill=NA)+
  geom_sf(data=crabstations)+
  theme_bw()

# get crab data
crab_swept <- read.csv("data/CrabSurvey/sabmpa_area_swept.csv")

spid <- read.csv("data/CrabSurvey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv") %>%
  rowwise() %>%
  mutate(sp_name = if_else(SPEC=="",
                           if_else(COMM=="",as.character(SPECCD_ID),COMM),
                           SPEC))

crab <- read.csv("data/CrabSurvey/SABMPA2023export.csv") %>%
  left_join(crab_swept,by = join_by(TRIP, "SET_NO"=="SET", STATION)) %>%
  mutate(CPUE = EST_NUM_CAUGHT/AREA_SWEPT,
         station_date=paste(STATION,BOARD_DATE,sep = ": "),
         LONGITUDE=-LONGITUDE) %>%
  left_join(spid,by = join_by(SPECCD_ID)) %>%
  mutate(sp_name=if_else(is.na(sp_name),
                         as.character(SPECCD_ID),
                         sp_name))

# get bathymetry data for SAB
bbsab <- st_bbox(sab)
#here removed getnoaa bathy and just put in the csv file 
noaabathy <- as.bathy(read.csv("marmap_coord_-60.6499968895566;45.4167096823514;-57.3666964444786;46.7833095612108_res_0.25.csv")) %>%
  fortify.bathy() %>%
  st_as_stars() %>%
  st_set_crs(4326)

# get benthoscape
benthoscape <- st_read("data/Shapefiles/benthoscape.shp")%>%
  st_make_valid() %>%
  st_transform(4326)

# set up crab_stations and join with benthoscape and bathy

crabstations <- crab %>%
  group_by(STATION) %>%
  reframe(LONGITUDE=mean(LONGITUDE),
          LATITUDE=mean(LATITUDE)) %>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"),
           crs = 4326) %>%
  st_join(noaabathy %>% st_as_sf()) %>%
  st_join(benthoscape%>%dplyr::select(Class_name,geometry)) %>%
  st_join(rv_strata%>%dplyr::select(STRAT,geoms)%>%rename(geometry=geoms))%>%
  mutate(Class_name=if_else(is.na(Class_name),
                            "Unclassified",
                            Class_name),
         strat_name=ifelse(is.na(STRAT),"outside rv",STRAT)) %>%
  dplyr::select(STATION,z,Class_name,strat_name) %>%
  mutate(STATION = as.character(STATION),
         z_class = if_else(z>=(-100),
                           "<100 m",
                           if_else(z>=(-150),
                                   "100 to 150 m",
                                   ">150 m")))

focal_species <- c("HIPPOGLOSSOIDES PLATESSOIDES",
                   "GLYPTOCEPHALUS CYNOGLOSSUS",
                   "SEBASTES",
                   "AMBLYRAJA RADIATA",
                   "GADUS MORHUA",
                   "CHIONOECETES OPILIO",
                   "ANARHICHAS LUPUS")


crab_df <- crab%>%
           mutate(STATION = as.character(STATION),# to match the crab stations
                  year=year(as.POSIXct(BOARD_DATE)))%>% 
           left_join(.,data.frame(crabstations))%>%
           data.frame()%>%
           filter(SPEC %in% focal_species)%>%
           rename(station = STATION,species=SPEC,abund=CPUE,bioclass_name = Class_name)%>%
           dplyr::select(station,year,species,abund,bioclass_name,strat_name,z_class)

crab_df_stratified <- crab_df%>%
                      group_by(station,year,z_class,bioclass_name,strat_name)%>%
                      pivot_wider(names_from = species,values_from = abund,values_fill = 0)%>%
                      mutate_at(focal_species,~replace_na(.,0))%>%
                      ungroup()%>%
                      mutate(sums = rowSums(dplyr::select(.,all_of(focal_species))))%>%
                      filter(sums>0)

species_abundance <- crab_df_stratified %>%
                     select(-station, -year, -bioclass_name, -strat_name, - z_class)
                         

perm_bioclassification <- adonis2(species_abundance ~ bioclass_name*year, data = crab_df_stratified,  method = "bray", permutations = 999)
perm_rvstrata <- adonis2(species_abundance ~ strat_name*year, data = crab_df_stratified, method = "bray", permutations = 999)
perm_depth_strata <- adonis2(species_abundance ~ z_class*year, data = crab_df_stratified, method = "bray", permutations = 999)


#compare explained variance
data.frame(method=c("Bioclassification strata","RV survey strata","Depth-based zonation"),
           r2 = c(perm_bioclassification$R2[1],perm_rvstrata$R2[1],perm_depth_strata$R2[1]))%>%
arrange(-r2)

print(perm_rvstrata)  
print(perm_bioclassification)  
print(perm_depth_strata)

#clustering

crab_df <- crab%>%
  mutate(STATION = as.character(STATION),# to match the crab stations
         year=year(as.POSIXct(BOARD_DATE)))%>% 
  left_join(.,data.frame(crabstations))%>%
  data.frame()%>%
  rename(station = STATION,species=SPEC,abund=CPUE,bioclass_name = Class_name)%>%
  dplyr::select(station,year,species,abund,bioclass_name,strat_name,z_class)

crab_df_stratified <- crab_df%>%
  filter(species != "")%>%
  group_by(station,year,z_class,bioclass_name,strat_name)%>%
  pivot_wider(names_from = species,values_from = abund,values_fill = 0)%>%
  mutate_at(crab_df%>%filter(species != "")%>%pull(species),~replace_na(.,0))%>%
  ungroup()%>%
  mutate(sums = rowSums(dplyr::select(.,all_of(focal_species))))%>%
  filter(sums>0)

species_abundance <- crab_df_stratified %>%
  select(-station, -year, -bioclass_name, -strat_name, - z_class)



#data for the clustering
nmds_df <- species_abundance %>% 
           mutate(site=1:nrow(crab_df_stratified))%>%
           column_to_rownames('site')%>%
           mutate(tot = rowSums(.,na.rm=T))%>%
            filter(tot>0)%>% #filter stations with 0 assignments -- WRK2.2 didn't have any
            select(-tot)%>%
            decostand(method = "total")

pa_data <- nmds_df%>%
  mutate_if(is.numeric, ~1 * (. > 0))%>%
  vegdist(., method = "jaccard", binary = TRUE)

nmds_pa <- pa_data%>%metaMDS(.,k=5)
nmds_bray <- metaMDS(nmds_df,k=5,distance="bray")

data.scores <- as.data.frame(scores(nmds_bray,"sites"))%>%
  mutate(station=rownames(.),
         method="bray",
         stress=round(nmds_bray$stress,3))%>%
  cbind(.,crab_df_stratified%>%dplyr::select(bioclass_name,strat_name,z_class))%>%
  rbind(.,
        as.data.frame(scores(nmds_pa,"sites"))%>%
          mutate(station=rownames(.),
                 method="pa",
                 stress=round(nmds_pa$stress,3))%>%
          cbind(.,crab_df_stratified%>%dplyr::select(bioclass_name,strat_name,z_class)))%>%
  select(NMDS1,NMDS2,station,bioclass_name,z_class,strat_name,method,stress)


#calculate contex hull
hull_data <- data.scores%>%
  gather("group","value",bioclass_name:strat_name)%>%
  group_by(group,value,method)%>%
  slice(chull(x=NMDS1,y=NMDS2))

data_scores_df <- data.scores%>%gather("group","value",bioclass_name:strat_name)


#extract the stress for each  
stresses <- hull_data%>%
  group_by(group,method)%>%
  summarise(stress=round(unique(stress),3))%>%
  ungroup()

#assemble plots

temp_list = list()

for(i in c("bray","pa")){
  
  for(j in c("bioclass_name","strat_name","z_class")){
    
    plot_df = hull_data%>%
      filter(method==i,group==j)%>%
      mutate(lab = ifelse(method=="bray","Bray-curtis","Jaccard"))
    
    p1 = ggplot()+
      geom_polygon(data=plot_df%>%filter(value !="Unclassified"),aes(x=NMDS1,y=NMDS2,fill=value),alpha=0.30,col="grey30",lwd=0.2)+
      geom_point(data=plot_df%>%filter(value !="Unclassified"),aes(x=NMDS1,y=NMDS2,shape=value,fill=value),col="black",size=3)+
      geom_point(data=plot_df%>%filter(value =="Unclassified"),aes(x=NMDS1,y=NMDS2),shape=19,size=0.25)+
      scale_shape_manual(values=c(21:25))+
      nmdstheme+
      facet_wrap(~lab,ncol=3)+
      theme(strip.background = element_rect(fill="white"))+
      annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5, 
               label = paste("Stress =", stresses%>%filter(method=="bray")%>%pull(stress), "k = 5")) 
    
    if(i == "bray" & j=="bioclass_name") {p1 = p1 + 
                                               theme(
                                                     axis.title = element_blank())+
                                                labs(fill="",shape="",title = "Bioclassifcation")}
    
    if(i == "pa" &  j == "bioclass_name"){p1 = p1 + 
                                               theme(
                                                     axis.title = element_blank())+
                                                labs(fill="",shape="")}
    
    if(i == "bray" & j=="strat_name") {p1 = p1 + 
                                              theme(strip.background = element_blank(),
                                                    strip.text = element_blank(),
                                                    axis.title.x = element_blank())+
                                              labs(title = "RV strata",fill="",shape="")}
    
    if(i == "pa" &  j == "strat_name"){p1 = p1 + 
                                              theme(strip.background = element_blank(),
                                                    strip.text = element_blank(),
                                                    axis.title = element_blank())+
                                              labs(fill="",shape="")}
    
    if(i == "bray" & j=="z_class") {p1 = p1 + 
                                          theme(strip.background = element_blank(),
                                                strip.text = element_blank(),
                                                axis.title.y = element_blank())+
                                          labs(title = "Depth-based stratification",fill="",shape="")}
    
    if(i == "pa" &  j == "z_class"){p1 = p1 + 
                                            theme(strip.background = element_blank(),
                                                  strip.text = element_blank(),
                                                  axis.title.y = element_blank())+
                                            labs(fill="",shape="")}
                                          
                                            
    temp_list[[paste(i,j,sep="-")]] = p1
    
  } # end of i loop
} # end of j loop


combo_plot = ((((temp_list[[1]] | temp_list[[4]]) + plot_layout(guides="collect")) & theme(legend.position="bottom",legend.text = element_text(size=12))))/
             ((((temp_list[[2]] | temp_list[[5]]) + plot_layout(guides="collect")) & theme(legend.position="bottom",legend.text = element_text(size=12))))/
             ((((temp_list[[3]] | temp_list[[6]]) + plot_layout(guides="collect",axis_titles = "collect")) & theme(legend.position="bottom",legend.text = element_text(size=12))))


ggsave("output/crab_clusters.png",combo_plot,height=8*2,width=5*2,units="in",dpi=300)




#conduct nmds 
nmds_bray <- metaMDS(dist_matrix,k=8,distance="bray")


cluster <- hclust(dist_matrix, method = "ward.D2")

plot(cluster)

## richness comparison
# Calculate species richness as the number of species with non-zero abundance at each station

crab_df_rich <- crab%>%
                mutate(STATION = as.character(STATION),# to match the crab stations
                       year=year(as.POSIXct(BOARD_DATE)))%>% 
                data.frame()%>%
                filter(CPUE>0)%>%
                group_by(year,STATION)%>%
                summarize(rich = length(unique(SPEC)))%>%
                ungroup()%>%
                left_join(.,data.frame(crabstations))


