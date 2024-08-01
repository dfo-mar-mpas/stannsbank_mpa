#Curdo Correction - the shapefile from NRCAN has it offset from the reported centriod and where we observed the bank in 2023

library(tidyr)
library(tidyverse)
library(sf)

#original banks
banks <- read_sf("data/Shapefiles/sab_banks_old.shp")%>%
          st_transform(latlong)

#get the offset Curdo
polygon <- banks%>%filter(name=="Curdo Bank")

#Where the bank should be centered
target_centroid <- data.frame(lat=45.893003, lon=-59.598049)
target_centroid_sf <- target_centroid%>%st_as_sf(coords=c("lon","lat"),crs=latlong)

#Where it the offset wrong location is
current_centroid <- st_centroid(polygon)
current_centroid_coords <- st_coordinates(current_centroid)

# Compute the translation vector
translation_vector <- data.frame(lon = target_centroid$lon - current_centroid_coords[1, "X"],
                                  lat = target_centroid$lat - current_centroid_coords[1, "Y"])

#fix it. 
curdo_new <- polygon%>%
                st_coordinates(.)%>%
                data.frame()%>%
                rename(lon=1,lat=2)%>%
                dplyr::select(lon,lat)%>%
                mutate(lon=lon + translation_vector$lon,
                       lat=lat + translation_vector$lat)%>%
                st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
                st_combine()%>%
                st_cast("POLYGON")%>%
                st_as_sf()%>%
                rename(geometry=1)%>%
                mutate(name="Curdo Bank")%>%
                dplyr::select(name,geometry)

banks_new <- rbind(curdo_new,banks%>%filter(name=="Scatarie Bank"))

st_write(banks_new, "data/Shapefiles/sab_banks.shp")

#show it worked

ggplot()+
  geom_sf(data=polygon,fill="cornflowerblue")+
  geom_sf(data=curdo_new,fill="coral2")+
  geom_sf(data=target_centroid_sf)
