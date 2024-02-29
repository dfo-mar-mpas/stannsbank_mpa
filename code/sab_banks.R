## code to make Scatarie Bank

#load libraries
library(dplyr)
library(sf)
library(ggplot2)
library(raster)
library(smoothr)

#function for converting raster to polygon
source("code/strat_function.R")

#load sab polygon
sab <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")

#load SAB dem - based on the Crabitat 'predSpace' model
dem_sab <- raster("data/Bathymetry/sab_dem.tif")

#get the banks
bank_df <- data.frame(bank=c("Curdo Bank","Scatarie Bank"),
                      lon=c(-59.598049,-59.213127),
                      lat=c(45.893003,45.988401))%>%
  st_as_sf(coords=c("lon","lat"),crs=latlong)

#range extent of Scattarie 
bank_range <- bank_df%>%
              filter(bank=="Scatarie Bank")%>%
              st_transform(utm)%>%
              st_buffer(10)%>%
              st_transform(latlong)%>%
              st_bbox()

#download curdo from the geojson link for the site - https://geonames.nrcan.gc.ca/search-place-names/unique?id=CAISM
curdo_json <- "https://geogratis.gc.ca/services/geoname/en/geonames/CAISM.geojson?expand=feature,concise,generic,province,status,language,provider&exclude=feature.geometry&_gl=1*zk41wn*_ga*Mjk0MjcwMjQ4LjE3MDkyMTk0Nzk.*_ga_C2N57Y7DX5*MTcwOTIxOTQ3OS4xLjEuMTcwOTIxOTQ5OS4wLjAuMA.."

curdo <- st_read(curdo_json)%>%
         st_transform(latlong)%>%
         mutate(name="Curdo Bank",
                geometry=st_union(geometry))%>%
         dplyr::select(name,geometry)

#scatarie is a point and not a polygon so it will have to be assembled manually

#Crop to the extent of Scatarie
dem_sab2 <- dem_sab%>%
            crop(.,bank_range%>%st_as_sfc()%>%st_as_sf()%>%st_transform(proj4string(dem_sab)))%>%
            projectRaster(.,crs=latlong)

#find the range of the plot
values(dem_sab2)[values(dem_sab2)>50] <- NA #this seems to make the most sense for the bank's shape 
values(dem_sab2)[values(dem_sab2)<10] <- NA
values(dem_sab2) <- values(dem_sab2)*-1

#assemble shape file based on polygon 
scatarie_bank <- stratFun(dem_sab2,min_depth = -10,max_depth = -50)%>%
                 smoothr::smooth(.,method="chaikin")%>%
                 slice(1)%>% #just the big piece
                 mutate(name="Scatarie Bank")%>%
                 dplyr::select(name,geometry)

#bind them together
sab_banks <- rbind(scatarie_bank,curdo)

#save the output for later use
write_sf(sab_banks,"data/Shapefiles/sab_banks.shp")
