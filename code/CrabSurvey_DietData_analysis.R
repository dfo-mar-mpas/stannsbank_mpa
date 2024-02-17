#Diet data exploration

##libraries
library(sp)
library(rgdal)
library(raster)
library(shape)
library(rgeos)
library(raster)
library(ggplot2)
library(mapdata)
library(marmap)
library(maptools)
library(dplyr)

#Read in diet data pulled from ISDB - note at the time of coding (Feb 2024) the 2023 data was not processed
diet <- read.csv("data/CrabSurvey/MPA.Diet.SnowCrabSurvey.Feb.2023.csv")
head(diet)

plot(diet$SLONGDD, diet$SLATDD)

sabdiet <- diet %>% filter(SLATDD > 45.1) %>% mutate(Year=as.Date(as.numeric(SDATE), tryFormats = c("%d-%m-%Y", "%d/%m/%Y")),"%Y")
dim(sabdiet)

## SAB RV survey diet database summary code
# Harri Pettitt-Wade
# harri.pettitt-wade@dfo-mpo.gc.ca





#Read in the St Anns Polygon
SAB <- readOGR("/SAB Gut contents/StAnnBank/StAnnsBank_MPA.shp")

diets_Oct2022 <- read.csv("SAB Gut contents/FH.DataRequest.OCT.2022.csv")

str(diets_Oct2022)

diets_splist <- read.csv("SAB Gut contents/FH.Species.List.OCT2022.csv")

diets$PRED_SEQ

diets_splist$SPECCD

#diets_ <- dets %>%
#  filter(station == "CAMP")

diets_splist_prey <- diets_splist
diets_splist_pred <- diets_splist

data.table::setnames(diets_splist_prey,  c("SPECCD"), 
                     c("PREYSPECCD")) 

diets_FULL <- merge(diets, diets_splist_prey, by.x = diets$PREYSPECCD, by.y = diets_splist$PREYSPECCD, all.x = TRUE, allow.cartesian = TRUE)

#dat
diets_FULL <- merge(diets,diets_splist_prey,by.x="PREYSPECCD",by.y="PREYSPECCD",all.x=F, allow.cartesian = TRUE)

diets_full_preds <- merge(diets,diets_splist_pred,by.x="SPEC",by.y="SPEC",all.x=F, allow.cartesian = TRUE)

data.table::setnames(diets_splist_pred,c("PREYSPECCD","COMMON","SPECIES","PHYLUM","CLASS","ORD","FAMILY","GENUS"),
                     c("SPEC","COMMON_pred","SPECIES_pred","PHYLUM_pred","CLASS_pred","ORD_pred","FAMILY_pred","GENUS_pred"))

data.table::setnames(diets_splist_pred,c("PREYSPECCD"),
                     c("SPEC"))

diets_FULL_predprey <- merge(diets_full_PREYSPECCD, diets_splist_pred, by.x = diets_full_PREYSPECCD$SPEC, by.y = diets_splist_pred$SPEC, all.x = TRUE, allow.cartesian = TRUE)

#Read in the diet data
data <- read.csv("RVdietdata4VN.csv")
data$slong <- data$slong*-1 # convert to decimal degrees if needed

diets_FULL <- merge(diets,diets_splist,by="ID")


ggplot(data = diets_full_preds, mapping = aes(x = diets_full_preds$DATASOURCE, y = diets_full_preds$RANK)) +
  geom_boxplot(mapping = aes(color=diets_full_preds$FAMILY))

library(stringr)
diets_full_preds$total_prey <- str_count(df1$Name,"a")
df1

my_data <- as_tibble(diets_full_preds)
my_data

summary <- my_data %>%
  group_by(my_data$ORD_pred) %>%
  summarise(
    count = n())

my_data_dt <- as.data.frame(my_data)

ggplot(data = my_data_dt, mapping = aes(x = ORD_pred, y = summarise(count=n())) +
         geom_histogram(mapping = aes(color=ORD_pred)
         )
       
       p <- ggplot(my_data_dt, aes(ORD_pred))
       p + theme_classic() + 
         geom_bar(mapping = aes(color=my_data_dt$COMMON_pred))
       
       total_common = n(my_data$COMMON_pred),
       total_genus = n(my_data$GENUS_pred),
       total_class = n(my_data$CLASS_pred),
       total_phylum = n(my_data$PHYLUM_pred),
       
       
       my_data %>%
         group_by(ORD_pred) %>%
         select() %>%
         summarise(
           mean.dispersal = mean(Release.Dispersal, na.rm = TRUE),
           max.dispersal = max(Release.Dispersal, na.rm = TRUE)
         )
       
       #take a subset of data and convert to a spatialPoints object 
       xy <- diets[,c("SLONGDD","SLATDD")]
       coordinates(xy) <- ~SLONGDD+SLATDD
       proj4string(xy) <- proj4string(SAB)
       
       sum(is.na(diets$SLATDD))
       sum(!is.na(diets$SLONGDD))
       
       #Projections ------------
       latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
       # <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
       
       
       # Make raster grid ----------------
       
       #grid <- xyz%>%
       #  st_as_sf(coords=c("Longitude","Latitude"),crs=latlong)%>%
       #  st_transform(planar)%>% #make in planar coordinates so you can set the resolution
       
       diets_NA <- subset(diets,is.na(diets$SLATDD)|is.na(diets$SLONGDD)) #subset NAs, have a look to see if I can figure why they don't have lat lon first
       diets_noNA <- subset(diets,!is.na(diets$SLATDD)|!is.na(diets$SLONGDD)) #subset the rest with lat lons
       
       diets_original <- diets #change names to keep it clear and backup
       diets <- diets_noNA #backup and simplify naming for quick coding
       
       #Create a projection of the polygon in UTM
       SAB_UTM <- spTransform(SAB,"+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")
       
       #Create buffered polygons at 5, 10 and 20 km from the boundaries and reproject back to decimal degrees, convert to SPatialPolygonsDataframe
       BufferFun <- function(x,buffer){
         buffpoly <- gBuffer(x,width=buffer)
         buffpoly <- spTransform(buffpoly,"+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
         p.df <- data.frame(ID=1) #need this to convert to SpatialPolygonsDataFrame
         buffpoly <- SpatialPolygonsDataFrame(buffpoly,p.df,match.ID = F)
         return(buffpoly)
       }
       
       
       SAB_buff_5 <- BufferFun(SAB_UTM,buffer=5)
       SAB_buff_15 <- BufferFun(SAB_UTM,buffer=15)
       SAB_buff_30 <- BufferFun(SAB_UTM,buffer=30)
       
       #Identify coordinates within the SAB with a logical vector
       data$inout <- complete.cases(over(xy,geometry(SAB)))
       data$inout5 <- complete.cases(over(xy,geometry(SAB_buff_5)))
       data$inout15 <- complete.cases(over(xy,geometry(SAB_buff_15)))
       data$inout30 <- complete.cases(over(xy,geometry(SAB_buff_30)))
       
       #Set plotting limits with a one degree buffer
       rextent <- extent(xy)
       Long.lim  <-  c(rextent[1]-0.5, rextent[2]+0.5)
       Lat.lim <-  c(rextent[3]-0.5, rextent[4]+0.5)
       
       #download bathymetry 
       bathy <- getNOAA.bathy(Long.lim[1],Long.lim[2],Lat.lim[1],Lat.lim[2],res=1,keep=T)
       
       #Get the mapping layers
       canada <- map_data("worldHires", "Canada")
       
       #Take the polygon objects and convert to tabular long-format for ggplot using the fortify function add a column with the buffer distance
       fortbathy <- fortify(bathy)
       fortSAB <- fortify(SAB)%>%mutate(buffer=0)
       fortSAB5 <- fortify(SAB_buff_5)%>%mutate(buffer=5)
       fortSAB15 <- fortify(SAB_buff_15)%>%mutate(buffer=15)
       fortSAB30 <- fortify(SAB_buff_30)%>%mutate(buffer=30)
       
       PolygonDB <- rbind(fortSAB,fortSAB5,fortSAB15,fortSAB30)
       PolygonDB$VAL <- paste(PolygonDB$group,PolygonDB$buffer,sep="-")
       PolygonDB$GROUP <- PolygonDB$VAL
       PolygonDB[grep("-0",PolygonDB$VAL),"GROUP"] <- "MPA"
       PolygonDB[grep("-5",PolygonDB$VAL),"GROUP"] <- "5 km"
       PolygonDB[grep("-15",PolygonDB$VAL),"GROUP"] <- "15 km"
       PolygonDB[grep("-30",PolygonDB$VAL),"GROUP"] <- "30 km"
       
       PolygonDB$GROUP <- factor(PolygonDB$GROUP,levels=c("MPA","5 km","15 km","30 km"))
       
       #labels
       data$lab <- "St. Anns Bank"
       data[!data$inout,"lab"] <- "Outside"
       
       #Plot it with ggplot - view plot
       p1 <- ggplot() +
         geom_polygon(data = fortSAB15, 
                      aes(x=long, y = lat, group = group), 
                      fill = "white", 
                      color="black") +
         geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
                      fill = "white", color="black") + 
         geom_point(data=data,aes(x=slong,y=slat,fill=lab),pch=21)+
         scale_fill_manual(values=c("white","black"))+
         geom_contour(data=bathy,aes(x=x,y=y,z=z),breaks=c(-250),lwd=0.1,colour="grey50")+
         coord_fixed(xlim = Long.lim,  ylim = Lat.lim, ratio = 1.2,expand = 0)+
         theme_bw()+
         theme(legend.position = "bottom",
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_rect(fill = "white", colour = "black"),
               plot.background = element_rect(colour = "white"),
               strip.background = element_rect(colour = "black", fill = "white"))+
         labs(x=expression(paste("Longitude ",degree,"W",sep="")),
              y=expression(paste("Latitude ",degree,"N",sep="")),
              fill="");p1
       
       #save plot
       ggsave("RVDiet_Plot.png",p1)
       
       
       #Plot it with ggplot - with buffers
       p2 <- ggplot() +
         geom_polygon(data = PolygonDB, 
                      aes(x=long, y = lat, group = VAL ,col=GROUP),
                      fill = NA) +
         geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
                      fill = "white", color="black") + 
         geom_point(data=data,aes(x=slong,y=slat),pch=19)+
         scale_fill_manual(values=c("white","black"))+
         geom_contour(data=bathy,aes(x=x,y=y,z=z),breaks=c(-250),lwd=0.1,colour="grey50")+
         coord_fixed(xlim = Long.lim,  ylim = Lat.lim, ratio = 1.2,expand = 0)+
         theme_bw()+
         theme(legend.position = "bottom",
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_rect(fill = "white", colour = "black"),
               plot.background = element_rect(colour = "white"),
               strip.background = element_rect(colour = "black", fill = "white"))+
         labs(x=expression(paste("Longitude ",degree,"W",sep="")),
              y=expression(paste("Latitude ",degree,"N",sep="")),
              fill="");p2
       
       #save plot
       ggsave("RVDiet_buffers.png",p2)
       
       write.csv(release_subset1, "SABMPA dets release subset.csv")
       
       # Sample data of the package
       data <- rd("Employee")
       
       # install.packages(lessR)
       library(lessR)
       # Donut chart
       PieChart(PREY_SEQ, data = diets,
                fill = "viridis",
                color = "black",
                lwd = 2,
                lty = 1,
                rows = (Gender == "F" & Salary > 45000),
                main = NULL)
       
