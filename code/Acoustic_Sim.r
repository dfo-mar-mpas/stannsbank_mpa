## Load required libraries
library(leaflet)
library(lwgeom)
library(readr)
library(gganimate)
library(viridis)
library(utils)
library(tidyterra)
library(sf)
library(nngeo)
library(tidyverse)
library(terra)
library(animation)
library(sfheaders)
library(rnaturalearth)
library(geosphere)

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#load polygons
sab <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
        st_transform(utm)%>%
        st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
        st_union()%>% #gets rid of the zones
        st_as_sf()

sab_banks <- read_sf("data/Shapefiles/sab_banks.shp")%>%st_transform(utm)

#load coastline and make basemap ------
coast_hr <- read_sf("data/shapefiles/NS_coastline_project_Erase1.shp")

basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
           filter(woe_name !="Nova Scotia")%>% #this is in HR with coast_hr
           st_as_sf()%>%
           st_union()%>%
           st_transform(utm)%>%
           st_as_sf()

plot_lims <- sab%>%
             st_buffer(150)%>%
             st_bbox()

#now get coordinates for a grid
sab_grid <- sab%>%
            st_make_grid(cellsize=10,what="centers")%>%
            st_as_sf()%>%
            st_intersection(sab)

#start location (Shelf -> Gulf) running centered on the MPA

sab_cent <- sab%>%
            st_transform(latlong)%>%
            st_centroid()

fish_bear <- -24.46022 #this is a bearing to go through roughly the centre of SAB from the SS to the Gulf

# start_locs <- data.frame(lon=c(-58.7,-59.7),
#                         lat=c(45.5,47),
#                         dir=c("start","end"))%>%
#                st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
#                st_transform(utm)

start_locs <- destPoint(as_Spatial(sab_cent),fish_bear+180,75*1000)%>%
              data.frame()%>%
              st_as_sf(coords=c("lon","lat"),crs=latlong)%>%
              rbind(.,
                    destPoint(as_Spatial(sab_cent),fish_bear,75*1000)%>%
                    data.frame()%>%
                    st_as_sf(coords=c("lon","lat"),crs=latlong))%>%
               st_transform(utm)%>%
               mutate(dir=c("start","end"))

lin_dr <- start_locs%>%st_combine()%>%st_cast("LINESTRING")

#set particle starting parameters
release_loc <- start_locs%>%
  filter(dir=="start")%>%
  st_buffer(20)%>%#20 km buffer
  st_bbox()

ggplot()+
  geom_sf(data=basemap)+
  geom_sf(data=coast_hr)+
  geom_sf(data=sab,fill=NA)+
  geom_sf(data=sab_grid)+
  geom_sf(data=sab_banks)+
  geom_sf(data=lin_dr,
          arrow = arrow(angle = 45, 
                        ends = "last", 
                        type = "closed", 
                        length = unit(0.25, "cm")))+
  geom_sf(data=start_locs,aes(col=dir),size=2)+
  geom_sf(data=release_loc%>%st_as_sfc(),fill=NA,lwd=0.5,lty=2)+
  theme_bw()+
  labs(col="")+
  coord_sf(expand=0,xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)])

# Number of desired particle to simulate
num_fish <- 1000
bdlng.sec <- 5  # swim speed
ping.rate <- 90

# particle attributes. Don't ask why they're called prey
virtual_fish <- data.frame(id = paste("B.prey", 1:num_fish, sep = ""),
                   lon = runif(num_fish, release_loc[1], release_loc[3]), #Start y location
                   lat = runif(num_fish, release_loc[2], release_loc[4]), #Start x location
                   size = rnorm(num_fish, mean = .15, sd = 0.015),#random size
                   time=1) %>% 
          mutate(speed =  size* #Speed as a function of size (bodylengths/sec)
                   bdlng.sec* # Body lengths per second
                   ping.rate) # Average ping rate) 

virtual_fish_sf <- virtual_fish%>%st_as_sf(coords=c("lon","lat"),crs=utm)

#set up model steps --- 
# current trajectory versus attractant movement ratio
cvar<-.95 #(eg. 94% of movement is based on current trajectory whereas 6% is biased towards point of attraction)

# sd of random direction (units in radians)
prey.noise_sd<-.22

# the general radian direct from seed location towards Gulf of Maine
fish_bear_rad<-fish_bear*pi/180

# Maximum number of steps within sim
num_steps<-6000 ## just a made up number


# Simulation progress bar
pb <- txtProgressBar(min = 0, max = num_steps, style = 3)

#Create empty list to populate for efficiency
sim_data <- rep(list(NULL),num_steps)
system.time(for (step in 1:num_steps) {
  ptm <- Sys.time()
  for (i in 1:num_fish) {
    
    cur_fish<-virtual_fish[i,]
    
    attr_bearing<-(pi*fish_bear_rad)%%(2*pi)
    
    if(step>2){cur_dat<-sim_data[[step-2]] %>% filter(id==cur_fish$id)}else{cur_dat<-data.frame()}
    
    if(nrow(cur_dat)>0){
      prev_dx <- cur_fish$lon -cur_dat$lon  
      prev_dy <- cur_fish$lat -cur_dat$lat
    }else{
      prev_dx <- cur_fish$lon 
      prev_dy <- cur_fish$lat 
    }
    
    bearing <- as.numeric(atan2(prev_dy, prev_dx))
    
    correlated_bearing <- (bearing  + rnorm(1, mean = 0, sd = prey.noise_sd))%% (2 * pi)
    
    new_x <- prey[i,]$x + (prey[i,]$speed*cos(correlated_bearing)*cvar)+
      (prey[i,]$speed*cos(attr_bearing)*(1-cvar))
    
    new_y <- prey[i,]$y + (prey[i,]$speed*sin(correlated_bearing)*cvar)+
      (prey[i,]$speed*sin(attr_bearing)*(1-cvar))
    new_point<-st_point(c(new_x, new_y))
    check<-terra::extract(bof.rst, cbind(new_x,new_y))[1,]
    while (check==0) {
      rnd.dir<-runif(1,0,6.28319)
      new_x <- prey$x[i] + prey[i,]$speed*cos(rnd.dir)
      new_y <- prey$y[i] + prey[i,]$speed*sin(rnd.dir)
      check<-terra::extract(bof.rst, cbind(new_x,new_y))[1,]
      new_point <- st_point(c(new_x, new_y))
      #print(paste(prey[i,"id"],"_",step,"stuck"))
    }
    
    prey[i,]$x <- new_x
    prey[i,]$y <- new_y
    
  }# End of Prey loop
  prey$time<-step
  
  for(i in 1:num_fish){
    Lefin<-terra::extract(finish, cbind(prey[i,]$x, prey[i,]$y))[1,]
    if (!is.na(Lefin)) {
      escaped <- rbind(escaped, prey[i, ])
      prey <- prey[-i, ]
      num_fish <- num_fish - 1  # Update the number of prey
    }
  }
  
  Sys.sleep(.01)
  current_data <- prey
  current_data$time <- rep(step, nrow(current_data))
  sim_data[[step]] <- current_data
  setTxtProgressBar(pb, step)
  if (nrow(prey)==0) break
  #setTxtProgressBar(pb, step)
  rtx <- tail(compact(sim_data),2*20*12) %>% # @90sec ping rate this is 12 hour tail
    bind_rows()
  xvals <- split(rtx$x, rtx$id)
  yvals <- split(rtx$y, rtx$id)
  
  # Create the scatter plot that follows A.prey
  plot(bof.rst)
  # Create lines for movement
  pst <- rtx %>% arrange(id, time) %>% group_by(id) %>% slice(n())
  mapply(lines, xvals, yvals, col = "black", lwd = 1)
  
  # Add points
  points(pst[, c('x', 'y')], col = "black",
         bg = 'grey', cex = 1, pch = 21)
  #plot(finish,add=T)
  #print(paste(step,"_",round(Sys.time() - ptm,1)))
  # Add "x" points for consumed prey
  #points(eaten[eaten$time == step, ]$x, eaten[eaten$time == step, ]$y, pch = "x")
  
}
)


escaped <- data.frame(id=NA,y=NA,x=NA,size=NA,time=NA,speed=NA) 

