## Load required libraries
library(leaflet)
library(lwgeom)
library(readr)
library(gganimate)
library(viridis)
library(utils)
library(tidyterra)
library(sf);
library(nngeo);
library(tidyverse)
library(terra);
library(animation)
library(sfheaders)


setwd("C:/Temp/Wilson/")

# UTM CRS Bay of Fundy
locations_projected <- st_read("BoFshp/iho.shp")  %>% 
  st_transform(crs = "EPSG:32620") 

# Lat/Long CRS Bay of Fundy
locations_projected.ll<-locations_projected %>% st_as_sf %>% 
  st_transform(crs = "EPSG:4326") 

# Aggregated rastre version of Bay of Fundy
r<-rast(stars::st_rasterize(locations_projected %>% dplyr::select(geometry)))

#bof.rst <- aggregate(r, factor = 1)
bof.rst<-ifel(is.na(r), 0, r)
plot(bof.rst)

# Finish region that defines successful migrants
finish <- data.frame(x=c(-67.24392,-66.08700,-66.17797,-67.43117),
                         y=c(44.66460,44.17642,44.02959,44.54224)) %>% 
  sfheaders::sf_polygon(
    #obj = finish
    x = 'x',
    y = 'y') %>%  
  st_set_crs( 4326 ) %>% 
  st_transform(finish,crs = "EPSG:32620") 
finish<-rast(stars::st_rasterize(finish%>% dplyr::select(geometry)))
#finish.ll <- st_intersection(locations_projected.ll, finish) 

# Number of receivers to deploy within defined region
n_points_wanted = 70

# Define region to be populated with n receivers
test_lines <- data.frame(x=c(-66.18511,-66.77016,-67.09989,-66.19178,-65.49455          ),
                         y=c(45.16801,45.05906,44.75639,44.40122,44.81564))

#################################################################
### where to find a way to convert the csv data coordinates #####
#################################################################

#newpoly <- sfheaders::sf_polygon(
#  obj = test_lines
#  , x = 'x'
#  , y = 'y') %>%  
#st_set_crs( 4326 ) 
#inter <- st_intersection(locations_projected.ll, newpoly)
#grid <- st_sample(inter,n_points_wanted,type="regular") 
#st_crs(grid)<-4326


## Import data: Longitude and Latitude are set as column “x” and “y” respectively in that file
sdata <- read.csv("C:/Temp/Wilson/square_grid_lsr_to_obof.csv", sep=",")

## Create grid as an sfc object from Long and Lat coordinates
## Solution provided by https://ryanpeek.org/2017-08-03-converting-xy-data-with-sf-package/

grid <- st_as_sf(sdata, coords = c("x", "y"), crs = 4326)

grid.utm<-grid %>% 
  st_transform(crs = "EPSG:32620")

#plotting:
ggplot() + 
  geom_spatraster(data=finish,na.rm = TRUE)+
  geom_sf(data = locations_projected,fill="transparent") +
  geom_sf(data = grid, color = 'red')+
  scale_fill_continuous(na.value = "transparent")

# Number of desired particle to simulate
num.prey<-1000
bdlng.sec<-5  # swim speed
ping.rate=90

# particle attributes. Don't ask why they're called prey
prey <- data.frame(id = paste("B.prey", 1:num.prey, sep = ""),
                     y = runif(num.prey, 45.37619, 45.47266), #Start y location
                     x = runif(num.prey, -65.20452, -64.97920), #Start x location
                     size = rnorm(num.prey, mean = .15, sd = 0.015),#random size
                     time=1) %>% 
  mutate(speed =  size*
           bdlng.sec* # Body lengths per second
           ping.rate # Average ping rate
         ) #Speed as a function of size (bodylengths/sec)

prey.sf<-prey%>% st_as_sf(coords = c("x","y"))%>% st_set_crs( 4326) %>% 
  st_transform(crs = "EPSG:32620")

prey<-prey.sf  %>%  as_tibble() %>% select(-geometry) %>% 
cbind(data.frame(st_coordinates(prey.sf[,1]))) %>% 
  rename(x=X,y=Y)

escaped <- data.frame(id=NA,y=NA,x=NA,size=NA,time=NA,speed=NA) 

# current trajectory versus attractant movement ratio
cvar<-.95 #(eg. 94% of movement is based on current trajectory whereas 6% is biased towards point of attraction)

# sd of random direction (units in radians)
prey.noise_sd<-.22

# the general radian direct from seed location towards Gulf of Maine
dir.rad<-1.185

# Maximum number of steps within sim
num_steps<-12500

# Simulation progress bar
pb <- txtProgressBar(min = 0, max = num_steps, style = 3)

#Create empty list to populate for efficiency
sim_data <- rep(list(NULL),num_steps)
system.time(for (step in 1:num_steps) {
  ptm <- Sys.time()
  for (i in 1:num.prey) {
          cur.prey<-prey[i,]
          attr<-pi*dir.rad
          attr_bearing <- attr%% (2 * pi)
          if(step>2){cur.dat<-sim_data[[step-2]] %>% filter(id==cur.prey$id)}else{cur.dat<-data.frame()}
          if(nrow(cur.dat)>0){
            prev_dx <- cur.prey$x -cur.dat$x  
            prev_dy <- cur.prey$y -cur.dat$y
          }else{
            prev_dx <- cur.prey$x 
            prev_dy <- cur.prey$y 
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
  
  for(i in 1:num.prey){
    Lefin<-terra::extract(finish, cbind(prey[i,]$x, prey[i,]$y))[1,]
    if (!is.na(Lefin)) {
      escaped <- rbind(escaped, prey[i, ])
      prey <- prey[-i, ]
      num.prey <- num.prey - 1  # Update the number of prey
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

close(pb)
sim_data.df <- do.call(rbind, sim_data)

dist.threshold<-650

grid.utm<-st_as_sf(grid)%>% 
  st_transform(crs = "EPSG:32620") 

rd_walk_full<-st_as_sf(sim_data.df,coords =c("x","y") ,crs="EPSG:32620" )

test<-st_distance(rd_walk_full,
                  grid.utm) %>%
  as.data.frame.table() %>% 
  pivot_wider(names_from = Var2,values_from = Freq)  

detects<-  cbind(sim_data.df,test) %>% 
  dplyr::select(-Var1) %>% 
  pivot_longer(cols = -c(id:y),names_to = 'Receiver',values_to = 'Distance.m') %>% 
  filter(Distance.m<=units::set_units(dist.threshold, 'm'))

grid.m<-grid%>%   st_transform(32620) %>% 
  st_coordinates() %>% as.data.frame() %>% 
  mutate(Receiver=names(test)[-1]) %>% 
  rename(deploy.x=X,deploy.y=Y)

detects<-detects %>% 
  left_join(grid.m,by="Receiver")

# ADD DISTANCE TO PREDICTED DETEION TO BE CONSERVATIVE AND COMPENSATE FOR DEPTH OF RECEIVER TO FISH AT SURFACE
predict.penalty<-100 # this penalty is in meters

mydata<-readRDS("mydata.rds")
# Run logistic regression model
model <- glm(pDet ~Dist, data=mydata, family=binomial(link="logit"))
detects<-cbind(detects,
               as.data.frame(predict(model, newdata = detects %>% 
                                       rename(Dist=Distance.m) %>% 
                                       mutate(Dist=as.numeric(Dist)+predict.penalty), 
                                             type="link", se=TRUE))
)

detects$fit <- model$family$linkinv(detects$fit)  # Rescale to 0-1

ggplot(detects, aes(x=as.numeric(Distance.m)))+ 
  geom_line(aes(y=fit)) + 
  scale_y_continuous(limits = c(0,1))+
  labs(x="Distance", y="Prob of Detection",title=paste("Penalized mean probability of detection for V7 at high power\nacross locations as a function of distance "),
       caption="Note: This represents the probability of a\nsingle transmission being detected")

Prob.detection<-detects %>% 
  group_by(id,Receiver) %>% 
  mutate(p.det.rec=1-prod(1-fit)) %>% #the probability of being detected at least once by a receiver
  filter(time == min(time)) %>%
  group_by(id) %>% 
  mutate(p.det.grid=1-(prod(1-p.det.rec))) %>% #the probability of being detected at least once within the grid
  arrange(id) %>% 
  ungroup()

ids.det=length(unique(Prob.detection$id))
ids.simed=length(unique(sim_data.df$id))

#=======================================
# MUY IMPORTANTE
# THE VALUE BELOW REPRESENTS THE OVERALL 
# DETECTION PROBABILITY OF THE GRID
#=======================================
Prob.detection %>% 
  distinct(id,p.det.grid) %>% 
  summarise(sum(p.det.grid)/200)

ggplot(Prob.detection)+
  geom_bar(aes(x=reorder(interaction(id,Receiver),-p.det.rec),y=p.det.rec,fill=interaction(id,Receiver)),stat="identity",position = 'dodge',show.legend = F)+
  theme(axis.text.x= element_blank(),
        axis.ticks.x= element_blank())+
  labs(subtitle = "Probabilty of being detected per adjacent receiver within the grid",
       x="",y="P(detection) at receiver",caption = paste("Detetions occured for",ids.det,"out of",ids.simed))+
  scale_fill_viridis_d()

Prob.detection %>% 
  dplyr::select(id,p.det.grid) %>% 
  mutate(id=reorder(id,-p.det.grid)) %>% 
  distinct() %>% 
  arrange(p.det.grid) %>% 
ggplot()+
  geom_bar(aes(x=(id),y=p.det.grid,fill=id),stat="identity",position = 'dodge',show.legend = F)+
  theme(axis.text.x= element_blank(),
        axis.ticks.x= element_blank())+
  labs(subtitle ="Probability of being detected within the grid",
       x="","P(detection) within grid",caption = paste("Detetions occured for",ids.det,"out of",ids.simed))+
  scale_fill_viridis_d()

ids.det=length(unique(Prob.detection$id))
ids.simed=length(unique(sim_data.df$id))

Prob.detection %>% 
  distinct(id,p.det.grid) %>% 
  summarise(sum(p.det.grid)/ids.simed)

#Count of receivers within dist.threshold per transmitter
detects %>% 
  group_by(id,Receiver) %>% 
  summarise(n=n()) %>% 
  group_by(id) %>% 
  summarise(n=n()) %>% 
  ggplot()+
  geom_bar(aes(x='zx',fill=factor(n)),position = "stack",color="grey20")+
  scale_y_continuous(name="receiver count by particle",limits = c(0,200),breaks=seq(0,200,25))

detects %>% 
  group_by(id,Receiver) %>% 
  summarise(n=n()) %>% 
  arrange(n) %>% 
  ggplot()+
  geom_bar(aes(x=id,y=n,fill=Receiver),stat="identity",show.legend = F)+
  #scale_y_continuous(name="receiver count by particle",limits = c(0,200),breaks=seq(0,200,25))+
  coord_flip()+
    theme(axis.text.y= element_blank(),
          axis.ticks.y= element_blank())+
    labs(caption ="Fill color reflects number \nof potential detections per receiver")

new_lines <- sim_data.df %>%
  st_as_sf(coords = c("x", "y"), crs = 32620) %>%
  group_by(id) %>%
  dplyr::summarize(do_union=FALSE) %>%  # do_union=FALSE doesn't work as well
  st_cast("LINESTRING") 
Prob.detects.sf<-Prob.detection %>% st_as_sf( coords = c("x","y")) %>% st_set_crs(32620)
detects.sf<-detects %>% st_as_sf( coords = c("x","y")) %>% st_set_crs(32620)

ggplot() +
  geom_sf(data = locations_projected) +
  geom_sf(data=new_lines,alpha=.05,linewidth=.05)+
  geom_sf(data = grid.utm, col ="red")+
  geom_sf(data = Prob.detects.sf,aes(col =p.det.grid))+
  scale_color_viridis()
 
  
leaflet(locations_projected.ll) %>% 
  addPolygons(weight=.2) %>% 
  addPolylines(data=new_lines %>% st_transform(4326),
               #fill=F,
               color = "black",
               #fillColor = 'grey30', 
               fillOpacity=.1,
               opacity=.1
               ) %>% 
  addCircleMarkers(radius=2,color="yellow",opacity = .4,data=grid%>% st_cast("POINT") %>% 
                     as("Spatial") %>%  st_as_sf(st_transform(4326))) %>% 
  addCircleMarkers(radius=2,color="red",opacity = 1,data=detects.sf%>% st_transform(4326))  %>%
  addCircleMarkers(radius=2,color="yellow",opacity = .4,data=grid%>% st_cast("POINT") %>% 
                     as("Spatial") %>%  st_as_sf(st_transform(4326))) %>% 
  addMeasure(
    position = "bottomleft",
    primaryLengthUnit = "meters") %>%
  leafem::addMouseCoordinates()




grid.m<-grid%>%   st_transform(32620) %>% 
  st_coordinates() %>% as.data.frame()
AnimTime=60

system.time(saveVideo({
  ani.options(interval = AnimTime/(max(sim_data.df$time)/1))
  for(i in seq(1,max(sim_data.df$time))) {
    rtx <- sim_data.df[between(sim_data.df$time, i - (2*20*12), i), ]
    xvals <- split(rtx$x, rtx$id)
    yvals <- split(rtx$y, rtx$id)
    # Create the scatter plot that follows A.prey
    plot(bof.rst,
         xlim = c(min(sim_data.df[sim_data.df$time==i,]$x-25000), 
                  max(sim_data.df[sim_data.df$time==i,]$x+25000)),
         ylim = c(min(sim_data.df[sim_data.df$time==i,]$y-25000), 
                  max(sim_data.df[sim_data.df$time==i,]$y+25000)),
         main = paste("Step = ", i)
    )
    # Create lines for movement
    pst <- rtx %>% arrange(id, time) %>% group_by(id) %>% slice(n())
    mapply(lines, xvals, yvals, col = "black", lwd = 1)
    # Add points
    points(pst[, c('x', 'y')], col = "black",
           bg = "red", cex = 1, pch = 21)
    points(grid.m[, c('X', 'Y')],
           col = "black",
           bg = "yellow", cex = 1, pch = 21)
    points(detects[detects$time==i, c('deploy.x', 'deploy.y')],
           col = "black",
           bg = "purple", cex = 2, pch = 21)
    print(i)
  }
  #ani.pause()  ## pause for a while ('interval')
  
},#End of Vid For Loop
video.name = paste("sim",'bodylng_sec',bdlng.sec,'avg_ping_rate',ping.rate,"cvar",cvar,"dir.rad",dir.rad,"sd.noise",prey.noise_sd,".mp4",sep="_"), 
ani.width = 1000, 
ani.height = 600
))
