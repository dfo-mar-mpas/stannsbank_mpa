library(sf)
library(ggplot2)
library(RColorBrewer)
library(arcpullr)
library(marmap)
library(stars)
library(lme4)
library(MASS)
library(parallel)
library(furrr)
library(future)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)

batch_n <- as.numeric(args[1])
batch_n <- ifelse(is.na(batch_n),15,batch_n)

b <- as.numeric(args[2])
b <- ifelse(is.na(b),1,b)

print(paste("batch_n =",batch_n))
print(paste("b =",b))
# Sys.sleep(b)
# saveRDS(data.frame(b,batch_n),paste0("power_crab_batch",b,"_sp_",".RDS"))

# get SAB
url <- "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Oceans_Act_Marine_Protected_Areas/MapServer/0"
sab <- get_spatial_layer(url,where="NAME_E='St. Anns Bank Marine Protected Area'")

# get crab
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




# get bathy
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
  st_join(benthoscape) %>%
  mutate(Class_name=if_else(is.na(Class_name),
                            "Unclassified",
                            Class_name)) %>%
  dplyr::select(STATION,z,Class_name) %>%
  mutate(STATION = as.character(STATION),
         z_class = if_else(z>=(-100),
                           "<100 m",
                           if_else(z>=(-150),
                                   "100 to 150 m",
                                   ">150 m")))

###################### richness #####################

richness <- crab %>%
  group_by(STATION,BOARD_DATE) %>%
  reframe(species=list(sp_name),
          Richness=length(unique(sp_name))) %>%
  mutate(BOARD_DATE = as.Date(BOARD_DATE),
         STATION = as.character(STATION),
         YEAR = as.numeric(format(BOARD_DATE,'%Y')),
         sp="Richness") %>%
  left_join(crabstations %>% as.data.frame() %>% dplyr::select(-geometry),by="STATION")


# set up richness for BACI
richness_BACI <- richness %>%
  as.data.frame() %>%
  filter(YEAR>2013) %>%
  mutate(PERIOD=if_else(YEAR>2018,
                        "After",
                        "Before"))

# set up crab for BACI
crab_BACI <- crab %>%
  complete(sp_name, nesting(STATION, BOARD_DATE),fill=list(CPUE=0)) %>%  # add in 0's
  mutate(station_date=paste(STATION,BOARD_DATE,sep = ": "),
         year=station_date %>%
           gsub(".*: ","",.) %>%

           gsub("-.*","",.)) %>%
  filter(year>2013) %>%
  mutate(PERIOD=if_else(year>2018,
                        "After",
                        "Before"),
         STATION = as.character(STATION)) %>%
  left_join(crabstations %>% as.data.frame() %>% dplyr::select(-geometry),by="STATION")


# combine crab and richness for power analysis
BACI <- data.frame(sp=c(richness_BACI$sp,crab_BACI$sp_name),
                   dependent=c(richness_BACI$Richness,crab_BACI$CPUE),
                   rand1=c(richness_BACI$Class_name,crab_BACI$Class_name),
                   rand2=c(richness_BACI$z_class,crab_BACI$z_class),
                   PERIOD=c(richness_BACI$PERIOD,crab_BACI$PERIOD)) #%>%

# prepare empty results df
results <- expand.grid(sp=c("Richness",
                            "HIPPOGLOSSOIDES PLATESSOIDES",
                            "GLYPTOCEPHALUS CYNOGLOSSUS",
                            "SEBASTES",
                            "AMBLYRAJA RADIATA",
                            "GADUS MORHUA",
                            "CHIONOECETES OPILIO",
                            "ANARHICHAS LUPUS"
),
replicate=1:1000,
sites=c(5,
        10,
        25,
        50,
        100),
trend=c("Increase","Decrease"),
effect_size=seq(0.1,1,0.2),
p=NA
) %>%
  arrange(sp,replicate,sites,trend,effect_size) %>%
  group_by(batch=(row_number()-1) %/% (n()/batch_n)) %>%
  ungroup


############# simulate! ####################

b_data <- BACI %>% filter(sp %in% results$sp[results$batch==b])
for(s in unique(b_data$sp)){
  # fit model to species data
  sp_data <- BACI %>% filter(sp==s)
  realdatamodel <- lme4::glmer.nb(formula= dependent ~ PERIOD+(1|rand1)+(1|rand2),
                                  data = sp_data)

  # get parameters from sp model
  theta <- lme4::getME(realdatamodel,"glmer.nb.theta")
  periodeffect <- ranef(realdatamodel) %>%
    as.data.frame()
  spresults <- results[results$batch==b&results$sp==s,]

  # simualtion all variables for sp
  spresults$p <- lapply(1:nrow(spresults),function(row){
    # spresults$p <- future_map(1:nrow(spresults),function(row){
    p <- try(log("a"),silent = TRUE)
    while(inherits(p,"try-error")){
      p <- try({
        # randomly select random effect variables for sites
        r1 <- rep(sample(unique(BACI$rand1),spresults$sites[row],replace=TRUE),2)
        r2 <- rep(sample(unique(BACI$rand2),spresults$sites[row],replace=TRUE),2)
        r1index <- r1 %>%
          lapply(function(x){
            which(x==periodeffect$grp)
          }) %>%
          unlist() %>%
          as.numeric()
        r2index <- r2 %>%
          lapply(function(x){
            which(x==periodeffect$grp)
          }) %>%
          unlist()

        PERIOD <- rep(c("Before","After"),each=spresults$sites[row])

        # simulate dependent variable
        simdata <- data.frame(PERIOD=PERIOD,
                              rand1=r1,
                              rand2=r2,
                              dependent=rnegbin(n=spresults$sites[row]*2,
                                                mu=exp(fixef(realdatamodel)[1]+
                                                         rnorm(spresults$sites[row],
                                                               periodeffect$condval[r1index],
                                                               periodeffect$condsd[r1index])+
                                                         rnorm(spresults$sites[row],
                                                               periodeffect$condval[r2index],
                                                               periodeffect$condsd[r2index])),
                                                theta = theta)*(rep(c(1,if_else(spresults$trend[row]=="Increase",
                                                                                1+spresults$effect_size[row],
                                                                                1-spresults$effect_size[row])),each=spresults$sites[row]))
        )
        # simdata %>% group_by(PERIOD) %>% reframe(cpue=mean(dependent,na.rm=T));spresults[row,]
        gc()
        # model fit with and without PERIOD as a fixed effect
        simdatawithperiod <- lme4::glmer.nb(formula= dependent ~ PERIOD+(1|rand1)+(1|rand2),
                                            data = simdata)
        simdatawithoutperiod <- lme4::glmer.nb(formula= dependent ~ (1|rand1)+(1|rand2),
                                               data = simdata)

        # Likelihood Ratio Test
        LRT <- anova(simdatawithperiod,simdatawithoutperiod)
        spresults$p[row] <- LRT$`Pr(>Chisq)`[2]
        p <- LRT$`Pr(>Chisq)`[2]
        rm(r1,r2,r1index,r2index,PERIOD,simdata,simdatawithperiod,simdatawithoutperiod,LRT)
        gc()
        return(p)

      }, silent = TRUE)
      # if(inherits(p,"try-error")) browser()
    }
    return(p)
    # }
  }) %>% unlist()
  results[results$batch==b&results$sp==s,] <- spresults
  saveRDS(spresults,paste0("power_crab_batch",b,"_sp_",s,".RDS"))


}


