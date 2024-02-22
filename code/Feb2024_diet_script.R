install.packages("readxl")
library(readxl)
library(arcpullr)
library(leaflet)
library(ggplot2)
library(dplyr)
library(sf)
library(taxize)
library(vegan)
library(rphylopic)

# FIGURE 1. Map
url <- "https://gisp.dfo-mpo.gc.ca/arcgis/rest/services/FGP/Oceans_Act_Marine_Protected_Areas/MapServer/0"
sab <- get_spatial_layer(url,where="NAME_E='St. Anns Bank Marine Protected Area'")
diet <- read.csv("data/CrabSurvey/MPA.Diet.SnowCrabSurvey.Feb.2024.csv")
sabdiet <- diet %>% filter(SLATDD > 45.1)
dim(sabdiet)

latitude <- sabdiet$SLATDD
longitude <- sabdiet$SLONGDD

myMap <- leaflet(sab) %>%
  addProviderTiles(providers$Esri.OceanBasemap) %>%
  addPolygons(popup = paste(sab$NAME_E,sab$ZONE_E))
myMap <- myMap %>%
  addCircleMarkers(data = diet, ~longitude, ~latitude, color="gray")
myMap

ggplot()+
  #geom_sf(data=novsco, fill=gray(.9),size=0)+
  #geom_sf(data=sab_benthoscape,aes(fill=class),lwd=0.25)+
  geom_sf(data=sab,colour="blue", fill=NA, linewidth=1)+
  #geom_sf(data=gully2, colour="red", fill=NA)+
  coord_sf(xlim=c(-60.5, -58.3), ylim=c(45.25,47.1), expand=F)+
  geom_point(data=sabdiet, aes(x=longitude, y=latitude))+
  labs(y="Latitude",x="Longitude")+
  theme_bw()+
  theme(panel.background = element_rect(fill="lightcyan1"), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title = element_text(size=16), legend.text = element_text(size=10), legend.position = "right")




# FIGURE 2: SPECIES ACCUMULATION
species <- diet$PREYSPECCD
codes <- read.csv("data/CrabSurvey/GROUNDFISH_GSSPECIES_ANDES_20230901.csv")

df <- diet[, c("SLATDD", "SLONGDD", "PREYSPECCD")]
df$species <- 0
for (i in seq_along(df$PREYSPECCD)) {
  if (length(which(codes$SPECCD_ID == df$PREYSPECCD[i])) > 0) {
    df$species[i] <- codes$SPEC[which(codes$SPECCD_ID == df$PREYSPECCD[i])]
  } else {
    df$species[i] <- "NA"
  }
} # PREDATOR AND PREY SPECIES

SPECIES <- unique(df$species)
sKey <- data.frame(species=SPECIES,
                   group=0)
for (i in seq_along(SPECIES)) {
  s <- SPECIES[i]
  if(is.na(s)|grepl("EGGS UNID",s)| s== "" | s == "NA"){
    s <- "NA"
  } else if(grepl("PORIFERA P.",s)){
    s <- "Porifera"
  } else if(grepl("ASTROTECTEN",s)){
    # typo in database as per https://publications.gc.ca/collections/collection_2015/mpo-dfo/Fs97-6-3149-eng.pdf
    s <- gsub("ASTROTECTEN","ASTROPECTEN",s)
  } else if(grepl("^ASTERIAS|SOLASTER PAPPOSUS
",s)){
    # Asteris genus alone returning a plant species, adding rubens correctly returns a starfish
    # Solaster returns a crustacean
    s <- "Asterias rubens"
  } else if(grepl("PURSE SMOOTH SKATE|RAJA EGGS",s)){
    # Could not get reasonable results automatically. Manually defining real scientific name
    s <- "Malacoraja senta"
  } else if(grepl("PURSE THORNY SKATE",s)){
    # Could not get reasonable results automatically. Manually defining real scientific name
    s <- "Amblyraja radiata"
  } else if(grepl("APHRODITA .*|POLYCHAETA .*",s)){
    # Could not get reasonable results automatically. Manually defining a scientific name
    s <- "Glycera dibranchiata"
  } else if(grepl("BOLTENIA .*|ASCIDIA .*",s)){
    # Could not get reasonable results automatically. Manually defining a scientific name
    s <- "Boltenia ovifera"
  } else if(grepl("LIPARIS SP.|LYCODES SP.|MACROZOARCES AMERICANUS",s)){
    # Could not get reasonable results automatically. Manually defining a scientific name
    s <- "Liparis liparis"
  } else if(grepl("CHLAMYS ISLANDICA",s)){
    # Could not get reasonable results automatically. Manually defining a scientific name
    s <- "PLACOPECTEN MAGELLANICUS"
  } else if(grepl("SCLEROCRANGON SP.",s)){
    # Could not get reasonable results automatically. Manually defining a scientific name
    s <- "SCLEROCRANGON boreas"
  }
  
  tN <- tax_name(s, get = 'subphylum', db = 'itis',ask=FALSE)$subphylum
  if (!(is.na(tN))) {
    if (!(tN == "Vertebrata")) {
      sKey$group[i] <- tN
    } else {
      tN <- tax_name(s, get = 'class', db = 'itis',ask=FALSE)$class
      sKey$group[i] <- tN
    }
  } else {
    sKey$group[i] <- "NA"
  }
}

sKey$group[which(is.na(sKey$group))] <- "NA"

# Add subphylum to df
df$subphylum <- 0
for (i in seq_along(sKey$species)) {
  keep <- which(df$species == sKey$species[i])
  df$subphylum[keep] <- sKey$group[i]
  
}
names(df) <- c("lat", "lon", "PRESPECCD", "species","subphylum")


DF<- df %>%
  mutate(lat_rounded = round(lat, 2), lon_rounded = round(lon, 2)) %>%
  group_by(lat_rounded, lon_rounded, subphylum) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = subphylum, values_from = count, values_fill = 0)
dietdata <- DF

# Diet Species (for accumulation)
df$species[which(df$species == "")] <- "NA"
DF2<- df %>%
  mutate(lat_rounded = round(lat, 2), lon_rounded = round(lon, 2)) %>%
  group_by(lat_rounded, lon_rounded, species)%>%
  summarise(count = n())%>%
  pivot_wider(names_from = species, values_from = count, values_fill = 0)


names(DF2)[1:2] <- c("lat", "lon")
dietdata2 <- DF2
# NEW
dietdata2 <-  dietdata2 %>%
  dplyr::select(-lat, -lon)
spacc2 <- specaccum(dietdata2,"random")

# prep data for ggplot
plotdata <- data.frame(Sites=spacc2$sites,
                       Richness=spacc2$richness,
                       SD=spacc2$sd,
                       Source = "Diet")
ggplot(plotdata) +
  geom_line(aes(x=Sites, y=Richness,col=Source), linewidth=1.5) +
  geom_ribbon(aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD),fill=Source),alpha=0.2)+
  ylab(label = "Species Richness")+
  coord_cartesian(ylim = c(0, 170), xlim = c(0, length(plotdata$Sites)+20))+
  theme_bw()+
  theme(legend.position = "none", text=element_text(size=21))

ggsave(filename = "DietData_Specaccum.png",plot = last_plot(), device = "png", path = "output/CrabSurvey/", width = 10, height=8, units = "in", dpi = 400)

# FIGURE 3: SPECIES ABUNDANCE

load(system.file("data/phylopic", "phylopic.RData", package = "SABapp"))

new_df <- df %>%
  group_by(species) %>%
  summarize(abundance = n()) %>%
  ungroup()


phylopic.sp <- get_uuid(name = new_df$species, n=1)

new_df$photo <- 0
for (i in seq_along(new_df$species)) {
  k <- which(phylopic$originalName == new_df$species[i])
  new_df$photo[i] <- unname(phylopic$photo[k])
}


ggplot(new_df %>% filter(abundance > 10)) +
  geom_col(aes(x=`species`,y=log(abundance)))+
  #geom_phylopic(aes(x=`species`,y=`abundance`+200,img=photo),size=0.5)+
  coord_flip()+
  theme_classic()+
  labs(title = "Taxonomic Group Abundance",
       x = "Taxonomic Group",
       y = "Abundance")


# Figure 4: Species diets per predator 
#Add predator species
diet2 <- merge(sabdiet, codes, by.x="SPEC", by.y="SPECCD_ID", all.x=T)
diet2 <- diet2 %>% rename(PRED.COMM= COMM)
diet2 <- diet2 %>% rename(PRED.SPECIES= SPEC.y)

diet3 <- merge(diet2, codes, by.x="PREYSPECCD", by.y="SPECCD_ID", all.x=T)
diet3 <- diet3 %>% rename(PREY.COMM=COMM)
diet3 <- diet3 %>% group_by(PRED.COMM)

#get stations to figure out which ones are inside/outside the MPA
stns <-read.csv("data/CrabSurvey/StationInfo.csv")
catch<- catch %>% select(c(TRIP_ID, TRIP, SET_NO, STATION))
diet4<-left_join(diet3, catch, by=c("MISSION", "SETNO")) %>% data.frame
diet5 <- merge(diet4, stns, by.x="STATION", by.y="Station")

prey.tab <- janitor::tabyl(diet3, PRED.COMM, PREY.COMM)
prey.tab <- reshape2::melt(prey.tab, na.rm=T)
prey.tab <- prey.tab[!grepl("NA_", prey.tab$variable),]
prey.tab <- prey.tab[!grepl("emptystring_", prey.tab$variable),]
#prey.tab <- prey.tab[!grepl("", prey.tab$variable),]


ggplot(prey.tab, aes(x=PRED.COMM, y=value, fill=variable))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(legend.position = "bottom", axis.text.x=element_text(angle=90, vjust=0.5))




#Plot up stomach fullness histogram
ggplot()+
  geom_histogram(data=diet5, aes(FULLNESS, colour=Inside, fill=Inside), binwidth = 1, position="dodge")+
  facet_wrap(vars(Year))+
  ylab(label = "Count")+
  xlab(label = "Stomach Fullness")+
  theme_bw()+
  theme(text=element_text(size=21))

ggsave(filename = "StomachFullness.png", plot = last_plot(), device = "png", path = "output/CrabSurvey/", width = 10, height=8, units="in", dpi=400)








