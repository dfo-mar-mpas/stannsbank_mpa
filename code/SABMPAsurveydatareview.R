#This script creates the required tables for the SABMPA post activity report
#---- set up-----

MPA.SAB.data = function (current_year) {
  
  ######NEED TO MODIFY THE REST SIMILAR TO MPA.gully.data############

getwd()
#setwd("C:/Users/GlassA/Desktop/covid copies/SABMPA")

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(viridis)


# ----get catch data for all tows inside the MPA ----
require(ROracle)
con= dbConnect(DBI::dbDriver("Oracle"), oracle.username, oracle.password, oracle.server)
#select data from the 14 tows in the SABMPA
sabdat=dbGetQuery(con, (" SELECT trip.trip_id, trip.trip,
trip.board_date, st.set_no, st.station, pr.latitude, pr.longitude,
 st.est_catch, ca.speccd_id, ca.est_num_caught, ca.est_discard_wt
FROM istrips trip, isgears gr, isfishsets st, iscatches ca, issetprofile pr
WHERE trip.tripcd_id = 7061
AND trip.trip_id = gr.trip_Id
AND pr.PNTCD_ID=2
AND (trip.trip_id = st.trip_Id
    AND gr.gear_id = st.gear_id)
AND (st.fishset_id = pr.fishset_Id
    AND st.set_no = pr.set_no)
AND (st.fishset_id = ca.fishset_id
AND st.haulccd_id = '1')
and st.station in ('016', '017', '018', '029', '031', '032', '034', '091', '093', '205', '511', '609', '611', '801')  
order by board_date, station, speccd_id"))

str(sab)
#separate date
sab$BOARD_DATE <- format(as.Date(sab$BOARD_DATE, tryFormats = c("%Y-%m-%d", "%Y/%m/%d"),"%Y"))
sab = sab %>% 
  mutate(BOARD_DATE = ymd(BOARD_DATE)) %>% 
  mutate_at(vars(BOARD_DATE), funs(year, month, day))
summary(sab)

#number of distinct species by year
sab %>%
  group_by(year) %>%
  summarise(unique = n_distinct(SPECCD_ID))

#change speccd_id to name and group certain species
data1=sab

data1$SPECCD_ID[data1$SPECCD_ID == 10] <- "Atlantic Cod"
data1$SPECCD_ID[data1$SPECCD_ID == 11] <- "Haddock"
data1$SPECCD_ID[data1$SPECCD_ID == 12] <- "White Hake"
data1$SPECCD_ID[data1$SPECCD_ID == 13] <- "Red Hake"
data1$SPECCD_ID[data1$SPECCD_ID == 14] <- "Silver Hake"
data1$SPECCD_ID[data1$SPECCD_ID == 15] <- "Cusk"
data1$SPECCD_ID[data1$SPECCD_ID == 16] <- "Pollock"
data1$SPECCD_ID[data1$SPECCD_ID == 17] <- "Tom Cod"
data1$SPECCD_ID[data1$SPECCD_ID == 23] <- "Redfish spp."
data1$SPECCD_ID[data1$SPECCD_ID == 30] <- "Halibut"
data1$SPECCD_ID[data1$SPECCD_ID == 31] <- "Turbot"
data1$SPECCD_ID[data1$SPECCD_ID == 40] <- "American Plaice"
data1$SPECCD_ID[data1$SPECCD_ID == 41] <- "Witch Flounder"
data1$SPECCD_ID[data1$SPECCD_ID == 42] <- "Yellowtail Flounder"
data1$SPECCD_ID[data1$SPECCD_ID == 43] <- "Winter Flounder"
data1$SPECCD_ID[data1$SPECCD_ID == 44] <- "Gulfstream Flounder"
data1$SPECCD_ID[data1$SPECCD_ID == 50] <- "Striped Atlantic Wolffish"
data1$SPECCD_ID[data1$SPECCD_ID == 51] <- "Spotted Wolffish"
data1$SPECCD_ID[data1$SPECCD_ID == 52] <- "Northern Wolffish"
data1$SPECCD_ID[data1$SPECCD_ID == 60] <- "Herring"
data1$SPECCD_ID[data1$SPECCD_ID == 61] <- "American Shad"
data1$SPECCD_ID[data1$SPECCD_ID == 62] <- "Alewife"
data1$SPECCD_ID[data1$SPECCD_ID == 64] <- "Capelin"
data1$SPECCD_ID[data1$SPECCD_ID == 70] <- "Mackerel"
data1$SPECCD_ID[data1$SPECCD_ID == 112] <- "Longfin Hake"
data1$SPECCD_ID[data1$SPECCD_ID == 114] <- "Fourbeard Rockling"
data1$SPECCD_ID[data1$SPECCD_ID == 118] <- "Greenland Cod"
data1$SPECCD_ID[data1$SPECCD_ID == 142] <- "Four-spot Flounder"
data1$SPECCD_ID[data1$SPECCD_ID == 143] <- "Windowpane Flounder"
data1$SPECCD_ID[data1$SPECCD_ID == 201] <- "Thorny Skate"
data1$SPECCD_ID[data1$SPECCD_ID == 202] <- "Smooth Skate"
data1$SPECCD_ID[data1$SPECCD_ID == 211] <- "Skates"
data1$SPECCD_ID[data1$SPECCD_ID == 241] <- "Northern Hagfish"
data1$SPECCD_ID[data1$SPECCD_ID == 300] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 301] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 302] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 303] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 304] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 306] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 311] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 317] <- "Sculpin NS"
data1$SPECCD_ID[data1$SPECCD_ID == 320] <- "Sea Raven"
data1$SPECCD_ID[data1$SPECCD_ID == 323] <- "Sculpins"
data1$SPECCD_ID[data1$SPECCD_ID == 350] <- "Atlantic Sea Poacher"
data1$SPECCD_ID[data1$SPECCD_ID == 351] <- "Alligatorfish"
data1$SPECCD_ID[data1$SPECCD_ID == 400] <- "Monkfish"
data1$SPECCD_ID[data1$SPECCD_ID == 410] <- "Marlinspike Grenadier"
data1$SPECCD_ID[data1$SPECCD_ID == 500] <- "Snailfish spp."
data1$SPECCD_ID[data1$SPECCD_ID == 501] <- "Lumpfish"
data1$SPECCD_ID[data1$SPECCD_ID == 502] <- "Spiny Lumpsucker"
data1$SPECCD_ID[data1$SPECCD_ID == 565] <- "Barracudina"
data1$SPECCD_ID[data1$SPECCD_ID == 590] <- "Sand Lances"
data1$SPECCD_ID[data1$SPECCD_ID == 619] <- "Eelpouts"
data1$SPECCD_ID[data1$SPECCD_ID == 622] <- "Snakeblenny"
data1$SPECCD_ID[data1$SPECCD_ID == 623] <- "Daubed Shanny"
data1$SPECCD_ID[data1$SPECCD_ID == 625] <- "Radiated Shanny"
data1$SPECCD_ID[data1$SPECCD_ID == 626] <- "Fourline Snakeblenny"
data1$SPECCD_ID[data1$SPECCD_ID == 627] <- "Eelpouts"
data1$SPECCD_ID[data1$SPECCD_ID == 630] <- "Wrymouth"
data1$SPECCD_ID[data1$SPECCD_ID == 640] <- "Ocean Pout"
data1$SPECCD_ID[data1$SPECCD_ID == 642] <- "Eelpouts"
data1$SPECCD_ID[data1$SPECCD_ID == 644] <- "Blennies,Shannies,Gunnels"
data1$SPECCD_ID[data1$SPECCD_ID == 647] <- "Eelpouts"
data1$SPECCD_ID[data1$SPECCD_ID == 712] <- "Barracudina"
data1$SPECCD_ID[data1$SPECCD_ID == 720] <- "Atlantic Saury"
data1$SPECCD_ID[data1$SPECCD_ID == 880] <- "Sculpin NS"
data1$SPECCD_ID[data1$SPECCD_ID == 1100] <- "Eggs Unid."
data1$SPECCD_ID[data1$SPECCD_ID == 1224] <- "Skate Eggs"
data1$SPECCD_ID[data1$SPECCD_ID == 1510] <- "Whelk Eggs"
data1$SPECCD_ID[data1$SPECCD_ID == 1821] <- "Sea Squirts"
data1$SPECCD_ID[data1$SPECCD_ID == 1823] <- "Sea Potato"
data1$SPECCD_ID[data1$SPECCD_ID == 2000] <- "Crustacea"
data1$SPECCD_ID[data1$SPECCD_ID == 2100] <- "Shrimp"
data1$SPECCD_ID[data1$SPECCD_ID == 2211] <- "Shrimp"
data1$SPECCD_ID[data1$SPECCD_ID == 2212] <- "Shrimp"
data1$SPECCD_ID[data1$SPECCD_ID == 2313] <- "Shrimp"
data1$SPECCD_ID[data1$SPECCD_ID == 2316] <- "Shrimp" 
data1$SPECCD_ID[data1$SPECCD_ID == 2411] <- "Shrimp" 
data1$SPECCD_ID[data1$SPECCD_ID == 2413] <- "Shrimp" 
data1$SPECCD_ID[data1$SPECCD_ID == 2415] <- "Shrimp" 
data1$SPECCD_ID[data1$SPECCD_ID == 2416] <- "Shrimp" 
data1$SPECCD_ID[data1$SPECCD_ID == 2511] <- "Jonah Crab"
data1$SPECCD_ID[data1$SPECCD_ID == 2513] <- "Atlantic Rock Crab"
data1$SPECCD_ID[data1$SPECCD_ID == 2521] <- "Lesser Toad Crab"
data1$SPECCD_ID[data1$SPECCD_ID == 2523] <- "Northern Stone Crab"
data1$SPECCD_ID[data1$SPECCD_ID == 2526] <- "Snow Crab"
data1$SPECCD_ID[data1$SPECCD_ID == 2527] <- "Toad Crab"
data1$SPECCD_ID[data1$SPECCD_ID == 2559] <- "Hermit Crab"
data1$SPECCD_ID[data1$SPECCD_ID == 2565] <- "Shrimp"
data1$SPECCD_ID[data1$SPECCD_ID == 2990] <- "Barnacles"
data1$SPECCD_ID[data1$SPECCD_ID == 3100] <- "Bristle Worm"
data1$SPECCD_ID[data1$SPECCD_ID == 3200] <- "Sea Mouse"
data1$SPECCD_ID[data1$SPECCD_ID == 4200] <- "Snails and Slugs"
data1$SPECCD_ID[data1$SPECCD_ID == 4210] <- "Whelks"
data1$SPECCD_ID[data1$SPECCD_ID == 4310] <- "Clams"
data1$SPECCD_ID[data1$SPECCD_ID == 4312] <- "Northern Propeller Clam"
data1$SPECCD_ID[data1$SPECCD_ID == 4321] <- "Sea Scallop"
data1$SPECCD_ID[data1$SPECCD_ID == 4322] <- "Icelandic Scallop"
data1$SPECCD_ID[data1$SPECCD_ID == 4340] <- "Cockles"
data1$SPECCD_ID[data1$SPECCD_ID == 4400] <- "Sea Slug"
data1$SPECCD_ID[data1$SPECCD_ID == 4430] <- "Canoe Shells"
data1$SPECCD_ID[data1$SPECCD_ID == 4511] <- "Shortfin Squid"
data1$SPECCD_ID[data1$SPECCD_ID == 4514] <- "Squid"
data1$SPECCD_ID[data1$SPECCD_ID == 4521] <- "Octopus"
data1$SPECCD_ID[data1$SPECCD_ID == 4522] <- "Bobtail Squid"
data1$SPECCD_ID[data1$SPECCD_ID == 5100] <- "Sea Spider"
data1$SPECCD_ID[data1$SPECCD_ID == 6100] <- "Sea Stars"
data1$SPECCD_ID[data1$SPECCD_ID == 6113] <- "Sea Stars"
data1$SPECCD_ID[data1$SPECCD_ID == 6115] <- "Mud Star"
data1$SPECCD_ID[data1$SPECCD_ID == 6117] <- "Horse Star"
data1$SPECCD_ID[data1$SPECCD_ID == 6119] <- "Blood Star"
data1$SPECCD_ID[data1$SPECCD_ID == 6121] <- "Purple Sunstar"
data1$SPECCD_ID[data1$SPECCD_ID == 6123] <- "Spiny Sunstar"
data1$SPECCD_ID[data1$SPECCD_ID == 6125] <- "Sea Stars"
data1$SPECCD_ID[data1$SPECCD_ID == 6128] <- "Sea Stars"
data1$SPECCD_ID[data1$SPECCD_ID == 6200] <- "Brittle Star"
data1$SPECCD_ID[data1$SPECCD_ID == 6300] <- "Basket Star"
data1$SPECCD_ID[data1$SPECCD_ID == 6400] <- "Sea Urchin"
data1$SPECCD_ID[data1$SPECCD_ID == 6411] <- "Sea Urchin"
data1$SPECCD_ID[data1$SPECCD_ID == 6413] <- "Sea Urchin"
data1$SPECCD_ID[data1$SPECCD_ID == 6500] <- "Sand Dollar"
data1$SPECCD_ID[data1$SPECCD_ID == 6511] <- "Sand Dollar"
data1$SPECCD_ID[data1$SPECCD_ID == 6600] <- "Sea Cucumbers"
data1$SPECCD_ID[data1$SPECCD_ID == 6717] <- "Sea Cucumbers"
data1$SPECCD_ID[data1$SPECCD_ID == 8300] <- "Anemones"
data1$SPECCD_ID[data1$SPECCD_ID == 8313] <- "Anemones"
data1$SPECCD_ID[data1$SPECCD_ID == 8318] <- "Sea Pen"
data1$SPECCD_ID[data1$SPECCD_ID == 8332] <- "Coral NS"
data1$SPECCD_ID[data1$SPECCD_ID == 8500] <- "Jellyfish"
data1$SPECCD_ID[data1$SPECCD_ID == 8520] <- "Jellyfish"
data1$SPECCD_ID[data1$SPECCD_ID == 8600] <- "Sponges"
data1$SPECCD_ID[data1$SPECCD_ID == 9300] <- "Seaweeds"
data1$SPECCD_ID[data1$SPECCD_ID == 641] <- "Arctic Eelpout"
data1$SPECCD_ID[data1$SPECCD_ID == 6127] <- "Sea Stars"
data1$SPECCD_ID[data1$SPECCD_ID == 6129] <- "Sea Stars"
data1$SPECCD_ID[data1$SPECCD_ID == 6213] <- "Ophiura sarsii"


data1 %>%  filter(SPECCD_ID == "Sea Urchin", year == 2016)

data1$Inside <- SABTOWS2$Inside

# separate date for new file
data1$BOARD_DATE <- as.Date(data1$BOARD_DATE, format = "%Y-%m-%d hh:mm:ss")
data1 = data1 %>% 
  mutate(BOARD_DATE = ymd(BOARD_DATE)) %>% 
  mutate_at(vars(BOARD_DATE), funs(year, month, day))
summary(data1)

# table of total species caught per year 
table1 <- data1 %>%
  group_by(year, SPECCD_ID, Inside, AREA_SWEPT) %>%
  summarise(numcaught = sum(EST_NUM_CAUGHT), wt = sum(EST_DISCARD_WT)) 
table1
write.table(table1,file="output/CrabSurvey/SABMPAspectableallyears.csv", sep=",")

table1a<-data1 %>% # individuals by spec by year
  group_by(SPEC, year) %>% 
  summarise(N = sum(EST_NUM_CAUGHT)) %>% 
  spread(year, N)
table1a
write.table(table1a,file="output/CrabSurvey/SABMPAspectableallyears_by_spec.csv", sep=",")

#this looks weird standardized by AREA_SWEPT, maybe change it up
specplot <- ggplot(table1 %>% filter(numcaught > 50)) +
  geom_point(aes(x = year, y = log(numcaught*(AREA_SWEPT/mean(AREA_SWEPT, na.rm=T))), col = SPECCD_ID), size = 3, shape = "circle") + 
  geom_line(aes(x = year, y = log(numcaught*(AREA_SWEPT/mean(AREA_SWEPT, na.rm=T))), col = SPECCD_ID)) +
  facet_wrap(vars(Inside), nrow=2)+
  xlab("Year") +
  ylab("Log(# Captured)") + 
  theme_bw(14) +
  scale_x_continuous(breaks = c(2015, 2016, 2017, 2018, 2019, 2021, 2022, 2023)) +
  scale_color_viridis(discrete = T, option ="D")+
  scale_fill_viridis(discrete = T)+
  #ggtitle("St Anns Bank MPA Captured Species (2015-2023)")+
  guides(col=guide_legend(title="Species"))
specplot

ggsave(filename = "AllSpeciesCaughtSAB.png", plot=specplot, device = "png", path = "output/CrabSurvey/", width = 10, height = 8, units = "in", dpi = 400)

# species caught per year, for current year
table2 <-data1 %>%
  filter(year == "2023") %>%
  group_by(SPEC) %>%
  summarise(total = sum(EST_NUM_CAUGHT), wt = sum(EST_DISCARD_WT))
table2

write.table(table2,file="output/SABMPAspectable2023.csv", sep=",")

#----stomach sample details----
require(ROracle)
con= dbConnect(DBI::dbDriver("Oracle"), oracle.username, oracle.password, oracle.server)
sabstom=dbGetQuery(con, ("SELECT stom.trip, stom.board_date, stom.set_no, sets.station, 
stom.speccd_id, stom.fish_no
FROM sncrabsets sets, snstomachdetails stom
WHERE sets.trip = stom.trip
AND sets.trip = stom.trip
And sets.set_no = stom.set_no
AND sets.setcd_id = '11'
order by trip, board_date, station, speccd_id"))



str(sabstom)
#separate date
sabstom$BOARD_DATE <- as.Date(sabstom$BOARD_DATE, format = "%Y-%m-%d hh:mm:ss")
sabstom = sabstom %>% 
  mutate(BOARD_DATE = ymd(BOARD_DATE)) %>% 
  mutate_at(vars(BOARD_DATE), funs(year, month, day))
summary(sabstom)

t1<-distinct(sabstom)# removes duplicates if there are any


# number of stomach samples per year, for current year, month is selected to only
#get the SABMPA samples
table3 <-t1 %>%
  filter(year == '2019', month == '9') %>%
  group_by(STATION) %>%
  count(SPECCD_ID)
table3
write.table(table3,file="SABMPAstombystation_spec_table.csv", sep=",")

# to get total SABMPA stomachs for current year
table4 <-t1 %>%
  filter(year == '2019', month == '9') %>%
  count(year) 
table4
write.table(table4,file="SABMPAstom_total.csv", sep=",")

#to get total number of stomach samples inside the SABMPA
table5in <-t1 %>%
  filter(year == '2019', month == '9', STATION %in% c('016', '018', '029', '032', '034', '091', '093', '511', '609')) %>%
  count(STATION)
table5in
write.table(table5in,file="SABMPAstom_in_by_tow_table.csv", sep=",")

table5in2 <-t1 %>%
  filter(year == '2019', month == '9', STATION %in% c('016', '018', '029', '032', '034', '091', '093', '511', '609')) %>%
  count(year)
table5in2
write.table(table5in2,file="SABMPAstom_in_total_table.csv", sep=",")

#to get total number of stomach samples taken adjacent to the SABMPA
table5out <-t1 %>%
  filter(year == '2019', STATION %in% c('036', '051', '204', '206', '507')) %>%
  count(STATION)
table5out
write.table(table5out,file="SABMPAstom_out_by_tow_table.csv", sep=",")

table5out2 <-t1 %>%
  filter(STATION %in% c('036', '051', '204', '206', '507'), year == '2019') %>%
  count(year)
table5out2
write.table(table5out2,file="SABMPAstom_out_total_table.csv", sep=",")

}
