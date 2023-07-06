### SABMPA_acoustic scripts ####################################################
# harri.pettitt-wade@dfo-mpo.gc.ca
# SABMPA acoustic telemetry data summaries
# temperature data from VR2AR receivers
################################################################################

rm(list = ls())
graphics.off()

##################################
# 0 - Package loading
##################################

# load packages
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(plotly)
library(tidyr)
library(stringr)
library(janitor)

#library(sp)
#library(sf)
#library(magrittr)
#library(scales)
#library(ggridges)

##################################
# 1 - Data loading
##################################

setwd(dir = "C:/Users/PETTITTWADEH/Documents/R/Projects/stannsbank_mpa")
getwd()

#### OTN line temp data - SABMPA, Halifax (hfx) and Cabot Straight (cbs) - this is the raw data from OTN
T20212022 <- read_csv("data/Temperature/cbs_hfx_sabmpa_temp_data_2021_to_2022.csv", 
                  col_types = cols(date = col_date(format = "%m/%d/%Y"), 
                                    data = col_number(),
                                    recover_ind = col_skip()))

# SAB crab lines data 2016-2020 & 2020-2021 (split by year from two separate raw data sets downloaded from OTN)

T2015 <- read_csv("data/Temperature/temp_data_2015.csv", 
                  col_types = cols(date = col_date(format = "%Y-%m-%d")))

T2016 <- read_csv("data/Temperature/temp_data_2016.csv", 
                  col_types = cols(date = col_date(format = "%Y-%m-%d")))

T2017 <- read_csv("data/Temperature/temp_data_2017.csv", 
                  col_types = cols(date = col_date(format = "%Y-%m-%d")))

T2018 <- read_csv("data/Temperature/temp_data_2018.csv", 
                  col_types = cols(date = col_date(format = "%Y-%m-%d")))

T2019 <- read_csv("data/Temperature/temp_data_2019.csv", 
                  col_types = cols(date = col_date(format = "%Y-%m-%d")))

T2020a <- read_csv("data/Temperature/temp_data_2020a.csv", 
                  col_types = cols(date = col_date(format = "%Y-%m-%d")))

#bind together so it represents the original data offloaded from OTN
T20152020a <- rbind(T2015,T2016,T2017,T2018,T2019,T2020a) 

summary(T20152020a) #check the data makes sense

T2020b <- read_csv("data/Temperature/temp_data_2020b.csv", 
                   col_types = cols(date = col_date(format = "%Y-%m-%d")))

T2021 <- read_csv("data/Temperature/temp_data_2021.csv", 
                   col_types = cols(date = col_date(format = "%Y-%m-%d")))

#bind together so it represents the original data offloaded from OTN
T2020b2021 <- rbind(T2020b,T2021) 

summary(T2020b2021) #check the data makes sense

#remove the objects that were merged to save memory
rm(list = c('T2015','T2016','T2017','T2018','T2019','T2020a',"T2020b",'T2021'))

# SAB crab lines data 2016-2020 - use this if you have the raw (not split) data from OTN
#SABMPA_temp_20162020_notclean <- read_csv("data/Temperature/V2LSCF_temp_data_2016_to_2020.csv", 
#                                                   col_types = cols(date = col_date(format = "%Y-%m-%d"), 
#                                                                    data = col_number(), units = col_skip(), 
#                                                                    recover_ind = col_skip()))

# SAB crab lines data 2020-2021 - use this if you have the raw (not split) data from OTN
#SABMPA_temp_20202021_notclean <- fread("data/Temperature/v2lscf_temps_2020-2021.csv", 
#                              sep = ",", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)



##################################
# 2 - Data cleaning + merging
##################################

# three separate SAB temperature data sets from OTN need cleaning prior to merge (V2LSCF 2015-2020, V2LSCF 2020-2021, OTN line 2021-2022)

### --------------------------------------------

# A. 2015-2020 SABMPA temperature data (North and South crab lines)
clean_me <- T20152020a

table(clean_me$description) # table counts of unique values in a column
table(clean_me$receiver) 
table(clean_me$station_name) 
class(clean_me)
str(clean_me)
head(clean_me)

clean_me <- clean_me[ ,c("station_name","date","receiver",
                         "description","data","rcv_serial_no",
                         "deploy_date","recover_date",
                         #"deploy_date_time","recover_date_time",
                          "dep_lat","dep_long","the_geom","catalognumber")] # pull out the columns to merge

setnames(clean_me, c("station_name","data",
                     "deploy_date","recover_date"),
         c("station","temperature_c",
           "deploy_date_time","recover_date_time"
         )) # change headers for consistency across datasets



## change station names so they are consistent
table(clean_me$station) 
clean_me$station <- gsub('V2LSCF-', '', clean_me$station)

clean_me$station <- clean_me$station %>%
  str_replace_all(c("101"="1","102"="2","103"="3","104"="4","105"="5","106"="6","107"="7","108"="8","109"="9","110"="10",
                    "111"="11","112"="12","113"="13","114"="14","115"="15","116"="16","117"="17","118"="18","119"="19",
                    "120"="20","121"="21","122"="22","123"="23"))

table(clean_me$station) # notable that some stations have many more rows of data than others (i.e., S1,2,4)
table(clean_me$station,clean_me$description) # check to see if the stats were collected consistently across stations
table(clean_me$description)

# adjust the date time class
clean_me <- as.data.table(clean_me %>% mutate(deploy_date_time = as.character(gsub("T"," ", deploy_date_time)))) # remove the 'T' if present for OTN format datetime
clean_me <- as.data.table(clean_me %>% mutate(recover_date_time = as.character(gsub("T"," ", recover_date_time))))

clean_me$deploy_date_time <- as.POSIXct(clean_me$deploy_date_time, tz = "UTC") # convert datetime to POSIXct UTC tz
clean_me$recover_date_time <- as.POSIXct(clean_me$recover_date_time, tz = "UTC") 

### identify and remove NAs ###
# may be an error here if some recover date times not provided - NAs for the most recent recoveries in the data set.
# identify and remove NAs
clean_me %>% summarise(across(everything(), ~ sum(is.na(.)))) # identify the NAs

#clean_me_recover_NAs <- clean_me %>%
#  filter(recover_date_time == "NA") # look at the data with NAs. #some recover date times are missing

#clean_me_NAs <- clean_me %>% filter_all(any_vars(is.na(.))) # filter out the rows with NAs

#str(clean_me_NAs) # look at the data with NAs (they are in recover datetime for 20162020 data)
#summary(clean_me_NAs) # most of the NAs are for Aug 2020 so most likely weren't updated for the final year

clean_me$date <- as.Date(clean_me$date,format="%Y-%m-%d")

#clean_me$date <- format(as.Date(clean_me$date,"%Y-%m-%d"))
class(clean_me$date)
head(clean_me)
summary(clean_me)

cleaning_dates <- as.data.frame(table(clean_me$date,clean_me$description))
# ambient_deg_c and ambient_temperature_deg_c appear duplicated (?)
# -> they both occur 1401 times on 2020-08-22
cleaning_dates2 <- as.data.frame(table(clean_me$date,clean_me$description,clean_me$station))

clean_me$Year <- format(as.Date(clean_me$date),"%Y")
clean_me$Month <- format(as.Date(clean_me$date),"%m")
clean_me$Day <- format(as.Date(clean_me$date),"%d")

table(clean_me$Year)
table(clean_me$Year,clean_me$description) # see the distribution of rows (data) by year and description
# 2015-2017 data had the 'Temperature' label and switched to 'ambient_deg_c' some time in 2017
# rename 'ambient_deg_c' = 'Temperature' 
# 'ambient_temperature_deg_c' occurs in 2019 and 2020 in addition to 'ambient_deg_c'. 
# Identical row numbers for both in 2020, indicating the data is a duplicate. 
# for 2020 - 'ambient_temperature_deg_c' could be average -> compare with 'ambient_deg_c' temp
# in general there is more data in 2018 on-wards and this could be due to what was extracted from the devices, 
# or the actual sample rate has increased due to a change in VEMCO programming. 
# This script currently uses the perviously merged 2016-2020 data from OTN
# future scripts could extract receiver event data for individual years and then merge.

str(clean_me)
head(clean_me)
tail(clean_me)

## change temperature description so they are consistent (and then I can extract just temperature data, exclude min/max/mean)
table(clean_me$description) 
library(stringr)
clean_me$description <- str_replace(clean_me$description,"ambient_mean_deg_c","Average temperature")
clean_me$description <- str_replace(clean_me$description,"ambient_max_deg_c","Maximum temperature")
clean_me$description <- str_replace(clean_me$description,"ambient_min_deg_c","Minimum temperature")
#clean_me$description <- str_replace(clean_me$description,"ambient_temperature_deg_c","Temperature_2") # could just leave this as it is as I am unsure if this is raw temperature or mean etc. Same number of rows as 'Temperature' in 2020.
clean_me$description <- str_replace(clean_me$description,"ambient_deg_c","Temperature")

table(clean_me$Year,clean_me$description)

# filter out the 'Temperature' data and leave min/max/average. 
# This should be daily temperatures on each receiver

## Extract only data for description = "Temperature"
clean_me2 <- clean_me %>%
  filter(description == "Temperature")
table(clean_me2$description)

head(clean_me2)
min(clean_me2$date) # check the date range matches the original data set date range
max(clean_me2$date)
summary(clean_me2$date)

SABMPA_temp_data_20152020_clean <- clean_me2

### --------------------------------------------
# B. 2020-2021 SABMPA temperature data (North and South crab lines)

clean_me2020 <- T2020b2021

clean_me2020 <- clean_me2020[ ,c("station_name","date","receiver",
                         "description","data","rcv_serial_no","deploy_date",
                         "recover_date","dep_lat","dep_long","the_geom","catalognumber")]

setnames(clean_me2020, c("station_name","data", "deploy_date","recover_date"),
         c("station","temperature_c","deploy_date_time","recover_date_time"))

str(clean_me2020)
head(clean_me2020)

## change station names so they are consistent
table(clean_me2020$station) 
clean_me2020$station <- clean_me2020$station %>%
  str_replace_all(c("101"="1","102"="2","103"="3","104"="4","105"="5","106"="6","107"="7","108"="8","109"="9","110"="10",
                    "111"="11","112"="12","113"="13","114"="14","115"="15","116"="16","117"="17","118"="18","119"="19",
                    "120"="20","121"="21","122"="22","123"="23"))
table(clean_me2020$station)

## Extract only data for description = "Temperature"
library(lubridate)
head(clean_me2020)
table(clean_me2020$description)


## change temperature description so they are consistent (and then I can extract just temperature data, exclude min/max/mean)
table(clean_me2020$description) # only ambient_deg_c, so I can just change that to Temperature to match the changes for the others. 
library(stringr)
clean_me2020$description <- str_replace(clean_me2020$description,"ambient_deg_c","Temperature")
table(clean_me2020$description)

# adjust the date time class
clean_me2020 <- as.data.table(clean_me2020 %>% mutate(deploy_date_time = as.character(gsub("T"," ", deploy_date_time)))) # remove the 'T' if present for OTN format datetime
clean_me2020 <- as.data.table(clean_me2020 %>% mutate(recover_date_time = as.character(gsub("T"," ", recover_date_time))))

clean_me2020$deploy_date_time <- as.POSIXct(clean_me2020$deploy_date_time, tz = "UTC") # convert datetime to POSIXct UTC tz
clean_me2020$recover_date_time <- as.POSIXct(clean_me2020$recover_date_time, tz = "UTC") 

head(clean_me2020)
table(clean_me2020$deploy_date_time)
table(clean_me2020$recover_date_time)

#make 20202021 data 'Date' format (not 'IDate to match other data sheets)
clean_me2020$date <- as.Date(clean_me2020$date,format="%Y-%m-%d")

# # check there are consistent number of data points for each day. Just first and last days should have less (day deployed/recovered).
min(clean_me2020$date)
max(clean_me2020$date)
table(clean_me2020$date) 
class(clean_me2020$date)

# remove data for each day deployed/recovered to remove the peaks at that time. Will result in 1 less day of data per year, consistent for all stations.
clean_me2020 <- clean_me2020 %>% 
  filter(!date %in% as.Date('2020-08-22') & !date %in% as.Date('2021-08-10'))
table(clean_me2020$date)

# add individual columns for Y/M/D to SABMPA_temp_data_20202021_clean to match 20162020 data sheet
clean_me2020$Year <- format(as.Date(clean_me2020$date),"%Y")
clean_me2020$Month <- format(as.Date(clean_me2020$date),"%m")
clean_me2020$Day <- format(as.Date(clean_me2020$date),"%d")

# check the temperature data makes sense
summary(clean_me2020$temperature_c)

SABMPA_temp_data_20202021_clean <- clean_me2020

### --------------------------------------------
# C. 2021-2022 SABMPA temperature data (OTN line)

clean_me_otn <- T20212022
head(clean_me_otn)
str(clean_me_otn)

# SABMPA_temp_20162021_clean_v4 <- SABMPA_temp_20162021_clean_v4 %>% ## remove / delete / drop a column
#  mutate(Date_deploy = NULL)

### --- add 'array' column to the OTN line filled with 'OTN2021' to identify this data from crab N and S lines
clean_me_otn <- clean_me_otn %>% # selecting only sabmpa project data (previous data set also had cabot straight (cbs) and halifax (hfx) data in there
  filter(project == "sabmpa")
table(clean_me_otn$project)

clean_me_otn <- clean_me_otn %>% mutate(array = "OTN2021") # add new column called 'array' and fill with 'OTN2021'
head(clean_me_otn)

clean_me_otn <- clean_me_otn[ ,c("project","station_name","date","receiver",
                                 "description","data","rcv_serial_no","deploy_date_time",
                                 "recover_date_time","dep_lat","dep_long","the_geom","catalognumber","array")]

setnames(clean_me_otn, c("station_name","data"),
         c("station","temperature_c"))


## change station names so they are consistent
table(clean_me_otn$station) 
clean_me_otn$station <- clean_me_otn$station %>%
  str_replace_all(c("001"="1","002"="2","003"="3","004"="4","005"="5","006"="6","007"="7","008"="8","009"="9","010"="10",
                    "011"="11","012"="12","013"="13","014"="14","015"="15","016"="16","017"="17","018"="18","019"="19",
                    "020"="20","021"="21","022"="22","023"="23","024"="24","025"="25","026"="26","027"="27","028"="28",
                    "029"="29","030"="30","031"="31","032"="32","033"="33","034"="34","035"="35","036"="36","037"="37",
                    "038"="38","039"="39","040"="40","041"="41","042"="42","043"="43","044"="44","045"="45","046"="46"))

clean_me_otn$station <- str_replace(clean_me_otn$station,"221","Y21") #previous code messed up the year '2021' in the stations, just change to Y21
table(clean_me_otn$station)

# # check there are consistent number of data points for each day. Just first and last days should have less (day deployed/recovered).
min(clean_me_otn$date)
max(clean_me_otn$date)
table(clean_me_otn$date) 
class(clean_me_otn$date)

# remove data for each day deployed/recovered, will result in 1 less day of data per year consistent for all stations.
#clean_me_otn2 <- clean_me_otn %>% 
#  filter(between(date, as.Date('2021-08-11'), as.Date('2022-08-05'))) # both of these methods work

clean_me_otn <- clean_me_otn %>% 
  filter(!date %in% as.Date('2021-10-10') & !date %in% as.Date('2022-08-06'))

# adjust the date time class
clean_me_otn <- as.data.table(clean_me_otn %>% mutate(deploy_date_time = as.character(gsub("T"," ", deploy_date_time)))) # remove the 'T' if present for OTN format datetime
clean_me_otn <- as.data.table(clean_me_otn %>% mutate(recover_date_time = as.character(gsub("T"," ", recover_date_time))))

clean_me_otn$deploy_date_time <- as.POSIXct(clean_me_otn$deploy_date_time, tz = "UTC") # convert datetime to POSIXct UTC tz
clean_me_otn$recover_date_time <- as.POSIXct(clean_me_otn$recover_date_time, tz = "UTC") 

head(clean_me_otn)
table(clean_me_otn$deploy_date_time)
table(clean_me_otn$recover_date_time)

# adding new columns for deploy and recover dates (no times)
clean_me_otn <- clean_me_otn %>% 
  mutate(date_deploy = as.Date(deploy_date_time)) %>%
  mutate(date_recover = as.Date(recover_date_time))
head(clean_me_otn)
summary(clean_me_otn)

# add individual columns for Y/M/D
clean_me_otn$Year <- format(as.Date(clean_me_otn$date),"%Y")
clean_me_otn$Month <- format(as.Date(clean_me_otn$date),"%m")
clean_me_otn$Day <- format(as.Date(clean_me_otn$date),"%d")

## change temperature description so they are consistent
table(clean_me_otn$description) # only Temperature, no need to change anything

head(clean_me_otn)
summary(clean_me_otn$temperature_c) # the max is high, but I have removed deploy recover days. Might be best to plot to see what's happening

SABMPA_temp_data_20212022_clean <- clean_me_otn

### --- NOTES ON OTN LINE TEMPERATURE DATA

head(clean_me_otn)
min(clean_me_otn$deploy_date_time)
max(clean_me_otn$deploy_date_time)
min(clean_me_otn$recover_date_time)
max(clean_me_otn$recover_date_time)
min(clean_me_otn$date)
max(clean_me_otn$date)

# not sure why but the data doesn't start till October despite deployment listed as August.

#######################################

### --------------------------------------------
# D. Merging the crab lines data (2015-2020, 2020-2021) and further cleaning to merge with OTN line data (2021-2022)

## this data is ready to begin merge. 
# Still needs some cleaning but easier to merge the two crab line data sets and clean together from here.
#SABMPA_temp_data_20152020_clean
#SABMPA_temp_data_20202021_clean
#SABMPA_temp_data_20212022_clean

### --- check min max dates for each dataset prior to merge 
# make sure no overlap and see how big the gaps are - make sure to use updated version names here
min(SABMPA_temp_data_20152020_clean$date)
max(SABMPA_temp_data_20152020_clean$date)
min(SABMPA_temp_data_20202021_clean$date)
max(SABMPA_temp_data_20202021_clean$date)
min(SABMPA_temp_data_20212022_clean$date)
max(SABMPA_temp_data_20212022_clean$date)

#compare column headers before merge
library(janitor)
compare_df_cols(SABMPA_temp_data_20152020_clean,SABMPA_temp_data_20202021_clean)
#simple bind rows together
SABMPA_temp_20162021_clean <- rbind(SABMPA_temp_data_20152020_clean,SABMPA_temp_data_20202021_clean) 

# check date range as an indication that merge worked
min(SABMPA_temp_20162021_clean$date) 
max(SABMPA_temp_20162021_clean$date)
head(SABMPA_temp_20162021_clean)
tail(SABMPA_temp_20162021_clean)

summary(SABMPA_temp_20162021_clean$temperature_c)
# we have not removed deploy/recover days yet, so still have some extreme temp values

### --- further data clean to merge with OTN line data

# add column for 'project' and fill with SAB to match the 20212022 OTN line data 
# can skip this step if working with only sabmpa data from 20212022, i.e., without bfx and cbs data
# but project column is useful once I bring in temp data from other MPAs
SABMPA_temp_20162021_clean_v2 <- SABMPA_temp_20162021_clean %>% mutate(project = "sabmpa")
head(SABMPA_temp_20162021_clean_v2)
SABMPA_temp_20162021_clean_v2[,c(13,1,2,3,4,5,6,7,8,9,10,11,12)] #rearrange to place 'project' in the first column

### --- adding array column - crab north and south lines before merge with OTN line data
SABMPA_temp_20162021_clean_v2 <- SABMPA_temp_20162021_clean_v2 %>% mutate(array = station) # copy station over to new column 'array'

# remove all the numeric values from 'array'
library(tidyverse)
SABMPA_temp_20162021_clean_v2 <- SABMPA_temp_20162021_clean_v2 %>% 
  mutate(array = trimws(str_remove(array, "(\\s+[A-Za-z]+)?[0-9-]+")))

# rename in array so N = North, S = South (array lines)
SABMPA_temp_20162021_clean_v2 <- SABMPA_temp_20162021_clean_v2 %>% mutate(array = as.character(gsub("N","North", array)))
SABMPA_temp_20162021_clean_v2 <- SABMPA_temp_20162021_clean_v2 %>% mutate(array = as.character(gsub("S","South", array)))

## NEED TO REMOVE THE EXTREME DATA THAT IS DUE TO RETRIEVAL AND DEPLOYMENT THROUGHOUT THE DATASET
summary(SABMPA_temp_20162021_clean_v2$temperature_c)
summary(SABMPA_temp_20162021_clean_v2)
# adding new columns for deploy and recover dates (no times)
SABMPA_temp_20162021_clean_v2 <- SABMPA_temp_20162021_clean_v2 %>% 
  mutate(date_deploy = as.Date(deploy_date_time)) %>%
  mutate(date_recover = as.Date(recover_date_time))
head(SABMPA_temp_20162021_clean_v2)

table(SABMPA_temp_20162021_clean_v2$date_deploy) # listing all the deploy dates
table(SABMPA_temp_20162021_clean_v2$date_recover) # listing all the recover dates

SABMPA_temp_20162021_cleanfilt <- SABMPA_temp_20162021_clean_v2 %>% # removing rows for data on days deployed
  filter(!date %in% as.Date('2015-05-11') & !date %in% as.Date('2016-04-17')
         & !date %in% as.Date('2017-05-05') & !date %in% as.Date('2018-06-10')
         & !date %in% as.Date('2019-08-07') & !date %in% as.Date('2020-08-22'))

SABMPA_temp_20162021_cleanfilt <- SABMPA_temp_20162021_cleanfilt %>% # removing rows for data on days recovered
  filter(!date %in% as.Date('2016-04-17') & !date %in% as.Date('2017-05-05')
         & !date %in% as.Date('2018-06-10') & !date %in% as.Date('2019-08-07')
         & !date %in% as.Date('2020-08-22') & !date %in% as.Date('2021-08-10'))

nrow(SABMPA_temp_20162021_clean_v2)-nrow(SABMPA_temp_20162021_cleanfilt) # checking to see an appropriate amount of rows have been removed

head(SABMPA_temp_20162021_cleanfilt)
summary(SABMPA_temp_20162021_cleanfilt)

### --- MERGING OTN LINE 2021-2022 DATA WITH CRAB LINES 2016-2021
library(janitor)
# compare column headers before merge
compare_df_cols(SABMPA_temp_20162021_cleanfilt,SABMPA_temp_data_20212022_clean) 

# simple bind rows together - note rbind should account for columns in different order, but check
SABMPA_temp_20162022_clean <- rbind(SABMPA_temp_20162021_cleanfilt,SABMPA_temp_data_20212022_clean) 

 # make sure I am merging the right versions of each data sheet, i.e., after filtering etc.
head(SABMPA_temp_20162022_clean)
summary(SABMPA_temp_20162022_clean)
summary(SABMPA_temp_20162022_clean$temperature_c)

table(SABMPA_temp_20162022_clean$project)
table(SABMPA_temp_20162022_clean$array)
table(SABMPA_temp_20162022_clean$description)
table(SABMPA_temp_20162022_clean$station)
      
##############################
# plotting
##############################

# A. plotting the merged N and S crab lines data (2016-2020 + 2020-2021)

p1 <- ggplot(data=SABMPA_temp_20162021_cleanfilt) +
  geom_smooth(aes(x=date, y=temperature_c, 
                  ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))),
                  group=as.factor(array),color=as.factor(array))) + # make it plot daily mean temperature with CIs, Colour by project with key.
  labs(x="", y=expression(paste("Temperature ",degree,"C")),color='Array') +
  scale_x_date(date_labels = "%b - %Y",date_breaks = "3 month") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        legend.position = "bottom")+
  guides(color=guide_legend(override.aes=list(fill=NA)))

### --------------------------------------------
# B. plotting the OTN line data

p2 <- ggplot(SABMPA_temp_data_20212022_clean) +
#  geom_point(size=1.3,aes(x=date, y=temperature_c)) +
  geom_smooth(aes(x=date, y=temperature_c, 
                  ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))),
                  group=as.factor(array),color=as.factor(array))) + # make it plot daily mean temperature with CIs, Colour by project with key.
  labs(x="", y=expression(paste("Temperature ",degree,"C")), color="Array") +
  scale_x_date(date_labels = "%b - %Y",date_breaks = "1 month") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        legend.position = "bottom")+
  guides(color=guide_legend(override.aes=list(fill=NA)))

### --------------------------------------------
# C. plotting the full merged data N and S crab lines (2016-2020 + 2020-2021) and OTN line (2021-2022)

p3 <- ggplot(SABMPA_temp_20162022_clean) +
  geom_smooth(aes(x=date,y=temperature_c,
                  ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))),
                  group=as.factor(array),color=as.factor(array))) + # make it plot daily mean temperature with CIs, Colour by project with key.
  labs(x="", y=expression(paste("Temperature ",degree,"C")),color='Array') +
  scale_x_date(date_labels = "%b - %Y",date_breaks = "3 month") +
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        legend.position = "bottom")+
  guides(color=guide_legend(override.aes=list(fill=NA)))

# getting monthly averages with sd for plotting

MonthSum <- SABMPA_temp_20162022_clean %>%
  group_by(Year, Month, array) %>%
  summarise(mn = mean(temperature_c, na.rm=T),
            sd = sd(temperature_c, na.rm=T)) %>%
  ungroup() %>%
  mutate(datetime=dmy(paste("01", Month, Year, sep="-"))) %>%
  data.frame()

p4 <- ggplot(MonthSum, aes(x=datetime, y=mn, col=array, group=array))+
  #geom_errorbar(aes(ymin=mn-sd,ymax=mn+sd))+
  geom_line() +
  geom_point(size=1.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom") +
  scale_x_date(date_labels = "%b - %y",date_breaks = "3 month") +
  labs(x="",
       y=expression(paste("Temperature ",degree,"C")),
       col="Array")

p5 <- ggplot(MonthSum, aes(x=datetime, y=mn, col=array, group=array))+
  geom_errorbar(aes(ymin=mn-sd,ymax=mn+sd))+
  geom_line() +
  geom_point(size=1.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom") +
  scale_x_date(date_labels = "%b - %y",date_breaks = "3 month") +
  labs(x="",
       y=expression(paste("Temperature ",degree,"C")),
       col="")

p6 <- ggplot(SABMPA_temp_20162022_clean, aes(x=Month, y=temperature_c, group=Month))+
  geom_boxplot() +
  stat_summary(fun = "mean",geom = "point", shape = 21, size = 2,col="blue")+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x="Month",
       y=expression(paste("Temperature ",degree,"C")),
       col="")

#saving the figures
ggsave("output/Temperature/SAB_temp_20152021_gam.png",p1,dpi=200,width=8,units="in")
ggsave("output/Temperature/SAB_temp_20212022_gam.png",p2,dpi=200,width=8,units="in")
ggsave("output/Temperature/SAB_temp_20152022_gam.png",p3,dpi=200,width=8,units="in")
ggsave("output/Temperature/SAB_temp_20152022_monthmean.png",p4,dpi=200,width=8,units="in")
ggsave("output/Temperature/SAB_temp_20152022_monthmeansd.png",p5,dpi=200,width=8,units="in")
ggsave("output/Temperature/SAB_temp_20152022_monthmean_boxp.png",p6,dpi=200,width=8,units="in")

# if you want to save the clean temperature data
#write.csv(SABMPA_temp_data_20152020_clean,"data/Temperature/SABMPA_temp_data_20152020_clean.csv")
#write.csv(SABMPA_temp_data_20202021_clean,"data/Temperature/SABMPA_temp_data_20202021_clean.csv")
#write.csv(clean_me_otn,"data/Temperature/SABMPA_temp_data_20212022_clean.csv")