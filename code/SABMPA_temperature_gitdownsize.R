#SABMPA_temperature_gitdownsize
# harri.pettitt-wade@dfo-mpo.gc.ca
# SABMPA acoustic telemetry data summaries
# temperature data from VR2AR receivers
################################################################################

## Temperature files are too large to upload to github. 
# This script is to downsize the files to be uploaded to github and used in the script 'SABMPA_temperature_v1.R'.

# Warning: due to the file sizes, this script links to data on a local drive. 
# The data is also available on the OTN data portal, so users could download the open data and change the directory in this script to the file download.


# to downsize: 
# 1.'V2LSCF_temp_data_2016_to_2020_notclean.csv'
# 2.'v2lscf_temps_2020-2021_notclean.csv'


rm(list = ls())
graphics.off()

library(dplyr)

setwd(dir = "C:/Users/PETTITTWADEH/Documents/R/projects/SABMPA_acoustic/data/Temperature")
getwd()

# SAB crab lines data 2016-2020 - split by year
TA <- read_csv("V2LSCF_temp_data_2016_to_2020_notclean.csv", 
                                          col_types = cols(date = col_date(format = "%Y-%m-%d"), 
                                                           data = col_number(), units = col_skip(), 
                                                           recover_ind = col_skip()))

splitdata <- split(TA, format(as.Date(TA$date), "T%Y"))
list2env(splitdata ,.GlobalEnv)

write.csv(T2015,"temp_data_2015.csv")
write.csv(T2016,"temp_data_2016.csv")
write.csv(T2017,"temp_data_2017.csv")
write.csv(T2018,"temp_data_2018.csv")
write.csv(T2019,"temp_data_2019.csv")
write.csv(T2020,"temp_data_2020a.csv")

# SAB crab lines data 2020-2021 - split by year
rm(list = ls())

setwd(dir = "C:/Users/PETTITTWADEH/Documents/R/projects/SABMPA_acoustic/data/Temperature")
getwd()

TB <- fread("v2lscf_temps_2020-2021_notclean.csv", 
            sep = ",", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

splitdata <- split(TB, format(as.Date(TB$date), "T%Y"))
list2env(splitdata ,.GlobalEnv)

write.csv(T2020,"temp_data_2020b.csv")
write.csv(T2021,"temp_data_2021.csv")
