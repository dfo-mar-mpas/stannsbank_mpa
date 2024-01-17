# README
## Description of data files in data/CrabSurvey

The data in this 2023 folder is for the 2024 SAB CSAS data review. 

1. SABMPA2023export.csv is from Amy Glass and contains all catch data for all SAB stations (enhanced and not) from 2015 to 2023. 
There is no real data from 2020 due to covid. The 2019 survey did extend into January 2020, so change the "2020" year to 2019. 

2. sabmpa_area_swept.csv is a small file of the surface area of sea floor swept per survey station and can be linked to other data by the TRIP and STATION columns

3. FishMorph.RData is from Brent Cameron and has all the fish and crab morphology data within it since 2004. This needs to be filtered by year and station to just be the SAB stations. 
   Nick saved this as "fishmorph2" in the RData analysis file
   
4. 2023CrabSurveyDat.RData is Nick's analysis objects from various R scripts also in the Github. It also has some shapefiles and analyses saved within it.

5. OldData/ folder has data from Amy and Brent that was missing some stations, so we don't really need this
