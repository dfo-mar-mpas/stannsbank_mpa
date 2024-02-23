#load libraries ----------

library(dplyr)
library(tidyr)
library(tibble)
library(taxize)
library(RCurl)
library(worrms)

#load functions -----------

#this is from the eDNA analysis from the offshore. 
script <- RCurl::getURL("https://raw.githubusercontent.com/rystanley/offshore_edna/main/code/ClassifyFunction.R", ssl.verifypeer = FALSE)
eval(parse(text = script),envir=.GlobalEnv)
rm(script)  

#load the snow crab species data -----
tax_df <- read.csv("output/taxonomic_raw.csv")%>%
          filter(!species_filter %in% c("SEAWEED, ALGAE ,KELP; THALLOPHYTA"),
                 !grepl("eggs",tolower(species_filter)))%>%
          mutate(species = ifelse(species_filter == "LEPTASTERIAS (HEXASTERIAS) POLARIS","LEPTASTERIAS POLARIS",species_filter),
                 species = ifelse(species_filter == "PORANIOMORPHA (PORANIOMORPHA) HISPIDA", "PORANIOMORPHA HISPIDA",species_filter))

PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#run the taxonomic classification ** note that there will be time that require input (to select amongst options. By default go with the 'accepted' and closest match)

tax_list <- data.frame()

for(i in 1:nrow(tax_df)){
  
  sp=tax_df[i,"species"]%>%tolower()
  
  message(paste0("Working on ",sp," ",i," of ",nrow(tax_df)))
  
  temp <- classification(sp,db="worms") #run classification
  
  temp2 <- temp[[1]]%>% #unpack classification
          data.frame()%>%
          select(rank,name)%>%
          spread(rank,name)%>%
          mutate(aphiaID = temp[[1]]%>%data.frame()%>%slice(n())%>%pull(id))
  
  temp3 <- temp2%>% #trim classification
           select(all_of(c(names(temp2)[names(temp2)%in%PhyloNames],"aphiaID")))
  
  #if data is missing or a certain taxonomic level isn't identified.
  missing_cols <- setdiff(c(PhyloNames,"aphiaID"),names(temp3))
  
  if(length(missing_cols)>0){
    
    temp3[,missing_cols] <- NA
    
    temp3 <- temp3[,c(PhyloNames,"aphiaID")]
    
  }
  
  temp3$species_filter <- tax_df[i,"species"] # this is for linking back in the original code
  
  tax_list <- rbind(tax_list,temp3)
  
}

#Now get the 'common' names for each

tax_list$common <- NA

for(i in 1:nrow(tax_list)){
  
  sp=tax_list[i,"Species"]%>%tolower()
  
  message(paste0("Working on ",sp," ",i," of ",nrow(tax_list)))

  if(!is.na(tax_list[i,"Species"])){
    
    comm_name <- sci2comm(tax_list[i,"Species"])%>%as.character()
    
    if(comm_name == "character(0)"){comm_name = NA}
    
    tax_list[i,"common"] <- comm_name
  
  }
}

#save the output 
write.csv(tax_list,"output/crab_taxa_clean.csv",row.names = FALSE)
