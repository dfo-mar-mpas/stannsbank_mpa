#eDNA summary analysis

#load libraries
library(worrms)
library(tidyverse)
library(taxize)

common_name <- function(x) {
  
  require(worrms)
  
  # Attempt to get common name records
  cn <- tryCatch(
    wm_common_id(x),
    error = function(e) NULL  # Return NULL if there's an error
  )
  
  # Check if there are any records
  if (is.null(cn) || nrow(cn) == 0) {
    return(NA)  # Return NA if no records
  }
  
  if(!"English" %in% cn$language){ #some species don't have any english name
    common <- cn%>%
    slice(1)%>%
      pull(vernacular) %>%
      sub("^(\\w)(\\w*)", "\\U\\1\\L\\2", ., perl = TRUE)
  }
  
  # Filter for English language and format the common name
  if("English" %in% cn$language){
    common <- cn %>%
    filter(language == "English") %>%
    pull(vernacular) %>%
    sub("^(\\w)(\\w*)", "\\U\\1\\L\\2", ., perl = TRUE)}
  
  # Return the first common name, or NA if none match the filter
  if (length(common) > 1) {
    return(common[1])
  } else {
    return(common)
  }
}

#for plotting and data collation 
PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#load the eDNA summary sheet
edna_raw <- read.csv("data/eDNA/2026_species_list.csv")

problem_sp <- data.frame(species=c("Balaenoptera musculus","Hemitripterus americanus", #bit of trial and error here. 
                                "Mesoplodon bidens","Acartia hudsonica","Pista maculata"),
                         aphiaid=c(137090,159518,
                                 137121,149751,868065))%>%
              left_join(.,edna_raw%>%
                          distinct(species,.keep_all=T)%>%
                          mutate(species = gsub("CMC0[1-4]", "", species), 
                                 species=trimws(species))%>%
                          dplyr::select(species,type))

edna_df <- edna_raw%>%
           mutate(species = gsub("CMC0[1-4]", "", species), 
                  species=trimws(species))%>%
           distinct(species,.keep_all=T)%>%
           filter(!species %in% problem_sp$species,
                  !is.na(species))%>%
           rowwise()%>%
           mutate(aphiaid=wm_name2id(species))%>%
           data.frame()%>%
           dplyr::select(names(problem_sp))%>%
           rbind(.,problem_sp)

#Get the full classification 

trad_species <- edna_df$species

tax_list_trad <- data.frame()

for(i in 1:length(trad_species)){
  
  sp <- edna_df[i,'aphiaid']
  
  message(paste0("Working on ",edna_df[i,'species']," ",i," of ",length(trad_species)))
  
  temp <- classification(sp,db="worms") #run classification
  
  temp2 <- temp[[1]]%>% #unpack classification
    data.frame()
  
  if(nrow(temp2[1])>1){
    
    temp2 <- temp2%>%
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
  }else{temp3 <- data.frame(matrix(NA, nrow = 1, ncol = length(PhyloNames)))
  names(temp3) = PhyloNames
  temp3$aphiaID = NA}
  
  temp3$species_filter <- sp # this is for linking back in the original code
  
  temp3$common <- common_name(sp)
  
  tax_list_trad <- rbind(tax_list_trad,temp3)
  
}                    


#clean up the final dataframe
edna_formatted <- tax_list_trad%>%
  dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter","common")))

          
write.csv(edna_formatted,file = "output/edna/2026_species_list_formatted.csv",row.names = FALSE)
