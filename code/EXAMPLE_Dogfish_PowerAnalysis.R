##code analyzing power analysis of Laurentian Channel trawl data
DFO.LC.CPUE<-function(){
  options(scipen=999)
  library(ggplot2)
  library(gridExtra)
  library(MASS)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(MASS)
  #install.packages("magrittr") 
  library(magrittr)
  #install.packages("dplyr") 
  #remove.packages("rlang")
  #install.packages("rlang")
  library(rlang)
  library(dplyr)
  library(imputeTS) ## replace NA to 0
  
   ##toggle whether to run power analysis - Y or N input
  power.analysis<-"Y"
  
  ##load up trawl data
  trawl<-read.table("Abundance and Biomass Data_MPA.csv", sep=",", header=T)
  ##change blanks to zeros
  trawl <- trawl %>% replace(is.na(.), 0)
  #trawl=na.replace(trawl, 0)
  head(trawl)
  ##pick out conservation objective species
  trawl<-trawl[,c(1:19)]
  
  # these code for species richness  
  # trawl<-read.table("Abundance and Biomass Data_MPA.csv", sep=",", header=T)
  
  ##split data to get row sums
  # trawl.info<-trawl[,2:16]
  #  trawl.data<-trawl[,27:229]
  #  trawl.PA<-trawl.data %>% mutate(replace(., .>0, 1))
  #  trawl.info$richness<-rowSums(trawl.PA, na.rm=T)
  
  
  with(trawl, plot(black_dogfish~surveyyear))
  plot1a<-ggplot(data=trawl, aes(surveyyear,black_dogfish, colour=strata))+geom_point() + geom_smooth() + ylab("Black Dogfish per trawl") + xlab("Survey year") 
  
  ##add artificial period category
  trawl$period<-rep(1, nrow(trawl))
  trawl[trawl$Yr<14,"period"]<-0
  trawl[trawl$Yr>19,"period"]<-2
  ##do mixed effects linear model using lme4 package glmer with Area, depth, season and unique.set(random) as model terms
  library(lme4)
  # model<-glmer(black_dogfish~as.character(period)+DEP_MEAN+(1|strata), data=trawl, family="poisson")
  # summary(model)
  
  ##first do the before-during comparison with full model
  ##model.nb<-glm.nb(black_dogfish~as.character(period)+DEP_MEAN, data=trawl)
  trawl$strata<-as.factor(trawl$strata)
  model.nb<-glmer.nb(black_dogfish~period+(1|strata), data=trawl)
  model.nb.reduced<-glmer.nb(black_dogfish~1+(1|strata), data=trawl)
  summary(model.nb)
  
  ##now compare the importance of the main effect removed using the proper likelihood ratio test
  LR.result<-anova(model.nb,model.nb.reduced, test="Chisq")
  LR.result
  ##note dispersion parameter with strata accounted for changes from 0.287 to 0.3294
  ##strata level variance listed at 5.046
  ##power analysis - need to extract dispersion parameter Family: Negative Binomial(0.3495)  ( log )
  neg.bin.dispersion<-0.3495
  
  ##extract strata-specific means for latest period
  strata.means<-trawl %>%
    filter(period==2) %>%
    group_by(strata) %>%
    summarize(mean.catch.strata=mean(black_dogfish))
  strata.means$strata<-as.character(strata.means$strata)
  
  ##toggle whether mixed effects model (strata as random effect) is used or just regular glm.nb.  Note data generated will always be based on stratified data (T/F). 
  mixed.effect<-"F"
  
  ##restrict analysis to core strata if toggled (T or F)
  core.areas<-"T"
  ##adjust core threshold as percentage of mean abundance considered to be core strata.  For example, if core.threshold is set to 0.5 and the mean abundance across all thresholds is 100, then only strata with a mean baseline abundance of 50 or greater will be included in the simulations.
  core.threshold<-0.25
  
  ##list strata areas (inside MPA) - which determines sampling effort
  strata<-as.character(c(305, 313, 705, 706, 711, 712, 713, 714, 715))
  # strata<-as.character(c(305,  711, 712, 713, 714)) # have to remove stra 313, 705,706, 715 in raw data (file)
  
  MPA.strata.areas<-c(545.89,332.14,665.08,1144.59,1502.87,2050.76,1909.24,3086.6,98.52)
  #MPA.strata.areas<-c(545.89,1502.87,2050.76,1909.24,3086.6)# have to remove stra 313, 705,706, 715 in raw data (file)
  
  
  strata.areas<-as.data.frame(cbind(strata, MPA.strata.areas))
  strata.areas$MPA.strata.areas<-as.numeric(as.character(strata.areas$MPA.strata.areas))
  ##determine proportionate area so we can allocate samples accordingly
  strata.areas$strata.prop.area<-strata.areas$MPA.strata.areas/sum(strata.areas$MPA.strata.areas)
  ##need to merge strata.areas with strata.means
  strata.list<-full_join(strata.areas, strata.means)
  
  ##define number of iterations
  ##set.seed(5)
  if(power.analysis=="Y"){
    iterations<-1000  ##kept low until machinery is working (take a day for each species)
    # it runs faster with lower iteration
    ##define biomass adjustment 
    biomass.adjustment<-1
    ##define effect size
    #effect.size<-c(0.1,0.2, 0.4,0.6)### org
    effect.size<-c(0.1,0.2, 0.4,0.6,0.8)### revised paper
    ##define number of years to run dataset over
    
    ##define trawl numbers to run power analysis on
    trawl.sample.size<-c(15,23,50,100) # number of tows
    
    ##clear the old simulation file if it exists
    if(exists("power.raw")) rm("power.raw", envir = globalenv())
    
    ##for every defined sample size...
    for(a in 1:length(trawl.sample.size)){
      ##for every defined effect size...
      for(b in 1:length(effect.size)){
        ##conduct the following functions according to the fixed number of iterations
        for(aa in 1:iterations){
          print(c(a,b,aa))
          treatment.data<-baseline.data<-NULL
          
          
          temp.rand.trawl.samp<-as.data.frame(table(sample(c(strata.list$strata), size = trawl.sample.size[a], replace = TRUE, prob =strata.list$strata.prop.area)))
          names(temp.rand.trawl.samp)<-c("strata","trawls")
          
          ##use those intercepts and use neg binomial to generate black dogfish abundance in trawl sets using the dispersion parameter
          
          ##determine baseline mean on to which we can add intercept adjustments - note these mean values are not on log scale now
          ##note biomass adjustment allows us to model lower levels of animals. default is no adjustment
          baseline.means<-strata.list$mean.catch.strata*biomass.adjustment
          treatment.mean<-baseline.means*(effect.size[b])
          
          ##generating black dogfish abundance specific to each strata.
          for(bb in 1:nrow(temp.rand.trawl.samp)){
            ##clear out temp data
            temp.strata.treatment<-temp.strata.baseline<-temp.strata<-temp.sample.size<-temp.effect.size<-temp.treatment.class<-temp.strata.baseline.data<-temp.strata.treatment.class<-temp.treatment.class<-temp.strata.treatment.data<-temp.strata.data<-NULL
            
            temp.strata.treatment<-rnbinom(n=temp.rand.trawl.samp$trawls[bb], size=neg.bin.dispersion, mu=exp(log(treatment.mean[bb])))
            
            temp.strata.baseline<-rnbinom(n=temp.rand.trawl.samp$trawls[bb], size=neg.bin.dispersion, mu=exp(log(baseline.means[bb])))
            
            
            
            ##add metadata and join with simulated catch numbers for both baseline and impact sites
            temp.strata<-rep(as.character(temp.rand.trawl.samp$strata[bb]), temp.rand.trawl.samp$trawls[bb])
            temp.sample.size<-rep(trawl.sample.size[a],temp.rand.trawl.samp$trawls[bb])
            temp.effect.size<-rep(effect.size[b],temp.rand.trawl.samp$trawls[bb])
            temp.strata.baseline.class<-rep("Baseline", length(temp.strata.baseline))
            temp.strata.baseline.data<-as.data.frame(cbind(temp.strata, temp.strata.baseline,temp.sample.size,temp.effect.size, temp.strata.baseline.class))
            temp.strata.treatment.class<-rep("After", length(temp.strata.treatment))
            temp.strata.treatment.data<-as.data.frame(cbind(temp.strata,temp.strata.treatment,temp.sample.size,temp.effect.size, temp.strata.treatment.class))
            
            ##add the treatment column and join both files
            names(temp.strata.treatment.data)<-c("strata","Abundance",'sample.size','effect.size','Period')
            names(temp.strata.baseline.data)<-c("strata","Abundance",'sample.size','effect.size','Period')
            temp.strata.data<-rbind(temp.strata.treatment.data, temp.strata.baseline.data)
            ifelse(bb==1, power.raw<-temp.strata.data, power.raw<-rbind(power.raw, temp.strata.data))
          }
          
          power.raw$Abundance<-as.numeric(as.character(power.raw$Abundance))
          power.raw$strata<-as.character(power.raw$strata)
          
          ##filter out strata that are not core strata (if core strata is toggled)
          
          if(core.areas=="T") {core.strata<-as.data.frame(strata.means[strata.means$mean.catch.strata>=mean(strata.means$mean.catch.strata)*core.threshold,])
          power.raw<-power.raw[power.raw$strata %in% core.strata$strata,]
          }
          
          ##do stats
          
          tryCatch({
            
            ##based on whether mixed.effects is T or F, we will use a negative binomial error structure within a mixed effects model or just a regular GLM
            
            ##ifelse(mixed.effect=="T", {  
            ##temp.model<-glmer.nb(Abundance~Period + (1|strata), data=power.raw)
            ##temp.model.reduced<-glmer.nb(Abundance~1+ (1|strata), data=power.raw)},
            ##{temp.model<-glm.nb(Abundance~Period, data=power.raw)
            ##temp.model.reduced<-glm.nb(Abundance~1, data=power.raw)
            ##})
            
            ##since ifelse statement causing errors, doing this with two if statements.  First for mixed model
            if(mixed.effect=="T") {  
              temp.model<-glmer.nb(Abundance~Period + (1|strata), data=power.raw)
              temp.model.reduced<-glmer.nb(Abundance~1+ (1|strata), data=power.raw)}
            ##then for non-mixed model option
            if(mixed.effect=="F") {  
              temp.model<-glm.nb(Abundance~Period, data=power.raw)
              temp.model.reduced<-glm.nb(Abundance~1, data=power.raw)}
            
            ##use maximum likelihood to get most accurate p value
            ML.result<-anova(temp.model,temp.model.reduced, test="Chisq")
            summary(temp.model)
            ##save p values for each scenario
            temp.p.value<-ML.result[[8]][[2]]
            ##add metadata and then join with all other effect size data sets and sample size data sets
            temp.scenario<-c(trawl.sample.size[a],effect.size[b],temp.p.value)
            ifelse(a==1 & b==1 & aa==1, scenario.values<-temp.scenario, scenario.values<-rbind(scenario.values,temp.scenario))
          }##end of if loop related to tryCatch
          , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }##end of aa loop - one run for each iteration
        
        
      }##end of b loop - one run for each effect size
    }##end of a loop - one run for each string sample size
    
    
    ##change scenario values to data frame and name columns
    scenario.values<-as.data.frame(scenario.values, row.names=F)
    names(scenario.values)<-c('sample.size','effect.size','p.value')
    print(scenario.values)
    
    
    ##now determine the power for each one
    prop.sign<-function(x){length(x[x<0.05])/length(x)}
    summary.stats<-aggregate(p.value~sample.size+effect.size, FUN=prop.sign, data=scenario.values)
    names(summary.stats)<-c("sample.size",'effect.size','prop.sign')
    write.table(summary.stats, "LC summary power strata.csv", row.names=F,sep=',')
  }
  ##do a plot of effect size, trawl number and power
  ggplot(summary.stats, aes(1-effect.size,prop.sign, color=as.character(sample.size)))+
    geom_line()+geom_point(size=1.5)+
    #xlab("Black Dodfish Abundance Decline")+
    xlab("Effect Size Applied")+
    ylab("Power")+
    geom_hline(yintercept=0.8, color="red", linetype=3) + 
    scale_color_discrete(name="Trawls")+
    scale_color_manual(name = "Trawls",
                       breaks = c("15", "23", "50","100"),
                       values = c("15" = "blue", "23" = "green", "50" = "purple","100"="red") )+
    theme_bw()+
    theme(axis.text.x = element_text(face="plain", color="black",size=10),
          axis.text.y = element_text(face="plain", color="black",size=10)) + 
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10,face="plain"))+
    
    #theme(axis.title.y = element_text(color="darkgreen", size=10, face="plain"))+
    theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+ 
    theme(legend.position="right")+
    #theme(legend.title = element_blank())+
    theme(legend.text = element_text(colour="black", size=10, face="plain"),
          legend.title=element_text(size=10))+
    scale_x_continuous(limits=c(0, 1),breaks=seq(0,1,0.2))+
    scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.2))
  
}

# write out the csv file for future use once we want to replot figure (don't need to run all power ana codes)
write.csv(summary.stats,"dogfish_power_data.csv",row.names=FALSE)
