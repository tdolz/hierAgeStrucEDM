### Created 8/3/22 ####
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(cowplot)
library(rEDM)
library(GPEDM)
library(Metrics)
library(parallel)

### Tanya request figure--- Simulation data with prediction
### 
#this is basically a 30 points data comp
#A figure comparing R2 values from the total abundance model with R2 for the hier model aggregated to total abundance, 
#for each simulation and empirical dataset. All would have 30 time points. For the simulations, there would only be one value per replicate, 
#so you could show the distribution among replicates 

#load preylists and split into lists. 
preylist1 <-read.csv("Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist1 <- preylist1 %>% split(f=preylist1$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist2 <-read.csv("Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist2 <- preylist2 %>% split(f=preylist2$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist3 <-read.csv("Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist3 <- preylist3 %>% split(f=preylist3$index)%>%lapply(function(x) x[!names(x) %in% c("index")])

#plist <-preylist1[[1]] #for troublesshooting


GP30 <-function(plist){
  plist<-as.data.frame(plist)%>%
    pivot_longer(2:21, names_to = "age_class", values_to = "value")%>%filter(age_class !="V21")
  #N_total list
  plistNT <-plist %>% group_by(time_step)%>%summarize(NewValue=sum(value))%>%mutate(age_class="all", value=log(NewValue))
  #back to the regular list, log it
  plist <-mutate(plist, value=log(value)) #LOGGED
  plist <-plist %>%group_by(age_class)%>%arrange(age_class, time_step)#make sure it is sorted by age_class then year
  
  #FORMAT THE DATA
  #30 years of 20 ages
  preyall <-filter(plist, time_step >=300 & time_step <=340)%>%as.data.frame()
  plags <-makelags(data=preyall, y="value",pop="age_class", E=round(sqrt(30)), tau=1)
  prey <-cbind(preyall, plags)%>%as.data.frame()
  prey.train = filter(prey, time_step <= (max(prey$time_step)-10))
  prey.test = filter(prey, time_step > (max(prey$time_step)-10))
  #total abundance
  preytotal <-filter(plistNT, time_step >=300 & time_step <=340)%>%as.data.frame()
  plagsNT <-makelags(data=preytotal, y="value",pop="age_class", E=round(sqrt(30)), tau=1)
  preyNT <-cbind(preytotal, plagsNT)%>%as.data.frame()
  NT.train = filter(preyNT, time_step <= (max(preyNT$time_step)-10))
  NT.test = filter(preyNT, time_step > (max(preyNT$time_step)-10))
  
  ##FIT THE GP
  #30 years of 20 ages
  fit_all <-fitGP(data = prey.train, y = "value",time="time_step", x=colnames(plags),newdata=prey.test,pop="age_class",scaling = "local",predictmethod = "loo")
  #total abundance
  fit_NT <-fitGP(data = NT.train, y = "value", x=colnames(plagsNT),newdata=NT.test,predictmethod = "loo")
  
  #EXTRACT THE RESULTS
  #30 years of 20 ages -overall result
  fit_all_stats <-fit_all$outsampfitstats
  #30 years of 20 ages - add results to get overall r2
  fitSum <-fit_all$outsampresults%>% mutate(newpred=exp(predmean), newobs=exp(obs))%>% filter(!is.na(newpred))%>%
    group_by(timestep)%>%summarise(predmean=log(sum(newpred)), Obs=log(sum(newobs))) #grouping by time steps, adding across pops. 
  fitsumStats <-c(getR2(obs=fitSum$Obs, pred=fitSum$predmean), Metrics::rmse(actual=fitSum$Obs, predicted=fitSum$predmean))
  #total abundance
  fit_NT_stats <-fit_NT$outsampfitstats
  #so you want to compare fitSum to fitNT
  Statsout <-cbind(fit_all_stats, fitsumStats, fit_NT_stats)
  #Statsout<-as.data.frame(t(Statsout))%>%rownames_to_column(var="approach")
  Statsout
}

#real differences in rmse, not so strong differences in R2 - maybe go to RMSE???

#TEST
shortlist <-preylist1[1:2]
numCores = detectCores()
system.time(foo <-mclapply(shortlist, GP30, mc.cores=numCores))

#APPLY - preylist1
numCores = detectCores()
system.time(parafoo <-mclapply(preylist1, GP30, mc.cores=numCores))
#corral results into data frame
df <-list()
for (i in 1:length(parafoo)){
  df[[i]] <-as.data.frame(parafoo[[i]])%>% mutate(iter=i)%>%rownames_to_column()
}
preystats1 <-df%>% bind_rows()%>%as.data.frame()%>%mutate(sim="I")

#APPLY - preylist2
numCores = detectCores()
system.time(parafoo <-mclapply(preylist2, GP30, mc.cores=numCores))
#corral results into data frame
df <-list()
for (i in 1:length(parafoo)){
  df[[i]] <-as.data.frame(parafoo[[i]])%>% mutate(iter=i)%>%rownames_to_column()
}
preystats2 <-df%>% bind_rows()%>%as.data.frame()%>%mutate(sim="II")

#APPLY - preylist3
numCores = detectCores()
system.time(parafoo <-mclapply(preylist3, GP30, mc.cores=numCores))
#corral results into data frame
df <-list()
for (i in 1:length(parafoo)){
  df[[i]] <-as.data.frame(parafoo[[i]])%>% mutate(iter=i)%>%rownames_to_column()
}
preystats3 <-df%>% bind_rows()%>%as.data.frame()%>%mutate(sim="III")

#put them all together.
preystats <-bind_rows(preystats1, preystats2, preystats3)
preystats <-pivot_longer(preystats, 2:4, names_to="approach", values_to="value")%>%as.data.frame()

#write csv
write.csv(preystats, "thirtypoints.csv")


