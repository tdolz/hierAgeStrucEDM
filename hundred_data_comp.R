#### this is the 100 data points analysis ############ 
#### we split the data into 5 ages 20 years, 4 ages 25 years, 10 ages 10 years and 20 ages 5 years. 
#### we do a hierarchical model of each
#### we do a total abundance model of each
#### we do a single age method of each.
#### Then we export the data to the Sim_data_vis.R script to visualize it. 
#### It produces Figure 2 in the manuscript 

#### 3/18/2022

source("analysis_functions.R")
devtools::install_github("tanyalrogers/GPEDM") #update frequently

library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(rEDM)
library(GPEDM)
library(Metrics)
library(corrplot)
library(pracma)
library(parallel)

## Read in the simulation data & make it long format ## BUT TURN IT INTO A LIST!!! 
preylist1 <-read.csv("Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
  pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()
preylist1 <- preylist1 %>% split(f=preylist1$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist2 <-read.csv("Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
  pivot_longer(3:23, names_to = "age_class")%>%as.data.frame()%>%split(f=preylist2$index)
preylist2 <- preylist2 %>% split(f=preylist2$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist3 <-read.csv("Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
  pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()%>%split(f=preylist3$index)
preylist3 <- preylist3 %>% split(f=preylist3$index)%>%lapply(function(x) x[!names(x) %in% c("index")])


#############THE 100 DATA POINTS FUNCTION############################################################
#############creates the data that are used in figure 2 and supplmentary figure s3
############# makes datasets of different combinations of 100 data points and also of various other sizes and lengths
############# fits the GP and saves only the fitstats and the hyperparameter (phi) values. 

GP100 <-function(plist){
 
 #inputs
 plist<-as.data.frame(plist)
 names(plist)<-c("time_step", "age_class","value")
 #create the Ntotal age class. 
 plistNT <-plist %>% group_by(time_step)%>%summarize(NewValue=sum(value))%>%mutate(age_class="all", value=log(NewValue))
 plist <-mutate(plist, value=log(value)) #REMEMBER TO TURN THIS ON AND OFF
 plist <-plist %>%group_by(age_class)%>%arrange(age_class, time_step)#make sure it is sorted by age_class then year
 plist<-filter(plist, age_class !="V21") # we're not doing the plus group ever 
 
 agelist <-unique(plist$age_class)
 agelist <-agelist[c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)]
 
 #############Format the data##############
 #5 AGE CLASSES, 20 YEARS
 prey5 <-filter(plist, time_step >=300 & time_step <= 330 & age_class %in% agelist[1:5])%>% as.data.frame()
 prey5Lags = makelags(data=prey5, yd="value", pop="age_class", E=round(sqrt(20)), tau=2)
 prey5 = cbind(prey5,prey5Lags)
 prey5.train = filter(prey5, time_step <= (max(prey5$time_step)-10))
 prey5.test = filter(prey5, time_step > (max(prey5$time_step)-10))
 
 #4 AGE CLASSES, 25 YEARS 
 prey4 <-filter(plist, time_step >=300 & time_step <= 335 & age_class %in% agelist[1:4])%>% as.data.frame()
 prey4Lags = makelags(data=prey4, yd="value", pop="age_class", E=round(sqrt(25)), tau=1)
 prey4 = cbind(prey4,prey4Lags)
 prey4.train = filter(prey4, time_step <= (max(prey4$time_step)-10))
 prey4.test = filter(prey4, time_step > (max(prey4$time_step)-10))
 
 # 10 AGE CLASSES, 10 YEARS ##
 prey10 <-filter(plist, time_step >=300 & time_step <= 320 & age_class %in% agelist[1:10])%>% as.data.frame()
 prey10Lags = makelags(data=prey10, yd="value", pop="age_class", E=round(sqrt(10)), tau=1)
 prey10 = cbind(prey10,prey10Lags)
 prey10.train = filter(prey10, time_step <= (max(prey10$time_step)-10))
 prey10.test = filter(prey10, time_step > (max(prey10$time_step)-10))
 
 # 20 AGE CLASSES, 5 YEARS ##
 prey205 <- filter(plist, time_step >=300 & time_step <= 315 & age_class %in% agelist[1:20])%>% as.data.frame()
 prey205Lags = makelags(data=prey205, yd="value", pop="age_class", E=round(sqrt(5)), tau=1)
 prey205 = cbind(prey205,prey205Lags)
 prey205.train = filter(prey205, time_step <= (max(prey205$time_step)-10))
 prey205.test = filter(prey205, time_step > (max(prey205$time_step)-10))
 
 # N_TOTAL GROUP***** collected from 20 age classes over 100 years.
 preyNT <-plistNT%>%filter(time_step >= 300 & time_step <= 410)%>% as.data.frame()
 preyNTLags = makelags(data=preyNT, yd="value", E=round(sqrt(100)), tau=1)
 preyNT = cbind(preyNT,preyNTLags)
 preyNT.train = filter(preyNT, time_step <= (max(preyNT$time_step)-10))
 preyNT.test = filter(preyNT, time_step > (max(preyNT$time_step)-10))
 
 
 #100 time points for individual ages model - do not include the Ntotal age class. 
 preypivall <- plist %>%filter(time_step > 299 & time_step < 410 & age_class !="all") %>% mutate(age_class=as.factor(age_class))%>% as.data.frame()
 
 #######################################################################################################################################
 #######################################################################################################################################
 # 20 YEARS 5 AGES 
 prey520_final <-fitGP(data = prey5.train, yd = "value", xd=colnames(prey5Lags),datanew=prey5.test,pop="age_class",scaling = "local",predictmethod = "loo")
 #Single age method
 preypivall520 <-filter(preypivall, age_class %in% c("V1","V2","V3","V4","V5"))%>%as.data.frame()
 prey520_ages <-indv_age(preypivall520)%>%mutate(model="prey520", tslength=100, approach="single_age")
 #aggregate the single age results "agg hier"
 prey520agg <-aggPred(preypivall520,20) # the aggregate r2 for 20 years. 
 #aggregate the hierarchical results "sum hier"
 prey520sumhier <-full_join(prey520_final$outsampresults, prey520_final$insampresults, by=c("timestep","pop","obs"))%>%
  mutate(newpred=exp(predmean.x), newobs=exp(obs))%>% filter(!is.na(newpred))%>%
  group_by(timestep)%>%summarise(predmean=log(sum(newpred)), Obs=log(sum(newobs)))%>%mutate(model="prey520")
 #Ntotal abundance, different years. 
 prey520NT <-ntotal_app(prey5[,1:3]) 
 prey520NT <-as.data.frame(t(prey520NT))%>%mutate(approach="TAindex",model="prey520")
 
 # 25 years of 4 age classes
 prey425_final <-fitGP(data = prey4.train, yd = "value", xd=colnames(prey4Lags),datanew=prey4.test,pop="age_class",scaling = "local",predictmethod = "loo")
 #now each age class individually
 preypivall425 <-filter(preypivall, age_class %in% c("V1","V2","V3","V4"))%>%as.data.frame()
 prey425_ages <-indv_age(preypivall425)%>%mutate(model="prey425",tslength=100, approach="single_age")
 prey425agg <-aggPred(preypivall425,25) # the aggregate r2 for 25 years. 
 prey425sumhier <-full_join(prey425_final$outsampresults, prey425_final$insampresults, by=c("timestep","pop","obs"))%>%
  mutate(newpred=exp(predmean.x), newobs=exp(obs))%>% filter(!is.na(newpred))%>%
  group_by(timestep)%>%summarise(predmean=log(sum(newpred)), Obs=log(sum(newobs)))%>%mutate(model="prey425")
 #Ntotal abundance, different years. 
 prey425NT <-ntotal_app(prey4[,1:3])
 prey425NT <-as.data.frame(t(prey425NT))%>%mutate(approach="TAindex",model="prey425")
 
 # 10 years of 10 age classes
 prey1010_final <-fitGP(data = prey10.train, yd = "value", xd=colnames(prey10Lags),datanew=prey10.test,pop="age_class",scaling = "local",predictmethod = "loo")
 #now each age class individually
 preypivall1010 <-filter(preypivall, age_class %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))%>%as.data.frame()
 prey1010_ages <-indv_age(preypivall1010)%>%mutate(model="prey1010",tslength=100, approach="single_age")
 prey1010agg <-aggPred(preypivall1010,10) # the aggregate r2 for 10 years. 
 prey1010sumhier <-full_join(prey1010_final$outsampresults, prey1010_final$insampresults, by=c("timestep","pop","obs"))%>%
  mutate(newpred=exp(predmean.x), newobs=exp(obs))%>% filter(!is.na(newpred))%>%
  group_by(timestep)%>%summarise(predmean=log(sum(newpred)), Obs=log(sum(newobs)))%>%mutate(model="prey1010")
 #Ntotal abundance, different years. 
 prey1010NT <-ntotal_app(prey10[,1:3])
 prey1010NT <-as.data.frame(t(prey1010NT))%>%mutate(approach="TAindex",model="prey1010")
 
 # 5 years of 20 age classes
 prey205_final <-fitGP(data = prey205.train, yd = "value", xd=colnames(prey205Lags),datanew=prey205.test,pop="age_class",scaling = "local",predictmethod = "loo")
 #now each age class individually
 preypivall205 <-filter(preypivall, age_class %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20"))%>%as.data.frame()#all the age classes. 
 prey205_ages <-indv_age(preypivall205)%>%mutate(model="prey205",tslength=100, approach="single_age")
 prey205agg <-aggPred(preypivall205,5) # the aggregate r2 for 5 years.
 prey205sumhier <-full_join(prey205_final$outsampresults, prey205_final$insampresults, by=c("timestep","pop","obs"))%>%
  mutate(newpred=exp(predmean.x), newobs=exp(obs))%>% filter(!is.na(newpred))%>%
  group_by(timestep)%>%summarise(predmean=log(sum(newpred)), Obs=log(sum(newobs)))%>%mutate(model="prey205")
 #Ntotal abundance, different years. 
 prey205NT <-ntotal_app(prey205[,1:3])
 prey205NT <-as.data.frame(t(prey205NT))%>%mutate(approach="TAindex",model="prey205")
 
 
 #################################################################################################################################  
 ################################# Extract Fit stats ###############################################################################
 outsamp <-bind_rows(prey205_final$outsampfitstats,prey1010_final$outsampfitstats,prey520_final$outsampfitstats,prey425_final$outsampfitstats)
 names(outsamp) <-c("OOS_R2","OOS_rmse")
 insamp <-bind_rows(prey205_final$insampfitstats,prey1010_final$insampfitstats,prey520_final$insampfitstats,prey425_final$insampfitstats)
 rhos <-c(tail(prey205_final$pars,1),tail(prey1010_final$pars,1),tail(prey520_final$pars,1),tail(prey425_final$pars,1))
 
 fitstats <-bind_cols(outsamp,insamp,rhos)%>%as.data.frame()
 rownames(fitstats) <-c("5yrs20ages","10yrs10ages","20yrs5ages","25yrs4ages")
 colnames(fitstats)[8]<-"rho"
 fitstats$tslength <-c(5,10,20,25)
 fitstats <- mutate(fitstats, approach="hier")
 
 #concatonate the hier results
 fitsumhier <-bind_rows(prey520sumhier,prey425sumhier,prey1010sumhier,prey205sumhier)%>%
  group_by(model)%>%summarize(OOS_R2=getR2(obs=Obs,pred=predmean))%>%mutate(approach="sumhier")%>%as.data.frame
 
 #concatonate apples to apples results
 Ntotsyrs <-bind_rows(prey205NT,prey1010NT,prey520NT,prey425NT)%>%mutate(`age class`=as.character(`age class`))
 
 #fitstats_agg
 ###### for now, we are saying it's OOS_R2, but it's really in sample R2
 fitagg <-data.frame(matrix(ncol=3, nrow=4))
 fitagg[,1] <-c(prey205agg, prey1010agg,prey520agg,prey425agg)
 fitagg[,2] <-c("prey205","prey1010","prey520","prey425")
 fitagg[,3] <-c(5,10,20,25)
 names(fitagg) <-c("OOS_R2","model","tslength")
 fitagg <-mutate(fitagg, approach="aggregate")
 
 ####### fitstats from individual ages. 
 fitstats_ages <-bind_rows(prey205_ages,prey1010_ages,prey520_ages,prey425_ages)%>%mutate(tslength=100)
 
 pars <-bind_rows(prey205_final$pars,prey1010_final$pars,prey520_final$pars,prey425_final$pars)%>%as.data.frame()
 rownames(pars) <-c("5yrs20ages","10yrs10ages","20yrs5ages","25yrs4ages")
 pars$tslength <-c(5,10,20,25)
 pars$approach <-c("hier","hier","hier","hier")
 pars <-mutate(pars, maxE=round(sqrt(tslength)))
 
 #outputs
 gplout <-list(pars,fitstats,fitstats_ages,fitagg,fitsumhier)
 names(gplout) <-c("pars","fitstats","fitstats_ages","fitagg","fitsumhier")
 gplout
}

#############################################################################################################################
#############################################################################################################################
############FUNCTION TO PROCESS OUTPUT FROM GP100 FUNCTION 
### corralling the elements from the list output##################### 
## the input is a list, which is the output of the GP100 function.

PROCESS100 <-function(parafoo){
#fitstats#
fstats <-lapply(parafoo, `[`, "fitstats")
for (i in 1:length(fstats)){
 fstats[[i]] <- mutate(fstats[[i]]$fitstats, iter=i)%>%rownames_to_column()
}
fstats <-fstats %>% bind_rows()%>%as.data.frame()%>%dplyr::rename(models=rowname)%>%mutate(age_class="all")
fstats <-mutate(fstats, model=ifelse(models=="5yrs20ages","prey205",ifelse(models=="20yrs5ages","prey520",ifelse(models=="10yrs10ages","prey1010",ifelse(models=="25yrs4ages","prey425","preyNT")))))

#fitages#
fitages <-lapply(parafoo, `[`, "fitstats_ages")
for (i in 1:length(fitages)){
 fitages[[i]]<- mutate(fitages[[i]]$fitstats_ages, iter=i)%>%rownames_to_column()
}
fitages <-fitages %>% bind_rows()%>%as.data.frame()
fitages <- mutate(fitages, models=ifelse(model=="prey205","5yrs20ages",ifelse(model=="prey520","20yrs5ages",ifelse(model=="prey1010","10yrs10ages",ifelse(model=="prey425","25yrs4ages","20agesdiffyrs")))))%>%
 dplyr::select(-rowname)%>%dplyr::rename(age_class="age class")

#aggfits#
fitagg <-lapply(parafoo, `[`, "fitagg")
for (i in 1:length(fitagg)){
 fitagg[[i]]<- mutate(fitagg[[i]]$fitagg, iter=i)%>%rownames_to_column()
}
fitagg <-fitagg %>% bind_rows()%>%as.data.frame()
fitagg <-dplyr::select(fitagg,-rowname)%>% mutate(age_class="all",models=ifelse(model=="prey205","5yrs20ages",ifelse(model=="prey520","20yrs5ages",ifelse(model=="prey1010","10yrs10ages","25yrs4ages"))))

#hierarchical sum fits#
fitsumhier <-lapply(parafoo, `[`, "fitsumhier")
for (i in 1:length(fitsumhier)){
 fitsumhier[[i]]<- mutate(fitsumhier[[i]]$fitsumhier, iter=i)%>%rownames_to_column()
}
fitsumhier <-fitsumhier %>% bind_rows()%>%as.data.frame()
fitsumhier <-dplyr::select(fitsumhier,-rowname)%>% mutate(age_class="all",models=ifelse(model=="prey205","5yrs20ages",ifelse(model=="prey520","20yrs5ages",ifelse(model=="prey1010","10yrs10ages","25yrs4ages"))))


#combine fitstats and fitages to newfits and return newfits#
newfits <-bind_rows(fstats,fitages,fitagg,fitsumhier)
newfits
}

########## For each simulation: Import sim data --> apply GPLOOP function --> format output --> save result. 
#simulation data was created in "make_simulation_data.R"

############# Simulation 1: Apply GP Loop and corral output #####################

######## Apply the GP100 function #########
## this takes a while, so we apply it in parallel.
numCores = detectCores()
system.time(parafoo <-mclapply(preylist1, GP100, mc.cores=numCores))

#process hyperparameters#
qpars <-lapply(parafoo, `[`, "pars")
for (i in 1:length(qpars)){
 qpars[[i]] <- mutate(qpars[[i]]$pars, iter=i)%>%rownames_to_column()
}
qpars <-qpars %>% bind_rows()%>%as.data.frame()

#process fitstats using process100 function
newfits <-PROCESS100(parafoo)

#######Save to csv######
write.csv(qpars, file="Sim1Phis_100data.csv")
write.csv(newfits, file="Sim1fitstats_100data.csv")

############# Simulation 2: Apply GP Loop and corral output #####################

######## Apply the GP100 function #########
## this takes a while, so we apply it in parallel.
numCores = detectCores()
system.time(parafoo <-mclapply(preylist2, GP100, mc.cores=numCores))

#process hyperparameters#
qpars <-lapply(parafoo, `[`, "pars")
for (i in 1:length(qpars)){
 qpars[[i]] <- mutate(qpars[[i]]$pars, iter=i)%>%rownames_to_column()
}
qpars <-qpars %>% bind_rows()%>%as.data.frame()

#process fitstats
newfits <-PROCESS100(parafoo)

#######Save to csv######
write.csv(qpars, file="Sim2Phis_100data.csv")
write.csv(newfits, file="Sim2fitstats_100data.csv")

############# Simulation 3: Apply GP Loop and corral output #####################

######## Apply the GP100 function #########
## this takes a while, so we apply it in parallel.
numCores = detectCores()
system.time(parafoo <-mclapply(preylist3, GP100, mc.cores=numCores))

#process hyperparameters#
qpars <-lapply(parafoo, `[`, "pars")
for (i in 1:length(qpars)){
 qpars[[i]] <- mutate(qpars[[i]]$pars, iter=i)%>%rownames_to_column()
}
qpars <-qpars %>% bind_rows()%>%as.data.frame()

#process fitstats
newfits <-PROCESS100(parafoo)

#######Save to csv######
write.csv(qpars, file="Sim3Phis_100data.csv")
write.csv(newfits, file="Sim3fitstats_100data.csv")
