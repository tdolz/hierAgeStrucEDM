### Simulation Sandbox ###
### Trying to add variation to the simulations #########
### 4/19/2022 ##

source("analysis_functions.R")
#devtools::install_github("tanyalrogers/GPEDM")
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





### SIMULATION I ################################################################################
# original parameter values# 

maxiter = 150. #try making 200 of them and deleting all the zero peak one
preylist<-list()
predlist<-list()
count0peaks <-list()
recnoise <-matrix(nrow=maxiter, ncol=5)
colnames(recnoise)<-c("recnoise","meanpeaks","varprey","meanPeriod","sdperiod")


for (m in 1:maxiter){
 
 ### The simulation (scenario 5) ########
 ### RECRUITMENT OPTIONS ####
 #ricker system param
 r = 1/1000
 sdrec= 0.1
 recnoise1 =  rnorm(1,0,0.2)
 ########
 Am = 20 # max age
 t = 500 # num time steps
 ages = 1:Am
 nsurv=0.8^(ages-1)
 N0 = 1000
 pop2=matrix(nrow=Am, ncol=t) #matrix of 10 ages (columns) and 500 time steps (rows)
 pop2[,1]=N0*nsurv #populate first time step (column)
 s = 0.8 #constant survival
 fec =  c(0,0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100) #constant fecundity after maturation at age 3.
 
 for (i in 1:t-1){ #columns
  if (i == 1) {
   next
  }
  for(j in 1:Am){ #rows
   if (j==1) {
    E <-sum(fec*pop2[,i-1])
    #E <-sum(pop2[,i-1])
    pop2[j,i]=E*exp(r*(1-E))*exp(sdrec*recnoise1)
   }
   else{
    pop2[j,i]=s*pop2[j-1,i-1]
   }
  }
 }
 
 prey <-as.data.frame(t(pop2)) %>% rownames_to_column(var="time_step") %>%
  mutate(time_step=as.numeric(time_step)); 
 preypiv <-pivot_longer(prey, 2:21, names_to = "age_class")
 
 #measure how many peaks in a ten year time period.
 countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
  dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
 count0peaks[[m]] <-sum(countpeaks$pks==0)
 
 #outputs
 groupmeans <-preypiv %>% group_by(age_class)%>%summarize(preymeans=mean(value))
 mgroups <-min(groupmeans$preymeans)
 preyvar <-preypiv %>%group_by(age_class)%>%summarize(varprey=var(value))
 
 #mean period
 meanper <-preypiv %>%filter(time_step > 299 & !is.na(value))%>%
  dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(periodt=period(value))
 
 recnoise[m,1]<- recnoise1
 recnoise[m,2]<-mean(countpeaks$pks)
 recnoise[m,3]<-mean(preyvar$varprey)
 recnoise[m,4]<-mean(meanper$periodt)
 recnoise[m,5]<-sd(meanper$periodt)
 
 
 if(sum(countpeaks$pks==0)< 15) { #we can't get them all in there unfortunately. 
  preylist[[m]] <- prey
 }else
  preylist[[m]] <- NA
 
 
 
}

##############################################################
####################
## View and save the preylist. 
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#how many peaks on average?
ntotalpeaks <-recnoises %>% group_by(index)%>%summarize(avpks=mean(meanpeaks))
mean(ntotalpeaks$avpks)
sd(ntotalpeaks$avpks)

#what is the period?
mean(recnoises$meanPeriod)
sd(recnoises$meanPeriod)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

##visualize the data to make sure it's ok. 
preylist.mean <-pivot_longer(preylists,3:22, names_to = 'age_class')%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),.groups='drop')

preylist.mean%>%
 filter(time_step >= 300 & time_step <=400)%>%
 ggplot(aes(time_step,mean.value))+
 #geom_line()+
 geom_point(size=0.5)+
 facet_wrap(~age_class,scales="free")+
 theme_classic()

### SIMULATION I ################################################################################
# NEW parameter values# 
set.seed(4649)
maxiter = 1000. #try making 200 of them and deleting all the zero peak one
preylist<-list()
predlist<-list()
count0peaks <-list()
recnoise <-matrix(nrow=maxiter, ncol=5)
colnames(recnoise)<-c("recnoise","meanpeaks","varprey","meanPeriod","sdperiod")


for (m in 1:maxiter){
  
  tryCatch({
 
 ### The simulation (scenario 5) ########
 ### RECRUITMENT OPTIONS ####
 #ricker system param
 r = 1/1000
 sdrec= 1 #increased
 recnoise1 =  rnorm(1,0,1) #increased SD
 ########
 Am = 20 # max age
 t = 500 # num time steps
 ages = 1:Am
 nsurv=0.8^(ages-1)
 N0 = 1000
 pop2=matrix(nrow=Am, ncol=t) #matrix of 10 ages (columns) and 500 time steps (rows)
 pop2[,1]=N0*nsurv #populate first time step (column)
 s = 0.8 #constant survival
 fec =  c(0,0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100) #constant fecundity after maturation at age 3.
 
 for (i in 1:t-1){ #columns
  if (i == 1) {
   next
  }
  for(j in 1:Am){ #rows
   if (j==1) {
    E <-sum(fec*pop2[,i-1])
    #E <-sum(pop2[,i-1])
    pop2[j,i]=E*exp(r*(1-E))*exp(sdrec*recnoise1)
   }
   else{
    pop2[j,i]=s*pop2[j-1,i-1]
   }
  }
 }
 
 prey <-as.data.frame(t(pop2)) %>% rownames_to_column(var="time_step") %>%
  mutate(time_step=as.numeric(time_step)); 
 preypiv <-pivot_longer(prey, 2:21, names_to = "age_class")
 
 #measure how many peaks in a ten year time period.
 countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
  dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
 count0peaks[[m]] <-sum(countpeaks$pks==0)
 
 #outputs
 groupmeans <-preypiv %>% group_by(age_class)%>%summarize(preymeans=mean(value))
 mgroups <-min(groupmeans$preymeans)
 preyvar <-preypiv %>%group_by(age_class)%>%summarize(varprey=var(value))
 
 #mean period
 meanper <-preypiv %>%filter(time_step > 299 & !is.na(value))%>%
  dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(periodt=period(value))
 
 recnoise[m,1]<- recnoise1
 recnoise[m,2]<-mean(countpeaks$pks)
 recnoise[m,3]<-mean(preyvar$varprey)
 recnoise[m,4]<-mean(meanper$periodt)
 recnoise[m,5]<-sd(meanper$periodt)
 
 
 if(mean(meanper$periodt)<12) { #we can't get them all in there unfortunately. 
  preylist[[m]] <- prey
 }else
  preylist[[m]] <- NA
 
  }, error=function(e){})
 
}

##############################################################
####################
## View and save the preylist. 
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[!sapply(preylist,is.null)] #remove all null elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#how many peaks on average?
ntotalpeaks <-recnoises %>% group_by(index)%>%summarize(avpks=mean(meanpeaks))
mean(ntotalpeaks$avpks, na.rm=T)
sd(ntotalpeaks$avpks, na.rm=T)

#what is the period?
mean(recnoises$meanPeriod, na.rm=T)
sd(recnoises$meanPeriod, na.rm=T)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

#just look at the first one, not the mean
preylists%>%pivot_longer(3:22, names_to = 'age_class')%>%
  filter(index==1)%>%
  filter(time_step >= 300 & time_step <=400)%>%
  #filter(time_step >= 300 & time_step <=320)%>%
  ggplot(aes(time_step,value))+
  #ggplot(aes(time_step,mean.value, color=age_class))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()


##visualize the data to make sure it's ok. 
preylist.mean <-pivot_longer(preylists,3:22, names_to = 'age_class')%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),.groups='drop')

preylist.mean%>%
 filter(time_step >= 300 & time_step <=400)%>%
 ggplot(aes(time_step,mean.value))+
 geom_line()+
 facet_wrap(~age_class,scales="free")+
 theme_classic()






### SIMULATION II ################################################################################
############# Original Parameter Values ################################

maxiter = 150. #try making 200 of them and deleting all the zero peak one.s 
preylist<-list()
predlist<-list()
count0peaks <-list()
recnoise <-matrix(nrow=maxiter, ncol=6)
colnames(recnoise)<-c("recnoise1","recnoise2","meanpeaks","varprey","meanPeriod","sdperiod")


for (m in 1:maxiter){
  
  ### The simulation (scenario 5) ########
  ### RECRUITMENT OPTIONS ####
  #Species 1 - Prey
  phi1 = 1/1000
  sdrec1 = .01  
  recnoise1 = rnorm(1,0,0.2)
  beta1 =0.002
  basefec1=c(0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100) # prey
  BM1 = 2 #basefec modifier
  
  #Species 2 - Predator
  phi2 = 1/100
  sdrec2 = .01 
  recnoise2 = rnorm(1,0,0.2)
  beta2 = 0.001
  basefec2=c(0,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
  BM2 = 1 #basefec modifier 2 is good
  s = 0.8 #constant survival
  
  # constants common to both species (for now)#
  t = 500 # num time steps
  N01 = 1000 #initial prey #1000 is good
  N02 = N01/10 #initial pred N0/10 is good. 
  Am = 21
  ages = seq(1,Am,1)
  nsurv=0.8^(ages-1) #both species start out with the same age structure. 
  
  #make two population matrices
  pop1=matrix(nrow=Am, ncol=t) #matrix of 11 ages (columns) and 500 time steps (rows)
  pop2=matrix(nrow=Am, ncol=t)
  pop1[,1]=N01*nsurv #populate first time step (column) for species 1
  pop2[,1]=N02*nsurv #populate the first time step for species 2
  
  
  ### TIME LOOP ####
  for (i in 1:t-1){ #columns
    if (i == 1) { #skip the first time step (column) because it is already populated.
      next
    }
    ####RECRUITMENT#####
    ####
    for(j in 1:Am){
      if (j == 1) {
        pop1[j,i]=sum(pop1[,i-1]*basefec1)*exp(-phi1*sum(pop1[,i-1]*basefec1)-beta1*sum(pop2[,i-1]))*exp(sdrec1*recnoise1)
        pop2[j,i]=sum(pop2[,i-1]*basefec2)*exp(-phi2*sum(pop2[,i-1]*basefec2)-beta2*sum(pop1[,i-1]))*exp(sdrec2*recnoise2) 
        ###SURVIVAL###
        ###plus group survival
      }else if (j == Am){  #survival in the plus group
        pop1[j,i]=pop1[j-1,i-1]*exp(-sum(pop2[,i-1])/(sum(pop2[,i-1])+sum(pop1[,i-1])))+
          pop1[j,i-1]*exp(-sum(pop2[,i-1])/(sum(pop2[,i-1])+sum(pop1[,i-1])))
        pop2[j,i]=pop2[j-1,i-1]*s
      }
      #regular survival   
      else{
        pop1[j,i]=pop1[j-1,i-1]*exp(-sum(pop2[,i-1])/(sum(pop2[,i-1])+sum(pop1[,i-1])))
        pop2[j,i]=pop2[j-1,i-1]*s
      }
    }
  }
  
  ### age specific ###
  prey <-as.data.frame(t(pop1)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); prey <-prey[-500,]
  preypiv <-pivot_longer(prey, 2:22, names_to = "age_class")
  pred <-as.data.frame(t(pop2)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); pred <-pred[-500,]
  predpiv <-pivot_longer(pred, 2:22, names_to = "age_class")
  
  #outputs
  groupmeans <-preypiv %>% group_by(age_class)%>%summarize(preymeans=mean(value))
  predmeans <-predpiv %>% group_by(age_class)%>%summarize(predmeans=mean(value))
  
  mgroups <-min(groupmeans$preymeans)
  mpreds <-min(predmeans$predmeans)
  
  #count peaks
  countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
  
  
  #mean period
  meanper <-preypiv %>%filter(time_step > 299)%>%
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(periodt=period(value))
  
  #prey variance
  preyvar <-preypiv %>%group_by(age_class)%>%summarize(varprey=var(value))
  
  #store recnoises and mean peaks
  recnoise[m,1]<-recnoise1
  recnoise[m,2]<-recnoise2
  recnoise[m,3]<-mean(countpeaks$pks)
  recnoise[m,4]<-mean(preyvar$varprey)
  recnoise[m,5]<-mean(meanper$periodt)
  recnoise[m,6]<-sd(meanper$periodt)
  
  
  #Option2#######################################
  if(sum(countpeaks$pks==0)< 8) { #we can't get them all in there unfortunately. 
    preylist[[m]] <- prey
  }else
    preylist[[m]] <- NA
  
  
}


##########################################
####################
## View and save the preylist. 
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#how many peaks on average?
ntotalpeaks <-recnoises %>% group_by(index)%>%summarize(avpks=mean(meanpeaks))
mean(ntotalpeaks$avpks)
sd(ntotalpeaks$avpks)

#what is the period?
mean(recnoises$meanPeriod)
sd(recnoises$meanPeriod)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

##visualize the data to make sure it's ok. 
preylist.mean <-pivot_longer(preylists,3:23, names_to = 'age_class')%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),.groups='drop')

preylist.mean%>%
  filter(time_step >= 300 & time_step <=400)%>%
  ggplot(aes(time_step,mean.value))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()


############# SIMULATION II New parameter values ################################
###########################################################################

maxiter = 200. #try making 200 of them and deleting all the zero peak one.s 
preylist<-list()
predlist<-list()
count0peaks <-list()
recnoise <-matrix(nrow=maxiter, ncol=6)
colnames(recnoise)<-c("recnoise1","recnoise2","meanpeaks","varprey","meanPeriod","sdperiod")


for (m in 1:maxiter){
  
  tryCatch({
  
  ### The simulation (scenario 5) ########
  ### RECRUITMENT OPTIONS ####
  #Species 1 - Prey
  phi1 = 1/1000
  sdrec1 = 1  
  recnoise1 = rnorm(1,0,1)
  beta1 =0.002
  basefec1=c(0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100) # prey
  BM1 = 2 #basefec modifier
  
  #Species 2 - Predator
  phi2 = 1/100
  sdrec2 = 1 
  recnoise2 = rnorm(1,0,1)
  beta2 = 0.001
  basefec2=c(0,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
  BM2 = 1 #basefec modifier 2 is good
  s = 0.8 #constant survival
  
  # constants common to both species (for now)#
  t = 500 # num time steps
  N01 = 1000 #initial prey #1000 is good
  N02 = N01/10 #initial pred N0/10 is good. 
  Am = 21
  ages = seq(1,Am,1)
  nsurv=0.8^(ages-1) #both species start out with the same age structure. 
  
  #make two population matrices
  pop1=matrix(nrow=Am, ncol=t) #matrix of 11 ages (columns) and 500 time steps (rows)
  pop2=matrix(nrow=Am, ncol=t)
  pop1[,1]=N01*nsurv #populate first time step (column) for species 1
  pop2[,1]=N02*nsurv #populate the first time step for species 2
  
  
  ### TIME LOOP ####
  for (i in 1:t-1){ #columns
    if (i == 1) { #skip the first time step (column) because it is already populated.
      next
    }
    ####RECRUITMENT#####
    ####
    for(j in 1:Am){
      if (j == 1) {
        pop1[j,i]=sum(pop1[,i-1]*basefec1)*exp(-phi1*sum(pop1[,i-1]*basefec1)-beta1*sum(pop2[,i-1]))*exp(sdrec1*recnoise1)
        pop2[j,i]=sum(pop2[,i-1]*basefec2)*exp(-phi2*sum(pop2[,i-1]*basefec2)-beta2*sum(pop1[,i-1]))*exp(sdrec2*recnoise2) 
        ###SURVIVAL###
        ###plus group survival
      }else if (j == Am){  #survival in the plus group
        pop1[j,i]=pop1[j-1,i-1]*exp(-sum(pop2[,i-1])/(sum(pop2[,i-1])+sum(pop1[,i-1])))+
          pop1[j,i-1]*exp(-sum(pop2[,i-1])/(sum(pop2[,i-1])+sum(pop1[,i-1])))
        pop2[j,i]=pop2[j-1,i-1]*s
      }
      #regular survival   
      else{
        pop1[j,i]=pop1[j-1,i-1]*exp(-sum(pop2[,i-1])/(sum(pop2[,i-1])+sum(pop1[,i-1])))
        pop2[j,i]=pop2[j-1,i-1]*s
      }
    }
  }
  
  ### age specific ###
  prey <-as.data.frame(t(pop1)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); prey <-prey[-500,]
  preypiv <-pivot_longer(prey, 2:22, names_to = "age_class")
  pred <-as.data.frame(t(pop2)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); pred <-pred[-500,]
  predpiv <-pivot_longer(pred, 2:22, names_to = "age_class")
  
  #outputs
  groupmeans <-preypiv %>% group_by(age_class)%>%summarize(preymeans=mean(value))
  predmeans <-predpiv %>% group_by(age_class)%>%summarize(predmeans=mean(value))
  
  mgroups <-min(groupmeans$preymeans)
  mpreds <-min(predmeans$predmeans)
  
  #count peaks
  countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
  
  
  #mean period
  meanper <-preypiv %>%filter(time_step > 299)%>%
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(periodt=period(value))
  
  #prey variance
  preyvar <-preypiv %>%group_by(age_class)%>%summarize(varprey=var(value))
  
  #store recnoises and mean peaks
  recnoise[m,1]<-recnoise1
  recnoise[m,2]<-recnoise2
  recnoise[m,3]<-mean(countpeaks$pks)
  recnoise[m,4]<-mean(preyvar$varprey)
  recnoise[m,5]<-mean(meanper$periodt)
  recnoise[m,6]<-sd(meanper$periodt)
  
  
  #Option2#######################################
  #if(sum(countpeaks$pks==0)< 8) { #we can't get them all in there unfortunately. 
   # preylist[[m]] <- prey
 # }else
   # preylist[[m]] <- NA
  
  #new governor which prevents population crash.
 # test <- filter(prey, time_step >=300)
 # if(mean(test$V12)> 10) { #age 12 cannot have crashed by step 300 
     #preylist[[m]] <- prey
   #  }else
   #  preylist[[m]] <- NA
  
  #another option to prevent population crash.
   # tolerance = 100 # the number of crashed values we will allow. 
  # test <- filter(preypiv, time_step >=300)
  # if(lessT1(preypiv$value)< tolerance) { #if there are any zeros 
  #preylist[[m]] <- prey
   # }else
    #preylist[[m]] <- NA
  #this is not a great option. 
  
  #mean period governor
  if(mean(meanper$periodt)<10) { #we can't get them all in there unfortunately. 
    preylist[[m]] <- prey
  }else
    preylist[[m]] <- NA
  
  #no governor
  #preylist[[m]] <- prey
  
  }, error=function(e){})
  
}


####################
####################
## View and save the preylist. 
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[!sapply(preylist,is.null)] #remove all null elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#how many peaks on average?
ntotalpeaks <-recnoises %>% group_by(index)%>%summarize(avpks=mean(meanpeaks))
mean(ntotalpeaks$avpks, na.rm=T)
sd(ntotalpeaks$avpks, na.rm=T)

#what is the period?
mean(recnoises$meanPeriod, na.rm=T)
sd(recnoises$meanPeriod, na.rm=T)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

#just look at the first one, not the mean
preylists%>%pivot_longer(3:23, names_to = 'age_class')%>%
  filter(index==1)%>%
  filter(time_step >= 300 & time_step <=400)%>%
  filter(time_step >= 300 & time_step <=320)%>%
  ggplot(aes(time_step,value))+
  #ggplot(aes(time_step,mean.value, color=age_class))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()



##visualize the data to make sure it's ok. 
preylist.mean <-pivot_longer(preylists,3:23, names_to = 'age_class')%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),.groups='drop')

preylist.mean%>%
  filter(time_step >= 300 & time_step <=400)%>%
  ggplot(aes(time_step,mean.value))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()

### SIMULATION III ################################################################################
########## Original parameters #######################################################

maxiter = 150. 
preylist<-list()
predlist<-list()
count0peaks <-list()
recnoise <-matrix(nrow=maxiter, ncol=6)
colnames(recnoise)<-c("recnoise1","recnoise2","meanpeaks","varprey","meanPeriod","sdperiod")


for (m in 1:maxiter){
  
  ##### RECRUITMENT OPTIONS ####
  #Species 1 - Prey
  phi1 = 1/1000
  sdrec1 = .02  
  recnoise1 =rnorm(1,0,0.2)
  qfec1 = c(0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1) #prey consumption conversion efficiency (< 1)
  qfec1=qfec1*8 
  CM1 = 0.01 #modifier of the pred_eat_prey matrix 0.01
  basefec1=c(0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100) # prey
  BM1 = 2 #basefec modifier
  
  #Species 2 - Predator
  phi2 = 1/100
  sdrec2 = .02 
  recnoise2=rnorm(1,0,0.2)
  qfec2 = c(0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1) #pred consumption conversion efficiency (< 1)
  qfec2=qfec2*8 #8 is good
  CM2 = .05 #prey eat pred 0.1, 0.05 for steves pred matrix. 
  basefec2=c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
  BM2 = 2 #basefec modifier 2 is good
  
  # constants common to both species (for now)#
  t = 500 # num time steps
  N01 = 1000 #initial prey #1000 is good
  N02 = N01/10 #initial pred N0/10 is good. 
  Am = 20
  ages = seq(1,Am,1)
  nsurv=0.8^(ages-1) #both species start out with the same age structure. 
  
  #make two population matrices
  pop1=matrix(nrow=Am, ncol=t) #matrix of 11 ages (columns) and 500 time steps (rows)
  pop2=matrix(nrow=Am, ncol=t)
  pop1[,1]=N01*nsurv #populate first time step (column) for species 1
  pop2[,1]=N02*nsurv #populate the first time step for species 2
  
  #fecundity matrices
  fec1=matrix(nrow=Am, ncol=t) #prey
  fec2=matrix(nrow=Am, ncol=t) #predator
  #Baseline fecundity: Populate the first column (time step) of the fecundity matrix
  fec1[,1]=basefec1
  fec2[,1]=basefec2
  
  #predators eat prey at a rate which increases with age.
  #prey never really escape predation but predation is diminished. 
  #Reverse the matrix order 
  #pred_eat_prey=read.csv("predatormatrix.csv", header=FALSE)%>%as.matrix()
  #pred_eat_prey=read.csv("stevespredatormatrix.csv", header=FALSE)%>%as.matrix()
  pred_eat_prey=read.csv("predmatrix_shortosc.csv", header=FALSE)%>%as.matrix() #the shorter oscillations. 
  colnames(pred_eat_prey)<- NULL
  #multiply by 0
  pred_eat_prey=pred_eat_prey*CM1 
  
  #prey only eat predators as eggs (for now)
  prey_eat_pred=read.csv("preymatrix.csv", header=FALSE)%>%as.matrix()
  colnames(prey_eat_pred)<- NULL
  #multiply by 0
  prey_eat_pred=prey_eat_pred*CM2
  
  ### TIME LOOP ####
  for (i in 1:t-1){ #columns
    if (i == 1) { #skip the first time step (column) because it is already populated.
      next
    }
    ### FECUNDITY
    for(j in 1:Am){  #rows
      #fecundity is the rate at which eater is eating*the thing being eaten across ages*
      fec1[j,i] = pop1[j,i-1]*(basefec1[j]+qfec1[j]*prey_eat_pred[j,]%*%pop2[,i-1]) #fecundity of prey
      fec2[j,i] = pop2[j,i-1]*(basefec2[j]+qfec2[j]*pred_eat_prey[j,]%*%pop1[,i-1]) #fecundity of pred
    }
    ####RECRUITMENT#####
    ####
    for(j in 1:Am){
      if (j == 1) {
        pop1[j,i]=sum(fec1[,i])*exp(-phi1*sum(fec1[,i]))*exp(sdrec1*recnoise1)
        pop2[j,i]=sum(fec2[,i])*exp(-phi2*sum(fec2[,i]))*exp(sdrec2*recnoise2) 
      }
      #regular survival
      else{
        pop1[j,i]=pop1[j-1,i-1]*exp(-sum(pred_eat_prey[,j-1]*pop2[,i-1]))
        pop2[j,i]=pop2[j-1,i-1]*exp(-sum(prey_eat_pred[,j-1]*pop1[,i-1]))
      }
    }
  }
  ### age specific ###
  prey <-as.data.frame(t(pop1)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); prey <-prey[-500,]
  preypiv <-pivot_longer(prey, 2:21, names_to = "age_class")
  pred <-as.data.frame(t(pop2)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); pred <-pred[-500,]
  predpiv <-pivot_longer(pred, 2:21, names_to = "age_class")
  
  #outputs
  groupmeans <-preypiv %>% group_by(age_class)%>%summarize(preymeans=mean(value))
  predmeans <-predpiv %>% group_by(age_class)%>%summarize(predmeans=mean(value))
  
  mgroups <-min(groupmeans$preymeans)
  mpreds <-min(predmeans$predmeans)
  
  #count peaks
  countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
  
  #mean period
  meanper <-preypiv %>%filter(time_step > 299)%>% #period after the burn-in period
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(periodt=period(value))
  
  #prey variance
  preyvar <-preypiv %>%group_by(age_class)%>%summarize(varprey=var(value))
  
  #store recnoises and mean peaks
  recnoise[m,1]<-recnoise1
  recnoise[m,2]<-recnoise2
  recnoise[m,3]<-mean(countpeaks$pks)
  recnoise[m,4]<-mean(preyvar$varprey)
  recnoise[m,5]<-mean(meanper$periodt)
  recnoise[m,6]<-sd(meanper$periodt)
  
  
  #if(mgroups > 0.1 & mpreds > 0.1 & mean(countpeaks$pks) > 1.2) {
  if(mean(countpeaks$pks) > 1.2) {
    preylist[[m]] <- prey
  }else
    preylist[[m]] <- NA
}

##############################################################
####################
## View and save the preylist. 
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#how many peaks on average?
ntotalpeaks <-recnoises %>% group_by(index)%>%summarize(avpks=mean(meanpeaks))
mean(ntotalpeaks$avpks)
sd(ntotalpeaks$avpks)

#what is the period?
mean(recnoises$meanPeriod)
sd(recnoises$meanPeriod)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

##visualize the data to make sure it's ok. 
preylist.mean <-pivot_longer(preylists,3:22, names_to = 'age_class')%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),.groups='drop')

preylist.mean%>%
  filter(time_step >= 300 & time_step <=400)%>%
  ggplot(aes(time_step,mean.value))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()

### SIMULATION III #################################
########## New Parameters #########################################

maxiter = 1000. 
preylist<-list()
predlist<-list()
count0peaks <-list()
recnoise <-matrix(nrow=maxiter, ncol=6)
colnames(recnoise)<-c("recnoise1","recnoise2","meanpeaks","varprey","meanPeriod","sdperiod")


for (m in 1:maxiter){
  
  tryCatch({
  
  ##### RECRUITMENT OPTIONS ####
  #Species 1 - Prey
  phi1 = 1/1000
  sdrec1 = 1  
  recnoise1 =rnorm(1,0,1)
  qfec1 = c(0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1) #prey consumption conversion efficiency (< 1)
  qfec1=qfec1*8 
  CM1 = 0.01 #modifier of the pred_eat_prey matrix 0.01
  basefec1=c(0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100) # prey
  BM1 = 2 #basefec modifier
  
  #Species 2 - Predator
  phi2 = 1/100
  sdrec2 = 1 
  recnoise2=rnorm(1,0,1)
  qfec2 = c(0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1) #pred consumption conversion efficiency (< 1)
  qfec2=qfec2*8 #8 is good
  CM2 = .05 #prey eat pred 0.1, 0.05 for steves pred matrix. 
  basefec2=c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
  BM2 = 2 #basefec modifier 2 is good
  
  # constants common to both species (for now)#
  t = 500 # num time steps
  N01 = 1000 #initial prey #1000 is good
  N02 = N01/10 #initial pred N0/10 is good. 
  Am = 20
  ages = seq(1,Am,1)
  nsurv=0.8^(ages-1) #both species start out with the same age structure. 
  
  #make two population matrices
  pop1=matrix(nrow=Am, ncol=t) #matrix of 11 ages (columns) and 500 time steps (rows)
  pop2=matrix(nrow=Am, ncol=t)
  pop1[,1]=N01*nsurv #populate first time step (column) for species 1
  pop2[,1]=N02*nsurv #populate the first time step for species 2
  
  #fecundity matrices
  fec1=matrix(nrow=Am, ncol=t) #prey
  fec2=matrix(nrow=Am, ncol=t) #predator
  #Baseline fecundity: Populate the first column (time step) of the fecundity matrix
  fec1[,1]=basefec1
  fec2[,1]=basefec2
  
  #predators eat prey at a rate which increases with age.
  #prey never really escape predation but predation is diminished. 
  #Reverse the matrix order 
  #pred_eat_prey=read.csv("predatormatrix.csv", header=FALSE)%>%as.matrix()
  #pred_eat_prey=read.csv("stevespredatormatrix.csv", header=FALSE)%>%as.matrix()
  pred_eat_prey=read.csv("predmatrix_shortosc.csv", header=FALSE)%>%as.matrix() #the shorter oscillations. 
  colnames(pred_eat_prey)<- NULL
  #multiply by 0
  pred_eat_prey=pred_eat_prey*CM1 
  
  #prey only eat predators as eggs (for now)
  prey_eat_pred=read.csv("preymatrix.csv", header=FALSE)%>%as.matrix()
  colnames(prey_eat_pred)<- NULL
  #multiply by 0
  prey_eat_pred=prey_eat_pred*CM2
  
  ### TIME LOOP ####
  for (i in 1:t-1){ #columns
    if (i == 1) { #skip the first time step (column) because it is already populated.
      next
    }
    ### FECUNDITY
    for(j in 1:Am){  #rows
      #fecundity is the rate at which eater is eating*the thing being eaten across ages*
      fec1[j,i] = pop1[j,i-1]*(basefec1[j]+qfec1[j]*prey_eat_pred[j,]%*%pop2[,i-1]) #fecundity of prey
      fec2[j,i] = pop2[j,i-1]*(basefec2[j]+qfec2[j]*pred_eat_prey[j,]%*%pop1[,i-1]) #fecundity of pred
    }
    ####RECRUITMENT#####
 
       ####
    for(j in 1:Am){
      if (j == 1) {
        pop1[j,i]=sum(fec1[,i])*exp(-phi1*sum(fec1[,i]))*exp(sdrec1*recnoise1)
        pop2[j,i]=sum(fec2[,i])*exp(-phi2*sum(fec2[,i]))*exp(sdrec2*recnoise2) 
      }
      #regular survival
      else{
        pop1[j,i]=pop1[j-1,i-1]*exp(-sum(pred_eat_prey[,j-1]*pop2[,i-1]))
        pop2[j,i]=pop2[j-1,i-1]*exp(-sum(prey_eat_pred[,j-1]*pop1[,i-1]))
      }
    }
  }
    
  ### age specific ###
  prey <-as.data.frame(t(pop1)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); prey <-prey[-500,]
  preypiv <-pivot_longer(prey, 2:21, names_to = "age_class")
  pred <-as.data.frame(t(pop2)) %>% rownames_to_column(var="time_step") %>%
    mutate(time_step=as.numeric(time_step)); pred <-pred[-500,]
  predpiv <-pivot_longer(pred, 2:21, names_to = "age_class")
  
  #outputs
  groupmeans <-preypiv %>% group_by(age_class)%>%summarize(preymeans=mean(value))
  predmeans <-predpiv %>% group_by(age_class)%>%summarize(predmeans=mean(value))
  
  mgroups <-min(groupmeans$preymeans)
  mpreds <-min(predmeans$predmeans)
  
  #count peaks
  countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
  
  #mean period
  meanper <-preypiv %>%filter(time_step > 299)%>% #period after the burn-in period
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(periodt=period(value))
  
  #prey variance
  preyvar <-preypiv %>%group_by(age_class)%>%summarize(varprey=var(value))
  
  #store recnoises and mean peaks
  recnoise[m,1]<-recnoise1
  recnoise[m,2]<-recnoise2
  recnoise[m,3]<-mean(countpeaks$pks)
  recnoise[m,4]<-mean(preyvar$varprey)
  recnoise[m,5]<-mean(meanper$periodt)
  recnoise[m,6]<-sd(meanper$periodt)
  
  #governor
  if(mean(meanper$periodt) < 10){
    preylist[[m]] <- prey
  }else
    preylist[[m]] <- NA
   
  }, error=function(e){})
     
}

##############################################################
####################
## View and save the preylist. 
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[!sapply(preylist,is.null)] #remove all null elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long - BUT ONLY IF IT IS MORE THAN 100

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#how many peaks on average?
ntotalpeaks <-recnoises %>% group_by(index)%>%summarize(avpks=mean(meanpeaks))
mean(ntotalpeaks$avpks, na.rm=T)
sd(ntotalpeaks$avpks, na.rm=T)

#what is the period?
mean(recnoises$meanPeriod, na.rm=T)
sd(recnoises$meanPeriod, na.rm=T)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

#just look at the first one, not the mean
preylists%>%pivot_longer(3:22, names_to = 'age_class')%>%
  filter(index==1)%>%
  filter(time_step >= 300 & time_step <=400)%>%
  filter(time_step >= 300 & time_step <=320)%>%
  ggplot(aes(time_step,value))+
  #ggplot(aes(time_step,mean.value, color=age_class))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()

##visualize the data to make sure it's ok. 
preylist.mean <-pivot_longer(preylists,3:22, names_to = 'age_class')%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),.groups='drop')

preylist.mean%>%
  filter(time_step >= 300 & time_step <=400)%>%
  filter(time_step >= 300 & time_step <=320)%>%
  ggplot(aes(time_step,mean.value))+
  #ggplot(aes(time_step,mean.value, color=age_class))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
