####### 1/20/21 ###############
####### Comparison of hierarchical, mixed age and individual age models ###########################
####### devtools::install_github("tanyalrogers/GPEDM") #update frequently

library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(rEDM)
devtools::install_github("tanyalrogers/GPEDM")
library(GPEDM)
library(Metrics)
library(corrplot)
library(pracma)
library(parallel)

######################################################################################################################
################FUNCTIONS ###########################################
### Function to fit GP to the age classes one by one### This is specific to this analysis. 
indv_age <-function(df, maxE){
 ages <-unique(df$age_class)
 fitstats <-list()
 for (m in 1:length(ages)){
  new_df <-filter(df, age_class==ages[m])
  newdfLags = makelags(data=new_df, y="value", E=maxE, tau=1)
  new_df = cbind(new_df,newdfLags)
  new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
  new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
  mod1 <-fitGP(data = new_df.train, y = "value", x=colnames(newdfLags),newdata=new_df.test,
               pop="age_class",scaling = "global",predictmethod = "loo")
  mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,as.character(ages[m]))
  names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","age class")
  fitstats[[m]] <- mod1_out
 }
 fitstats <-bind_rows(fitstats)%>%as.data.frame()%>%mutate(across(OOS_R2:df, as.numeric))
}
######################################################################################################################
################# Mixed age function ###########

##**function for one before one after
#one before, one after#
#### Function to fit GP for one model based on it's own lags and lags of age classes above and below###
mixed_age <-function(df, maxE){
 ages <-unique(df$age_class)
 #maxE <-round(sqrt(dim(df)[1]/length(ages)))
 fitstats <-list()
 for (m in 1:length(ages)){
  #we might not need this if statement. 
  if (m==1){
   new_df <-filter(df, age_class==ages[m] | age_class == ages[m+1])
   
  } else if (m==length(ages)) {
   new_df <-filter(df, age_class==ages[m] | age_class == ages[m+1])
  } else {
   new_df <-filter(df, age_class==ages[m] | age_class==ages[m+1] | age_class==ages[m-1])
  }
  #new_df <-dplyr::select(new_df, -age)%>%
  new_df <-pivot_wider(new_df, names_from = "age_class",values_from = "value")
  #manually scale:
  lastcol <-as.numeric(dim(new_df)[2])
  new_df <- new_df %>% mutate(across(2:all_of(lastcol),scale)) 
  dfLags = makelags(y=new_df[,2:all_of(lastcol)], E=maxE, tau=1)
  dfdata = as.data.frame(cbind(new_df,dfLags))
  #testing and training (not for empirical data)
  #df.train = filter(dfdata, Year < (max(Year)-10))
  #df.test = filter(dfdata, Year >= (max(Year)-10))
  #mod1 = fitGP(data=df.train, y=ages[m], x=colnames(dfLags), newdata=df.test, predictmethod = "loo")
  mod1 <-fitGP(data = dfdata, y = paste(ages[m]), x=colnames(dfLags), predictmethod = "loo")
  mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,as.character(ages[m]))
  names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","age class")
  fitstats[[m]] <- mod1_out
 }
 fitstats <-bind_rows(fitstats)%>%as.data.frame()%>%mutate(across(OOS_R2:df, as.numeric))
}
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
############### FUNCTION TO COMPARE #######################

MIXM30 <- function(plist,maxE){
 plist<-as.data.frame(plist)
 #create the Ntotal age class. 
 plistNT <-plist %>% group_by(time_step)%>%summarize(value=sum(value))%>%mutate(age_class="all") #sum total
 plist <-full_join(plist, plistNT)
 plist <-mutate(plist, value=log(value)) #REMEMBER TO TURN THIS ON AND OFF
 plist <-plist %>%group_by(age_class)%>%arrange(age_class, time_step)#make sure it is sorted by age_class then year
 plist<-filter(plist, age_class !="V21") # we're not doing the plus group ever 

 ### instead of doing 100 years, we're going to do 30. 
 plist2 <-filter(plist, time_step >=300 & time_step < 330)%>% as.data.frame()
 plist3 <-filter(plist2, age_class != "all")%>% as.data.frame()
 #maxE = round(sqrt(30))
 
 #Ntotal model
 plistNT <-filter(plist2, age_class=="all")
 modNT <-fitGP(data = plistNT, y = "value", E=maxE, tau=1, predictmethod = "loo")
 
 #Hierarchical model
 modHier <-fitGP(data = plist3, y = "value", pop="age_class",scaling = "local", E=maxE, tau=1, predictmethod = "loo")
 
 #individual ages
 modIndAge <-indv_age(plist3,maxE)%>%mutate(model="IndvAge")
 
 #mixed ages
 modMixed <-mixed_age(plist3,maxE)
 
 ### Extract fitstats for each age from Hier ####
 outsamp <-bind_cols(modHier$outsampfitstatspop$R2pop,modHier$outsampfitstatspop$rmsepop)
 names(outsamp) <-c("OOS_R2","OOS_rmse")
 insamp <-bind_cols(modHier$insampfitstatspop$R2pop,modHier$insampfitstatspop$rmsepop)
 names(insamp) <-c("R2","rmse")
 fitstats_Hier <-cbind(outsamp, insamp)
 fitstats_Hier <-mutate(fitstats_Hier, model="Hier")%>%rownames_to_column("age")
 
 #extract individual age fitstats from individual age model
 fitstats_IndvAge <- mutate(modIndAge, across(OOS_R2:df, as.numeric))%>%rownames_to_column("age")
 
 #extract mixed age fitstats from mixed age model
 fitstats_Mixed <- mutate(modMixed, across(OOS_R2:df, as.numeric))%>%rownames_to_column("age")%>%mutate(model="Mixed")
 
 fitstats <-bind_rows(fitstats_Hier, fitstats_IndvAge,fitstats_Mixed)
 
 #extract fitstats for Ntotal and Hierarchical_all
 outsamp <-bind_rows(modHier$outsampfitstats,modNT$outsampfitstats)
 names(outsamp) <-c("OOS_R2","OOS_rmse")
 insamp <-bind_rows(modHier$insampfitstats,modNT$insampfitstats)
 rhos <-c(tail(modHier$pars,1), tail(modNT$pars,1))
 fitstats2 <-bind_cols(outsamp,insamp,rhos)%>%as.data.frame()
 rownames(fitstats2) <-c("modHier","modNT")
 colnames(fitstats2)[8]<-"rho"
 fitstats2 <-rownames_to_column(fitstats2, var="model")%>%mutate(`age class`="all", age="all")%>% mutate(across(OOS_R2:df, as.numeric))
 fitstats <-bind_rows(fitstats,fitstats2)
}

################################## IMPORT THE DATA #############################################
preylist1 <-read.csv("Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
  pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()
preylist1 <-split(preylist1, f=preylist1$index)

preylist2 <-read.csv("Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
  pivot_longer(3:23, names_to = "age_class")%>%as.data.frame()%>%filter(age_class !="V21")
preylist2 <-split(preylist2, f=preylist2$index)

preylist3 <-read.csv("Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
  pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()
preylist3 <-split(preylist3, f=preylist3$index)



################################# SIMULATION I ############################################

## This grid does not work for simulation II for some reason, so I am going to use max E =10 and tau = 1. 
## Determine E and tau with a grid
pgrid <-as.data.frame(preylist1[[1]]) %>%filter(age_class !="V21" & time_step >=300 & time_step < 330)
pgrid <-mutate(pgrid, value=log(value))

Ees <-seq(2,10,1)
taus <-seq(1,3,1)
var_pairs = expand.grid(Ees, taus) # Combinations of vars, 2 at a time
ETdf <-matrix(nrow=dim(var_pairs)[1],ncol=4)
ETdf[,1]<-var_pairs[,1]
ETdf[,2]<-var_pairs[,2]
r2matrix1 = array(NA, dim = c(length(Ees), length(taus)), dimnames = list(Ees,taus)) 
rmsematrix1 = array(NA, dim = c(length(Ees), length(taus)), dimnames = list(Ees,taus)) 
for (i in 1:nrow(var_pairs)) {
  try({
    fit1 <-fitGP(data = pgrid, y = "value", pop="age_class",scaling = "local", E=var_pairs[i,1], tau=var_pairs[i,2], predictmethod = "loo")
    fit1_r2 <-fit1$outsampfitstats[[1]]
    fit1_rmse <-fit1$outsampfitstats[[2]]
    r2matrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_r2
    ETdf[i,3] <-fit1_r2
    ETdf[i,4] <-fit1_rmse
    rmsematrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_rmse
  },silent=F)
}
r2matrix1
rmsematrix1
#grab the position of the best E and tau from the matrix. 
bestET <-which(rmsematrix1==min(rmsematrix1,na.rm=T),arr.ind=T)
bestE <-as.numeric(noquote(rownames(bestET)))
bestTau <-as.numeric(bestET[2])
#E=8, Tau=1

############## Test the function
foo <-MIXM30(preylist1[[1]], maxE=10) 

######## Test with the mclapply function
numCores = detectCores()
system.time(parafoo1 <-mclapply(preylist1, MIXM30, maxE=10, mc.cores=numCores))

#### check outputs
length(purrr::keep(parafoo1, is.list)) #how much will it reduce the length by?
parafoo1 <-purrr::keep(parafoo1, is.list) # permanently remove all NA elements

##########bind the list into a DF
Output30_1 <-parafoo1 %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

#############################################################################################
######################## SIMULATION II ###########################################################

## Determine E and tau with a grid
pgrid <-as.data.frame(preylist2[[1]]) %>%filter(age_class !="V21" & time_step >=300 & time_step < 330)
pgrid <-mutate(pgrid, value=log(value))

Ees <-seq(2,10,1)
taus <-seq(1,3,1)
var_pairs = expand.grid(Ees, taus) # Combinations of vars, 2 at a time
ETdf <-matrix(nrow=dim(var_pairs)[1],ncol=4)
ETdf[,1]<-var_pairs[,1]
ETdf[,2]<-var_pairs[,2]
r2matrix2 = array(NA, dim = c(length(Ees), length(taus)), dimnames = list(Ees,taus)) 
rmsematrix2 = array(NA, dim = c(length(Ees), length(taus)), dimnames = list(Ees,taus)) 
for (i in 1:nrow(var_pairs)) {
  try({
    fit1 <-fitGP(data = pgrid, y = "value", pop="age_class",scaling = "local", E=var_pairs[i,1], tau=var_pairs[i,2], predictmethod = "loo")
    fit1_r2 <-fit1$outsampfitstats[[1]]
    fit1_rmse <-fit1$outsampfitstats[[2]]
    r2matrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_r2
    ETdf[i,3] <-fit1_r2
    ETdf[i,4] <-fit1_rmse
    rmsematrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_rmse
  },silent=F)
}
r2matrix2
rmsematrix2


############## Test the function
foo <-MIXM30(preylist2[[1]], maxE=10) 

######## Test with the mclapply function
numCores = detectCores()
system.time(parafoo2 <-mclapply(preylist2, MIXM30,maxE=10, mc.cores=numCores))

#### check outputs
length(purrr::keep(parafoo2, is.list)) #how much will it reduce the length by?
parafoo2 <-purrr::keep(parafoo2, is.list) # permanently remove all NA elements

##########bind the list into a DF
Output30_2 <-parafoo2 %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()%>%
  mutate(Sim="II")

#############################################################################################
######################## SIMULATION III ###########################################################

## Determine E and tau with a grid
pgrid <-as.data.frame(preylist3[[1]]) %>%filter(age_class !="V21" & time_step >=300 & time_step < 330)
pgrid <-mutate(pgrid, value=log(value))

Ees <-seq(2,10,1)
taus <-seq(1,3,1)
var_pairs = expand.grid(Ees, taus) # Combinations of vars, 2 at a time
ETdf <-matrix(nrow=dim(var_pairs)[1],ncol=4)
ETdf[,1]<-var_pairs[,1]
ETdf[,2]<-var_pairs[,2]
r2matrix3 = array(NA, dim = c(length(Ees), length(taus)), dimnames = list(Ees,taus)) 
rmsematrix3 = array(NA, dim = c(length(Ees), length(taus)), dimnames = list(Ees,taus)) 
for (i in 1:nrow(var_pairs)) {
  try({
    fit1 <-fitGP(data = pgrid, y = "value", pop="age_class",scaling = "local", E=var_pairs[i,1], tau=var_pairs[i,2], predictmethod = "loo")
    fit1_r2 <-fit1$outsampfitstats[[1]]
    fit1_rmse <-fit1$outsampfitstats[[2]]
    r2matrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_r2
    ETdf[i,3] <-fit1_r2
    ETdf[i,4] <-fit1_rmse
    rmsematrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_rmse
  },silent=F)
}
r2matrix3
rmsematrix3


############## Test the function
foo <-MIXM30(preylist3[[1]], maxE=10) 

######## Test with the mclapply function
numCores = detectCores()
system.time(parafoo3 <-mclapply(preylist3, MIXM30,maxE=10, mc.cores=numCores))

#### check outputs
length(purrr::keep(parafoo3, is.list)) #how much will it reduce the length by?
parafoo3 <-purrr::keep(parafoo3, is.list) # permanently remove all NA elements

##########bind the list into a DF
Output30_3 <-parafoo3 %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()%>%
  mutate(Sim="III")

MixedAgeOUT <-bind_rows(Output30_1,Output30_2,Output30_3)

write.csv(MixedAgeOUT,"mixedageoutnew.csv")
