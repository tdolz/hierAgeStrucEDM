### New (8/10/22) Parallel Crossfits script ### 
### 
### updated from the version modified on 12/12/21 in the aws folder. but also checked for dates against the aws server. 

# There are three levels of complexity of the simulation. 
# Here we are writing a version of the script to take advantage of parallel cores. 
# The practice parallel script can be found in "practiceparallel.R".
# The original crossfits are found in "SpeciesCrossfitsxxl.R"
# the purpose of this is to run the crossfits 100x each to create confidence intervals on the rsquared and rmse scores.

### 8/20/22 implementing a 5 year testing set (as a test)

#packages
library(MASS)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(rEDM)
#devtools::install_github("tanyalrogers/GPEDM") #, force=T) #update frequently
library(GPEDM)
library(forcats)
library(parallel)
library(pracma)
#############################################################################################################################################################################################################


## Read in the simulation data & make it long format ## BUT TURN IT INTO A LIST!!! 
preylist1 <-read.csv("simulated_data/Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
 pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()
preylist1 <- preylist1 %>% split(f=preylist1$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist2 <-read.csv("simulated_data/Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
 pivot_longer(3:23, names_to = "age_class")%>%as.data.frame()
preylist2 <- preylist2 %>% split(f=preylist2$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist3 <-read.csv("simulated_data/Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
 pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()
preylist3 <- preylist3 %>% split(f=preylist3$index)%>%lapply(function(x) x[!names(x) %in% c("index")])



################################################## THE XFIT FUNCTION. ######################################################################
# should honestly be the same for all of them. 
# updated 11.28.2021 after discussion with Tanya 

Xfit <-function(plist){

   #test set length
   test_length =5
   #test_length =10
   
 #inputs
 plist<-as.data.frame(plist)
 plist <-mutate(plist, value=log(value)) #REMEMBER TO TURN THIS ON AND OFF
 plist <-plist %>%group_by(age_class)%>%arrange(age_class, time_step)#make sure it is sorted by age_class then year
 plist<-filter(plist, age_class !="V21") # we're not doing the plus group ever 
 agelist <-unique(plist$age_class)
 agelist <-agelist[c(1,12,14,15,16,17,18,19,20,2,3,4,5,6,7,8,9,10,11,13)] #NEW! sort this 
 tslengthlist <-seq(test_length+5,(30+test_length),1) # test set length
 var_pairs=expand.grid(agelist, tslengthlist)
 var_pairs = as.matrix(t(var_pairs))
 var_pairs[2,]=gsub(" ", "", var_pairs[2,], fixed = TRUE)
 
 #outputs
 rho_matrix = array(NA, dim = c(length(agelist), length(tslengthlist)), dimnames = list(agelist,tslengthlist)) 
 E_matrix = array(NA, dim = c(length(agelist), length(tslengthlist)), dimnames = list(agelist,tslengthlist)) 
 phi_matrix = array(NA, dim = c(length(agelist), length(tslengthlist)), dimnames = list(agelist,tslengthlist)) 
 maxphi_matrix = array(NA, dim = c(length(agelist), length(tslengthlist)), dimnames = list(agelist,tslengthlist)) 
 rmse_matrix = array(NA, dim = c(length(agelist), length(tslengthlist)), dimnames = list(agelist,tslengthlist)) 
 dyrho_matrix = array(NA, dim = c(length(agelist), length(tslengthlist)), dimnames = list(agelist,tslengthlist))
 obs_sd = array(NA, dim = c(length(agelist), length(tslengthlist)), dimnames = list(agelist,tslengthlist))
 
 for (m in 1:length(var_pairs[1,])) {
  df <-filter(plist, time_step > 299)%>% #WHEN THE SIMULATION HAS STABILIZED!
   mutate(time_step=time_step-299)%>%
   filter(time_step <= as.numeric(var_pairs[2,m]))
  f <-which(agelist==var_pairs[1,m])
  df <-filter(df, age_class %in% agelist[1:f])
  df <-mutate(df, value=as.numeric(as.character(value)))%>%as.data.frame()
  
  #### MAX E RULES #########
  #OPTION 1 (new)
  realts <-length(unique(df$time_step))-test_length
  maxE=round(sqrt(realts))
  
  #OPTION 2 (new)
  #realts <-length(unique(df$time_step))-10
  #if (round(sqrt(realts)) < 3) {
  #maxE = 3
  #}else{maxE = round(sqrt(realts))}
  #maxE =3 
  
  try({
   #fit1 <-fitGP(data = df, yd = "value", pop="age_class",scaling = "local", E=maxE, tau=1, predictmethod = "loo")
   #df.train <-filter(df, time_step < (max(time_step)-10)) #train up until the last 10 years.
   #df.test <-filter(df, time_step >= (max(time_step)-10)) #test on last 10 years
   #fit1 <-fitGP(data=df.train, yd="value", E=maxE, tau=1, pop="age_class", datanew=df.test)
   
   #NEW with MAKE LAGS 
   dfLags = makelags(data=df, y="value", pop="age_class", E=maxE, tau=1)
   dfdata = cbind(df,dfLags)
   df.train = filter(dfdata, time_step < (max(time_step)-test_length))
   df.test = filter(dfdata, time_step >= (max(time_step)-test_length))
   fit1 = fitGP(data=df.train, y="value", x=colnames(dfLags), pop="age_class", newdata=df.test, predictmethod = "loo", scaling="local")
   
   #extract rmse and r-squared
   fitr <-round(fit1$outsampfitstats[[1]],4)
   fitrmse <-round(fit1$outsampfitstats[[2]],4)
   obsSD <-sd(df$value)
   
   #find the number of embedding dimenstions used by ARD
   phis <-fit1$pars[1:(length(fit1$pars)-3)]
   relevant_phis <-phis[which(phis>=0.01)]
   #total number of lags used by ARD
   lengthphi <-length(relevant_phis)
   #the maximum lag used by ARD. 
   maxphi <-max(as.numeric(as.character(gsub("phi","",names(relevant_phis)))))
   #the dynamic correlation
   dyrho <-as.numeric(tail(fit1$pars,1))
   
   rmse_matrix[var_pairs[1,m], var_pairs[2,m]]<-fitrmse
   maxphi_matrix[var_pairs[1,m], var_pairs[2,m]]<-maxphi
   phi_matrix[var_pairs[1,m], var_pairs[2,m]]<-lengthphi
   rho_matrix[var_pairs[1,m], var_pairs[2,m]] <-fitr
   E_matrix[var_pairs[1,m], var_pairs[2,m]] <-maxE
   dyrho_matrix[var_pairs[1,m], var_pairs[2,m]] <-dyrho
   obs_sd[var_pairs[1,m], var_pairs[2,m]] <-obsSD
  },silent=T)
 }
 
 
 #unmatrix
 rhoD <-data.frame(col=colnames(rho_matrix)[col(rho_matrix)], row=rownames(rho_matrix)[row(rho_matrix)], dist=c(rho_matrix))
 names(rhoD)<-c("time_step","age_class","Rsq")
 eD <-data.frame(col=colnames(E_matrix)[col(E_matrix)], row=rownames(E_matrix)[row(E_matrix)], dist=c(E_matrix))
 names(eD)<-c("time_step","age_class","E")
 maxphim <-data.frame(col=colnames(maxphi_matrix)[col(maxphi_matrix)], row=rownames(maxphi_matrix)[row(maxphi_matrix)], dist=c(maxphi_matrix))
 names(maxphim)<-c("time_step","age_class","max_phi")
 rmsem <-data.frame(col=colnames(rmse_matrix)[col(rmse_matrix)], row=rownames(rmse_matrix)[row(rmse_matrix)], dist=c(rmse_matrix))
 names(rmsem)<-c("time_step","age_class","rmse")
 dyrhom <-data.frame(col=colnames(dyrho_matrix)[col(dyrho_matrix)], row=rownames(dyrho_matrix)[row(dyrho_matrix)], dist=c(dyrho_matrix))
 names(dyrhom)<-c("time_step","age_class","dy_rho")
 phism <-data.frame(col=colnames(phi_matrix)[col(phi_matrix)], row=rownames(phi_matrix)[row(phi_matrix)], dist=c(phi_matrix))
 names(phism)<-c("time_step","age_class","num_phis")
 obssd <-data.frame(col=colnames(obs_sd)[col(obs_sd)], row=rownames(obs_sd)[row(obs_sd)], dist=c(obs_sd))
 names(obssd)<-c("time_step","age_class","obs_SD")
 rhoD <-full_join(rhoD, eD)%>%full_join(maxphim)%>%full_join(rmsem)%>%full_join(phism)%>%full_join(dyrhom)%>%full_join(obssd)
 
 #output
 rhoD
}

####################################################################################################################################################################################################

#############TEST TEST TEST TEST ####################

#Test Function
shortlist <-list(preylist1[[1]],preylist1[[2]])
numCores = detectCores()
system.time(foo <-mclapply(shortlist, Xfit, mc.cores=numCores))
head(foo)

foo = purrr::keep(foo, is.data.frame)
for(t in 1:length(foo)){
 foo[[t]] <-mutate(foo[[t]], iter=t)
}
testfoo<-bind_rows(foo)%>%as.data.frame()
#write.csv(testfoo,"testfoo.csv")

############################# SIM I ############################################

##Apply##
numCores = detectCores()
system.time(
 paraXfit <-mclapply(preylist1, Xfit, mc.cores=numCores) 
)

##Process##
#remove ones that didn't work from the list
# first view and write down the ones that didn't work.
#paraXfit = paraXfit[-c(19,45,67,93)]
paraXfit = purrr::keep(paraXfit, is.data.frame)
for(t in 1:length(paraXfit)){
 paraXfit[[t]] <-mutate(paraXfit[[t]], iter=t)
}
crossfits1<-bind_rows(paraXfit)%>%as.data.frame()

##Write to CSV##
write.csv(crossfits1, "tslength_ageclass_outputs/crossfits_sim1.csv") # 10 year testing set
#write.csv(crossfits1, "crossfits_sim1_5test.csv") # 5 year testing set

############################# SIM II ############################################

##Apply##
numCores = detectCores()
system.time(
 paraXfit <-mclapply(preylist2, Xfit, mc.cores=numCores) 
)

##Process##
#remove ones that didn't work from the list
# first view and write down the ones that didn't work.
#paraXfit = paraXfit[-c(19,45,67,93)]
paraXfit = purrr::keep(paraXfit, is.data.frame)
for(t in 1:length(paraXfit)){
 paraXfit[[t]] <-mutate(paraXfit[[t]], iter=t)
}
crossfits2<-bind_rows(paraXfit)%>%as.data.frame()

##Write to CSV##
write.csv(crossfits2, "tslength_ageclass_outputs/crossfits_sim2.csv") 

############################# SIM III ############################################

##Apply##
numCores = detectCores()
system.time(
 paraXfit <-mclapply(preylist3, Xfit, mc.cores=numCores) 
)

##Process##
#remove ones that didn't work from the list
# first view and write down the ones that didn't work.
#paraXfit = paraXfit[-c(19,45,67,93)]
paraXfit = purrr::keep(paraXfit, is.data.frame)
for(t in 1:length(paraXfit)){
 paraXfit[[t]] <-mutate(paraXfit[[t]], iter=t)
}
crossfits3<-bind_rows(paraXfit)%>%as.data.frame()

##Write to CSV##
write.csv(crossfits3, "tslength_ageclass_outputs/crossfits_sim3.csv") 




