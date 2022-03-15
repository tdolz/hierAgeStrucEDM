####### 1/20/21 ###############
####### Comparison of hierarchical, mixed age and individual age models ###########################
####### devtools::install_github("tanyalrogers/GPEDM") #update frequently

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

######################################################################################################################
################FUNCTIONS ###########################################
### Function to fit GP to the age classes one by one### This is specific to this analysis. 
indv_age <-function(df, maxE){
 ages <-unique(df$age_class)
 fitstats <-list()
 for (m in 1:length(ages)){
  new_df <-filter(df, age_class==ages[m])
  newdfLags = makelags(data=new_df, yd="value", E=maxE, tau=1)
  new_df = cbind(new_df,newdfLags)
  new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
  new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
  mod1 <-fitGP(data = new_df.train, yd = "value", xd=colnames(newdfLags),datanew=new_df.test,
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
  dfLags = makelags(yd=new_df[,2:all_of(lastcol)], E=maxE, tau=1)
  dfdata = as.data.frame(cbind(new_df,dfLags))
  #testing and training (not for empirical data)
  #df.train = filter(dfdata, Year < (max(Year)-10))
  #df.test = filter(dfdata, Year >= (max(Year)-10))
  #mod1 = fitGP(data=df.train, yd=ages[m], xd=colnames(dfLags), datanew=df.test, predictmethod = "loo")
  mod1 <-fitGP(data = dfdata, yd = paste(ages[m]), xd=colnames(dfLags), predictmethod = "loo")
  mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,as.character(ages[m]))
  names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","age class")
  fitstats[[m]] <- mod1_out
 }
 fitstats <-bind_rows(fitstats)%>%as.data.frame()%>%mutate(across(OOS_R2:df, as.numeric))
}
######################################################################################################################
######## count peaks function ###############
count_peaks <-function(x){
 mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak. 
 x <-findpeaks(x, minpeakheight = mph)
 x <-dim(x)[1]
 x[is.null(x)] <-0 #if there are no peaks return 0. 
 x
}
######################################################################################################################
######################################################################################################################
######################################################################################################################
############### FUNCTION TO COMPARE #######################
plist <-preylist[[1]]

MIXM30 <- function(plist,maxE){
 plist<-as.data.frame(plist)
 #create the Ntotal age class. 
 plistNT <-plist %>% group_by(time_step)%>%summarize(value=sum(value))%>%mutate(age_class="all")
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
 modNT <-fitGP(data = plistNT, yd = "value", E=maxE, tau=1, predictmethod = "loo")
 
 #Hierarchical model
 modHier <-fitGP(data = plist3, yd = "value", pop="age_class",scaling = "local", E=maxE, tau=1, predictmethod = "loo")
 
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

################################ TWO SPP XXL LOOP ####################################################################
########################### CREATE THE PREYLISTS ###################################################################################################

maxiter = 2 #try making 200 of them and deleting all the zero peak one.s 
preylist<-list()
predlist<-list()
count0peaks <-list()
recnoise <-matrix(nrow=maxiter, ncol=4)
colnames(recnoise)<-c("recnoise1","recnoise2","meanpeaks","varprey")


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
 
 #prey variance
 preyvar <-preypiv %>%group_by(age_class)%>%summarize(varprey=var(value))
 
 #store recnoises and mean peaks
 recnoise[m,1]<-recnoise1
 recnoise[m,2]<-recnoise2
 recnoise[m,3]<-mean(countpeaks$pks)
 recnoise[m,4]<-mean(preyvar$varprey)
 
 
 #if(mgroups > 0.1 & mpreds > 0.1 & mean(countpeaks$pks) > 1.2) {
 if(mean(countpeaks$pks) > 1.2) {
  preylist[[m]] <- preypiv
 }else
  preylist[[m]] <- NA
}



##################################################################################################################################

## View and save the preylist. 
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)
preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()%>%full_join(recnoises)

#how many peaks on average?
ntotalpeaks <-preylists %>% group_by(index)%>%summarize(avpks=mean(meanpeaks))
mean(ntotalpeaks$avpks)
sd(ntotalpeaks$avpks)

## Determine E and tau with a grid
pgrid <-as.data.frame(preylist[[1]]) %>%filter(age_class !="V21" & time_step >=300 & time_step < 330)

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
  fit1 <-fitGP(data = pgrid, yd = "value", pop="age_class",scaling = "local", E=var_pairs[i,1], tau=var_pairs[i,2], predictmethod = "loo")
  fit1_r2 <-fit1$outsampfitstats[[1]]
  fit1_rmse <-fit1$outsampfitstats[[2]]
  r2matrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_r2
  ETdf[i,3] <-fit1_r2
  ETdf[i,4] <-fit1_rmse
  rmsematrix1[var_pairs[i,1], var_pairs[i,2]] = fit1_rmse
 },silent=T)
}
r2matrix1
rmsematrix1
# best is Tau = 1, E =6, but the real best is Tau =3, E=9 (Throws out a lot of data)
#using an E=5 is not substantially worse. 

#write.csv(preylists,"preylists2sppXXLSHORTOSC_mixed100.csv")

##################################################################################################################################
##### Apply the function ########

############## Test the function
foo <-MIXM30(preylist[[1]], maxE=5) 

######## Test with the mclapply function #########
numCores = detectCores()
system.time(parafoo <-mclapply(preylist, MIXM30, mc.cores=numCores))

#### check outputs
length(purrr::keep(parafoo, is.list)) #how much will it reduce the length by?
parafoo <-purrr::keep(parafoo, is.list) # permanently remove all NA elements

################# bind the list into a DF #################################

Output30 <-parafoo %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

