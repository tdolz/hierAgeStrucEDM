#two species XXL loop R. 
#install.packages("devtools") #if required
devtools::install_github("tanyalrogers/GPEDM") #update frequently

### updated 3/14/22

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

################FUNCTIONS ###########################################
### Function to fit GP to the age classes one by one### This is specific to this analysis. 
indv_age <-function(df){
  ages <-unique(df$age_class)
  maxE <-10 # in all cases it's going to be 10 because the ts is always 100 points long for this one. 
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

####### Function to create the aggregate individual age prediction #######################
## basically, make a prediction for each single age and add up all the single age and 

aggPred <-function(df,ts){ #basically the indv_age function, but you want to be able to specify time series length. 
  max_time <-min(df$time_step)+ts+10
  df <-filter(df, time_step <=max_time)
  ages <-unique(df$age_class)
  maxE <-round(sqrt(length(unique(df$time_step))-10)) 
  mod1res <-list()
  for (m in 1:length(ages)){
    new_df <-filter(df, age_class==ages[m])
    newdfLags = makelags(data=new_df, yd="value", E=maxE, tau=1)
    new_df = cbind(new_df,newdfLags)
    new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
    new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
    mod1 <-fitGP(data = new_df.train, yd = "value", xd=colnames(newdfLags),datanew=new_df.test,
                 pop="age_class",scaling = "global",predictmethod = "loo")
    res <-as.data.frame(mod1$outsampresults)%>%mutate(age=ages[m])
    mod1res[[m]] <-res
  }
  modres <-bind_rows(mod1res)%>%as.data.frame() #bind all the predicted and observed
  modsum <-modres %>% group_by(timestep) %>%summarize(aggpred =sum(predmean), aggobs=sum(obs)) #summarize by time step to get an aggregate predicted and observed
  aggR2 <-mutate(modsum, ssr=(aggobs - aggpred)^2, sst=(aggobs-mean(modsum$aggobs))^2)
  r2 <-1-sum(aggR2$ssr)/sum(aggR2$sst)
  r2
}


################################ ntotal years. 
ntotal_yrs <-function(df){
  ntotalstats <-list()
  yrs <-c(15,20,30,35)
  startime <-min(df$time_step)
  for (m in 1:length(yrs)){
    new_df <-filter(df, time_step <= startime+yrs[m])
    maxE=round(sqrt(yrs[m]-10))
    newdfLags = makelags(data=new_df, yd="value", E=maxE, tau=1)
    new_df = cbind(new_df,newdfLags)
    new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
    new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
    mod1 <-fitGP(data = new_df.train, yd = "value", xd=colnames(newdfLags),datanew=new_df.test,
                 pop="age_class",scaling = "local",predictmethod = "loo")
    mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,as.character(yrs[m]))
    names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","tslength")
    ntotalstats[[m]] <- mod1_out
  }
  ntotalstats <-bind_rows(ntotalstats)%>%as.data.frame()%>%mutate(across(OOS_R2:df, as.numeric))
}

################################ Apples to Apples###############################
ntotal_app<-function(df){
  newdf <-df %>% group_by(time_step)%>%summarize(NewValue=sum(value))
  tslength <-(max(newdf$time_step)-min(newdf$time_step))-10
  maxE=round(sqrt(tslength))
  ages=n_distinct(df$age_class)
  newdfLags = makelags(data=newdf, yd="NewValue", E=maxE, tau=1)
  new_df = cbind(newdf,newdfLags)
  new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
  new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
  mod1 <-fitGP(data = new_df.train, yd = "NewValue", xd=colnames(newdfLags),datanew=new_df.test,
               scaling = "global",predictmethod = "loo")
  mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,tslength,ages)
  names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","tslength","age class")
  mod1_out
  
}

######## count peaks function ###############
count_peaks <-function(x){
  mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak. 
  x <-findpeaks(x, minpeakheight = mph)
  x <-dim(x)[1]
  x[is.null(x)] <-0 #if there are no peaks return 0. 
  x
}

################ find period function ###############
period <-function(x){
  mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak. 
  x <-findpeaks(x, minpeakheight = mph)
  x <-x[,1:2]
  x <-as.data.frame(x)
  names(x) <-c("ypeak","xpeak")
  x <-mutate(x, dist = xpeak-lag(xpeak,1))
  mean.period = mean(x$dist,na.rm=T)
  mean.period
}
###################################################################
######################################################################################################################


#############THE GP FUNCTION############################################################

GPLOOP <-function(plist){
  
  #inputs
  plist<-as.data.frame(plist)
  #create the Ntotal age class. 
  plistNT <-plist %>% group_by(time_step)%>%summarize(NewValue=sum(value))%>%mutate(age_class="all", value=log(NewValue))
  plist <-full_join(plist, plistNT)
  plist <-mutate(plist, value=log(value)) #REMEMBER TO TURN THIS ON AND OFF
  plist <-plist %>%group_by(age_class)%>%arrange(age_class, time_step)#make sure it is sorted by age_class then year
  plist<-filter(plist, age_class !="V21") # we're not doing the plus group ever 
  
  agelist <-unique(plist$age_class)
  agelist <-agelist[c(2,13,15,16,17,18,19,20,21,3,4,5,6,7,8,9,10,11,12,14,1)]
  
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
  # Total abundance (100 years)
  N_total_final <-fitGP(data = preyNT.train, yd = "value", xd=colnames(preyNTLags),datanew=preyNT.test,predictmethod = "loo")
  
  # 20ages (5,10,20,25 years)
  twentyagesyrs <-ntotal_yrs(preypivall)%>%mutate(model="20agesdiffyrs", `age class`="all", tslength=as.numeric(tslength)-10, approach="hier")
  
  #Ntotal abundance, different years -DISCONTINUED
  #preypivNT <-plistNT %>%filter(time_step > 299 & time_step < 410)
  #Ntotsyrs <-ntotal_yrs(preypivNT)%>%mutate(model="NT_diffyears", `age class`="N_total", tslength=as.numeric(tslength)-10, approach="TAindex")
  
  # 20 years of 5 age classes 
  prey520_final <-fitGP(data = prey5.train, yd = "value", xd=colnames(prey5Lags),datanew=prey5.test,pop="age_class",scaling = "local",predictmethod = "loo")
  #now each age class individually (make sure ages are in the correct order. )
  preypivall520 <-filter(preypivall, age_class %in% c("V1","V2","V3","V4","V5"))%>%as.data.frame()
  prey520_ages <-indv_age(preypivall520)%>%mutate(model="prey520", tslength=100, approach="single_age")
  prey520agg <-aggPred(preypivall520,20) # the aggregate r2 for 20 years. 
  #aggregate the hierarchical results "agg hier"
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
  outsamp <-bind_rows(N_total_final$outsampfitstats,prey205_final$outsampfitstats,prey1010_final$outsampfitstats,prey520_final$outsampfitstats,prey425_final$outsampfitstats)
  names(outsamp) <-c("OOS_R2","OOS_rmse")
  insamp <-bind_rows(N_total_final$insampfitstats,prey205_final$insampfitstats,prey1010_final$insampfitstats,prey520_final$insampfitstats,prey425_final$insampfitstats)
  rhos <-c(tail(N_total_final$pars,1),tail(prey205_final$pars,1),tail(prey1010_final$pars,1),tail(prey520_final$pars,1),tail(prey425_final$pars,1))
  
  fitstats <-bind_cols(outsamp,insamp,rhos)%>%as.data.frame()
  rownames(fitstats) <-c("N_total","5yrs20ages","10yrs10ages","20yrs5ages","25yrs4ages")
  colnames(fitstats)[8]<-"rho"
  fitstats$tslength <-c(100,5,10,20,25)
  fitstats$approach <-c("TAindex","hier","hier","hier","hier")
  
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
  fitstats_ages <-bind_rows(prey205_ages,prey1010_ages,prey520_ages,prey425_ages)%>%mutate(tslength=100)%>%bind_rows(twentyagesyrs)%>%bind_rows(Ntotsyrs)
  
  pars <-bind_rows(N_total_final$pars,prey205_final$pars,prey1010_final$pars,prey520_final$pars,prey425_final$pars)%>%as.data.frame()
  rownames(pars) <-c("N_total","5yrs20ages","10yrs10ages","20yrs5ages","25yrs4ages")
  pars$tslength <-c(100,5,10,20,25)
  pars$approach <-c("TAindex","hier","hier","hier","hier")
  pars <-mutate(pars, maxE=round(sqrt(tslength)))
  
  #outputs
  gplout <-list(pars,fitstats,fitstats_ages,fitagg,fitsumhier)
  names(gplout) <-c("pars","fitstats","fitstats_ages","fitagg","fitsumhier")
  gplout
}
################################ TWO SPP XXL SHORT OSC PREYLISTS #####################################################

########################### CREATE THE PREYLISTS ###################################################################################################

maxiter = 150. #try making 200 of them and deleting all the zero peak one.s 
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
  
  
  #if(mgroups > 0.1 & mpreds > 0.1 & mean(countpeaks$pks) > 1.2) {
  if(mean(countpeaks$pks) > 1.2) {
    preylist[[m]] <- preypiv
  }else
    preylist[[m]] <- NA
}



##################################################################################################################################
##################################################################################################################################
##################################################################################################################
####################
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

#what is the period?
mean(preylists$meanPeriod)
sd(preylists$meanPeriod)


write.csv(preylists,"preylists2sppXXLSHORTOSC_GPLOOP100.csv")


##############################################################
############## Test the function
foo <-GPLOOP(preylist[[1]]) #it works.

#### shortlist test ###########
shortlist <-preylist[c(1,2)]
numCores = detectCores()
system.time(parafoo <-mclapply(shortlist, GPLOOP, mc.cores=numCores))


######## Test with the mclapply function #########
numCores = detectCores()
system.time(parafoo <-mclapply(preylist, GPLOOP, mc.cores=numCores))

#### some of them did not deliver... I wonder why?
length(purrr::keep(parafoo, is.list)) #how much will it reduce the length by?
parafoo <-purrr::keep(parafoo, is.list) # permanently remove all NA elements


### corralling the elements from the list output##################### 

#pars#
qpars <-lapply(parafoo, `[`, "pars")
for (i in 1:length(qpars)){
  qpars[[i]] <- mutate(qpars[[i]]$pars, iter=i)%>%rownames_to_column()
}
qpars <-qpars %>% bind_rows()%>%as.data.frame()

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


#combine fitstats and fitages to newfits#
newfits <-bind_rows(fstats,fitages,fitagg,fitsumhier)

#######Save to csv############################### 
write.csv(qpars, file="TwoSppXLPhisLN100.csv")
write.csv(newfits, file="TwoSppXLFitstatsLN100.csv")

###########################################################################################################################
############################### TWO SPECIES MED ############################################################################


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
    preylist[[m]] <- preypiv
  }else
    preylist[[m]] <- NA
  
  
}


##################################################################################################################################
##################################################################################################################################
##################################################################################################################
####################
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

#what is the period?
mean(preylists$meanPeriod)
sd(preylists$meanPeriod)


write.csv(preylists,"preylists2sppMEDGPLOOP100.csv")


##############################################################
############## Test the function
system.time(foo <-GPLOOP(preylist[[1]])) #it works.

######## Test with the mclapply function #########
numCores = detectCores()
system.time(parafoo <-mclapply(preylist, GPLOOP, mc.cores=numCores))


### corralling the elements from the list output##################### 
#pars#
qpars <-lapply(parafoo, `[`, "pars")
for (i in 1:length(qpars)){
  qpars[[i]] <- mutate(qpars[[i]]$pars, iter=i)%>%rownames_to_column()
}
qpars <-qpars %>% bind_rows()%>%as.data.frame()

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


#combine fitstats and fitages to newfits#
newfits <-bind_rows(fstats,fitages,fitagg,fitsumhier)


#######Save to csv############################### 
write.csv(qpars, file="TwoSppMEDPhisLN100.csv")
write.csv(newfits, file="TwoSppMEDFitstatsLN100.csv")

##################################################### ONE SPECIES #####################################################################################
#######################################################################################################################################################


maxiter = 150. #try making 200 of them and deleting all the zero peak one.s 
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
  
  #measure how many peaks in a ten year time period.
  countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
    dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
  count0peaks[[m]] <-sum(countpeaks$pks==0)
  
  
  if(sum(countpeaks$pks==0)< 15) { #we can't get them all in there unfortunately. 
    preylist[[m]] <- preypiv
  }else
    preylist[[m]] <- NA
  
 
  
}


##################################################################################################################
####################
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

#what is the period?
mean(preylists$meanPeriod)
sd(preylists$meanPeriod)


write.csv(preylists,"preylistsONESPPGPLOOP100.csv")


##############################################################
############## Test the function
system.time(foo <-GPLOOP(preylist[[1]])) #it works.

######## Test with the mclapply function #########
numCores = detectCores()
system.time(parafoo <-mclapply(preylist, GPLOOP, mc.cores=numCores))


### corralling the elements from the list output##################### 

#pars#
qpars <-lapply(parafoo, `[`, "pars")
for (i in 1:length(qpars)){
  qpars[[i]] <- mutate(qpars[[i]]$pars, iter=i)%>%rownames_to_column()
}
qpars <-qpars %>% bind_rows()%>%as.data.frame()

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


#combine fitstats and fitages to newfits#
newfits <-bind_rows(fstats,fitages,fitagg,fitsumhier)

#######Save to csv############################### 
write.csv(qpars, file="ONEsppPhisLN100.csv")
write.csv(newfits, file="OneSppFitstatsLN100.csv")


