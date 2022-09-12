#### MAKE SIMULATION DATA ######


### created 3/18/22
## updated 3/24/22 Great job!
## UPDATED 4/28/22. Original values can be found in "simulation_sandbox.R"
## updated 5/20/22 from "comparison of simulations for steve.Rmd" 


source("analysis_functions.R")
devtools::install_github("tanyalrogers/GPEDM")
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
library(grid)
library(gridExtra)
library(cowplot)

##################################################### Simulation I #####################################################################################
##################################################### ONE SPECIES ######################################################################################

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

#some data formatting
preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()
preypiv <- preylists%>%pivot_longer(3:22, names_to = 'age_class')
preypiv$age_class = fct_relevel(preypiv$age_class, c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20"))
preypiv$age_class = fct_recode(preypiv$age_class, "age 1"="V1","age 2"="V2","age 3"="V3","age 4"="V4","age 5"="V5","age 6"="V6", "age 7"="V7", "age 8"="V8","age 9"="V9",
                               "age 10"="V10","age 11"="V11","age 12"="V12","age 13"="V13", "age 14"="V14",'age 15'='V15',"age 16"='V16','age 17'= "V17",'age 18'="V18","age 19"="V19","age 20"="V20")

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#what is the period?
mean(recnoises$meanPeriod, na.rm=T)
sd(recnoises$meanPeriod, na.rm=T)
median(recnoises$meanPeriod, na.rm=T)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

##summarize across runs, mean and sd.. 
preylist.mean <-preypiv%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),sd.value=sd(value),.groups='drop')

#just look at the first one, not the mean
preypiv%>%
  filter(index==1)%>%
  filter(time_step >= 300 & time_step <=400)%>%
  #filter(time_step >= 300 & time_step <=320)%>%
  ggplot(aes(time_step,value))+
  #ggplot(aes(time_step,mean.value, color=age_class))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim1_run1.png")
#ggsave("sim1_run1.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

preylist.mean%>%
 ggplot(aes(time_step,mean.value))+
 geom_line()+
  geom_ribbon(aes(ymin=mean.value-sd.value, ymax=mean.value+sd.value), fill="grey", alpha=0.5)+
 facet_wrap(~age_class,scales="free")+
 theme_classic()
ggsave("sim_diagnostics/sim1_pre_burnin.png")
#ggsave("sim1_pre_burnin.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

preylist.mean%>%
 filter(time_step >= 300 & time_step <= 400)%>%
 ggplot(aes(time_step,mean.value))+
 geom_line()+
  geom_ribbon(aes(ymin=mean.value-sd.value, ymax=mean.value+sd.value), fill="grey", alpha=0.5)+
 facet_wrap(~age_class,scales="free")+
 theme_classic()
ggsave("sim_diagnostics/sim1_post_burnin.png")
#ggsave("sim1_post_burnin.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()


#grob plot of each age by itself the previous year.  
agelist <-unique(preypiv$age_class)
plotlist <-list()
for (i in 2:length(agelist)){
  
  df <-preypiv%>%
    filter(time_step >= 299 & time_step <=401)%>%
    filter(age_class==agelist[i] | age_class==agelist[i-1])%>%
    pivot_wider(id_cols=1:2, names_from = age_class, values_from = value)
  
  prev_yrcls_lag=lag(df[,3],1)
  df<-as.data.frame(cbind(df,prev_yrcls_lag))
  names(df)<-c("index","time_step","yrcls","yrPlus1","yrlag")
  df <-dplyr::select(df, -yrcls)%>%filter(time_step >= 300 & time_step <=400)
  
  p <-ggplot(df,aes(yrlag,yrPlus1, color=as.factor(index)))+
    geom_line()+
    ylab(paste(agelist[i],"(t)"))+xlab(paste(agelist[i-1], "(t-1)"))+
    guides(color="none")+
    theme_classic()
  
  plotlist[[i]]<-p
}
plotlist<-plotlist[2:20]

row1 <-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],ncol=5)
row2 <-plot_grid(plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],ncol=5)
row3 <-plot_grid(plotlist[[11]],plotlist[[12]],plotlist[[13]],plotlist[[14]],plotlist[[15]],ncol=5)
row4 <-plot_grid(plotlist[[16]],plotlist[[17]],plotlist[[18]],plotlist[[19]],ncol=5)

y.grob <-textGrob("Year Class (t)", gp=gpar(fontface="bold",col="black",fontsize=12),rot=90)
x.grob <-textGrob("Prev Year Class (t-1)", gp=gpar(fontface="bold",col="black",fontsize=12))
top.grob <-textGrob("Year class vs. itself", gp=gpar(fontface="bold",col="black",fontsize=12))

mainplot <-plot_grid(row1,row2,row3,row4, nrow=4)
grid.arrange(arrangeGrob(mainplot, left=y.grob, bottom=x.grob,top=top.grob, padding=unit(0,"line")))
ggsave("sim_diagnostics/sim1_ycvsitself.png")

####NO LAG###
#grob plot of each age by itself the previous year.  
agelist <-unique(preypiv$age_class)
plotlist <-list()
for (i in 2:length(agelist)){
  
  df <-preypiv%>%
    filter(time_step >= 300 & time_step <=400)%>%
    filter(age_class==agelist[i] | age_class==agelist[i-1])%>%
    pivot_wider(id_cols=1:2, names_from = age_class, values_from = value)
  names(df)<-c("index","time_step","yrcls","yrPlus1")
  
  p <-ggplot(df,aes(yrcls,yrPlus1, color=as.factor(index)))+
    geom_line()+
    ylab(paste(agelist[i],"(t)"))+xlab(paste(agelist[i-1], "(t)"))+
    guides(color="none")+
    theme_classic()
  
  plotlist[[i]]<-p
}
plotlist<-plotlist[2:20]
#plotlist

row1 <-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],ncol=5)
row2 <-plot_grid(plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],ncol=5)
row3 <-plot_grid(plotlist[[11]],plotlist[[12]],plotlist[[13]],plotlist[[14]],plotlist[[15]],ncol=5)
row4 <-plot_grid(plotlist[[16]],plotlist[[17]],plotlist[[18]],plotlist[[19]],ncol=5)

y.grob <-textGrob("age i (t)", gp=gpar(fontface="bold",col="black",fontsize=12),rot=90)
x.grob <-textGrob("age i-1 (t)", gp=gpar(fontface="bold",col="black",fontsize=12))
top.grob <-textGrob("age vs previous age", gp=gpar(fontface="bold",col="black",fontsize=12))

mainplot <-plot_grid(row1,row2,row3,row4, nrow=4)
grid.arrange(arrangeGrob(mainplot, left=y.grob, bottom=x.grob,top=top.grob, padding=unit(0,"line")))
ggsave("sim_diagnostics/sim1_ycvsprevyr.png")



preylists <-preylists %>%
 filter(time_step > 299)

write.csv(preylists,"simulated_data/Simulation1_data.csv")
write.csv(recnoises,"simulated_data/sim1_info.csv")


###############################Simulation II############################################################################
############################### formerly known as TWO SPECIES MED ############################################################################
set.seed(4649)
maxiter = 1000. 
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
    
    
    if(mean(meanper$periodt)<12) { #we can't get them all in there unfortunately. 
      preylist[[m]] <- prey
    }else
      preylist[[m]] <- NA
    
    #no governor
    #preylist[[m]] <- prey
    
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

#some data formatting
preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()
preypiv <- preylists%>%pivot_longer(3:23, names_to = 'age_class')
preypiv$age_class = fct_relevel(preypiv$age_class, c("V1","V2","V3","V4","V5","V6","V7","V8",                                                     "V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20"))
preypiv$age_class = fct_recode(preypiv$age_class, "age 1"="V1","age 2"="V2","age 3"="V3","age 4"="V4","age 5"="V5","age 6"="V6", "age 7"="V7", "age 8"="V8","age 9"="V9","age 10"="V10","age 11"="V11","age 12"="V12","age 13"="V13", "age 14"="V14",'age 15'='V15',"age 16"='V16','age 17'= "V17",'age 18'="V18","age 19"="V19","age 20"="V20", "age 21"="V21")


## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#what is the period?
mean(recnoises$meanPeriod, na.rm=T)
sd(recnoises$meanPeriod, na.rm=T)
median(recnoises$meanPeriod, na.rm=T)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

##visualize the data to make sure it's ok. 
preylist.mean <-preypiv%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),sd.value=sd(value),.groups='drop')

#just look at the first one, not the mean
preypiv%>%
  filter(index==1)%>%
  filter(time_step >= 300 & time_step <=400)%>%
  #filter(time_step >= 300 & time_step <=320)%>%
  ggplot(aes(time_step,value))+
  #ggplot(aes(time_step,mean.value, color=age_class))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim_diagnostics/sim2_run1.png")
#ggsave("sim2_run1.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

preylist.mean%>%
  ggplot(aes(time_step,mean.value))+
  geom_line()+
  geom_ribbon(aes(ymin=mean.value-sd.value, ymax=mean.value+sd.value), fill="grey", alpha=0.5)+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim_diagnostics/sim2_pre_burnin.png")
#ggsave("sim2_pre_burnin.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

preylist.mean%>%
  filter(time_step >= 300 & time_step <= 400)%>%
  ggplot(aes(time_step,mean.value))+
  geom_line()+
  geom_ribbon(aes(ymin=mean.value-sd.value, ymax=mean.value+sd.value), fill="grey", alpha=0.5)+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim_diagnostics/sim2_post_burnin.png")
#ggsave("sim2_post_burnin.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

#grob plot of each age by itself the previous year.  
agelist <-unique(preypiv$age_class)
plotlist <-list()
for (i in 2:length(agelist)){
  
  df <-preypiv%>%
    filter(time_step >= 299 & time_step <=401)%>%
    filter(age_class==agelist[i] | age_class==agelist[i-1])%>%
    pivot_wider(id_cols=1:2, names_from = age_class, values_from = value)
  
  prev_yrcls_lag=lag(df[,3],1)
  df<-as.data.frame(cbind(df,prev_yrcls_lag))
  names(df)<-c("index","time_step","yrcls","yrPlus1","yrlag")
  df <-dplyr::select(df, -yrcls)%>%filter(time_step >= 300 & time_step <=400)
  
  p <-ggplot(df,aes(yrlag,yrPlus1, color=as.factor(index)))+
    geom_line()+
    ylab(paste(agelist[i],"(t)"))+xlab(paste(agelist[i-1], "(t-1)"))+
    guides(color="none")+
    theme_classic()
  
  plotlist[[i]]<-p
}
plotlist<-plotlist[2:20]
#plotlist

row1 <-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],ncol=5)
row2 <-plot_grid(plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],ncol=5)
row3 <-plot_grid(plotlist[[11]],plotlist[[12]],plotlist[[13]],plotlist[[14]],plotlist[[15]],ncol=5)
row4 <-plot_grid(plotlist[[16]],plotlist[[17]],plotlist[[18]],plotlist[[19]],ncol=5)

y.grob <-textGrob("Year Class (t)", gp=gpar(fontface="bold",col="black",fontsize=12),rot=90)
x.grob <-textGrob("Prev Year Class (t-1)", gp=gpar(fontface="bold",col="black",fontsize=12))
top.grob <-textGrob("Year class vs. itself", gp=gpar(fontface="bold",col="black",fontsize=12))

mainplot <-plot_grid(row1,row2,row3,row4, nrow=4)
grid.arrange(arrangeGrob(mainplot, left=y.grob, bottom=x.grob,top=top.grob, padding=unit(0,"line")))
ggsave("sim_diagnostics/sim2_ycvsitself.png")

####NO LAG###
#grob plot of each age by itself the previous year.  
agelist <-unique(preypiv$age_class)
plotlist <-list()
for (i in 2:length(agelist)){
  
  df <-preypiv%>%
    filter(time_step >= 300 & time_step <=400)%>%
    filter(age_class==agelist[i] | age_class==agelist[i-1])%>%
    pivot_wider(id_cols=1:2, names_from = age_class, values_from = value)
  names(df)<-c("index","time_step","yrcls","yrPlus1")
  
  p <-ggplot(df,aes(yrcls,yrPlus1, color=as.factor(index)))+
    geom_line()+
    ylab(paste(agelist[i],"(t)"))+xlab(paste(agelist[i-1], "(t)"))+
    guides(color="none")+
    theme_classic()
  
  plotlist[[i]]<-p
}
plotlist<-plotlist[2:20]
#plotlist

row1 <-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],ncol=5)
row2 <-plot_grid(plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],ncol=5)
row3 <-plot_grid(plotlist[[11]],plotlist[[12]],plotlist[[13]],plotlist[[14]],plotlist[[15]],ncol=5)
row4 <-plot_grid(plotlist[[16]],plotlist[[17]],plotlist[[18]],plotlist[[19]],ncol=5)

y.grob <-textGrob("age i (t)", gp=gpar(fontface="bold",col="black",fontsize=12),rot=90)
x.grob <-textGrob("age i-1 (t)", gp=gpar(fontface="bold",col="black",fontsize=12))
top.grob <-textGrob("age vs previous age", gp=gpar(fontface="bold",col="black",fontsize=12))

mainplot <-plot_grid(row1,row2,row3,row4, nrow=4)
grid.arrange(arrangeGrob(mainplot, left=y.grob, bottom=x.grob,top=top.grob, padding=unit(0,"line")))
ggsave("sim_diagnostics/sim2_ycvsprevyr.png")


preylists <-preylists %>%
  filter(time_step > 299)

write.csv(preylists,"simulated_data/Simulation2_data.csv")
write.csv(recnoises,"simulated_data/sim2_info.csv")

################################ Simulation III ##################################################
#### formerly known as TWO SPP XXL SHORT OSC PREYLISTS

set.seed(4649)
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
    if(mean(meanper$periodt)<8) { #we can't get them all in there unfortunately. 
      preylist[[m]] <- prey
    }else
      preylist[[m]] <- NA
    
    
  }, error=function(e){})
  
}


##############################################################
####################
preylist <-preylist[!is.na(preylist)] #remove all NA elements
length(preylist)
preylist <-preylist[!sapply(preylist,is.null)] #remove all null elements
length(preylist)
preylist <-preylist[1:100] #make sure it's 100 units long

#some data formatting
preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()
preypiv <- preylists%>%pivot_longer(3:22, names_to = 'age_class')
preypiv$age_class = fct_relevel(preypiv$age_class, c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20"))
preypiv$age_class = fct_recode(preypiv$age_class, "age 1"="V1","age 2"="V2","age 3"="V3","age 4"="V4","age 5"="V5","age 6"="V6", "age 7"="V7", "age 8"="V8","age 9"="V9",
                               "age 10"="V10","age 11"="V11","age 12"="V12","age 13"="V13", "age 14"="V14",'age 15'='V15',"age 16"='V16','age 17'= "V17",'age 18'="V18","age 19"="V19","age 20"="V20")
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
recnoises$count0peaks <-unlist(count0peaks)

#what is the period?
mean(recnoises$meanPeriod, na.rm=T)
sd(recnoises$meanPeriod, na.rm=T)
median(recnoises$meanPeriod, na.rm=T)

preylists <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()

##summarize sim runs to a mean. 
preylist.mean <-preypiv%>%group_by(time_step,age_class)%>%summarize(mean.value=mean(value),sd.value=sd(value),.groups='drop')

#just look at the first one, not the mean
preypiv%>%
  filter(index==1)%>%
  filter(time_step >= 300 & time_step <=400)%>%
  #filter(time_step >= 300 & time_step <=320)%>%
  ggplot(aes(time_step,value))+
  #ggplot(aes(time_step,mean.value, color=age_class))+
  geom_line()+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim_diagnostics/sim3_run1.png")
#ggsave("sim3_run1.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

#stacked age classes of run 1 with colors
preypiv%>%
  filter(index==1)%>%
  filter(time_step >= 300 & time_step <=400)%>%
  #filter(time_step >= 300 & time_step <=320)%>%
  #ggplot(aes(time_step,value))+
  ggplot(aes(time_step,value, color=age_class))+
  geom_line()+
  #facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim_diagnostics/sim3_run1colors.png")
#ggsave("sim3_run1.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

preylist.mean%>%
  ggplot(aes(time_step,mean.value))+
  geom_line()+
  geom_ribbon(aes(ymin=mean.value-sd.value, ymax=mean.value+sd.value), fill="grey", alpha=0.5)+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim_diagnostics/sim3_pre_burnin.png")
#ggsave("sim3_pre_burnin.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

preylist.mean%>%
  filter(time_step >= 300 & time_step <= 400)%>%
  ggplot(aes(time_step,mean.value))+
  geom_line()+
  geom_ribbon(aes(ymin=mean.value-sd.value, ymax=mean.value+sd.value), fill="grey", alpha=0.5)+
  facet_wrap(~age_class,scales="free")+
  theme_classic()
ggsave("sim_diagnostics/sim3_post_burnin.png")
#ggsave("sim3_post_burnin.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()


#grob plot of each age by itself the previous year.  
agelist <-unique(preypiv$age_class)
plotlist <-list()
for (i in 2:length(agelist)){
  
  df <-preypiv%>%
    filter(time_step >= 299 & time_step <=401)%>%
    filter(age_class==agelist[i] | age_class==agelist[i-1])%>%
    pivot_wider(id_cols=1:2, names_from = age_class, values_from = value)
  
  prev_yrcls_lag=lag(df[,3],1)
  df<-as.data.frame(cbind(df,prev_yrcls_lag))
  names(df)<-c("index","time_step","yrcls","yrPlus1","yrlag")
  df <-dplyr::select(df, -yrcls)%>%filter(time_step >= 300 & time_step <=400)
  
  p <-ggplot(df,aes(yrlag,yrPlus1, color=as.factor(index)))+
    geom_line()+
    ylab(paste(agelist[i],"(t)"))+xlab(paste(agelist[i-1], "(t-1)"))+
    guides(color="none")+
    theme_classic()
  
  plotlist[[i]]<-p
}
plotlist<-plotlist[2:20]
#plotlist

row1 <-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],ncol=5)
row2 <-plot_grid(plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],ncol=5)
row3 <-plot_grid(plotlist[[11]],plotlist[[12]],plotlist[[13]],plotlist[[14]],plotlist[[15]],ncol=5)
row4 <-plot_grid(plotlist[[16]],plotlist[[17]],plotlist[[18]],plotlist[[19]],ncol=5)

y.grob <-textGrob("Year Class (t)", gp=gpar(fontface="bold",col="black",fontsize=12),rot=90)
x.grob <-textGrob("Prev Year Class (t-1)", gp=gpar(fontface="bold",col="black",fontsize=12))
top.grob <-textGrob("Year class vs. itself", gp=gpar(fontface="bold",col="black",fontsize=12))

mainplot <-plot_grid(row1,row2,row3,row4, nrow=4)
grid.arrange(arrangeGrob(mainplot, left=y.grob, bottom=x.grob,top=top.grob, padding=unit(0,"line")))
ggsave("sim_diagnostics/sim3_ycvsitself.png")

####NO LAG###
#grob plot of each age by itself the previous year.  
agelist <-unique(preypiv$age_class)
plotlist <-list()
for (i in 2:length(agelist)){
  
  df <-preypiv%>%
    filter(time_step >= 300 & time_step <=400)%>%
    filter(age_class==agelist[i] | age_class==agelist[i-1])%>%
    pivot_wider(id_cols=1:2, names_from = age_class, values_from = value)
  names(df)<-c("index","time_step","yrcls","yrPlus1")
  
  p <-ggplot(df,aes(yrcls,yrPlus1, color=as.factor(index)))+
    geom_line()+
    ylab(paste(agelist[i],"(t)"))+xlab(paste(agelist[i-1], "(t)"))+
    guides(color="none")+
    theme_classic()
  
  plotlist[[i]]<-p
}
plotlist<-plotlist[2:20]
#plotlist

row1 <-plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],ncol=5)
row2 <-plot_grid(plotlist[[6]],plotlist[[7]],plotlist[[8]],plotlist[[9]],plotlist[[10]],ncol=5)
row3 <-plot_grid(plotlist[[11]],plotlist[[12]],plotlist[[13]],plotlist[[14]],plotlist[[15]],ncol=5)
row4 <-plot_grid(plotlist[[16]],plotlist[[17]],plotlist[[18]],plotlist[[19]],ncol=5)

y.grob <-textGrob("age i (t)", gp=gpar(fontface="bold",col="black",fontsize=12),rot=90)
x.grob <-textGrob("age i-1 (t)", gp=gpar(fontface="bold",col="black",fontsize=12))
top.grob <-textGrob("age vs previous age", gp=gpar(fontface="bold",col="black",fontsize=12))

mainplot <-plot_grid(row1,row2,row3,row4, nrow=4)
grid.arrange(arrangeGrob(mainplot, left=y.grob, bottom=x.grob,top=top.grob, padding=unit(0,"line")))
ggsave("sim_diagnostics/sim3_ycvsprevyr.png")

## SAVE 
preylists <-preylists %>%
  filter(time_step > 299)

#write.csv(preylists,"Simulation3_dataNEW.csv")
write.csv(preylists,"simulated_data/Simulation3_data.csv")
#write.csv(recnoises,"sim3_infoNEW.csv")
write.csv(recnoises,"simulated_data/sim3_info.csv")

