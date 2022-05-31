### Simulation data plotted over predictions##
### 4/5/22, updated 5/30/22
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(cowplot)
library(rEDM)
library(GPEDM)

### Tanya request figure--- Simulation data with prediction
### 
### After playing with this, I think the best course of action is to put all the age specific
### facet plots as supplementary data. 

#load preylists and take only the first simulation run. 
preylist1 <-read.csv("Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
 filter(index=="1")%>%
 pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()%>%mutate(Simulation="I")
preylist2 <-read.csv("Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
 filter(index=="1")%>%
 pivot_longer(3:23, names_to = "age_class")%>%as.data.frame()%>%mutate(Simulation="II")%>%filter(age_class !="V21")
preylist3 <-read.csv("Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
 filter(index=="1")%>%
 pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()%>%mutate(Simulation="III")

### Plot of all data
preyall <-bind_rows(preylist1,preylist2,preylist3)

#recode the age classes and organize for the graph. 
preyall$age_class <-fct_recode(preyall$age_class, "Age 1"="V1", "Age 2"="V2", "Age 3"="V3", "Age 4"="V4",
                               "Age 5"="V5","Age 6"="V6","Age 7"="V7","Age 8"="V8","Age 9"="V9", "Age 10"="V10", "Age 11"="V11",
                               "Age 12" = "V12", "Age 13"="V13", "Age 14"= "V14","Age 15"="V15","Age 16"="V16","Age 17"="V17",
                               "Age 18"="V18", "Age 19"="V19","Age 20"="V20")
preyall$age_class <-fct_relevel(preyall$age_class, "Age 1","Age 2","Age 3","Age 4","Age 5","Age 6", "Age 7",
                                "Age 8","Age 9","Age 10","Age 11",'Age 12',"Age 13","Age 14","Age 15","Age 16","Age 17","Age 18",
                                "Age 19","Age 20")

#simulation 1 plot
preyall %>%
        filter(time_step <=400)%>%
        filter(Simulation == "I")%>%
        ggplot()+
        xlab("Time")+ylab("Abundance")+
        geom_point(aes(x=time_step, y=value), size=0.25)+ #value is raw observations
        facet_wrap(age_class~., scales="free_y")+
        #facet_grid(age_class~Simulation, scales="free")+
        theme_classic()



### SIM 1 ###################################################################################################################
#Make 30 year datasets + ten year test dataset to fit GPs#

#preylist1 <-preylists%>%pivot_longer(3:22, names_to = 'age_class')%>%
      #  filter(index==1)

#find the best E and tau - 
pgrid <-as.data.frame(preylist1) %>%filter(age_class !="V21" & time_step >=300 & time_step < 330)
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

#E=8, tau=1

#########Make 30 year datasets + ten year test dataset to fit GPs###
### Remember to log transform the value ###
p1 <-filter(preylist1, time_step >=300 & time_step <= 350)%>% as.data.frame()%>%mutate(value=log(value))

p1Lags = makelags(data=p1, y="value", pop="age_class", E=bestE, tau=bestTau)
#p1Lags = makelags(data=p1, y="value", pop="age_class", E=round(sqrt(30)), tau=1)
p1 = cbind(p1,p1Lags)
p1.train = filter(p1, time_step <= (max(p1$time_step)-20))
p1.test = filter(p1, time_step > (max(p1$time_step)-20))

p1gp <-fitGP(data = p1.train, y = "value", x=colnames(p1Lags),newdata=p1.test,pop="age_class",scaling = "local",predictmethod = "loo")


#### extract results ####
#### remember since the data is logged you will have to add it up differently###
sumtotal <-p1gp$outsampresults %>% 
 mutate(pred_mean=exp(predmean), OBS=exp(obs),Predfsd=exp(predfsd), Predsd=exp(predsd))%>%
 group_by(timestep)%>%
 summarise(pop="aggregate", predmean=sum(pred_mean), obs=sum(OBS), predfsd=sqrt(sum((Predfsd^2))),predsd=sqrt(sum((Predsd^2))))%>%
 mutate(time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%select(-timestep)%>%
 mutate(ymax=predmean+predfsd, ymin=predmean-predfsd)
 
#REMEMBER TO LOG AND UNLOG
 #sumtotal <-mutate(sumtotal,predmean=log(predmean), obs=log(obs), predfsd=log(predfsd), predsd=log(predsd), ymax=log(ymax),ymin=log(ymin))

# I am having a hard time trying to understand whether, because we Log transformed the data, we need
# to then exponentiate the predicted mean in order to add the standard deviation to show the
# uncertainty around the prediction. But I can't do that and then re-log it because of issues log transforming
#tried to just exponentiate (unlog) everything. So that is what we are doing. Presenting the unlogged data. 
p1res <-mutate(p1gp$outsampresults, time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%
 mutate(ymin=exp(predmean)-exp(predfsd), ymax=exp(predmean)+exp(predfsd))%>%
 mutate(predmean=exp(predmean))
p1.5 <-mutate(p1, value=exp(value))%>%dplyr::select(time_step, age_class, value) 

#REMEMBER TO LOG AND UNLOG
#p1res <-mutate(p1res,ymin=log(ymin), ymax=log(ymax), predmean=log(predmean))
#p1.5 <-mutate(p1.5, value=log(value))

p1res <-full_join(p1.5,p1res)%>%select(-timestep)

p1res <-bind_rows(sumtotal, p1res)%>%arrange(time_step, age_class)

#in order to get data for the "value" column of p1res for the aggregate prediction, 
#we will have to add the raw data - note that it is unlogged at present. 
rawta <-filter(preylist1, time_step >=300 & time_step <= 350)%>% 
        group_by(time_step)%>%summarize(TA=sum(value))%>%
        mutate(age_class="aggregate")

##REMEMBER TO LOG AND UNLOG
#rawta <-mutate(rawta, TA=log(TA))

p1.5res <-full_join(p1res,rawta)
p1.5res <-mutate(p1.5res, value=ifelse(age_class=="aggregate" & is.na(value), TA, value))

#fit a new GP to the total abundance dataset
ta1 <-rawta%>% as.data.frame()%>%mutate(value=log(TA))
ta1Lags = makelags(data=ta1, y="value", pop="age_class", E=round(sqrt(30)), tau=1)
ta1 = cbind(ta1,ta1Lags)
ta1.train = filter(ta1, time_step <= (max(ta1$time_step)-20))
ta1.test = filter(ta1, time_step > (max(ta1$time_step)-20))
#fit the GP
ta1gp <-fitGP(data = ta1.train, y = "value", x=colnames(ta1Lags),newdata=ta1.test,pop="age_class",scaling = "local",predictmethod = "loo")
#extract the results and convert out of log scale
ta1res <-mutate(ta1gp$outsampresults, time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%
        mutate(ymin=exp(predmean)-exp(predfsd), ymax=exp(predmean)+exp(predfsd))%>%
        mutate(predmean=exp(predmean))
ta1res$age_class <-fct_recode(ta1res$age_class,"Aggregate"="aggregate")

##REMEMBER TO LOG AND UNLOG
#ta1res <-mutate(ta1res, ymin=log(ymin), ymax=log(ymax), predmean=log(predmean))

#recode the age classes and organize for the graph. 
p1.5res$age_class <-fct_recode(p1.5res$age_class, "Age 1"="V1", "Age 2"="V2", "Age 3"="V3", "Age 4"="V4",
        "Age 5"="V5","Age 6"="V6","Age 7"="V7","Age 8"="V8","Age 9"="V9", "Age 10"="V10", "Age 11"="V11",
        "Age 12" = "V12", "Age 13"="V13", "Age 14"= "V14","Age 15"="V15","Age 16"="V16","Age 17"="V17",
        "Age 18"="V18", "Age 19"="V19","Age 20"="V20","Aggregate"="aggregate")
p1.5res$age_class <-fct_relevel(p1.5res$age_class, "Age 1","Age 2","Age 3","Age 4","Age 5","Age 6", "Age 7",
                "Age 8","Age 9","Age 10","Age 11",'Age 12',"Age 13","Age 14","Age 15","Age 16","Age 17","Age 18",
                "Age 19","Age 20","Aggregate")

 p1.5res %>%
         filter(age_class !="Aggregate")%>%
 ggplot()+
         xlab("Time")+ylab("Abundance")+
 geom_point(aes(x=time_step, y=value), size=0.75)+ #value is raw observations
 geom_ribbon(aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="blue") + #out of sample
 geom_line(aes(x=time_step, y=predmean), color="blue")+
 facet_wrap(age_class~., scales="free_y")+
 theme_classic()
 #save
 #ggsave("Sim1obspreds_age.png",height=6, width=8,dpi=300,  path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
 ggsave("Sim1obspreds_age.png",height=6, width=8,dpi=300)
 dev.off()
 
 p1.5res %>%
         filter(age_class=="Aggregate")%>%
         ggplot()+
         geom_point(aes(x=time_step, y=value))+ #value is raw observations
         xlab("Time")+ylab("Abundance")+
         geom_ribbon(aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="blue") + #out of sample
         geom_line(aes(x=time_step, y=predmean), color="blue")+
         geom_ribbon(data=ta1res, aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="red") + #out of sample
         geom_line(data=ta1res,aes(x=time_step, y=predmean), color="red")+
         facet_wrap(age_class~., scales="free")+
         theme_classic()
 #save
 #ggsave("Sim1obspreds_agg.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
 ggsave("Sim1obspreds_agg.png")
 dev.off()
 #using the total abundance index produces an identical fit which is super annoying. 
 #so maybe don't include it (the red)

 
#### SIM 2 ###################################################################################################################
 #find the best E and tau - 
 pgrid <-as.data.frame(preylist2) %>%filter(age_class !="V21" & time_step >=300 & time_step < 330)
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
 r2matrix2
 rmsematrix2
 #grab the position of the best E and tau from the matrix. 
 bestET <-which(rmsematrix2==min(rmsematrix2,na.rm=T),arr.ind=T)
 bestE <-as.numeric(noquote(rownames(bestET)))
 bestTau <-as.numeric(bestET[2])
 #previously best E=9 best tau=2
 
 
 #Make 30 year datasets + ten year test dataset to fit GPs#
 ### Remember to log transform the value ###
 p2 <-filter(preylist2, time_step >=300 & time_step <= 340, age_class != "V21")%>% as.data.frame()%>%mutate(value=log(value))
 
 p2Lags = makelags(data=p2, y="value", pop="age_class", E=bestE, tau=bestTau)
 #p2Lags = makelags(data=p2, y="value", pop="age_class", E=round(sqrt(30)), tau=1)
 p2 = cbind(p2,p2Lags)
 p2.train = filter(p2, time_step <= (max(p2$time_step)-10))
 p2.test = filter(p2, time_step > (max(p2$time_step)-10))
 
 p2gp <-fitGP(data = p2.train, y = "value", x=colnames(p2Lags),newdata=p2.test,pop="age_class",scaling = "local",predictmethod = "loo")
 
 
 #### extract results ####
 #### remember since the data is logged you will have to add it up differently###
 sumtotal <-p2gp$outsampresults %>% 
         mutate(pred_mean=exp(predmean), OBS=exp(obs),Predfsd=exp(predfsd), Predsd=exp(predsd))%>%
         group_by(timestep)%>%
         summarise(pop="aggregate", predmean=sum(pred_mean), obs=sum(OBS), predfsd=sqrt(sum((Predfsd^2))),predsd=sqrt(sum((Predsd^2))))%>%
         mutate(time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%select(-timestep)%>%
         mutate(ymax=predmean+predfsd, ymin=predmean-predfsd)
 #sumtotal <-mutate(sumtotal,predmean=log(pred_mean), obs=log(OBS), predfsd=log(Predfsd), predsd=log(Predsd))
 
 # I am having a hard time trying to understand whether, because we Log transformed the data, we need
 # to then exponentiate the predicted mean in order to add the standard deviation to show the
 # uncertainty around the prediction. But I can't do that and then re-log it because of issues log transforming
 #tried to just exponentiate (unlog) everything. So that is what we are doing. Presenting the unlogged data. 
 p2res <-mutate(p2gp$outsampresults, time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%
         mutate(ymin=exp(predmean)-exp(predfsd), ymax=exp(predmean)+exp(predfsd))%>%
         mutate(predmean=exp(predmean))
 #mutate(ymin=log(ymin), ymax=log(ymax))
 p2.5 <-mutate(p2, value=exp(value))%>%dplyr::select(time_step, age_class, value) 
 p2res <-full_join(p2.5,p2res)%>%select(-timestep)
 
 p2res <-bind_rows(sumtotal, p2res)%>%arrange(time_step, age_class)
 
 #in order to get data for the "value" column of p2res for the aggregate prediction, 
 #we will have to add the raw data - note that it is unlogged at present. 
 rawta <-filter(preylist2, time_step >=300 & time_step <= 340, age_class!="V21")%>% 
         group_by(time_step)%>%summarize(TA=sum(value))%>%
         mutate(age_class="aggregate")
         #mutate(time_step=time_step +1) ### IS THIS FIXING MISALIGN OR MAKING IT WORSE?
 p2.5res <-full_join(p2res,rawta)
 p2.5res <-mutate(p2.5res, value=ifelse(age_class=="aggregate" & is.na(value), TA, value))
 
 #fit a new GP to the total abundance dataset
 ta2 <-rawta%>% as.data.frame()%>%mutate(value=log(TA))
 ta2Lags = makelags(data=ta2, y="value", pop="age_class", E=round(sqrt(30)), tau=1)
 ta2 = cbind(ta2,ta2Lags)
 ta2.train = filter(ta2, time_step <= (max(ta2$time_step)-10))
 ta2.test = filter(ta2, time_step > (max(ta2$time_step)-10))
 #fit the GP
 ta2gp <-fitGP(data = ta2.train, y = "value", x=colnames(ta2Lags),newdata=ta2.test,pop="age_class",scaling = "local",predictmethod = "loo")
 #extract the results and convert out of log scale
 ta2res <-mutate(ta2gp$outsampresults, time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%
         mutate(ymin=exp(predmean)-exp(predfsd), ymax=exp(predmean)+exp(predfsd))%>%
         mutate(predmean=exp(predmean))
 ta2res$age_class <-fct_recode(ta2res$age_class,"Aggregate"="aggregate")
 
 #recode the age classes and organize for the graph. 
 p2.5res$age_class <-fct_recode(p2.5res$age_class, "Age 1"="V1", "Age 2"="V2", "Age 3"="V3", "Age 4"="V4",
                                "Age 5"="V5","Age 6"="V6","Age 7"="V7","Age 8"="V8","Age 9"="V9", "Age 10"="V10", "Age 11"="V11",
                                "Age 12" = "V12", "Age 13"="V13", "Age 14"= "V14","Age 15"="V15","Age 16"="V16","Age 17"="V17",
                                "Age 18"="V18", "Age 19"="V19","Age 20"="V20","Aggregate"="aggregate")
 p2.5res$age_class <-fct_relevel(p2.5res$age_class, "Age 1","Age 2","Age 3","Age 4","Age 5","Age 6", "Age 7",
                                 "Age 8","Age 9","Age 10","Age 11",'Age 12',"Age 13","Age 14","Age 15","Age 16","Age 17","Age 18",
                                 "Age 19","Age 20","Aggregate")
 
 p2.5res %>%
         filter(age_class !="Aggregate")%>%
         ggplot()+
         xlab("Time")+ylab("Abundance")+
         geom_point(aes(x=time_step, y=value), size=0.75)+ #value is raw observations
         geom_ribbon(aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="blue") + #out of sample
         geom_line(aes(x=time_step, y=predmean), color="blue")+
         facet_wrap(age_class~., scales="free_y")+
         theme_classic()
 #save
 #ggsave("Sim2obspreds_age.png",height=6, width=8,dpi=300, path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
 ggsave("Sim2obspreds_age.png",height=6, width=8,dpi=300)
 dev.off()
 
 p2.5res %>%
         filter(age_class=="Aggregate")%>%
         ggplot()+
         geom_point(aes(x=time_step, y=value))+ #value is raw observations
         xlab("Time")+ylab("Abundance")+
         geom_ribbon(aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="blue") + #out of sample
         geom_line(aes(x=time_step, y=predmean), color="blue")+
         geom_ribbon(data=ta2res, aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="red") + #out of sample
         geom_line(data=ta2res,aes(x=time_step, y=predmean), color="red")+
         facet_wrap(age_class~., scales="free")+
         theme_classic()
 #save
 ggsave("Sim2obspreds_agg.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
 ggsave("Sim2obspreds_agg.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
 dev.off()
 #using the total abundance index produces an identical fit which is super annoying. 
 #so maybe don't include it (the red)
 
 ### SIM 3 ###################################################################################################################
 
 #find the best E and tau - 
 pgrid <-as.data.frame(preylist3) %>%filter(age_class !="V21" & time_step >=300 & time_step < 330)
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
     r2matrix3[var_pairs[i,1], var_pairs[i,2]] = fit1_r2
     ETdf[i,3] <-fit1_r2
     ETdf[i,4] <-fit1_rmse
     rmsematrix3[var_pairs[i,1], var_pairs[i,2]] = fit1_rmse
   },silent=F)
 }
 r2matrix3
 rmsematrix3
 #grab the position of the best E and tau from the matrix. 
 bestET <-which(rmsematrix3==min(rmsematrix3,na.rm=T),arr.ind=T)
 bestE <-as.numeric(noquote(rownames(bestET)))
 bestTau <-as.numeric(bestET[2])
 #previously best E=8 best tau=1
 
 #Make 30 year datasets + ten year test dataset to fit GPs#
 ### Remember to log transform the value ###
 p3 <-filter(preylist3, time_step >=300 & time_step <= 340)%>% as.data.frame()%>%mutate(value=log(value))
 
 p3Lags = makelags(data=p3, y="value", pop="age_class", E=bestE, tau=bestTau)
 #p3Lags = makelags(data=p3, y="value", pop="age_class", E=round(sqrt(30)), tau=1)
 p3 = cbind(p3,p3Lags)
 p3.train = filter(p3, time_step <= (max(p3$time_step)-10))
 p3.test = filter(p3, time_step > (max(p3$time_step)-10))
 
 p3gp <-fitGP(data = p3.train, y = "value", x=colnames(p3Lags),newdata=p3.test,pop="age_class",scaling = "local",predictmethod = "loo")
 
 
 #### extract results ####
 #### remember since the data is logged you will have to add it up differently###
 sumtotal <-p3gp$outsampresults %>% 
         mutate(pred_mean=exp(predmean), OBS=exp(obs),Predfsd=exp(predfsd), Predsd=exp(predsd))%>%
         group_by(timestep)%>%
         summarise(pop="aggregate", predmean=sum(pred_mean), obs=sum(OBS), predfsd=sqrt(sum((Predfsd^2))),predsd=sqrt(sum((Predsd^2))))%>%
         mutate(time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%select(-timestep)%>%
         mutate(ymax=predmean+predfsd, ymin=predmean-predfsd)
 #sumtotal <-mutate(sumtotal,predmean=log(pred_mean), obs=log(OBS), predfsd=log(Predfsd), predsd=log(Predsd))
 
 # I am having a hard time trying to understand whether, because we Log transformed the data, we need
 # to then exponentiate the predicted mean in order to add the standard deviation to show the
 # uncertainty around the prediction. But I can't do that and then re-log it because of issues log transforming
 #tried to just exponentiate (unlog) everything. So that is what we are doing. Presenting the unlogged data. 
 p3res <-mutate(p3gp$outsampresults, time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%
         mutate(ymin=exp(predmean)-exp(predfsd), ymax=exp(predmean)+exp(predfsd))%>%
         mutate(predmean=exp(predmean))
 #mutate(ymin=log(ymin), ymax=log(ymax))
 p3.5 <-mutate(p3, value=exp(value))%>%dplyr::select(time_step, age_class, value) 
 p3res <-full_join(p3.5,p3res)%>%select(-timestep)
 
 p3res <-bind_rows(sumtotal, p3res)%>%arrange(time_step, age_class)
 
 #in order to get data for the "value" column of p3res for the aggregate prediction, 
 #we will have to add the raw data - note that it is unlogged at present. 
 rawta <-filter(preylist3, time_step >=300 & time_step <= 340)%>% 
         group_by(time_step)%>%summarize(TA=sum(value))%>%
         mutate(age_class="aggregate")
 p3.5res <-full_join(p3res,rawta)
 p3.5res <-mutate(p3.5res, value=ifelse(age_class=="aggregate" & is.na(value), TA, value))
 
 #fit a new GP to the total abundance dataset
 ta3 <-rawta%>% as.data.frame()%>%mutate(value=log(TA))
 ta3Lags = makelags(data=ta3, y="value", pop="age_class", E=round(sqrt(30)), tau=1)
 ta3 = cbind(ta3,ta3Lags)
 ta3.train = filter(ta3, time_step <= (max(ta3$time_step)-10))
 ta3.test = filter(ta3, time_step > (max(ta3$time_step)-10))
 #fit the GP
 ta3gp <-fitGP(data = ta3.train, y = "value", x=colnames(ta3Lags),newdata=ta3.test,pop="age_class",scaling = "local",predictmethod = "loo")
 #extract the results and convert out of log scale
 ta3res <-mutate(ta3gp$outsampresults, time_step=timestep+330)%>%dplyr::rename(age_class=pop)%>%
         mutate(ymin=exp(predmean)-exp(predfsd), ymax=exp(predmean)+exp(predfsd))%>%
         mutate(predmean=exp(predmean))
 ta3res$age_class <-fct_recode(ta3res$age_class,"Aggregate"="aggregate")
 
 #recode the age classes and organize for the graph. 
 p3.5res$age_class <-fct_recode(p3.5res$age_class, "Age 1"="V1", "Age 2"="V2", "Age 3"="V3", "Age 4"="V4",
                                "Age 5"="V5","Age 6"="V6","Age 7"="V7","Age 8"="V8","Age 9"="V9", "Age 10"="V10", "Age 11"="V11",
                                "Age 12" = "V12", "Age 13"="V13", "Age 14"= "V14","Age 15"="V15","Age 16"="V16","Age 17"="V17",
                                "Age 18"="V18", "Age 19"="V19","Age 20"="V20","Aggregate"="aggregate")
 p3.5res$age_class <-fct_relevel(p3.5res$age_class, "Age 1","Age 2","Age 3","Age 4","Age 5","Age 6", "Age 7",
                                 "Age 8","Age 9","Age 10","Age 11",'Age 12',"Age 13","Age 14","Age 15","Age 16","Age 17","Age 18",
                                 "Age 19","Age 20","Aggregate")
 
 p3.5res %>%
         filter(age_class !="Aggregate")%>%
         ggplot()+
         xlab("Time")+ylab("Abundance")+
         geom_point(aes(x=time_step, y=value), size=0.75)+ #value is raw observations
         geom_ribbon(aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="blue") + #out of sample
         geom_line(aes(x=time_step, y=predmean), color="blue")+
         facet_wrap(age_class~., scales="free_y")+
         theme_classic()
 #save
 ggsave("Sim3obspreds_age.png", height= 6, width=8, dpi=300)
 #ggsave("Sim3obspreds_age.png", height= 6, width=8, dpi=300, path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
 dev.off()
 
 p3.5res %>%
         filter(age_class=="Aggregate")%>%
         ggplot()+
         geom_point(aes(x=time_step, y=value))+ #value is raw observations
         xlab("Time")+ylab("Abundance")+
         geom_ribbon(aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="blue") + #out of sample
         geom_line(aes(x=time_step, y=predmean), color="blue")+
         geom_ribbon(data=ta3res, aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="red") + #out of sample
         geom_line(data=ta3res,aes(x=time_step, y=predmean), color="red")+
         facet_wrap(age_class~., scales="free")+
         theme_classic()
 #save
 ggsave("Sim3obspreds_agg.png")
 #ggsave("Sim3obspreds_agg.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
 dev.off()
 #using the total abundance index produces an identical fit which is super annoying. 
 #so maybe don't include it (the red)
 
 ############### Combine simulations into one panel plot #########################
 
 p1.5resag <-mutate(p1.5res, Simulation ="I")%>%filter(age_class=="Aggregate")
 p2.5resag <-mutate(p2.5res, Simulation ="II")%>%filter(age_class=="Aggregate")
 p3.5resag <-mutate(p3.5res, Simulation ="III")%>%filter(age_class=="Aggregate")
 aggres <-bind_rows(p1.5resag, p2.5resag, p3.5resag)
 ta1res <-mutate(ta1res, Simulation="I")
 ta2res <-mutate(ta2res, Simulation="II")
 ta3res <-mutate(ta3res, Simulation="III")
 aggta <-bind_rows(ta1res,ta2res,ta3res)
 
 ##REMEMBER TO LOG AND UNLOG
 #aggta <-mutate(aggta, ymin=log(ymin), ymax=log(ymax), predmean=log(predmean))
 #aggres <-mutate(aggres, predmean=log(predmean), ymax=log(ymax), ymin=log(ymin), value=log(value))
 
aggres %>%
         ggplot()+
         geom_point(aes(x=time_step, y=value))+ #value is raw observations
         xlab("Time")+ylab("Abundance")+
         geom_ribbon(aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="blue") + #out of sample
         geom_line(aes(x=time_step, y=predmean), color="blue")+
         geom_ribbon(data=aggta, aes(x=time_step,y=predmean,ymin=ymin,ymax=ymax),alpha=0.4,fill="red") + #out of sample
         geom_line(data=aggta,aes(x=time_step, y=predmean), color="red")+
         facet_wrap(Simulation~., scales="free")+
         theme_classic()
ggsave("aggregatepreds20ages.png")
#ggsave("aggregatepreds20ages.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

#################### Compare the different ages and years ##############