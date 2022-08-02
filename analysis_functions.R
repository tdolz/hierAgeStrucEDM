################ FUNCTIONS YOU WILL NEED FOR THE ANALYSIS ##############################
################ 
################ 3/18/22

################# Fit GP to one age class at a time for 100 years of data, export fitstats ######
##### This is used in the "hundred_data_comp.R" script
indv_age <-function(df){
 ages <-unique(df$age_class)
 maxE <-10 # in all cases it's going to be 10 because the ts is always 100 points long for this one. 
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


################# ####### Function to create the aggregate individual age prediction #######################
## basically, make a prediction for each single age and add up all the single age predictions to create an aggregate N-total. 
## #basically the indv_age function, but you want to be able to specify time series length. 
## This is mostly deprecated, but was used in the  "hundred_data_comp.R"

aggPred <-function(df,ts){ 
  max_time <-min(df$time_step)+ts+10
  df <-filter(df, time_step <=max_time)
  ages <-unique(df$age_class)
  maxE <-round(sqrt(length(unique(df$time_step))-10)) 
  mod1res <-list()
  for (m in 1:length(ages)){
    new_df <-filter(df, age_class==ages[m])
    newdfLags = makelags(data=new_df, y="value", E=maxE, tau=1)
    new_df = cbind(new_df,newdfLags)
    new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
    new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
    mod1 <-fitGP(data = new_df.train, y = "value", x=colnames(newdfLags),newdata=new_df.test,
                 pop="age_class",scaling = "global",predictmethod = "loo")
    res <-as.data.frame(mod1$outsampresults)%>%mutate(age=ages[m])
    mod1res[[m]] <-res
  }
  modres <-bind_rows(mod1res)%>%as.data.frame() #bind all the predicted and observed
  modsum <-modres %>% group_by(timestep) %>%mutate(predmean=exp(predmean), obs=exp(obs))%>% #exponentiate
    summarize(aggpred =sum(predmean), aggobs=sum(obs))%>% #summarize by time step to get an aggregate predicted and observed
    mutate(aggpred=log(aggpred), aggobs=log(aggobs)) #log again
  aggR2dos <-getR2(modsum$aggobs,modsum$aggpred)
  aggR2dos
  # aggR2 <-mutate(modsum, ssr=(aggobs - aggpred)^2, sst=(aggobs-mean(modsum$aggobs))^2)
  # r2 <-1-sum(aggR2$ssr)/sum(aggR2$sst)
  # r2
}

########################################Create NTOTAL predictions for specified time series length  ########################### 
############# the input has to be an n_total index #####################
# used in the hundred_data_comp.R
ntotal_yrs <-function(df){
 ntotalstats <-list()
 yrs <-c(15,20,30,35)
 startime <-min(df$time_step)
 for (m in 1:length(yrs)){
  new_df <-filter(df, time_step <= startime+yrs[m])
  maxE=round(sqrt(yrs[m]-10))
  newdfLags = makelags(data=new_df, y="value", E=maxE, tau=1)
  new_df = cbind(new_df,newdfLags)
  new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
  new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
  mod1 <-fitGP(data = new_df.train, y = "value", x=colnames(newdfLags),newdata=new_df.test,
               pop="age_class",scaling = "local",predictmethod = "loo")
  mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,as.character(yrs[m]))
  names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","tslength")
  ntotalstats[[m]] <- mod1_out
 }
 ntotalstats <-bind_rows(ntotalstats)%>%as.data.frame()%>%mutate(across(OOS_R2:df, as.numeric))
}


########################################Create NTOTAL predictions for specified time series length  ########################### 
############# the input is a long form df with age structure where the values are already log transformed #####################
############# this has been edited since the last draft was sent around (new as of 3/18/22)
# used in the hundred_data_comp.R
ntotal_app<-function(df){
 newdf <-df %>% mutate(newvalue=exp(value))%>%group_by(time_step)%>%summarize(NewValue=sum(newvalue))%>%mutate(value=log(NewValue))
 tslength <-(max(newdf$time_step)-min(newdf$time_step))-10
 maxE=round(sqrt(tslength))
 ages=n_distinct(df$age_class)
 newdfLags = makelags(data=newdf, y="value", E=maxE, tau=1)
 new_df = cbind(newdf,newdfLags)
 new_df.train = filter(new_df, time_step <= (max(new_df$time_step)-10))
 new_df.test = filter(new_df, time_step > (max(new_df$time_step)-10))
 mod1 <-fitGP(data = new_df.train, y = "value", x=colnames(newdfLags),newdata=new_df.test,
              scaling = "local",predictmethod = "loo")
 phis <-mod1$pars[substr(names(mod1$pars),1,3)=="phi"]
 mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,tslength,ages)
 names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","tslength","age class")
 mod1_out <-list(mod1_out,phis)
 mod1_out
}

########################################### count peaks function ###############
########### This counts the number of peaks in a cycle of the simulated data #
########### used so we can limit period of the simulated time series in script "make_simulation_data.R"
count_peaks <-function(x){
 mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak. 
 x <-findpeaks(x, minpeakheight = mph)
 x <-dim(x)[1]
 x[is.null(x)] <-0 #if there are no peaks return 0. 
 x
}

##################################### positive integers ###############
####### use this to count if there are zeros, decimals, anything less than one
### to detect population crashes

lessT1 <-function(x){ #x is a vector of any length
  v <-x < 1
  length(v[v==T])
}

##################################### find period function ###############
###### this measures the period of the oscillating simulated time series 
###### used in "make_simulation_data.R"
period <-function(x){
 #mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak.
 mph <-mean(x)+(sd(x)) #the minimum peak height has to be 1 standard deviation from the mean
 x <-findpeaks(x, minpeakheight = mph)
 x <-x[,1:2]
 x <-as.data.frame(x)
 names(x) <-c("ypeak","xpeak")
 x <-mutate(x, dist = xpeak-lag(xpeak,1))
 mean.period = mean(x$dist,na.rm=T)
 mean.period
}