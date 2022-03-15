##Best Embed functions###
##
##July 24, 2021
##This is a script for final versions of these functions. 
##The sandbox for this script is BestEmbedGP.Rmd and BestE.Rmd
##
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(rEDM)
library(GPEDM)
library(Metrics)
library(MLmetrics)
library(mltools)
library(pracma)

######## count peaks function ###############
count_peaks <-function(x){
  mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak. 
  x <-findpeaks(x, minpeakheight = mph)
  x <-dim(x)[1]
  x[is.null(x)] <-0 #if there are no peaks return 0. 
  x
}




######r-squared function##########
rsq <-function(obs,preds,time_step){
 df <-cbind(obs,preds)%>%cbind(time_step)%>%as.data.frame()
 df<-na.omit(df)
 RSS <-c()
 TSS<-c()
 for (i in 1:length(df$obs)){
  RSS[i] <-(df$obs[i]-df$preds[i])^2
  TSS[i]<-((df$obs[i]-mean(df$obs))^2)}
 1-(mean(RSS)/mean(TSS))}
# 1-(sum(RSS)/sum(TSS))} # if you want to do sums instead. 

######RMSE function#########
RMSE = function(m, o){sqrt(mean((m - o)^2))}
########################################


######## GET BEST E USING SIMPLEX ######## 
GetE <-function(Elength, #the number of Es from 1 to x you want to try. Writing 10 will try E 1:10.
                df, #the data frame
                liblength, #the length of the library from the data frame you want to try.
                targ, #target column in df, must be in quotes
                cols) #predictor columns in df, must be in quotes
{
 
 #get params for simplex
 all = c((NROW(df)-100), (NROW(df)-0)) #the length of the library
 Ees_out <-matrix(ncol=6,nrow=Elength)
 Ees=seq(1,Elength,1)
 
 #run the simplex
 for (i in 1:10){
  simplex_out <- Simplex(dataFrame = df, lib = all, pred = all, columns = cols, E=Ees[i],
                         target = targ)
  Ees_out[i,1] <-Ees[i]
  Ees_out[i,2] <-rsq(simplex_out$Observations, simplex_out$Predictions, simplex_out$time_step) #rsq
  Ees_out[i,3] <-mse(actuals =simplex_out$Observations, preds= simplex_out$Predictions, na.rm=TRUE)#MSE
  ers <-(ComputeError(simplex_out$Observations, simplex_out$Predictions))
  Ees_out[i,4] <-ers$MAE #MAE
  Ees_out[i,5] <-ers$rho #rho
  Ees_out[i,6] <-ers$RMSE #RMSE
 }
 Ees_out <-as.data.frame(Ees_out)
 names(Ees_out)<-c("E","RSq","MSE","MAE","rho","RMSE")
 
 #diagnostic plots
 plot(Ees_out$E, Ees_out$MSE, type="l", xlab="Embedding Dimension (E)",ylab="MSE")
 plot(Ees_out$E, Ees_out$MAE, type="l", xlab="Embedding Dimension (E)",ylab="MAE")
 plot(Ees_out$E, Ees_out$rho, type="l", xlab="Embedding Dimension (E)",ylab="rho")
 plot(Ees_out$E, Ees_out$RSq, type="l", xlab="Embedding Dimension (E)",ylab="Rsq")
 
 #bestE
 bestE_rho <-Ees_out[which.max(Ees_out$rho),'E']
 bestE_MAE <-Ees_out[which.min(Ees_out$MAE),'E']
 bestE_MSE <-Ees_out[which.min(Ees_out$MSE),'E']
 bestE_RSq <-Ees_out[which.max(Ees_out$RSq),'E'] #Tanya used a rounded version of this one. 
 bestE <-c(bestE_rho,bestE_MAE,bestE_MSE, bestE_RSq)
 names(bestE) <- c("rho","MAE","MSE","RSq")
 Ees_out <<-Ees_out
 print(bestE)
}
###############################################


#########Best Embed with GP function##################################################################################################################

### Find best embedding dimension for GP 
GetE_GP <-function(maxE, #the number of Es from 1 to x you want to try. Writing 10 will try E 1:10.
                   df, #the data frame
                   response_yd, ##response (target) columns in df, must be in quotes
                   pop, ## population (if applicable)
                   scaling,# "local" or "global"
                   plot) #boolean
{
  
  #outputs
  Ees_out <-matrix(ncol=7,nrow=maxE)
  Ees=seq(1,maxE,1)
  
  #run the GP
  for (i in 1:maxE){
    
    gp <-fitGP(data = df, yd = response_yd, pop=pop,scaling = scaling, E=Ees[i], tau=1, predictmethod = "loo")
    
    Ees_out[i,1] <-Ees[i]
    Ees_out[i,2] <- gp$outsampfitstats[1] #out of sample rsq
    Ees_out[i,3] <- gp$insampfitstats[1] #in sample rsq
    Ees_out[i,4] <-mse(actuals =gp$outsampresults$obs, preds= gp$outsampresults$predmean, na.rm=TRUE)# OOS MSE
    Ees_out[i,5] <-mse(actuals =gp$insampresults$obs, preds= gp$insampresults$predmean, na.rm=TRUE)# in samp MSE
    Ees_out[i,6] <-gp$insampfitstats[4] #in sample NLL
    Ees_out[i,7]<-gp$insampfitstats[5] # degrees of freedom
    
  }
  Ees_out <-as.data.frame(Ees_out)
  names(Ees_out)<-c("E","OOS_R2","InSamp_R2","OOS_MSE","InSampMSE","InSamp_NLL","df")
  
  if(plot==TRUE){
    #diagnostic plots
    plot(Ees_out$E, Ees_out$OOS_R2, type="l", xlab="Embedding Dimension (E)",ylab="Out of Sample Rsq")
    plot(Ees_out$E, Ees_out$InSamp_R2, type="l", xlab="Embedding Dimension (E)",ylab="In Sample Rsq")
    plot(Ees_out$E, Ees_out$rho, type="l", xlab="Embedding Dimension (E)",ylab="In Sample Negative Log Likelihood")
  }
  
  #bestE
  bestE = c(Ees_out[which.max(Ees_out$OOS_R2),'E'], Ees_out[which.min(Ees_out$InSamp_R2),'E'])
  Ees_out2 <<-Ees_out
}
#############################################################################

#############determine the best E using Tanya's ratio and simplex############
showdiffS <- function(df){
df <-mutate(df, MSEratio = MSE/lead(MSE,1))
p1 <- ggplot(df, aes(E,MSEratio)) + geom_line()+
  geom_hline(yintercept=1.1, linetype="dashed")+
  scale_x_continuous(breaks=seq(0,max(df$E),1), labels=seq(0,max(df$E),1))+
    theme_classic()
df <-mutate(df,sigdif=ifelse(MSEratio < 1.1, "keep lower E", "use higher E"))
bestE <- min(df$E[which(df$sigdif=="keep lower E")])
Ediffs <<-df

multi_return <- function() {
my_list <- list("bestE" = bestE, "df" = Ediffs, "plot" = p1)
return(my_list) 
}
multi_return()
}
#############################################################################

#############determine the best E using Tanya's ratio and GP#################
showdiffGP <- function(df){
  df <-mutate(df, MSEratio = OOS_MSE/lead(OOS_MSE,1))
  p1 <-ggplot(df, aes(E,MSEratio)) + geom_line()+
    geom_hline(yintercept=1.1, linetype="dashed")+
    scale_x_continuous(breaks=seq(0,max(df$E),1), labels=seq(0,max(df$E),1))+
    theme_classic()
  df <-mutate(df,sigdif=ifelse(MSEratio < 1.1 | is.na(MSEratio), "keep lower E", "use higher E"))
  bestE <- min(df$E[which(df$sigdif=="keep lower E")])
  Ediffs <<-df
  multi_return <- function() {
    my_list <- list("bestE" = bestE, "df" = Ediffs, "plot" = p1)
    return(my_list) 
  }
  multi_return()
}
##############################################################################
#################### LIKELIHOOD RATIO TEST FOR COMPARING BESTE_GP AND ARD MODELS #######################

LRT <-function(modBestE, modARD){
  teststat = -2*(modBestE$insampfitstats[4]-modARD$insampfitstats[4])
  degdiff = abs(modBestE$insampfitstats[5]-modARD$insampfitstats[5])
  p.val=pchisq(teststat, df=degdiff, lower.tail = FALSE)
  p.val
}

#############################################################################
########################## function to parse important lags from ARD ################################
#e is the max lag supplied to the model using ARD
#model is the model. 
bestlagARD <-function(e,model){
  phis <-model$pars[1:e]
  phis <-phis[phis>=0.01]
  ees <-names(phis)
  ees <- as.numeric(str_replace(ees,"phi",""))
  bestlags <<-ees
}

####################################################################################################################### 
############
### Function to fit GP to the age classes one by one###
indv_age <-function(df,scaling){
  ages <-unique(df$age_class)
  maxE <-round(sqrt(dim(df)[1]/length(ages)))
  fitstats <-list()
  for (m in 1:length(ages)){
    new_df <-filter(df, age_class==ages[m])
    mod1 <-fitGP(data = new_df, yd = "value", scaling = scaling, E=maxE, tau=1, predictmethod = "loo")
    mod1_out<-c(mod1$outsampfitstats, mod1$insampfitstats,as.character(ages[m]))
    names(mod1_out)<-c("OOS_R2","OOS_rmse","R2","rmse", "ln_post", "lnL_LOO","df","age class")
    fitstats[[m]] <- mod1_out
  }
  fitstats <-bind_rows(fitstats)%>%as.data.frame()%>%mutate(across(OOS_R2:df, as.numeric))
}


##############################################################################################################################
##############################################################################################################################
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
    if ("age" %in% colnames(new_df)){
      new_df <-dplyr::select(new_df, -age)
    }
    
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




### VERSION WHERE MAX E IS SPECIFIED DIFFERENTLY
############################################################################################