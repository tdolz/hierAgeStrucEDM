##dynamic rho parallel##

#12/31/21
#to calculate lagged correlation and pairwise dynamic rho for each simulation. 

#### UPDATED 8/17/22 ####

library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(rEDM)
#devtools::install_github("tanyalrogers/GPEDM") #update frequently
library(GPEDM)
library(forcats)
library(parallel)
library(pracma)
library(plotly)
library(metR)
library(corrplot)


############################### 1. IMPORT AND FORMAT SIM DATA #############################

## Read in the simulation data & make it long format ## BUT TURN IT INTO A LIST!!! 
preylist1 <-read.csv("simulated_data/Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist1 <- preylist1 %>% split(f=preylist1$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist2 <-read.csv("simulated_data/Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist2 <- preylist2 %>% split(f=preylist2$index)%>%lapply(function(x) x[!names(x) %in% c("index")])
preylist3 <-read.csv("simulated_data/Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist3 <- preylist3 %>% split(f=preylist3$index)%>%lapply(function(x) x[!names(x) %in% c("index")])


############################### 2. CREATE THE FUNCTION FOR LAGGED CORRELATION 30 YEARS OF DATA #############################

LagCor <- function(plist){
 blockP <-plist[,1:21]
 blockP <-filter(blockP, time_step > 330) # 30 years of data only. 
 Pblock <-make_block(blockP,max_lag = 20, tau=1)
 Pblock2 <-dplyr::select(Pblock, "V1(t+0)", "V2(t+1)","V3(t+2)","V4(t+3)","V5(t+4)","V6(t+5)","V7(t+6)","V8(t+7)","V9(t+8)",
                         "V10(t+9)","V11(t+10)","V12(t+11)","V13(t+12)","V14(t+13)","V15(t+14)","V16(t+15)","V17(t+16)",
                         "V18(t+17)","V19(t+18)","V20(t+19)")
 names(Pblock2) <-c("age1 (t+0)", "age2 (t+1)","age3(t+2)","age4 (t+3)","age5 (t+4)","age6 (t+5)","age7 (t+6)","age8 (t+7)",
                    "age9 (t+8)","age10 (t+9)","age11 (t+10)","age12 (t+11)","age13 (t+12)","age14 (t+13)","age15 (t+14)",
                    "age16 (t+15)","age17 (t+16)","age18 (t+17)","age19 (t+18)","age20 (t+19)")   
 M <-cor(Pblock2, use="pairwise.complete.obs")
 #corrplot(M, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'black', diag=F,tl.col = "black")
 
 M2 =data.frame(col=colnames(M)[col(M)], row=rownames(M)[row(M)], dist=c(M))
 M2 <-na.omit(M2)
 length(M2$dist) == length(unique(M2$dist))
 vec <-c(mean(M2$dist),sd(M2$dist))
 #outputs
 gplout <-list(vec,M)
 names(gplout) <-c("summary","lagcorr_matrix")
 gplout
}

############################### 3. FUNCTION TO ORGANIZE THE LAGGED CORR OUTPUT #############################. 
#summary stats#

SumStats_corr <-function(lagfoo){
qsum <-lapply(lagfoo, `[`, "summary")  #the summary is the mean and sd of paired distance
qstats<-matrix(nrow=100,ncol=2)
for (i in 1:length(qsum)){
  qstats[i,1]<-qsum[[i]]$summary[1]
  qstats[i,2]<-qsum[[i]]$summary[2]
}
qstats <-as.data.frame(qstats)
names(qstats)<-c("mean_pairwisedist","sd_pairwisedist")
qstats <- rownames_to_column(qstats,var="iter")

#correlation matrix --> collapse into one summary correlation matrix. 
LC <-lapply(lagfoo, `[`, "lagcorr_matrix")

#collapse the matrices
lcmeans <-c()
lcsd <-c()
for (i in 1:length(LC)){
  m2 <-LC[[i]]$lagcorr_matrix
  m2[lower.tri(m2)]<-NA
  lcmeans[i] <-mean(!is.na(m2))
  m2vec <-as.vector(m2)
  m2vec <-m2vec[!is.na(m2vec)]
  lcsd[i]<-sd(m2vec)
}
rhomat_sum <-as.data.frame(cbind(lcmeans,lcsd))
names(rhomat_sum)<-c('mean_rho',"sd_rho")  #this is the mean and sd correlations for the whole matrix. 
qstats <-bind_cols(qstats,rhomat_sum)
qstats
}

Supermat_corr <-function(lagfoo){
#correlation matrix --> collapse into one summary correlation matrix. 
LC <-lapply(lagfoo, `[`, "lagcorr_matrix")

#combine the matrices
LCD <-list()
for(i in 1:length(LC)){
  LCD[[i]]<-LC[[i]]$lagcorr_matrix
}
supermat <-Reduce('+',LCD)/length(LCD)  #this is the summary matrix from all the iterations
supermat
}

############################### 4. FUNCTION FOR DYNAMIC RHO 30 YEARS OF DATA ONLY #############################

DYRHO <- function(plist){
  df <-pivot_longer(plist, 2:21, names_to = "age_class")%>%filter(time_step > 299 & time_step < 330) %>%na.omit()%>% as.data.frame()
  maxE =5 #because sqrt of 30
  vars = colnames(plist[2:21])
  var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time
  rho_matrix = array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,vars)) 
  for (i in 1:ncol(var_pairs)) {
    df2 = filter(df, age_class %in% c(var_pairs[1,i], var_pairs[2,i]))
    #try({
    fit1 <-fitGP(data = df2, y = "value", pop="age_class",scaling = "local", E=maxE, tau=1, predictmethod = "loo")
    fit1_rho <-tail(fit1$pars,1)
    #},silent=T)
    rho_matrix[var_pairs[1,i], var_pairs[2,i]] = fit1_rho
  }
  #mean and SD of dynamic correlation
  M2 =data.frame(col=colnames(rho_matrix)[col(rho_matrix)], row=rownames(rho_matrix)[row(rho_matrix)], dist=c(rho_matrix))
  M2 <-na.omit(M2)
  vec <-c(mean(M2$dist),sd(M2$dist))
  #outputs
  gplout <-list(vec,rho_matrix)
  names(gplout) <-c("summary","rho_matrix")
  gplout
}

############################### 5. FUNCTION TO ORGANIZE THE DYNAMIC RHO OUTPUT #############################. 
#summary stats#
SumStats_dy <-function(dyfoo){
dsum <-lapply(dyfoo, `[`, "summary")
dstats<-matrix(nrow=100,ncol=2)
for (i in 1:length(dsum)){
  dstats[i,1]<-dsum[[i]]$summary[1]
  dstats[i,2]<-dsum[[i]]$summary[2]
}
dstats <-as.data.frame(dstats)
names(dstats)<-c("mean","sd")
dstats <- rownames_to_column(dstats,var="iter")
dstats


#correlation matrix --> collapse into one summary correlation matrix. 
DC <-lapply(dyfoo, `[`, "rho_matrix")

#collapse the matrices
dcmeans <-c()
dcsd <-c()
for (i in 1:length(DC)){
  m2 <-DC[[i]]$rho_matrix
  m2[lower.tri(m2)]<-NA
  dcmeans[i] <-mean(!is.na(m2))
  m2vec <-as.vector(m2)
  m2vec <-m2vec[!is.na(m2vec)]
  dcsd[i]<-sd(m2vec)
}
dyrhomat_sum <-as.data.frame(cbind(dcmeans,dcsd))
names(dyrhomat_sum)<-c('mean_dyrho',"sd_dyrho")
dyrhomat_sum
dystats <-bind_cols(dstats,dyrhomat_sum)
}

Supermat_dyrho <-function(dyfoo){
  #correlation matrix --> collapse into one summary correlation matrix. 
  DC <-lapply(dyfoo, `[`, "rho_matrix")
  
  #combine the matrices
  LCD <-list()
  for(i in 1:length(DC)){
    LCD[[i]]<-DC[[i]]$rho_matrix
  }
  supermat <-Reduce('+',LCD)/length(LCD)  #this is the summary matrix from all the iterations
  supermat
}


############################################## 6. APPLY ##################################################################################

############################### SIM 1 #############################
### Apply Corr function ###
numCores = detectCores()
system.time(lagfoo <-mclapply(preylist1, LagCor, mc.cores=numCores))
### organize function output ####
lagcorstats1 <-SumStats_corr(lagfoo)
lagcormatrix1 <-Supermat_corr(lagfoo)

### Apply Dyrho function ####
numCores = detectCores()
system.time(dyfoo <-mclapply(preylist1, DYRHO, mc.cores=numCores))
### organize function output ###
dyrhostats1 <-SumStats_dy(dyfoo)
dyrhomatrix1 <-Supermat_dyrho(dyfoo)
#######################################################################################


############################### SIM 2 #############################
### Apply Corr function ###
numCores = detectCores()
system.time(lagfoo <-mclapply(preylist2, LagCor, mc.cores=numCores))
### organize function output ####
lagcorstats2 <-SumStats_corr(lagfoo)
lagcormatrix2 <-Supermat_corr(lagfoo)

### Apply Dyrho function ####
numCores = detectCores()
system.time(dyfoo <-mclapply(preylist2, DYRHO, mc.cores=numCores))
### organize function output ###
dyrhostats2 <-SumStats_dy(dyfoo)
dyrhomatrix2 <-Supermat_dyrho(dyfoo)
#######################################################################################


############################### SIM 3 #############################
### Apply Corr function ###
numCores = detectCores()
system.time(lagfoo <-mclapply(preylist3, LagCor, mc.cores=numCores))
### organize function output ####
lagcorstats3 <-SumStats_corr(lagfoo)
lagcormatrix3 <-Supermat_corr(lagfoo)

### Apply Dyrho function ####
numCores = detectCores()
system.time(dyfoo <-mclapply(preylist3, DYRHO, mc.cores=numCores))
### organize function output ###
dyrhostats3 <-SumStats_dy(dyfoo)
dyrhomatrix3 <-Supermat_dyrho(dyfoo)
#######################################################################################

############################################## 7. ORGANIZE & SAVE THE CSVS #######################################

#combine lag stats into one csv
lagcorstats1 <-mutate(lagcorstats1, sim="I")
lagcorstats2 <-mutate(lagcorstats2, sim="II")
lagcorstats3 <-mutate(lagcorstats3, sim="III")
lagcorstats <-bind_rows(lagcorstats1,lagcorstats2,lagcorstats3)
#write csv
write.csv(lagcorstats, "pairwise_rho_outputs/lagged_correlation_sim_stats.csv")

#combine lagged cor supermatrix into one csv
lcm1 <-as.data.frame(lagcormatrix1)%>%mutate(sim="I")
lcm2 <-as.data.frame(lagcormatrix2)%>%mutate(sim="II")
lcm3 <-as.data.frame(lagcormatrix3)%>%mutate(sim="III")
lcm <-bind_rows(lcm1,lcm2,lcm3)
#write csv
write.csv(lcm, "pairwise_rho_outputs/lagged_correlation_supermatrices.csv")

#combine dyrho stats into one csv
dyrhostats1 <-mutate(dyrhostats1, sim="I")
dyrhostats2 <-mutate(dyrhostats2, sim="II")
dyrhostats3 <-mutate(dyrhostats3, sim="III")
dyrhostats <-bind_rows(dyrhostats1,dyrhostats2,dyrhostats3)
#write csv
write.csv(dyrhostats, "pairwise_rho_outputs/dyrho_sim_stats.csv")

#combine lagged cor supermatrix into one csv
drm1 <-as.data.frame(dyrhomatrix1)%>%mutate(sim="I")
drm2 <-as.data.frame(dyrhomatrix2)%>%mutate(sim="II")
drm3 <-as.data.frame(dyrhomatrix3)%>%mutate(sim="III")
drm <-bind_rows(drm1,drm2,drm3)
#write csv
write.csv(drm, "pairwise_rho_outputs/dyrho_supermatrices.csv")

