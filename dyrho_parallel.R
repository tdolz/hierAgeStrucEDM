##dynamic rho parallel##

#12/31/21
#to calculate lagged correlation and pairwise dynamic rho for each simulation. 
#

library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(rEDM)
devtools::install_github("tanyalrogers/GPEDM") #update frequently
library(GPEDM)
library(forcats)
library(parallel)
library(pracma)
library(plotly)
library(metR)

#
#
#
#
### CREATE PREYLIST ######
######## count peaks function ###############
count_peaks <-function(x){
 mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak. 
 x <-findpeaks(x, minpeakheight = mph)
 x <-dim(x)[1]
 x[is.null(x)] <-0 #if there are no peaks return 0. 
 x
}


#run the simulation and populate a list with the outputs. 
maxiter = 2
preylist<-list()
predlist<-list()
count0peaks<-list()
recnoise <-matrix(nrow=maxiter, ncol=2)
colnames(recnoise)<-c("recnoise1","recnoise2")

for (m in 1:maxiter){
 ### RECRUITMENT OPTIONS ####
 #Species 1 - Prey
 phi1 = 1/1000
 sdrec1 = .02  
 recnoise1 =rnorm(1,0,0.5)
 #recnoise1 = 1
 qfec1 = c(0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1) #prey consumption conversion efficiency (< 1)
 qfec1=qfec1*8 
 CM1 = 0.01 #modifier of the pred_eat_prey matrix 0.01
 #CM1 = 0
 basefec1=c(0,0,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100) # prey
 BM1 = 2 #basefec modifier
 
 #Species 2 - Predator
 phi2 = 1/100
 sdrec2 = .02 
 recnoise2=rnorm(1,0,0.5)
 recnoise2= 1
 qfec2 = c(0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1) #pred consumption conversion efficiency (< 1)
 qfec2=qfec2*8 #8 is good
 CM2 = 0.1 #prey eat pred 0.1
 #CM2 =0
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
 pred_eat_prey=read.csv("predatormatrix.csv", header=FALSE)%>%as.matrix()
 #pred_eat_prey=read.csv("stevespredatormatrix.csv", header=FALSE)%>%as.matrix()
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
 #outputs
 groupmeans <-preypiv %>% group_by(age_class)%>%summarize(preymeans=mean(value))
 predmeans <-predpiv %>% group_by(age_class)%>%summarize(predmeans=mean(value))
 
 mgroups <-min(groupmeans$preymeans)
 mpreds <-min(predmeans$predmeans)
 
 #measure how many peaks in a ten year time period.
 countpeaks <-preypiv %>%filter(time_step > 299 & time_step < 310)%>%
  dplyr::select(-time_step)%>% group_by(age_class)%>%summarize(pks = count_peaks(value))
 count0peaks[[m]] <-sum(countpeaks$pks==0)
 
 #save recruitment noise
 recnoise[m,1]<- recnoise1
 recnoise[m,2]<- recnoise2
 
 if(mgroups > 0.1 & mpreds > 0.1 & sum(countpeaks$pks==0)< 2) {
  preylist[[m]] <- prey
 }else
  preylist[[m]] <- NA
}

## View and save the preylist. 

## create separate dataframe of preylist and the recruitment noise parameters and save separately. 
recnoises <-as.data.frame(recnoise)%>%rownames_to_column(var="index")
preylistswide <-preylist %>% map(~as_tibble(.)) %>% bind_rows(.id="index")%>%as.data.frame()%>%full_join(recnoises)


############################ VISUALLY INSPECT PREYLIST ################################################
preylistswide%>%
 pivot_longer(3:22, names_to = "age_class")%>%
 filter(time_step > 299 & time_step < 310)%>%
 #filter(age_class=="V2")%>%
 #filter(index==1)%>%
 #filter(index <= 25)%>%  #look at 25 plots at a time or it crashes.. 
 #filter(index > 25 & index <= 50)%>%
 #filter(index > 50 & index <= 75)%>%
 #filter(index > 75)%>%
 ggplot(aes(time_step,value))+
 geom_line()+
 facet_wrap(~age_class, scales="free")+
 #facet_wrap(~index)+
 theme_classic()
#it has to show a full cycle. And it does so it's fine.  

############################### 2. CREATE THE FUNCTION FOR LAGGED CORRELATION #############################
plist <-preylist[[1]]

LagCor <- function(plist){
 blockP <-plist[1:21]
 Pblock <-make_block(blockP,max_lag = 20, tau=1)
 Pblock2 <-dplyr::select(Pblock, "V1(t+0)", "V2(t+1)","V3(t+2)","V4(t+3)","V5(t+4)","V6(t+5)","V7(t+6)","V8(t+7)","V9(t+8)",
                         "V10(t+9)","V11(t+10)","V12(t+11)","V13(t+12)","V14(t+13)","V15(t+14)","V16(t+15)","V17(t+16)",
                         "V18(t+17)","V19(t+18)","V20(t+19)")
 names(Pblock2) <-c("age1 (t+0)", "age2 (t+1)","age3(t+2)","age4 (t+3)","age5 (t+4)","age6 (t+5)","age7 (t+6)","age8 (t+7)",
                    "age9 (t+8)","age10 (t+9)","age11 (t+10)","age12 (t+11)","age13 (t+12)","age14 (t+13)","age15 (t+14)",
                    "age16 (t+15)","age17 (t+16)","age18 (t+17)","age19 (t+18)","age20 (t+19)")   
 M <-cor(Pblock2, use="pairwise.complete.obs")
 #corrplot(M, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'black', diag=F,tl.col = "black")
 
 corlist <-t(combn(colnames(M),2))
 M2 =data.frame(corlist, dist=M[corlist])
 length(M2$dist) == length(unique(M2$dist))
 vec <-c(mean(M2$dist),sd(M2$dist))
 #outputs
 gplout <-list(vec,M)
 names(gplout) <-c("summary","lagcorr_matrix")
 gplout
 
}

############################### 3. APPLY THE FUNCTION  #############################
ptest <-LagCor(preylist[[1]])
# test mclapply function
numCores = detectCores()
system.time(lagfoo <-mclapply(preylist, LagCor, mc.cores=numCores))

############################### 4. CORRAL THE LIST OUTPUT #############################. 
#summary stats#
qsum <-lapply(lagfoo, `[`, "summary")
 qsum <-as.data.frame(qsum)
 names(qsum)<-c("mean","sd")
 qsum <- rownames_to_column(qsum,var="iter")

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
names(rhomat_sum)<-c('mean_rho',"sd_rho")

#combine the matrices
LCD <-list()
for(i in 1:length(LC)){
 LCD[[i]]<-LC[[i]]$lagcorr_matrix
}
supermat <-Reduce('+',LCD)/length(LCD)


############################### 5. CREATE THE FUNCTION FOR DYNAMIC RHO #############################
plist <-preylist[[1]]

DYRHO <- function(plist){
 df <-pivot_longer(plist, 2:21, names_to = "age_class")%>%filter(time_step > 299 & time_step <=400) %>%na.omit()%>% as.data.frame()
 maxE =10
 vars = colnames(plist[2:21])
 var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time
 rho_matrix = array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,vars)) 
 for (i in 1:ncol(var_pairs)) {
  df2 = filter(df, age_class %in% c(var_pairs[1,i], var_pairs[2,i]))
  #try({
  fit1 <-fitGP(data = df2, yd = "value", pop="age_class",scaling = "local", E=maxE, tau=1, predictmethod = "loo")
  fit1_rho <-tail(fit1$pars,1)
  #},silent=T)
  rho_matrix[var_pairs[1,i], var_pairs[2,i]] = fit1_rho
 }
 #mean and SD of dynamic correlation
 vec <-c(mean(!is.na(rho_matrix)),sd(as.vector(!is.na(rho_matrix))))
 #outputs
 gplout <-list(vec,rho_matrix)
 names(gplout) <-c("summary","rho_matrix")
 gplout
}

############################### 6. APPLY THE FUNCTION  #############################
ptest2 <-DYRHO(preylist[[1]])
# test mclapply function
numCores = detectCores()
system.time(dyfoo <-mclapply(preylist, DYRHO, mc.cores=numCores))

############################### 7. CORRAL THE LIST OUTPUT #############################. 
#summary stats#
dsum <-lapply(dyfoo, `[`, "summary")
dsum <-as.data.frame(dsum)
names(dsum)<-c("mean","sd")
dsum <- rownames_to_column(dsum,var="iter")

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
