### dynamic rho visualization ####
### 
### created 8/17/22 to visualize the outputs from dyrho_parallel.R
### this is only the simulated corrplots for lagged correlation and dynamic rho correlation
### the empirical plots are in the empirical_analysis.Rmd script
### 

# load packages
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(corrplot)
library(RColorBrewer)
library(viridisLite)

######################import the data#####################################################################################
dystats <-read.csv("dyrho_sim_stats.csv",header=T)
lagcorstats<-read.csv("lagged_correlation_sim_stats.csv",header=T)
dyrho_matrix <-read.csv("dyrho_supermatrices.csv",header=T)
lagcor_matrix <-read.csv("lagged_correlation_supermatrices.csv",header=T)


########################## SUMMARY STATISTICS #################################################################################
sumstats_dy <-dystats %>%group_by(sim)%>% summarize(mn=mean(mean), sem =sd(mean), msd=mean(sd_dyrho))
sumstats_dy 

sumstats_lagcor <-lagcorstats %>%group_by(sim)%>% summarize(mn=mean(mean_pairwisedist), sem =sd(mean_pairwisedist), msd=mean(sd_pairwisedist))
sumstats_lagcor 

######################## LAGGED CORRELATION CORRPLOT #########################################################################
lagged_ages <-c("age1 (t+0)", "age2 (t+1)","age3(t+2)","age4 (t+3)","age5 (t+4)","age6 (t+5)","age7 (t+6)","age8 (t+7)",
                "age9 (t+8)","age10 (t+9)","age11 (t+10)","age12 (t+11)","age13 (t+12)","age14 (t+13)","age15 (t+14)",
                "age16 (t+15)","age17 (t+16)","age18 (t+17)","age19 (t+18)","age20 (t+19)")

##export as .png with height =1000 and width = 1000##

### SIM I #############
lcm1 <-filter(lagcor_matrix, sim=="I")%>%select(-sim,-X)%>%as.matrix()
colnames(lcm1)<-lagged_ages
rownames(lcm1)<-lagged_ages
corrplot(lcm1, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'grey', diag=F,tl.col = "black")
#rescale
corrplot(lcm1, type="upper", method="color", title="", mar=c(0,0,2,0), 
         addCoef.col = 'grey', 
         diag=F,tl.col = "black",is.corr = T, col.lim = c(0.9, 1.0),
         #col=colorRampPalette(c("white","navy"),alpha=FALSE,space = "Lab")(50))
         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(50),cl.pos="n")


### SIM II ##############
lcm2 <-filter(lagcor_matrix, sim=="II")%>%select(-sim,-X)%>%as.matrix()
colnames(lcm2)<-lagged_ages
rownames(lcm2)<-lagged_ages
corrplot(lcm2, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'light grey', diag=F,tl.col = "black")
#rescale
corrplot(lcm2, type="upper", method="color", title="", mar=c(0,0,2,0), 
         addCoef.col = 'black', 
         diag=F,tl.col = "black",is.corr = F, col.lim = c(0.9, 1.0),
         #col=colorRampPalette(c("white","navy"),alpha=FALSE,space = "Lab")(50))
         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(50),cl.pos="n")



### SIM III ################
lcm3 <-filter(lagcor_matrix, sim=="III")%>%select(-sim,-X)%>%as.matrix()
colnames(lcm3)<-lagged_ages
rownames(lcm3)<-lagged_ages
corrplot(lcm3, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'grey', diag=F,tl.col = "black")
#rescale
corrplot(lcm3, type="upper", method="color", title="", mar=c(0,0,2,0), 
         addCoef.col = 'grey', 
         diag=F,tl.col = "black",is.corr = T, col.lim = c(0.9, 1.0),
         #col=colorRampPalette(c("white","navy"),alpha=FALSE,space = "Lab")(50))
         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(50),cl.pos="n")


######################## DYNAMIC RHO CORRPLOT #########################################################################
rho_ages <-c("age1","age2","age3","age4","age5","age6","age7","age8","age9","age10","age11","age12","age13","age14","age15",
             "age16","age17","age18","age19","age20")

##export as .png with height =1000 and width = 1000##


### SIM I #############
drm1 <-filter(dyrho_matrix, sim=="I")%>%select(-sim,-X)%>%as.matrix()
colnames(drm1)<-rho_ages
rownames(drm1)<-rho_ages
corrplot(drm1, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'black', diag=F,tl.col = "black",is.corr = TRUE)
#rescaled
corrplot(drm1, type="upper", method="shade", title="", mar=c(0,0,2,0), 
         addCoef.col = 'grey', 
         diag=F,tl.col = "black",is.corr = T, 
         col.lim = c(0.9, 1.0),
        #col=colorRampPalette(c("white","navy"),alpha=FALSE)(50))
col = colorRampPalette(c("white", "deepskyblue", "blue4"), alpha=F)(50),cl.pos="n")


### SIM II #############
drm2 <-filter(dyrho_matrix, sim=="II")%>%select(-sim,-X)%>%as.matrix()
colnames(drm2)<-rho_ages
rownames(drm2)<-rho_ages
corrplot(drm2, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'black', diag=F,tl.col = "black",is.corr = TRUE)
#rescaled
corrplot(drm2, type="upper", method="color", title="", mar=c(0,0,2,0), 
         addCoef.col = 'black', 
         diag=F,tl.col = "black",is.corr = F, col.lim = c(0.9, 1.0),
         #col=colorRampPalette(c("white","navy"),alpha=FALSE,space = "Lab")(50))
         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(50),cl.pos="n")



### SIM III #############
drm3 <-filter(dyrho_matrix, sim=="III")%>%select(-sim,-X)%>%as.matrix()
colnames(drm3)<-rho_ages
rownames(drm3)<-rho_ages
corrplot(drm3, type="upper", method="color", title="", mar=c(0,0,2,0), addCoef.col = 'black', diag=F,tl.col = "black",is.corr = TRUE)
#rescaled
corrplot(drm3, type="upper", method="color", title="", mar=c(0,0,2,0), 
         addCoef.col = 'black', 
         diag=F,tl.col = "black",is.corr = F, col.lim = c(0.9, 1.0),
         #col=colorRampPalette(c("white","navy"),alpha=FALSE,space = "Lab")(50))
         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(50),cl.pos="n")

