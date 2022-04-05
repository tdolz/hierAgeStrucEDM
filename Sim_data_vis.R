#### Simulation Data Vis ####
#### 3/18/22
#### Produces most of the visualizations for the manuscript, excluding some of the empirical data visualizations. 
#### 

library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(Metrics)
library(corrplot)
library("plotly")
library(kableExtra)
library(metR)
library(cowplot)
library(grid)
library(gridExtra)
library(ggrepel)
#geom_flat_violin (if needed)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

##############################################################################################################################
############# Tileplot of ARD values for hyperparameters: Figure S1 ############

phi_Sim1 <-read.csv("Sim1Phis_100data.csv", header=T)%>%mutate(sim="I")
phi_Sim2 <-read.csv("Sim2Phis_100data.csv", header=T)%>%mutate(sim="II")
phi_Sim3 <-read.csv("Sim3Phis_100data.csv", header=T)%>%mutate(sim="III")

#create giant dataframe from all simulations
phis <-bind_rows(phi_Sim1, phi_Sim2, phi_Sim3)%>%dplyr::select(-X)
phis <-phis[,c(1,2,3,7,8,9,4,5,6,10,12,13,14,11)]
pivphi <-pivot_longer(phis, 2:6, names_to = "model")%>%as.data.frame()

#create pivot version where it's the mean only. 
pivphimean <-pivphi %>% group_by(rowname,model,sim)%>% summarize(phimean=mean(value),phisd=sd(value),.groups="drop")
pivphimean[is.na(pivphimean)] <-0
#pivphimean$model <-fct_relevel(pivphimean$model,"phi10", after = 9)
pivphimean <-mutate(pivphimean, tslength=ifelse(rowname=="N_total",100,substr(rowname,1,2)))%>%mutate(tslength=ifelse(tslength=="5y",5,tslength))%>%mutate(maxE=round(sqrt(as.numeric(tslength))))

#discretize the color scale
cols <- c("(0,0.2]"="#f1eef6", "(0.2,0.4]" = "#d0d1e6", "(0.4,0.6]" = "#a6bddb", "(0.6,0.8]" = "#74a9cf", "(0.8,1.0]"  = "#2b8cbe", "(1.0,1.2]"="#045a8d", "(1.2,1.4]" ="#034e7b", "(1.4,1.6]" ="#011f4b")
pivphimean<-mutate(pivphimean, phimean_cut = cut(phimean, breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6)))

# need to fix the legend on this plot. Also maybe the low value color shouldn't be white because it's hard to see. 
ggplot(pivphimean, aes(rowname,model))+
 geom_tile(aes(fill=phimean_cut), show.legend = TRUE)+
 #scale_fill_gradient(low="white",high="blue")+
 scale_fill_manual(values=cols, na.value = "white")+
 xlab("")+ylab("")+coord_flip()+
 guides(fill=guide_legend(title="mean value of phi"))+
 facet_wrap(~sim)+
 theme_bw()+
 theme(axis.text.x = element_text(angle = 90),panel.background = element_rect(fill = "transparent", colour = "black"),
       strip.background =element_rect(fill="white"),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank())
ggsave("phiplot.png", path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
dev.off()

##############################################################################################################################


### Tanya request figure--- Simulation data with prediction

#load preylists and take only the first simulation run. 
preylist1 <-read.csv("Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
        filter(index=="1")%>%
        pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()
preylist2 <-read.csv("Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
        filter(index=="1")%>%
        pivot_longer(3:23, names_to = "age_class")%>%as.data.frame()
preylist3 <-read.csv("Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)%>%
        filter(index=="1")%>%
        pivot_longer(3:22, names_to = "age_class")%>%as.data.frame()

#Function to make 30 year datasets + ten year test dataset to fit GPs#
### Remember to log transform the value ###
p1 <-filter(preylist1, time_step >=300 & time_step <= 340)%>% as.data.frame()%>%mutate(value=log(value))
p1Lags = makelags(data=p1, yd="value", pop="age_class", E=round(sqrt(30)), tau=1)
p1 = cbind(p1,p1Lags)
p1.train = filter(p1, time_step <= (max(p1$time_step)-10))
p1.test = filter(p1, time_step > (max(p1$time_step)-10))

p1gp <-fitGP(data = p1.train, yd = "value", xd=colnames(p1Lags),datanew=p1.test,pop="age_class",scaling = "local",predictmethod = "loo")


