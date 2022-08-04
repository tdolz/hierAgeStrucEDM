### Created 8/3/22 ####
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(cowplot)
library(rEDM)
library(GPEDM)
library(Metrics)
library(corrplot)
library(pracma)
library(parallel)

### Tanya request figure--- Simulation data with prediction
### 
#this is basically a 30 points data comp
#A figure comparing R2 values from the total abundance model with R2 for the hier model aggregated to total abundance, 
#for each simulation and empirical dataset. All would have 30 time points. For the simulations, there would only be one value per replicate, 
#so you could show the distribution among replicates 

#load preylists and take only the first simulation run. 
preylist1 <-read.csv("Simulation1_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist2 <-read.csv("Simulation2_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)
preylist3 <-read.csv("Simulation3_data.csv", header=T,row.names=NULL)%>%dplyr::select(-X)


GP30 <-function(plist){
  
}