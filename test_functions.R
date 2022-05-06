#test period function. 
library(GPEDM)
source("analysis_functions.R")

###### this measures the period of the oscillating simulated time series 
###### used in "make_simulation_data.R"
period <-function(x){
 #mph <-max(x)/10 #the minimum peak height has to be 1/10 the height of the peak.
 mph <-mean(x)+sd(x) #the minimum peak height has to be 1 standard deviation from the mean
 x <-findpeaks(x, minpeakheight = mph)
 x <-x[,1:2]
 x <-as.data.frame(x)
 names(x) <-c("ypeak","xpeak")
 x <-mutate(x, dist = xpeak-lag(xpeak,1))
 mean.period = mean(x$dist,na.rm=T)
 mean.period
}


data("thetalog2pop")

pA=subset(thetalog2pop,Population=="PopA")
N=nrow(pA)
plot(Abundance~Time,data=pA,type="l",main="PopA")

meanper <-pA %>%
 dplyr::select(-Time)%>% group_by(Population)%>%summarize(periodt=period(Abundance))
meanper

cycle(pA$Abundance)
frequency(pA$Abundance)

