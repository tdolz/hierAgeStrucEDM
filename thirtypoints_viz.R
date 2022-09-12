### Thirtypoints data viz ############

### copied from Sim_Data_Vis100.Rmd on 9/12/22

library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(Metrics)
library(corrplot)
library(kableExtra)
library(metR)
library(cowplot)
library(grid)
library(gridExtra)
library(ggrepel)
library(ggbreak)


points30 <-read.csv("modelcomparison_outputs/thirtypoints.csv", header=T)

points30$approach <-fct_recode(points30$approach, "TA"="fit_NT_stats","HA"="fitsumStats", "Hierarchical overall"="fit_all_stats")


apocolr =c("Total Abundance"="#b2df8a","Hierarchical aggregate"="#1f78b4", "Hierarchical overall"="#a6cee3") 

#if you want to change to hierarchical all instead of hierarchical summed
#points30$approach <-fct_recode(points30$approach, "Total Abundance"="fit_NT_stats","Hierarchical"="fit_all_Stats")


points30%>%
 filter(rowname=="R2")%>%
 #filter(approach !="fitsumStats")%>%
 ggplot(aes(approach, value, fill=approach))+
 geom_point(aes(y = value, color =approach), 
            position = position_jitterdodge(), 
            size = 1, alpha = 0.2) +
 geom_boxplot(outlier.shape = NA, alpha = 0.8,
              position=position_dodge2(padding=0.5)) + 
 scale_fill_manual(values = apocolr)+
 scale_color_manual(values = apocolr)+
 ylab("Out-of-sample r-squared")+xlab("Approach")+
 facet_wrap(~sim, scales="free",ncol=3)+
 guides(fill=guide_legend(title="Analysis"), color=guide_legend(title="Analysis"))+
 theme(
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  axis.title.x=element_text(size=12,face="bold",margin=margin(t=8,unit="pt")),
  axis.title.y=element_text(size=12,face="bold",margin=margin(r=8,unit="pt")),
  strip.text = element_text(size=12),
  legend.key = element_rect(fill = "transparent"),
  legend.text= element_text(size=12),
  legend.title=element_text(size=12, face="bold"),
  legend.position = "bottom",
  strip.background =element_rect(fill="white"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = "dark grey"),
  plot.background = element_rect(fill = "white",colour = "white", size=0.5),
  panel.border = element_rect(colour="black",fill=NA, size=0.3))

#ggsave("thirtypointsfig_3.png",path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
#dev.off()
#
#

#Thirty points figure plotted separately and reassembled in InDesign**  
 #with ggbreak for axis breaks  
#https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html#feature-12-compatible-with-patchwork 

# there is probably a way to do scale-break with facets but I am just going to do it manually for now. 
#apocolr =c("TA"="#b2df8a","Hier"="#1f78b4") 

#points30$approach <-fct_recode(points30$approach, "TA"="fit_NT_stats","HA"="fitsumStats", "Hierarchical overall"="fit_all_stats")

points30$approach <-fct_recode(points30$approach, "TA"="Total Abundance","HO"="Hierarchical overall", "HA"="Hierarchical aggregate")
apocolr =c("TA"="#b2df8a","HA"="#1f78b4", "HO"="#a6cee3") 
# 
# 
### I
points30%>%
 filter(rowname=="R2")%>%
 #filter(approach !="fitsumStats")%>%
 filter(approach !="fit_all_stats")%>%
 filter(sim=="I")%>%
 ggplot(aes(approach, value, fill=approach))+
 geom_point(aes(y = value, color =approach), 
            position = position_jitterdodge(), 
            size = 1, alpha = 0.8) +
 geom_boxplot(width=0.3, outlier.shape = NA, alpha = 0.8,
              position=position_dodge2(padding=0.5)) + 
 scale_fill_manual(values = apocolr)+
 scale_color_manual(values = apocolr)+
 ylab("")+xlab("")+
 #facet_wrap(~sim, scales="free",ncol=3)+
 #guides(fill=guide_legend(title="Analysis"), color=guide_legend(title="Analysis"))+
 theme(
  axis.text.x = element_text(size=16),
  axis.text.y = element_text(size=16),
  axis.title.x=element_text(size=16,face="bold",margin=margin(t=8,unit="pt")),
  axis.title.y=element_text(size=16,face="bold",margin=margin(r=8,unit="pt")),
  strip.text = element_text(size=16),
  legend.key = element_rect(fill = "transparent"),
  legend.text= element_text(size=16),
  legend.title=element_text(size=16, face="bold"),
  legend.position = "none",
  strip.background =element_rect(fill="white"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = "dark grey"),
  plot.background = element_rect(fill = "white",colour = "white", size=0.5),
  panel.border = element_rect(colour="black",fill=NA, size=0.3))
#ggsave("thirtypointsfigI_3.png",height=6, width=4,dpi=300, path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
ggsave("figures/thirtypointsfigI_3.png",height=6, width=4,dpi=300)

### II
points30%>%
 filter(rowname=="R2")%>%
 #filter(approach !="fitsumStats")%>%
 filter(approach !="fit_all_stats")%>%
 filter(sim=="II")%>%
 ggplot(aes(approach, value, fill=approach))+
 geom_point(aes(y = value, color =approach), 
            position = position_jitterdodge(), 
            size = 1, alpha = 0.8) +
 geom_boxplot(width=0.3, outlier.shape = NA, alpha = 0.8,
              position=position_dodge2(padding=0.5)) + 
 scale_fill_manual(values = apocolr)+
 scale_color_manual(values = apocolr)+
 #scale_y_break(c(0.5,0.7), scales=5)+
 scale_y_break(c(0.8,0.925), scales=2)+
 ylab("")+xlab("")+
 #facet_wrap(~sim, scales="free",ncol=3)+
 #guides(fill=guide_legend(title="Analysis"), color=guide_legend(title="Analysis"))+
 theme(
  axis.text.x = element_text(size=16),
  axis.text.y = element_text(size=16),
  axis.title.x=element_text(size=16,face="bold",margin=margin(t=8,unit="pt")),
  axis.title.y=element_text(size=16,face="bold",margin=margin(r=8,unit="pt")),
  strip.text = element_text(size=16),
  legend.key = element_rect(fill = "transparent"),
  legend.text= element_text(size=16),
  legend.title=element_text(size=16, face="bold"),
  legend.position = "none",
  strip.background =element_rect(fill="white"),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = "dark grey"),
  plot.background = element_rect(fill = "white",colour = "white", size=0.5),
  panel.border = element_rect(colour="black",fill=NA, size=0.3))
#ggsave("thirtypointsfigII_3.png",height=6, width=4,dpi=300, path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
ggsave("figures/thirtypointsfigII_3.png",height=6, width=4,dpi=300)

### III
points30%>%
 filter(rowname=="R2")%>%
 filter(approach !="fit_all_stats")%>%
 #filter(approach !="fitsumStats")%>%
 filter(sim=="III")%>%
 ggplot(aes(approach, value, fill=approach))+
 geom_point(aes(y = value, color =approach), 
            position = position_jitterdodge(), 
            size = 1, alpha = 0.8) +
 geom_boxplot(width=0.3, outlier.shape = NA, alpha = 0.8,
              position=position_dodge2(padding=0.5)) + 
 scale_fill_manual(values = apocolr)+
 scale_color_manual(values = apocolr)+
 #scale_y_break(c(0.94,0.96), scales=2)+
 #scale_y_break(c(0.36,0.92), scales=5)+
 scale_y_break(c(0.94,0.975), scales=4)+
 #ylim(0.9,1)+
 ylab("")+xlab("")+
 #facet_wrap(~sim, scales="free",ncol=3)+
 #guides(fill=guide_legend(title="Analysis"), color=guide_legend(title="Analysis"))+
 theme(
  axis.text.x = element_text(size=16),
  axis.text.y = element_text(size=16),
  axis.title.x=element_text(size=16,face="bold",margin=margin(t=8,unit="pt")),
  axis.title.y=element_text(size=16,face="bold",margin=margin(r=8,unit="pt")),
  strip.text = element_text(size=16),
  legend.key = element_rect(fill = "transparent"),
  legend.text= element_text(size=16),
  legend.title=element_text(size=16, face="bold"),
  legend.position = "none",
  strip.background =element_rect(fill="white"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = "dark grey"),
  plot.background = element_rect(fill = "white",colour = "white", size=0.5),
  panel.border = element_rect(colour="black",fill=NA, size=0.3))
#ggsave("thirtypointsfigIII_3.png",height=6, width=4,dpi=300, path="/Users/tdolan/documents/postdoc/age structure/agestructfigs")
ggsave("figures/thirtypointsfigIII_3.png",height=6, width=4,dpi=300)
