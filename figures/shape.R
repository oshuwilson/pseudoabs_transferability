#script to produce plots that compare extrapolation and transferability
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(tidyverse)
  library(lme4)
  library(jtools)
  library(cowplot)
}


###Temporal###
#read in Boyce scores and mixed-effects model
boyce <- readRDS("output/leave-year-out/boyce_filtered.RDS")

#plot - log(median)
p1 <- ggplot(boyce, aes(x=Median, y=score)) + geom_point(colour="lightgrey") + theme_classic() +
  geom_smooth(method="gam", aes(col=pseudo, linetype=algorithm), se=F, linewidth = 1.5) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"), labels=c("Background", "Buffer", "CRW")) + 
  xlab("Median Shape Extrapolation Value") + 
  scale_x_log10(breaks = c(10, 20, 50, 100, 200)) + ylab("Continuous Boyce Index Score") +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
p1
  
###Spatial###
rm(list=setdiff(ls(), "p1"))

#read in Boyce scores and mixed-effects model
boyce <- readRDS("output/spatial/boyce_final.RDS")

#plot - log(median)
p2 <- ggplot(boyce, aes(x=Median, y=score)) + geom_point(colour="lightgrey") + theme_classic() +
  geom_smooth(method="gam", aes(col=pseudo, linetype=algorithm), se=F, linewidth = 1.5) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"), labels=c("Background", "Buffer", "CRW")) + 
  xlab("Median Shape Extrapolation Value") + 
  scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500)) + ylab("Continuous Boyce Index Score") +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
p2

###BOTH###
#modify p1 + p2
p1 <- p1 + guides(linetype = "none", colour = "none")
p2 <- p2 + ylab("")

#add titles
p1 <- p1 + ggtitle("Temporal Transfer")
p2 <- p2 + ggtitle("Spatial Transfer") +
  theme(legend.key.width = unit(3, "line"))

#create multiplot
plots <- plot_grid(p1, p2 + theme(legend.position="none"), 
          align="v", ncol=2)

#extract legend
legend <- get_legend(p2)

#plots and legend together
plot_grid(plots, legend, rel_widths = c(4, .5))

