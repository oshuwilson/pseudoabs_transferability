#plot sample sizes vs boyce and shape
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
boyce <- readRDS("output/leave-year-out/boyce_final.RDS")

#redefine levels
levels(boyce$pseudo) <- c("Background", "Buffer", "CRW")

#plot - sample size
p1 <- ggplot(boyce, aes(x=n, y=score)) + geom_point(colour="lightgrey") + theme_classic() +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=3), 
              aes(col=pseudo, linetype=algorithm), se=F) +
  scale_color_manual(values=c("red3", "steelblue4", "thistle")) + 
  xlab("Training Data Sample Size") + 
  ylab("Continuous Boyce Index Score") +
  scale_x_log10(breaks = c(2000, 20000, 200000), labels = scales::comma) +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
p1

#plot - sample size vs extrapolation
p2 <- ggplot(boyce, aes(x=n, y=Median)) + geom_point(colour="lightgrey") + theme_classic() +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=3), 
              aes(col=pseudo), se=F) +
  scale_color_manual(values=c("red3", "steelblue4", "thistle")) + 
  xlab("Training Data Sample Size") + 
  ylab("Median Shape Value") +
  scale_x_log10(breaks = c(2000, 20000, 200000), labels = scales::comma) +
  guides(colour = guide_legend(title = "Pseudo-Absence"))
p2

###Spatial###
rm(list=setdiff(ls(), c("p1", "p2")))

#read in Boyce scores and mixed-effects model
boyce <- readRDS("output/spatial/boyce_final.RDS")

#redefine levels
levels(boyce$pseudo) <- c("Background", "Buffer", "CRW")

#plot - sample size
p3 <- ggplot(boyce, aes(x=n, y=score)) + geom_point(colour="lightgrey") + theme_classic() +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=3), 
              aes(col=pseudo, linetype=algorithm), se=F) +
  scale_color_manual(values=c("red3", "steelblue4", "thistle"), labels = c("Background", "Buffer", "CRW")) + 
  xlab("Training Data Sample Size") + 
  ylab("Continuous Boyce Index Score") +
  scale_x_log10(breaks = c(3000, 20000, 200000), labels = scales::comma) +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
p3

#plot - sample size vs extrapolation
p4 <- ggplot(boyce, aes(x=n, y=Median)) + geom_point(colour="lightgrey") + theme_classic() +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=3), 
              aes(col=pseudo), se=F) +
  scale_color_manual(values=c("red3", "steelblue4", "thistle"), labels = c("Background", "Buffer", "CRW")) + 
  xlab("Training Data Sample Size") + 
  ylab("Median Shape Value") +
  scale_x_log10(breaks = c(3000, 20000, 200000), labels = scales::comma) +
  guides(colour = guide_legend(title = "Pseudo-Absence"))
p4

###BOTH###

#sample size
#modify p1 and p3
p1 <- p1 + guides(linetype = "none", colour = "none")
p3 <- p3 + ylab("")

#create multiplot
plots <- plot_grid(p1, p3 + theme(legend.position="none"),
                   ncol=2, align="v", labels = "AUTO")
plots

#extract legend
legend <- get_legend(p3)

#plots and legend together
plot_grid(plots, legend, rel_widths = c(4,0.5))

#sample size vs extrapolation
#modify p2 and p4
p2 <- p2 + guides(linetype = "none", colour = "none")
p4 <- p4 + ylab("")

#create multiplot
plots <- plot_grid(p2, p4 + theme(legend.position="none"),
                   ncol=2, align="v", labels = "AUTO")
plots

#extract legend
legend <- get_legend(p4)

#plots and legend together
plot_grid(plots, legend, rel_widths = c(4,0.5))
