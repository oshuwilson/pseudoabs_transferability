#script to produce ridge plots
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(tidyverse)
  library(ggridges)
  library(lme4)
  library(jtools)
  library(cowplot)
}


###Temporal###
#read in Boyce scores and mixed-effects model
boyce <- readRDS("output/leave-year-out/boyce_filtered.RDS")
levels(boyce$algopseudo)
levels(boyce$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW",
                                    "GAM Background", "GAM Buffer", "GAM CRW",
                                    "RF Background", "RF Buffer", "RF CRW")

#run mixed effect model 
boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "RF Background"))
boyce <- boyce %>% mutate_if(is.character, as.factor)
m1 <- lmer(score ~ algopseudo + (1|species) + (1|site) + (1|stage), data = boyce)
summ(m1)

#plot
p1 <- ggplot(boyce, (aes(x=score, y=fct_reorder(algopseudo, score, .fun=median), fill=fct_reorder(algopseudo, score, .fun=median)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_brewer(guide = "none") + scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", size = 0.2) +
  theme(axis.text.y = element_text(vjust = -1.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
p1 

###Spatial###
rm(list=setdiff(ls(), "p1"))

#read in Boyce scores and mixed-effects model
boyce <- readRDS("output/spatial/boyce_final.RDS")
levels(boyce$algopseudo)
levels(boyce$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW",
                              "GAM Background", "GAM Buffer", "GAM CRW",
                              "RF Background", "RF Buffer", "RF CRW")

#run mixed effect model 
boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "GAM Buffer"))
boyce <- boyce %>% mutate_if(is.character, as.factor)
m1 <- lmer(score ~ algopseudo + (1|species) + (1|site) + (1|stage), data = boyce)
summ(m1)

#plot
p2 <- ggplot(boyce, (aes(x=score, y=fct_reorder(algopseudo, score, .fun=median), fill=fct_reorder(algopseudo, score, .fun=median)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + theme_minimal() +
  xlab("Continuous Boyce Index Score") + ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_brewer(guide = "none") + scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", size = 0.2) +
  theme(axis.text.y = element_text(vjust = -1.5),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
p2

###BOTH###
#edit p2
p2 <- p2 + ylab("")

#add titles
p1 <- p1 + ggtitle("Temporal Transfer")
p2 <- p2 + ggtitle("Spatial Transfer")

#change font sizes
p1 <- p1 + theme_minimal(base_size = 7)
p2 <- p2 + theme_minimal(base_size = 7)

#plot together
p3 <- plot_grid(p1, p2, ncol=2, align="h")
p3

#export png
ggsave("text/figures/final figs/03_ridge_plot.png", p3,
       width = 180, height = 80, dpi = 300, units = "mm")

#export eps 
ggsave("text/figures/final figs/03_ridge_plot.eps", p3,
       width = 180, height = 80, dpi = 300, units = "mm")
