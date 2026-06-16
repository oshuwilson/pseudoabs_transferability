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
#read in transferability scores
boyce <- readRDS("output/temporal/transferability_scores.RDS")

#run mixed effect model 
boyce <- boyce %>% mutate_if(is.character, as.factor)
boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "BRT Background"))
boyce <- boyce %>% mutate(code = as.factor(paste(species, site, stage)))
m1 <- lmer(transferability ~ algopseudo + (1|code), data = boyce)
summ(m1)

# repeat mixed effects models for GLMMs
boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "GLMM Background"))
m2 <- lmer(transferability ~ algopseudo + (1|code), data = boyce)
summ(m2)

boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "GLMM CRWs"))
m3 <- lmer(transferability ~ algopseudo + (1|code), data = boyce)
summ(m3)

boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "GLMM Buffer"))
m4 <- lmer(transferability ~ algopseudo + (1|code), data = boyce)
summ(m4)

# mean score per algopseudo
mean_scores <- boyce %>% 
  group_by(algopseudo) %>% 
  summarise(mean_score = mean(transferability, na.rm = T),
            n = n())

# how many for each algopseudo score below 0
subzero <- boyce %>% 
  group_by(algopseudo) %>% 
  summarise(n_below_zero = sum(transferability < 0, na.rm = T),
            n_total = n(),
            prop_below_zero = n_below_zero/n_total,
            prop_below_.5 = sum(transferability < 0.5, na.rm = T)/n_total)

#plot
p1 <- ggplot(boyce, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Temporal Transfer") +
  theme(plot.title = element_text(face = "bold", size = 16))
p1

###Spatial###
rm(list=setdiff(ls(), "p1"))

#read in Boyce scores
boyce <- readRDS("output/spatial/transferability_scores.RDS")

#run mixed effect model 
boyce <- boyce %>% mutate_if(is.character, as.factor)
boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "mSTPP Buffer"))
boyce <- boyce %>% mutate(code = as.factor(paste(species, site, stage)))
m1 <- lmer(transferability ~ algopseudo + (1|code), data = boyce)
summ(m1)

#plot
p2 <- ggplot(boyce, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean, na.rm = T), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Spatial Transfer") +
  theme(plot.title = element_text(face = "bold", size = 16))
p2

###BOTH###
#remove xlab
p1 <- p1 + xlab("")

#plot together
p3 <- plot_grid(p1, p2, ncol=1, align = "hv")
p3 + ggview::canvas(8, 16)

#export png
ggsave("text/figures/ecography/04_ridge_plot.png", p3,
       width = 8, height = 16, dpi = 300)

#export eps
ggsave("text/figures/ecography/04_ridge_plot.eps", p3,
       width = 8, height = 16, dpi = 300)

