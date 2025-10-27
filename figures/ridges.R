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
boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "BART Background"))
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

#plot
ggplot(boyce, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = pseudo))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_manual(values = c("red3", "steelblue4", "orange2"),
                    guide = "none") +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Temporal Transfer") +
  theme(plot.title = element_text(face = "bold", size = 16))

###Spatial###
rm(list=setdiff(ls(), "p1"))

#read in Boyce scores and mixed-effects model
#boyce <- readRDS("output/spatial/boyce_final.RDS")
boyce <- readRDS("output/thinned_spat2/boyce_final.rds")
#if score is NA, revalue to -1
boyce <- boyce %>% mutate(score = ifelse(is.na(score), -1, score),
                          algopseudo = as.factor(paste(algorithm, pseudo, sep = " ")))

levels(boyce$algopseudo)
levels(boyce$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW",
                              "GAMM Background", "GAMM Buffer", "GAMM CRW",
                              "GLMM Background", "GLMM Buffer", "GLMM CRW",
                              "RF Background", "RF Buffer", "RF CRW")

#run mixed effect model 
boyce <- within(boyce, algopseudo <- relevel(algopseudo, ref = "BRT Background"))
boyce <- boyce %>% mutate_if(is.character, as.factor)
boyce <- boyce %>% mutate(code = as.factor(paste(species, site, stage)))
m1 <- lmer(score ~ algopseudo + (1|code), data = boyce)
summ(m1)

#plot
p2 <- ggplot(boyce, (aes(x=score, y=fct_reorder(algopseudo, score, .fun=median), fill=fct_reorder(algopseudo, score, .fun=median)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_manual(values = c("#F7FBFF", "#E1E9F2", "#CCD6E4",
                               "#B6C4D7", "#A0B1C9", "#8A9FBC",
                               "#758CAE", "#5F7AA1", "#496793",
                               "#335586", "#1E4278", "#08306B"),
                    guide = "none") + 
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", size = 0.2) +
  ggtitle("Spatial Transfer") +
  theme(plot.title = element_text(face = "bold", size = 16))
p2

###BOTH###
#remove xlab
# p2 <- p2 + ylab("")
p1 <- p1 + xlab("")

#plot together
p3 <- plot_grid(p1, p2, ncol=1, align = "hv")
p3

#export png
ggsave("text/figures/ecography/04_ridge_plot.png", p3,
       width = 8, height = 12, dpi = 300)

# #export eps 
# ggsave("text/figures/ecography/04_ridge_plot.eps", p3,
#        width = 8, height = 12, dpi = 300)

