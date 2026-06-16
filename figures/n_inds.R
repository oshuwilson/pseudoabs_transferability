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
#read in Boyce scores 
boyce <- readRDS("output/temporal/transferability_scores.RDS")

# reorder algorithm to group machine learning models, regression models, and spatial correlation models
boyce$algorithm <- factor(boyce$algorithm, levels = c("GLMM", "GAMM", "MaxEnt", "BART", "BRT", "RF", "mSTPP", "INLA-SPDE"))


#plot - sample size
p1 <- ggplot(boyce, aes(x=n_ind, y=transferability, group = algopseudo)) + theme_classic() +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm", 
              aes(col=pseudo, linetype=algorithm), se=F, linewidth=0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"),
                     labels = c("Background", "Buffer", "CRW")) + 
  xlab("Number of Individuals in Training Data") + 
  scale_x_log10(breaks = c(3, 5, 8, 10, 20, 30, 50)) +
  ylab("Continuous Boyce Index Score") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
  theme(legend.key.width = unit(1, "cm"))
p1 + ggview::canvas(8, 6)

p1 <- ggplot(boyce, aes(x=n_ind, y=transferability, group = algorithm)) +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm",
              aes(col=algorithm), se=F, linewidth=0.6) + 
  scale_color_manual(values = c("red1", "red3", "red4", "steelblue1", "steelblue3", "steelblue4", "orange1", "orange3")) +
  xlab("Number of Individuals in Training Data") + 
  ylab("Continuous Boyce Index Score") +
  scale_x_log10(breaks = c(3, 5, 10, 20, 30, 50), labels = scales::comma, expand = c(0.01,0.01)) +
  facet_wrap(~pseudo) +
  theme_minimal() +
  geom_hline(yintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
  theme(legend.key.width = unit(1, "cm"),
        panel.border = element_rect(color = "grey40", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Algorithm")) 
p1 + ggview::canvas(12,7) 

n_inds <- boyce %>% 
  group_by(n_ind, algopseudo) %>%
  summarise(mean_perf = mean(transferability, na.rm = T))

# export plot
ggsave("text/Ecological Modelling Submission/figs/supp/n_inds.png",
       width = 12, height = 7, units = "in", dpi = 300, p1)
