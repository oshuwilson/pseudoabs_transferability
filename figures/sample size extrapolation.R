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

#plot - sample size
ggplot(boyce, aes(x=n, y=transferability, group = algopseudo)) + theme_classic() +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm", 
              aes(col=pseudo, linetype=algorithm), se=F, linewidth=0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"),
                     labels = c("Background", "Buffer", "CRW")) + 
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Continuous Boyce Index Score") +
  scale_x_log10(breaks = c(500, 2000, 10000), labels = scales::comma) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash", "longdash")) +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))

# plot - extrapolation
ggplot(boyce, aes(x=Median, y=transferability)) + 
  geom_point(colour="lightgrey", size = 0.1) + 
  theme_classic() +
  geom_smooth(method="lm", aes(col=pseudo, linetype=algorithm), se=F, linewidth = 0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"), labels=c("Background", "Buffer", "CRW")) + 
  xlab("Median Shape Extrapolation Value") + 
  scale_x_log10(breaks = c(10, 20, 50, 100, 200)) + ylab("Continuous Boyce Index Score") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash", "longdash")) +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
p2

#plot - sample size vs extrapolation
p3 <- ggplot(boyce, aes(x=n, y=Median)) + 
  geom_point(colour="lightgrey", size = 0.1) + theme_classic() +
  geom_smooth(method="lm", 
              aes(col=pseudo), se=F, linewidth = 0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"),
                     labels = c("Background", "Buffer", "CRW")) + 
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Median Shape Extrapolation Value") +
  scale_x_log10(breaks = c(500, 2000, 10000), labels = scales::comma) +
  guides(colour = guide_legend(title = "Pseudo-Absence")) 
p3

###Spatial###
rm(list=setdiff(ls(), c("p1", "p2", "p3")))

#read in Boyce scores and mixed-effects model
#boyce <- readRDS("output/spatial/boyce_final.RDS")
boyce <- readRDS("output/thinned_spat2/boyce_final.rds")
#if score is NA, revalue to -1
boyce <- boyce %>% mutate(score = ifelse(is.na(score), -1, score))

#redefine levels
levels(boyce$pseudo) <- c("Background", "Buffer", "CRW")
levels(boyce$algorithm) <- c("BRT", "GAMM", "GLMM", "RF")

#create algopseudo
boyce <- boyce %>%
  mutate(algopseudo = paste(algorithm, pseudo, sep = " "))

#plot - sample size
p4 <- ggplot(boyce, aes(x=n, y=score)) + theme_classic() +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm",
              aes(col=pseudo, linetype=algorithm), se=F, linewidth=0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"), 
                     labels = c("Background", "Buffer", "CRW")) + 
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Continuous Boyce Index Score") +
  scale_x_log10(breaks = c(500, 2000, 10000), labels = scales::comma) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash"),
                        labels=c("BRT", "GAMM", "GLMM", "RF")) +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
  theme(legend.key.width = unit(3, "line"))
p4

# plot - extrapolation
p5 <- ggplot(boyce, aes(x=Median, y=score)) + 
  geom_point(colour="lightgrey", size = 0.1) + 
  theme_classic() +
  geom_smooth(method="lm", aes(col=pseudo, linetype=algorithm), se=F, linewidth = 0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"), labels=c("Background", "Buffer", "CRW")) + 
  xlab("Median Shape Extrapolation Value") + 
  scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash"),
                        labels=c("BRT", "GAMM", "GLMM", "RF")) +
  scale_x_log10(breaks = c(10, 20, 50, 100, 200, 500)) + ylab("Continuous Boyce Index Score") +
  guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
         colour = guide_legend(title = "Pseudo-Absence")) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
p5

#plot - sample size vs extrapolation
p6 <- ggplot(boyce, aes(x=n, y=Median)) + 
  geom_point(colour="lightgrey", size = 0.1) + theme_classic() +
  geom_smooth(method="lm", 
              aes(col=pseudo), se=F, linewidth = 0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"), 
                     labels = c("Background", "Buffer", "CRW")) + 
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Median Shape Extrapolation Value") +
  scale_x_log10(breaks = c(500, 2000, 10000), labels = scales::comma) +
  guides(colour = guide_legend(title = "Pseudo-Absence"))+
  theme(legend.key.width = unit(3, "line")) 
p6

###PLOT SAMPLE SIZE AND EXTRAPOLATION VS CBI###

# get the legend of a plot
legend <- get_legend(p1)

# remove the legend from all plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")

# remove ylabs from shape plots
p2 <- p2 + ylab("")
p5 <- p5 + ylab("")

# temporal plots with title
p_temp <- plot_grid(p1 + ggtitle("Temporal Transfer") + theme(plot.title = element_text(face = "bold", size = 16)),
                    p2, align = "hv", labels = c("a", "b"), label_y = 0.95)
p_temp


# spatial plots with title
p_spat <- plot_grid(p4 + ggtitle("Spatial Transfer")  + theme(plot.title = element_text(face = "bold", size = 16)), 
                    p5, align = "hv", labels = c("c", "d"), label_y = 0.95)
p_spat

# both together 
p_both <- plot_grid(p_temp, p_spat, ncol = 1, scale = 0.95)
p_both

# add legend
p_fin <- plot_grid(p_both, legend, rel_widths = c(0.8, 0.2))
p_fin


###PLOT SAMPLE SIZE VS EXTRAPOLATION - SUPP FIGURE###

# remove the legend from all plots
p3 <- p3 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")

# remove ylab from spatial plot
p6 <- p6 + ylab("")

# add titles
p3 <- p3 + ggtitle("Temporal Transfer")
p6 <- p6 + ggtitle("Spatial Transfer")

# plot together
p_supp <- plot_grid(p3, p6, legend, nrow = 1, rel_widths = c(0.4, 0.4, 0.2))

# export pngs
ggsave("text/figures/ecography/05_sample_size_extrapolation.png", p_fin,
       width=13, height=9, dpi = 300)
ggsave("text/figures/ecography/s02_sample_size_vs_extrapolation.png", p_supp,
       width = 13, height = 4.5, dpi = 300)

# export eps
# ggsave("text/figures/ecography/05_sample_size_extrapolation.eps", p_fin,
#        width=13, height=9, dpi = 300)

