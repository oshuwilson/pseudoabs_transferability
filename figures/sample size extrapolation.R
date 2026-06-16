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

# #plot - sample size
# p1 <- ggplot(boyce, aes(x=n, y=transferability, group = algopseudo)) + theme_classic() +
#   geom_point(colour="lightgrey", size = 0.1) +
#   geom_smooth(method="lm", 
#               aes(col=pseudo, linetype=algorithm), se=F, linewidth=0.6) +
#   scale_color_manual(values=c("red3", "steelblue4", "orange2"),
#                      labels = c("Background", "Buffer", "CRW")) + 
#   xlab("Number of Thinned Locations in Training Data") + 
#   ylab("Continuous Boyce Index Score") +
#   scale_x_log10(breaks = c(100, 500, 2000, 5000), labels = scales::comma) +
#   scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
#   guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
#          colour = guide_legend(title = "Pseudo-Absence")) +
#   scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
#   theme(legend.key.width = unit(1, "cm"))
# p1 

# reorder algorithm to group machine learning models, regression models, and spatial correlation models
boyce$algorithm <- factor(boyce$algorithm, levels = c("GLMM", "GAMM", "MaxEnt", "BART", "BRT", "RF", "mSTPP", "INLA-SPDE"))

p1 <- ggplot(boyce, aes(x=n, y=transferability, group = algorithm)) +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm",
              aes(col=algorithm), se=F, linewidth=0.6) + 
  scale_color_manual(values = c("red1", "red3", "red4", "steelblue1", "steelblue3", "steelblue4", "orange1", "orange3")) +
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Continuous Boyce Index Score") +
  scale_x_log10(breaks = c(100, 500, 2000, 5000), labels = scales::comma, expand = c(0.01,0.01)) +
  facet_wrap(~pseudo) +
  theme_minimal() +
  geom_hline(yintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
  theme(legend.key.width = unit(1, "cm"),
        panel.border = element_rect(color = "grey40", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Algorithm")) 
p1 + ggview::canvas(15,7) 

# # plot - extrapolation
# p2 <- ggplot(boyce, aes(x=Median, y=transferability)) + 
#   geom_point(colour="lightgrey", size = 0.1) + 
#   theme_classic() +
#   geom_smooth(method="lm", aes(col=pseudo, linetype=algorithm), se=F, linewidth = 0.6) +
#   scale_color_manual(values=c("red3", "steelblue4", "orange2"), labels=c("Background", "Buffer", "CRW")) + 
#   xlab("Median Shape Extrapolation Value") + 
#   scale_x_log10(breaks = c(10, 20, 50, 100, 200)) + ylab("Continuous Boyce Index Score") +
#   scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
#   guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
#          colour = guide_legend(title = "Pseudo-Absence")) +
#   scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
#   theme(legend.key.width = unit(1, "cm"))
# p2

p2 <- ggplot(boyce, aes(x=Median, y=transferability, group = algorithm)) +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm",
              aes(col=algorithm), se=F, linewidth=0.6) + 
  scale_color_manual(values = c("red1", "red3", "red4", "steelblue1", "steelblue3", "steelblue4", "orange1", "orange3")) +
  xlab("Median Shape Extrapolation Value") + 
  scale_x_log10(breaks = c(10, 20, 50, 100, 200), expand = c(0.01, 0.01)) + ylab("Continuous Boyce Index Score") +
  ylab("Continuous Boyce Index Score") +
  facet_wrap(~pseudo) +
  theme_minimal() +
  geom_hline(yintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
  theme(legend.key.width = unit(1, "cm"),
        panel.border = element_rect(color = "grey40", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Algorithm"))
p2 + ggview::canvas(15,7) 

#plot - sample size vs extrapolation
p3 <- ggplot(boyce, aes(x=n, y=Median)) + 
  geom_point(colour="lightgrey", size = 0.1) + theme_classic() +
  geom_smooth(method="lm", 
              aes(col=pseudo), se=F, linewidth = 0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"),
                     labels = c("Background", "Buffer", "CRW")) + 
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Median Shape Extrapolation Value") +
  scale_x_log10(breaks = c(100, 500, 2000, 5000), labels = scales::comma) +
  guides(colour = guide_legend(title = "Pseudo-Absence")) 
p3

###Spatial###
rm(list=setdiff(ls(), c("p1", "p2", "p3")))

#read in Boyce scores 
boyce <- readRDS("output/spatial/transferability_scores.RDS")

# reorder algorithm to group machine learning models, regression models, and spatial correlation models
boyce$algorithm <- factor(boyce$algorithm, levels = c("GLMM", "GAMM", "MaxEnt", "BART", "BRT", "RF", "mSTPP", "INLA-SPDE"))


#plot - sample size
# p4 <- ggplot(boyce, aes(x=n, y=transferability, group = algopseudo)) + theme_classic() +
#   geom_point(colour="lightgrey", size = 0.1) +
#   geom_smooth(method="lm", 
#               aes(col=pseudo, linetype=algorithm), se=F, linewidth=0.6) +
#   scale_color_manual(values=c("red3", "steelblue4", "orange2"),
#                      labels = c("Background", "Buffer", "CRW")) + 
#   xlab("Number of Thinned Locations in Training Data") + 
#   ylab("Continuous Boyce Index Score") +
#   scale_x_log10(breaks = c(200, 500, 2000, 5000), labels = scales::comma) +
#   scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
#   guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
#          colour = guide_legend(title = "Pseudo-Absence")) +
#   scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
# p4

p4 <- ggplot(boyce, aes(x=n, y=transferability, group = algorithm)) +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm",
              aes(col=algorithm), se=F, linewidth=0.6) + 
  scale_color_manual(values = c("red1", "red3", "red4", "steelblue1", "steelblue3", "steelblue4", "orange1", "orange3")) +
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Continuous Boyce Index Score") +
  scale_x_log10(breaks = c(100, 500, 2000, 5000), labels = scales::comma, expand = c(0.01,0.01)) +
  facet_wrap(~pseudo) +
  theme_minimal() +
  geom_hline(yintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
  theme(legend.key.width = unit(1, "cm"),
        panel.border = element_rect(color = "grey40", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Algorithm")) 
p4 + ggview::canvas(15,7)

# plot - extrapolation
# p5 <- ggplot(boyce, aes(x=Median, y=transferability)) + 
#   geom_point(colour="lightgrey", size = 0.1) + 
#   theme_classic() +
#   geom_smooth(method="lm", aes(col=pseudo, linetype=algorithm), se=F, linewidth = 0.6) +
#   scale_color_manual(values=c("red3", "steelblue4", "orange2"), labels=c("Background", "Buffer", "CRW")) + 
#   xlab("Median Shape Extrapolation Value") + 
#   scale_x_log10(breaks = c(10, 20, 50, 100, 200)) + ylab("Continuous Boyce Index Score") +
#   scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
#   guides(linetype = guide_legend(override.aes = list(colour = "black"), title = "Algorithm"),
#          colour = guide_legend(title = "Pseudo-Absence")) +
#   scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01))
# p5

p5 <- ggplot(boyce, aes(x=Median, y=transferability, group = algorithm)) +
  geom_point(colour="lightgrey", size = 0.1) +
  geom_smooth(method="lm",
              aes(col=algorithm), se=F, linewidth=0.6) + 
  scale_color_manual(values = c("red1", "red3", "red4", "steelblue1", "steelblue3", "steelblue4", "orange1", "orange3")) +
  xlab("Median Shape Extrapolation Value") + 
  scale_x_log10(breaks = c(10, 20, 50, 100, 200), expand = c(0.01,0.01)) + ylab("Continuous Boyce Index Score") +
  ylab("Continuous Boyce Index Score") +
  facet_wrap(~pseudo) +
  theme_minimal() +
  geom_hline(yintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  scale_y_continuous(limits=c(-1, 1), expand=c(0.01,0.01)) +
  theme(legend.key.width = unit(1, "cm"),
        panel.border = element_rect(color = "grey40", linewidth = 0.5),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Algorithm"))
p5 + ggview::canvas(15,7)

#plot - sample size vs extrapolation
p6 <- ggplot(boyce, aes(x=n, y=Median)) + 
  geom_point(colour="lightgrey", size = 0.1) + theme_classic() +
  geom_smooth(method="lm", 
              aes(col=pseudo), se=F, linewidth = 0.6) +
  scale_color_manual(values=c("red3", "steelblue4", "orange2"),
                     labels = c("Background", "Buffer", "CRW")) + 
  xlab("Number of Thinned Locations in Training Data") + 
  ylab("Median Shape Extrapolation Value") +
  scale_x_log10(breaks = c(200, 500, 2000, 5000), labels = scales::comma) +
  guides(colour = guide_legend(title = "Pseudo-Absence")) 
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
                    p2, align = "v", labels = c("a", "b"), label_y = 0.95, ncol = 1)
p_temp + ggview::canvas(13,11)


# spatial plots with title
p_spat <- plot_grid(p4 + ggtitle("Spatial Transfer") + theme(plot.title = element_text(face = "bold", size = 16)),
                    p5, align = "v", labels = c("c", "d"), label_y = 0.95, ncol = 1)
p_spat + ggview::canvas(13,11)

# both together 
p_both <- plot_grid(p_temp, p_spat, ncol = 1, scale = 0.95)
p_both + ggview::canvas(13, 22)

# add legend
p_fin <- plot_grid(p_both, legend, rel_widths = c(0.9, 0.1))
p_fin + ggview::canvas(13,20)


###PLOT SAMPLE SIZE VS EXTRAPOLATION - SUPP FIGURE###

# get legend
legend2 <- get_legend(p3)

# remove the legend from all plots
p3 <- p3 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")

# remove ylab from spatial plot
p6 <- p6 + ylab("")

# add titles
p3 <- p3 + ggtitle("Temporal Transfer")
p6 <- p6 + ggtitle("Spatial Transfer")

# plot together
p_supp <- plot_grid(p3, p6, legend2, nrow = 1, rel_widths = c(0.4, 0.4, 0.2))
p_supp

# export pngs
ggsave("text/figures/ecography/05_sample_size_extrapolation.png", p_fin,
       width=13, height=20, dpi = 300)
ggsave("text/figures/ecography/s02_sample_size_vs_extrapolation.png", p_supp,
       width = 13, height = 4.5, dpi = 300)

# export eps
ggsave("text/figures/ecography/05_sample_size_extrapolation.eps", p_fin,
       width=13, height=20, dpi = 300)



#-------------------------------------------------------------------------------
# Get R-Squared values of lines in graphs
#-------------------------------------------------------------------------------

rm(list=ls())

###Temporal###
#read in Boyce scores 
boyce <- readRDS("output/temporal/transferability_scores.RDS")

# for each algopseudo combo, fit a linear model and get R2 for sample size vs boyce
r2_results <- boyce %>%
  group_by(algopseudo) %>%
  do(model = lm(transferability ~ log10(n), data = .)) %>%
  summarise(algopseudo = unique(algopseudo),
            r_squared = summary(model)$r.squared)
print(r2_results)

# for each algopseudo combo, fit a linear model and get R2 for extrapolation vs boyce
r2_results2 <- boyce %>%
  group_by(algopseudo) %>%
  do(model = lm(transferability ~ log10(Median), data = .)) %>%
  summarise(algopseudo = unique(algopseudo),
            r_squared = summary(model)$r.squared)
print(r2_results2)

###Spatial###
rm(list=ls())

#read in Boyce scores
boyce <- readRDS("output/spatial/transferability_scores.RDS")

# for each algopseudo combo, fit a linear model and get R2 for sample size vs boyce
r2_results3 <- boyce %>%
  group_by(algopseudo) %>%
  do(model = lm(transferability ~ log10(n), data = .)) %>%
  summarise(algopseudo = unique(algopseudo),
            r_squared = summary(model)$r.squared)
print(r2_results3)

# for each algopseudo combo, fit a linear model and get R2 for
r2_results4 <- boyce %>%
  mutate(Median = ifelse(Median == 0, 0.000001, Median)) %>%
  group_by(algopseudo) %>%
  do(model = lm(transferability ~ log10(Median), data = .)) %>%
  summarise(algopseudo = unique(algopseudo),
            r_squared = summary(model)$r.squared)
print(r2_results4)
