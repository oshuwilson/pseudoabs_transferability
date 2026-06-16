# plot dropoff from self-tests to transferability tests

###TEMPORAL###
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(tidyverse)
  library(cowplot)
}

#read in Boyce scores 
boyce <- readRDS("output/temporal/transferability_scores.RDS")

#create dropoff column
boyce <- boyce %>%
  mutate(dropoff = transferability - self_test)

# calculate mean dropoff per configuration
boyce_mean <- boyce %>%
  group_by(algopseudo) %>%
  summarise(mean_dropoff = mean(dropoff, na.rm = TRUE))
boyce <- boyce %>%
  left_join(boyce_mean, by = "algopseudo")

# plot
p1 <- ggplot(boyce_mean, aes(x = fct_reorder(algopseudo, mean_dropoff, .fun=mean), y = mean_dropoff)) +
  geom_point(data = boyce, aes(y = dropoff), col = "grey80") +
  geom_bar(stat = "identity", alpha = 0.7, aes(fill = mean_dropoff), col = "black") +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 1) +
  coord_flip() +
  scale_fill_gradient2(limits = c(-0.8, 0.8), breaks = c(-0.8, -0.4, 0, 0.4, 0.8), name = "Mean Difference") +
  theme_minimal() +
  ylim(-2, 2) + 
  xlab("") +
  ylab("Performance Dropoff") + 
  ggtitle("Temporal Transfer")
p1

###SPATIAL###

rm(list = setdiff(ls(), "p1"))

# read in Boyce scores
boyce <- readRDS("output/spatial/transferability_scores.RDS")

# create dropoff column
boyce <- boyce %>%
  mutate(dropoff = transferability - self_test)

# calculate mean dropoff per configuration
boyce_mean <- boyce %>%
  group_by(algopseudo) %>%
  summarise(mean_dropoff = mean(dropoff, na.rm = TRUE))
boyce <- boyce %>%
  left_join(boyce_mean, by = "algopseudo")

# plot
p2 <- ggplot(boyce_mean, aes(x = fct_reorder(algopseudo, mean_dropoff, .fun=mean), y = mean_dropoff)) +
  geom_point(data = boyce, aes(y = dropoff), col = "grey80") +
  geom_bar(stat = "identity", alpha = 0.7, aes(fill = mean_dropoff), col = "black") +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 1) +
  coord_flip() +
  scale_fill_gradient2(limits = c(-0.8, 0.8), name = "Mean Difference") +
  theme_minimal() +
  ylim(-2, 2) + 
  xlab("") +
  ylab("Performance Dropoff") + 
  ggtitle("Spatial Transfer") 
p2

# combine plots
legend <- get_legend(p1)
plot(legend)

# remove legends from individual plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

# plot all together
drops <- plot_grid(p1, p2, ncol = 1) 
all <- plot_grid(drops, legend, ncol = 2, rel_widths = c(2, 0.4)) 
all + ggview::canvas(12, 16)

# export
ggsave("text/figures/ecography/s05_dropoff.png", all,
       width = 12, height = 16, units = "in", dpi = 300)
