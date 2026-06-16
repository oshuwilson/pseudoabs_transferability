#-------------------------------------------------------------------------------
# Plot temporal distance effects
#-------------------------------------------------------------------------------

rm(list = ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(tidyverse)

# read in temporal distance results
rf <- readRDS("output/temporal/distance/rf.rds")
bart <- readRDS("output/temporal/distance/bart.rds")
brt <- readRDS("output/temporal/distance/brt.rds")

# combine together
data <- bind_rows(rf, brt, bart)

# change site name WAP to AP
data <- data %>%
  mutate(site = ifelse(site == "WAP", "AP", site))

# plot
p1 <- ggplot(data, aes(x = gap, y = boyce, col = paste(site, stage),
                 fill = paste(site, stage))) +
  geom_point() +
  geom_smooth(method = "lm", se = F, linewidth = 1) +
  facet_wrap(~algo) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.border = element_rect(colour = "grey40", fill = NA),
        panel.grid.minor.x = element_blank()) +
  xlim(0, NA) +
  scale_color_viridis_d(end = 0.8) +
  labs(x = "Temporal Interval (years)",
       y = "Continuous Boyce Index",
       col = "Site and Stage",
       fill = "Site and Stage") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_continuous(limits = c(0, NA), breaks = seq(0, 14, 2))
p1 + ggview::canvas(12, 6)

# export 
ggsave("text/figures/ecography/sxx_temporal_distance.png", p1,
       width = 12, height = 6, dpi = 300)
