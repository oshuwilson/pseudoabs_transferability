#-------------------------------------------------------------------------------
# Plot transferability as a function of different species traits
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(tidyverse)
library(ggridges)
library(cowplot)

# read in data
temp <- readRDS("output/temporal/transferability_scores.RDS")
spat <- readRDS("output/spatial/transferability_scores.RDS")


#-------------------------------------------------------------------------------
# 1. Pagophilic vs non-pagophilic species
#-------------------------------------------------------------------------------

# define pagophilic species
pago <- c("EMPE", "CRAS", "WESE", "ADPE", "ANPE", "HUWH")
temp <- temp %>%
  mutate(pagophilic = ifelse(species %in% pago, "Ice-Associated", "Ice-Free"))
spat <- spat %>%
  mutate(pagophilic = ifelse(species %in% pago, "Ice-Associated", "Ice-Free"))

p1 <- ggplot(temp, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Temporal Transfer") +
  facet_wrap(~ pagophilic) +
  theme(plot.title = element_text(face = "bold", size = 16))

p2 <- ggplot(spat, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Spatial Transfer") +
  facet_wrap(~ pagophilic) +
  theme(plot.title = element_text(face = "bold", size = 16))

# plot together
pago_grid <- plot_grid(p1, p2) 
pago_grid + ggview::canvas(12, 10)

# export
ggsave("text/Ecological Modelling Submission/figs/supp/traits_pago.png",
       width = 12, height = 10, units = "in", dpi = 300, pago_grid)


#-------------------------------------------------------------------------------
# Wide-ranging vs narrow-ranging species/stages
#-------------------------------------------------------------------------------

# define cases as wide-ranging or narrow-ranging
temp <- temp %>%
  mutate(range = case_when(
    species == "ADPE" ~ "narrow-ranging",
    species == "ANFS" & stage == "breeding" ~ "narrow-ranging",
    species == "ANFS" & stage == "post-moult" ~ "wide-ranging",
    species == "ANPE" ~ "narrow-ranging",
    species == "CRAS" ~ "wide-ranging",
    species == "EMPE" ~ "wide-ranging",
    species == "GHAL" ~ "narrow-ranging",
    species == "HUWH" ~ "wide-ranging",
    species == "MAPE" & stage != "post-breeding" ~ "narrow-ranging",
    species == "MAPE" & stage == "post-breeding" ~ "wide-ranging",
    species == "SOES" & stage == "post-breeding" ~ "narrow-ranging",
    species == "SOES" & stage == "post-moult" ~ "wide-ranging",
    species == "SUFS" ~ "narrow-ranging"
  )) %>%
  mutate(range = ifelse(range == "wide-ranging", "Unrestricted", "Central-Place Forager"))

spat <- spat %>%
  mutate(range = case_when(
    species == "ADPE" ~ "narrow-ranging",
    species == "ANFS" & stage == "breeding" ~ "narrow-ranging",
    species == "ANFS" & stage == "post-moult" ~ "wide-ranging",
    species == "ANPE" ~ "narrow-ranging",
    species == "CRAS" ~ "wide-ranging",
    species == "EMPE" ~ "wide-ranging",
    species == "GHAL" ~ "narrow-ranging",
    species == "HUWH" ~ "wide-ranging",
    species == "MAPE" & stage != "post-breeding" ~ "narrow-ranging",
    species == "MAPE" & stage == "post-breeding" ~ "wide-ranging",
    species == "SOES" & stage == "post-breeding" ~ "narrow-ranging",
    species == "SOES" & stage == "post-moult" ~ "wide-ranging",
    species == "SUFS" ~ "narrow-ranging"
  )) %>%
  mutate(range = ifelse(range == "wide-ranging", "Unrestricted", "Central-Place Forager"))


p3 <- ggplot(temp, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Temporal Transfer") +
  facet_wrap(~ range) +
  theme(plot.title = element_text(face = "bold", size = 16))

p4 <- ggplot(spat, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Spatial Transfer") +
  facet_wrap(~ range) +
  theme(plot.title = element_text(face = "bold", size = 16))

# plot together
range_grid <- plot_grid(p3, p4)
range_grid + ggview::canvas(12, 10)

# export
ggsave("text/Ecological Modelling Submission/figs/supp/traits_range.png",
       width = 12, height = 10, units = "in", dpi = 300, range_grid)


#-------------------------------------------------------------------------------
# Seabirds vs Marine Mammals
#-------------------------------------------------------------------------------

unique(temp$species)

temp <- temp %>%
  mutate(taxa = case_when(
    species %in% c("ADPE", "MAPE", "EMPE", "ANPE", "GHAL") ~ "Seabirds",
    species %in% c("ANFS", "SUFS", "SOES", "CRAS", "HUWH") ~ "Marine Mammals"
  ))

spat <- spat %>%
  mutate(taxa = case_when(
    species %in% c("ADPE", "MAPE", "EMPE", "ANPE", "GHAL") ~ "Seabirds",
    species %in% c("ANFS", "SUFS", "SOES", "CRAS", "HUWH") ~ "Marine Mammals"
  ))

p5 <- ggplot(temp, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Temporal Transfer") +
  facet_wrap(~ taxa) +
  theme(plot.title = element_text(face = "bold", size = 16))

p6 <- ggplot(spat, (aes(x=transferability, y=fct_reorder(algopseudo, transferability, .fun=mean), fill = fct_reorder(algopseudo, transferability, .fun=mean)))) + 
  geom_density_ridges2(scale=2, linewidth = 0.1) + 
  theme_minimal() +
  xlab("Continuous Boyce Index Score") + 
  ylab("Algorithm and Pseudo-Absence Technique") +
  scale_fill_viridis_d(guide = "none", option = "mako", begin = 0.2) +
  scale_x_continuous(expand=c(0,0), limits=c(-1,1)) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed", linewidth = 0.2) +
  ggtitle("Spatial Transfer") +
  facet_wrap(~ taxa) +
  theme(plot.title = element_text(face = "bold", size = 16))

# plot together
taxa_grid <- plot_grid(p5, p6)
taxa_grid + ggview::canvas(12, 10)

# export
ggsave("text/Ecological Modelling Submission/figs/supp/traits_taxa.png",
       width = 12, height = 10, units = "in", dpi = 300, taxa_grid)


#-------------------------------------------------------------------------------
# Linear mixed-effects model
#-------------------------------------------------------------------------------

library(lme4)

# temporal
temp_lmm <- lmer(transferability ~ pagophilic + (1|species/stage) + (1|algopseudo), data = temp)
summary(temp_lmm)
