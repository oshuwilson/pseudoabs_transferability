#rankings

###TEMPORAL###
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(tidyverse)
  library(gt)
  library(cowplot)
}

#read in Boyce scores and mixed-effects model
#boyce <- readRDS("output/leave-year-out/boyce_filtered.RDS")
boyce <- readRDS("output/thinned_temp2/boyce_final.RDS")

# create algopseudo
boyce <- boyce %>%
  mutate(algopseudo = as.factor(paste(algorithm, pseudo, sep = " ")))
levels(boyce$algopseudo)
levels(boyce$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW",
                              "GAMM Background", "GAMM Buffer", "GAMM CRW",
                              "GLMM Background", "GLMM Buffer", "GLMM CRW",
                              "RF Background", "RF Buffer", "RF CRW")

#average rank (12 Highest 1 Lowest)
avg_rank <- boyce %>% arrange(score) %>% group_by(species, site, stage, season) %>%
  mutate(rank=rank(score, na.last = F)) %>% group_by(algopseudo) %>% summarise(average_rank = mean(rank)) %>%
  arrange(average_rank)

#individual ranks for each row
ranked <- boyce %>% arrange(score) %>% group_by(species, site, stage, season) %>%
  mutate(rank=round(rank(score, na.last = F)))
ggplot(ranked, aes(x=rank, fill=algopseudo)) + geom_bar(position="dodge", stat="count") +
  scale_fill_viridis_d() 

# % of tests where model is top ranking
ranked$top <- if_else(ranked$rank == 12, 1, 0)

# number of case studies
n_case <- nrow(boyce)/12

# summarise top proportions
top <- ranked %>% group_by(algopseudo) %>% summarise(prop = sum(top)/n_case) 


# % where models are within 0.05 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, season) %>% 
  nest()

within_5 <- function(data) {
  data <- data %>% mutate(lim = max(score, na.rm=T)-0.05)
  data <- data %>% mutate(within = if_else(score > lim, 1, 0))
  print(data)
}

test <- comp[1,] %>% unnest(cols = c(data))
test <- within_5(data = test)

all <- NULL

for(i in 1:n_case){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_5(data = test)
  all <- rbind(all, test)
}

within_0.05 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/n_case) 


# % where models are within 0.1 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, season) %>% 
  nest()

within_10 <- function(data) {
  data <- data %>% mutate(lim = max(score, na.rm=T)-0.1)
  data <- data %>% mutate(within = if_else(score > lim, 1, 0))
  print(data)
}

test <- comp[1,] %>% unnest(cols = c(data))
test <- within_10(data = test)

all <- NULL

for(i in 1:n_case){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_10(data = test)
  all <- rbind(all, test)
}

within_0.1 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/n_case) 

#check examples where RF background isn't within 0.1
rfback <- all %>% filter(algopseudo == "RF Background", within == 0)
topscores <- all %>% filter(top == 1) %>% 
  select(species, stage, site, season, algopseudo, score) %>%
  rename(topscore = score, topconf = algopseudo)
rfback <- rfback %>% left_join(topscores)

brtback <- all %>% filter(algopseudo == "BRT Background") %>% 
  select(species, site, stage, season, score) %>%
  rename(brtback = score)
rfback <- rfback %>% left_join(brtback)
fail <- rfback %>% filter(score < 0.5 & brtback < 0.5)

#merge table
top <- top %>% rename(top = prop)
within_0.05 <- within_0.05 %>% rename(top5 = prop)
within_0.1 <- within_0.1 %>% rename(top10 = prop)
props <- top %>% left_join(within_0.05) %>% left_join(within_0.1)
props$top <- format(round(props$top, 3), nsmall = 3)
props$top5 <- format(round(props$top5, 3), nsmall = 3)
props$top10 <- format(round(props$top10, 3), nsmall = 3)
props <- props %>% mutate_if(is.character, as.numeric)
props <- props %>% mutate(top = top * 100, top5 = top5 * 100, top10 = top10 * 100)
levels(props$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW", 
                              "GAMM Background", "GAMM Buffer", "GAMM CRW",
                              "GLMM Background", "GLMM Buffer", "GLMM CRW",
                              "RF Background", "RF Buffer", "RF CRW")

#create gt table for paper
props <- props %>% rownames_to_column() %>% mutate_all(as.character) %>% 
  pivot_longer(-rowname) %>% pivot_wider(names_from=rowname)
colnames(props) <- props[1,]
props <- props[-1,]
props$algopseudo <- as.factor(props$algopseudo)
props <- props %>% mutate_if(is.character, as.numeric)
levels(props$algopseudo) <- c("Top Model (%)", "Within 0.1 of Top Model (%)", "Within 0.05 of Top Model (%)")
props <- props %>% select(algopseudo,
                 `BRT Background`, `RF Background`, `GAMM Background`, `GLMM Background`,
                 `BRT CRW`, `RF CRW`, `GAMM CRW`, `GLMM CRW`,
                 `BRT Buffer`, `RF Buffer`, `GAMM Buffer`, `GLMM Buffer`)
tab <- props %>% gt() %>% cols_align(align = "center") %>%
  data_color(direction = "row", palette = "Purples", method = "numeric") 
tab <- tab %>%
  tab_spanner(label = "Background", columns = 2:5) %>%
  tab_spanner(label = "CRW", columns = 6:9) %>%
  tab_spanner(label = "Buffer", columns = 10:13) %>%
  tab_spanner(label = "Pseudo-Absence", columns = 1) %>%
  cols_label(algopseudo = "Algorithm",
    `RF Background` = "RF", `GAMM Background` = "GAMM", `GLMM Background` = "GLMM", `BRT Background` = "BRT",
    `RF CRW` = "RF", `GAMM CRW` = "GAMM", `GLMM CRW` = "GLMM", `BRT CRW` = "BRT",
    `RF Buffer` = "RF", `GAMM Buffer` = "GAMM", `GLMM Buffer` = "GLMM", `BRT Buffer` = "BRT") %>%
  cols_width(algopseudo ~ px(140), everything() ~ px(70)) 
tab

#gtsave(tab, filename = "text/figures/final figs/t02_temp_table_thinned_spatiotemp.png")


#top ranked model scores
best <- filter(ranked, top == 1) 
min(best$score, na.rm=T)
levels(best$algopseudo)

# relevel to put crw before buffer
best$algopseudo <- fct_relevel(best$algopseudo, 
             "BRT Background", "RF Background", "GAMM Background", "GLMM Background",
             "BRT CRW", "RF CRW", "GAMM CRW", "GLMM CRW",
             "BRT Buffer", "RF Buffer", "GAMM Buffer", "GLMM Buffer")

p1 <- ggplot(best, aes(x=score)) + 
  geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) + 
  theme_classic() +
  scale_fill_manual(values = c("red1", "red2", "red3", "red4",
                                "steelblue1", "steelblue2", "steelblue3", "steelblue4",
                               "darkorange1", "darkorange4")) +
  guides(fill = guide_legend(ncol=3, title = "Algorithm and Pseudo-Absence Technique")) +
  theme(legend.position = c(0.3, 0.9)) + 
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) + 
  ylab("Count") + xlab("Continuous Boyce Index Score")
p1

###SPATIAL###
rm(list=setdiff(ls(), "p1"))
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")
library(tidyverse)
library(gt)

#read in Boyce scores and mixed-effects model
#boyce <- readRDS("output/spatial/boyce_final.RDS")
boyce <- readRDS("output/thinned_spat2/boyce_final.rds")
#if score is NA, revalue to -1
boyce <- boyce %>% mutate(score = ifelse(is.na(score), -1, score))

# rename site and test_site
boyce <- boyce %>%
  rename(train_site = site) %>%
  rename(site = test_site)

#average rank (12 Highest 1 Lowest)
avg_rank <- boyce %>% arrange(score) %>% group_by(species, site, stage, train_site) %>%
  mutate(rank=rank(score, na.last = F)) %>% group_by(algopseudo) %>% summarise(average_rank = mean(rank)) %>%
  arrange(average_rank)

#individual ranks for each row
ranked <- boyce %>% arrange(score) %>% group_by(species, site, stage, train_site) %>%
  mutate(rank=round(rank(score, na.last = F))) 
ggplot(ranked, aes(x=rank, fill=algopseudo)) + geom_bar(position="dodge", stat="count") +
  scale_fill_viridis_d() 

#number of case studies
n_case <- nrow(boyce)/12

# % of tests where model is top ranking
ranked$top <- if_else(ranked$rank == 12, 1, 0)
top <- ranked %>% group_by(algopseudo) %>% summarise(prop = sum(top)/n_case) 


# % where models are within 0.05 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, train_site) %>% 
  nest()

within_5 <- function(data) {
  data <- data %>% mutate(lim = max(score)-0.05)
  data <- data %>% mutate(within = if_else(score > lim, 1, 0))
  print(data)
}

test <- comp[1,] %>% unnest(cols = c(data))
test <- within_5(data = test)

all <- NULL

for(i in 1:n_case){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_5(data = test)
  all <- rbind(all, test)
}

within_0.05 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/n_case) 


# % where models are within 0.1 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, train_site) %>% 
  nest()

within_10 <- function(data) {
  data <- data %>% mutate(lim = max(score)-0.1)
  data <- data %>% mutate(within = if_else(score > lim, 1, 0))
  print(data)
}

test <- comp[1,] %>% unnest(cols = c(data))
test <- within_10(data = test)

all <- NULL

for(i in 1:n_case){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_10(data = test)
  all <- rbind(all, test)
}

within_0.1 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/n_case) 


#merge table
top <- top %>% rename(top = prop)
within_0.05 <- within_0.05 %>% rename(top5 = prop)
within_0.1 <- within_0.1 %>% rename(top10 = prop)
props <- top %>% left_join(within_0.05) %>% left_join(within_0.1)
props$top <- format(round(props$top, 3), nsmall = 3)
props$top5 <- format(round(props$top5, 3), nsmall = 3)
props$top10 <- format(round(props$top10, 3), nsmall = 3)
props <- props %>% mutate(algopseudo = as.factor(algopseudo))
props <- props %>% mutate_if(is.character, as.numeric)
props <- props %>% mutate(top = top * 100, top5 = top5 * 100, top10 = top10 * 100)
levels(props$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW", 
                              "GAMM Background", "GAMM Buffer", "GAMM CRW",
                              "GLMM Background", "GLMM Buffer", "GLMM CRW",
                              "RF Background", "RF Buffer", "RF CRW")

#create gt table for paper
props <- props %>% rownames_to_column() %>% mutate_all(as.character) %>% 
  pivot_longer(-rowname) %>% pivot_wider(names_from=rowname)
colnames(props) <- props[1,]
props <- props[-1,]
props$algopseudo <- as.factor(props$algopseudo)
props <- props %>% mutate_if(is.character, as.numeric)
levels(props$algopseudo) <- c("Top Model (%)", "Within 0.1 of Top Model (%)", "Within 0.05 of Top Model (%)")
props <- props %>% select(algopseudo,
                          `BRT Background`, `RF Background`, `GAMM Background`, `GLMM Background`,
                          `BRT CRW`, `RF CRW`, `GAMM CRW`, `GLMM CRW`,
                          `BRT Buffer`, `RF Buffer`, `GAMM Buffer`, `GLMM Buffer`)
tab <- props %>% gt() %>% cols_align(align = "center") %>%
  data_color(direction = "row", palette = "Purples", method = "numeric") 
tab <- tab %>%
  tab_spanner(label = "Background", columns = 2:5) %>%
  tab_spanner(label = "CRW", columns = 6:9) %>%
  tab_spanner(label = "Buffer", columns = 10:13) %>%
  tab_spanner(label = "Pseudo-Absence", columns = 1) %>%
  cols_label(algopseudo = "Algorithm",
             `RF Background` = "RF", `GAMM Background` = "GAMM", `GLMM Background` = "GLMM", `BRT Background` = "BRT",
             `RF CRW` = "RF", `GAMM CRW` = "GAMM", `GLMM CRW` = "GLMM", `BRT CRW` = "BRT",
             `RF Buffer` = "RF", `GAMM Buffer` = "GAMM", `GLMM Buffer` = "GLMM", `BRT Buffer` = "BRT") %>%
  cols_width(algopseudo ~ px(140), everything() ~ px(70)) 
tab

#gtsave(tab, filename = "text/figures/final figs/t03_spat_table_spatiotemp.png")


#top ranked model scores
best <- filter(ranked, top == 1) 
min(best$score, na.rm=T)
best$algopseudo <- as.factor(best$algopseudo)
levels(best$algopseudo)
levels(best$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW", 
                              "GAMM Background",
                              "GLMM Background", "GLMM Buffer", "GLMM CRW",
                              "RF Background", "RF Buffer", "RF CRW")
best$algopseudo <- fct_relevel(best$algopseudo, 
            "BRT Background", "RF Background", "GAMM Background", "GLMM Background",
            "BRT CRW", "RF CRW", "GLMM CRW",
            "BRT Buffer", "RF Buffer", "GLMM Buffer")


p2 <- ggplot(best, aes(x=score)) + 
  geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) + 
  theme_classic() +
  scale_fill_manual(values = c("red1", "red2", "red3", "red4",
                               "steelblue1", "steelblue2", "steelblue4",
                               "darkorange1", "darkorange2", "darkorange4")) +
  guides(fill = guide_legend(ncol=3, title = "Algorithm and Pseudo-Absence Technique")) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) + 
  ylab("Count") + xlab("Continuous Boyce Index Score")
p2

# make version with full legend
mod <- data.frame(
  site = "test",
  algorithm = "test",
  train_site = "test",
  species = "test",
  stage = "test",
  algopseudo = as.factor(c("GLMM Background", "GAMM Background", "RF Background", "BRT Background",
                          "GLMM CRW", "GAMM CRW", "RF CRW", "BRT CRW",
                          "GLMM Buffer", "GAMM Buffer", "RF Buffer", "BRT Buffer")),
  n = 1,
  pseudo = "test",
  score = 1,
  shape = 1,
  rank = 12,
  top = 1)

# make algopseudo a factor
mod$algopseudo <- as.factor(mod$algopseudo)
mod$algopseudo <- fct_relevel(mod$algopseudo, 
                               "BRT Background", "RF Background", "GAMM Background", "GLMM Background",
                               "BRT CRW", "RF CRW", "GAMM CRW", "GLMM CRW",
                               "BRT Buffer", "RF Buffer", "GAMM Buffer", "GLMM Buffer")
unique(mod$algopseudo)

p3 <- ggplot(mod, aes(x=score)) + 
  geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) + 
  theme_classic() +
  scale_fill_manual(values = c("red1", "red2", "red3", "red4",
                               "steelblue1", "steelblue2", "steelblue3", "steelblue4",
                               "darkorange1", "darkorange2", "darkorange3", "darkorange4")) +
  guides(fill = guide_legend(ncol=3, title = "Algorithm and Pseudo-Absence Technique")) +
  scale_x_continuous(limits = c(0,1.01), breaks = seq(0, 1, 0.2)) + 
  ylab("Count") + xlab("Continuous Boyce Index Score")
p3


#plot both together
rm(list=setdiff(ls(), c("p1", "p2", "p3")))

#remove legends
p1 <- p1 + guides(fill="none")
p2 <- p2 + guides(fill="none")

#add titles
p1 <- p1 + ggtitle("Temporal Transfer")
p2 <- p2 + ggtitle("Spatial Transfer")

#customise font sizes
p1 <- p1 + theme_classic(base_size = 7)
p2 <- p2 + theme_classic(base_size = 7)

#get legend
legend <- get_legend(p3)
legend

#plot without legend
plots <- plot_grid(p1, p2 + theme(legend.position="none", ncol=2, align="v"))
plots

#plots and legend together
p3 <- plot_grid(legend, plots, ncol=1, rel_heights = c(1.5,6))
p3

#export
# ggsave("text/figures/final figs/s04_spat_temp_transfer_spatiotemp.png", p3, 
#        width = 220, height = 160, units = "mm", dpi = 300)
