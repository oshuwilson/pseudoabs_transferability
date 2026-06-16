#rankings

###TEMPORAL###
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(tidyverse)
  library(gt)
  library(cowplot)
}

#read in Boyce scores 
boyce <- readRDS("output/temporal/transferability_scores.RDS")

#average rank (24 Highest 1 Lowest)
avg_rank <- boyce %>% arrange(transferability) %>% group_by(species, site, stage, season) %>%
  mutate(rank=rank(transferability, na.last = F)) %>% group_by(algopseudo) %>% summarise(average_rank = mean(rank)) %>%
  arrange(average_rank)

#individual ranks for each row
ranked <- boyce %>% arrange(transferability) %>% group_by(species, site, stage, season) %>%
  mutate(rank=round(rank(transferability, na.last = F)))
ggplot(ranked, aes(x=rank, fill=algopseudo)) + geom_bar(position="dodge", stat="count") +
  scale_fill_viridis_d() 

# % of tests where model is top ranking
ranked$top <- if_else(ranked$rank == 24, 1, 0)

# number of case studies
n_case <- nrow(boyce)/24

# summarise top proportions
top <- ranked %>% group_by(algopseudo) %>% summarise(prop = sum(top)/n_case) 

# % where models are within 0.05 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, season) %>% 
  nest()

within_5 <- function(data) {
  data <- data %>% mutate(lim = max(transferability, na.rm=T)-0.05)
  data <- data %>% mutate(within = if_else(transferability > lim, 1, 0))
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
  data <- data %>% mutate(lim = max(transferability, na.rm=T)-0.1)
  data <- data %>% mutate(within = if_else(transferability > lim, 1, 0))
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
props <- props %>%
  mutate(top = as.numeric(top), top5 = as.numeric(top5), top10 = as.numeric(top10))
props <- props %>% mutate(top = top * 100, top5 = top5 * 100, top10 = top10 * 100)

#create gt table for paper
props <- props %>% rownames_to_column() %>% mutate_all(as.character) %>% 
  pivot_longer(-rowname) %>% pivot_wider(names_from=rowname)
colnames(props) <- props[1,]
props <- props[-1,]
props$algopseudo <- as.factor(props$algopseudo)
props <- props %>% mutate_if(is.character, as.numeric)
levels(props$algopseudo) <- c("Top Model (%)", "Within 0.1 of Top Model (%)", "Within 0.05 of Top Model (%)")
props <- props %>% select(algopseudo,
                 `BART Background`, `BRT Background`, `RF Background`, `MaxEnt Background`, `GAMM Background`, `GLMM Background`,
                 `BART CRWs`, `BRT CRWs`, `RF CRWs`, `MaxEnt CRWs`, `GAMM CRWs`, `GLMM CRWs`,
                 `BART Buffer`, `BRT Buffer`, `RF Buffer`, `MaxEnt Buffer`, `GAMM Buffer`, `GLMM Buffer`)
tab <- props %>% gt() %>% cols_align(align = "center") %>%
  data_color(direction = "row", palette = "Purples", method = "numeric") 
tab <- tab %>%
  tab_spanner(label = "Background", columns = 2:7) %>%
  tab_spanner(label = "CRW", columns = 8:13) %>%
  tab_spanner(label = "Buffer", columns = 14:19) %>%
  tab_spanner(label = "Pseudo-Absence", columns = 1) %>%
  cols_label(algopseudo = "Algorithm",
    `RF Background` = "RF", `GAMM Background` = "GAMM", `GLMM Background` = "GLMM", `BRT Background` = "BRT", `BART Background` = "BART", `MaxEnt Background` = "MaxEnt",
    `RF CRWs` = "RF", `GAMM CRWs` = "GAMM", `GLMM CRWs` = "GLMM", `BRT CRWs` = "BRT", `BART CRWs` = "BART", `MaxEnt CRWs` = "MaxEnt",
    `RF Buffer` = "RF", `GAMM Buffer` = "GAMM", `GLMM Buffer` = "GLMM", `BRT Buffer` = "BRT", `BART Buffer` = "BART", `MaxEnt Buffer` = "MaxEnt") %>%
  cols_width(algopseudo ~ px(140), everything() ~ px(70)) 
tab

# gtsave(tab, filename = "text/figures/ecography/ts05_temp_table_thinned_spatiotemp.png",
#        vwidth = 1600)


#top ranked model scores
best <- filter(ranked, top == 1) 
min(best$transferability, na.rm=T)
best$algopseudo <- as.factor(best$algopseudo)

# relevel to put crw before buffer
# best$algopseudo <- fct_relevel(best$algopseudo, 
#              "BART Background", "BRT Background", "RF Background", "MaxEnt Background", "GAMM Background", "GLMM Background",
#              "BART CRWs", "BRT CRWs", "RF CRWs", "MaxEnt CRWs", "GAMM CRWs", "GLMM CRWs",
#              "BART Buffer", "BRT Buffer", "RF Buffer", "MaxEnt Buffer", "GAMM Buffer", "GLMM Buffer")
best$algopseudo <- fct_relevel(best$algopseudo, 
                               "BART Background", "BRT Background", "RF Background", "MaxEnt Background", "GAMM Background", "GLMM Background",
                               "BART CRWs", "BRT CRWs", "RF CRWs", "MaxEnt CRWs", "GAMM CRWs", "GLMM CRWs",
                               "BART Buffer", "RF Buffer", "GAMM Buffer", "GLMM Buffer")

# plot - CHANGE LEGEND
# p1 <- ggplot(best, aes(x=transferability)) + 
#   geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) + 
#   theme_classic() +
#   scale_fill_manual(values = c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "#820D0F",
#                                 "#e3f1ff", "#89bff8", "#408ee0", "#1b75be", "#0051a2", "#003366",
#                                "#ffe5ab", "#b8622f", "#883906", "#622A04")) +
#   guides(fill = guide_legend(ncol=3, title = "Algorithm and Pseudo-Absence Technique")) +
#   theme(legend.position = c(0.3, 0.9)) + 
#   scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) + 
#   ylab("Count") + xlab("Continuous Boyce Index Score")
p1 <- ggplot(best, aes(x=transferability)) + 
  geom_histogram(binwidth = 0.05, boundary = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) + 
  ylab("Count") + xlab("Continuous Boyce Index Score")
p1

###SPATIAL###
rm(list=setdiff(ls(), "p1"))
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")
library(tidyverse)
library(gt)

#read in Boyce scores and mixed-effects model
boyce <- readRDS("output/spatial/transferability_scores.RDS")

# rename site 
boyce <- boyce %>%
  rename(train_site = site) %>%
  rename(site = test_site)

#average rank (24 Highest 1 Lowest)
avg_rank <- boyce %>% arrange(transferability) %>% group_by(species, train_site, site, stage) %>%
  mutate(rank=rank(transferability, na.last = F)) %>% group_by(algopseudo) %>% summarise(average_rank = mean(rank)) %>%
  arrange(average_rank)

#individual ranks for each row
ranked <- boyce %>% arrange(transferability) %>% group_by(species, site, stage, train_site) %>%
  mutate(rank=round(rank(transferability, na.last = F)))
ggplot(ranked, aes(x=rank, fill=algopseudo)) + geom_bar(position="dodge", stat="count") +
  scale_fill_viridis_d() 

#number of case studies
n_case <- nrow(boyce)/24

# % of tests where model is top ranking
ranked$top <- if_else(ranked$rank == 24, 1, 0)
top <- ranked %>% group_by(algopseudo) %>% summarise(prop = sum(top)/n_case) 

# % where models are within 0.05 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, train_site) %>% 
  nest()

within_5 <- function(data) {
  data <- data %>% mutate(lim = max(transferability, na.rm=T)-0.05)
  data <- data %>% mutate(within = if_else(transferability > lim, 1, 0))
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
  data <- data %>% mutate(lim = max(transferability, na.rm=T)-0.1)
  data <- data %>% mutate(within = if_else(transferability > lim, 1, 0))
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

# #check examples where RF background isn't within 0.1
# rfback <- all %>% filter(algopseudo == "RF Background", within == 0)
# toptransferabilitys <- all %>% filter(top == 1) %>% 
#   select(species, stage, site, season, algopseudo, transferability) %>%
#   rename(topscore = transferability, topconf = algopseudo)
# rfback <- rfback %>% left_join(topscores)
# 
# brtback <- all %>% filter(algopseudo == "BRT Background") %>% 
#   select(species, site, stage, season, transferability) %>%
#   rename(brtback = transferability)
# rfback <- rfback %>% left_join(brtback)
# fail <- rfback %>% filter(transferability < 0.5 & brtback < 0.5)

#merge table
top <- top %>% rename(top = prop)
within_0.05 <- within_0.05 %>% rename(top5 = prop)
within_0.1 <- within_0.1 %>% rename(top10 = prop)
props <- top %>% left_join(within_0.05) %>% left_join(within_0.1)
props$top <- format(round(props$top, 3), nsmall = 3)
props$top5 <- format(round(props$top5, 3), nsmall = 3)
props$top10 <- format(round(props$top10, 3), nsmall = 3)
props <- props %>%
  mutate(top = as.numeric(top), top5 = as.numeric(top5), top10 = as.numeric(top10))
props <- props %>% mutate(top = top * 100, top5 = top5 * 100, top10 = top10 * 100)

#create gt table for paper
props <- props %>% rownames_to_column() %>% mutate_all(as.character) %>% 
  pivot_longer(-rowname) %>% pivot_wider(names_from=rowname)
colnames(props) <- props[1,]
props <- props[-1,]
props$algopseudo <- as.factor(props$algopseudo)
props <- props %>% mutate_if(is.character, as.numeric)
levels(props$algopseudo) <- c("Top Model (%)", "Within 0.1 of Top Model (%)", "Within 0.05 of Top Model (%)")
props <- props %>% select(algopseudo,
                          `BART Background`, `BRT Background`, `RF Background`, `MaxEnt Background`, `GAMM Background`, `GLMM Background`,
                          `BART CRWs`, `BRT CRWs`, `RF CRWs`, `MaxEnt CRWs`, `GAMM CRWs`, `GLMM CRWs`,
                          `BART Buffer`, `BRT Buffer`, `RF Buffer`, `MaxEnt Buffer`, `GAMM Buffer`, `GLMM Buffer`)
tab <- props %>% gt() %>% cols_align(align = "center") %>%
  data_color(direction = "row", palette = "Purples", method = "numeric") 
tab <- tab %>%
  tab_spanner(label = "Background", columns = 2:7) %>%
  tab_spanner(label = "CRW", columns = 8:13) %>%
  tab_spanner(label = "Buffer", columns = 14:19) %>%
  tab_spanner(label = "Pseudo-Absence", columns = 1) %>%
  cols_label(algopseudo = "Algorithm",
             `RF Background` = "RF", `GAMM Background` = "GAMM", `GLMM Background` = "GLMM", `BRT Background` = "BRT", `BART Background` = "BART", `MaxEnt Background` = "MaxEnt",
             `RF CRWs` = "RF", `GAMM CRWs` = "GAMM", `GLMM CRWs` = "GLMM", `BRT CRWs` = "BRT", `BART CRWs` = "BART", `MaxEnt CRWs` = "MaxEnt",
             `RF Buffer` = "RF", `GAMM Buffer` = "GAMM", `GLMM Buffer` = "GLMM", `BRT Buffer` = "BRT", `BART Buffer` = "BART", `MaxEnt Buffer` = "MaxEnt") %>%
  cols_width(algopseudo ~ px(140), everything() ~ px(70)) 
tab

# gtsave(tab, filename = "text/figures/ecography/ts06_spatial_table_thinned_spatiotemp.png",
#        vwidth = 1600)

#top ranked model scores
best <- filter(ranked, top == 1) 
min(best$transferability, na.rm=T)
best$algopseudo <- as.factor(best$algopseudo)

# relevel to put crw before buffer
best$algopseudo <- fct_relevel(best$algopseudo, 
                               "BART Background", "BRT Background", "RF Background", "MaxEnt Background", "GAMM Background", "GLMM Background",
                               "BART CRWs", "BRT CRWs", "RF CRWs", "MaxEnt CRWs", "GAMM CRWs", "GLMM CRWs",
                               "BART Buffer", "BRT Buffer", "RF Buffer", "MaxEnt Buffer", "GAMM Buffer", "GLMM Buffer")

# p2 <- ggplot(best, aes(x=transferability)) + 
#   geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) + 
#   theme_classic() +
#   scale_fill_manual(values = c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "#820D0F",
#                                "#e3f1ff", "#89bff8", "#408ee0", "#1b75be", "#0051a2", "#003366",
#                                "#ffe5ab", "#ffb87f", "#e88c56", "#b8622f", "#883906", "#622A04")) +
#   guides(fill = guide_legend(ncol=3, title = "Algorithm and Pseudo-Absence Technique")) +
#   theme(legend.position = c(0.3, 0.9)) + 
#   scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) + 
#   ylab("Count") + xlab("Continuous Boyce Index Score")
p2 <- ggplot(best, aes(x=transferability)) + 
  geom_histogram(binwidth = 0.05, boundary = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) + 
  ylab("Count") + xlab("Continuous Boyce Index Score")
p2

# # make version with full legend
# mod <- data.frame(
#   site = "test",
#   algorithm = "test",
#   train_site = "test",
#   species = "test",
#   stage = "test",
#   algopseudo = as.factor(c("GLMM Background", "GAMM Background", "RF Background", "BRT Background", "BART Background",
#                           "GLMM CRW", "GAMM CRW", "RF CRW", "BRT CRW", "BART CRW",
#                           "GLMM Buffer", "GAMM Buffer", "RF Buffer", "BRT Buffer", "BART Buffer")),
#   n = 1,
#   pseudo = "test",
#   score = 1,
#   shape = 1,
#   rank = 24,
#   top = 1)
# 
# # make algopseudo a factor
# mod$algopseudo <- as.factor(mod$algopseudo)
# mod$algopseudo <- fct_relevel(mod$algopseudo,
#                                "BART Background", "BRT Background", "RF Background", "GAMM Background", "GLMM Background",
#                                "BART CRW", "BRT CRW", "RF CRW", "GAMM CRW", "GLMM CRW",
#                                "BART Buffer", "BRT Buffer", "RF Buffer", "GAMM Buffer", "GLMM Buffer")
# unique(mod$algopseudo)
# 
# p3 <- ggplot(mod, aes(x=score)) +
#   geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) +
#   theme_classic() +
#   scale_fill_manual(values = c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15",
#                                "#e3f1ff", "#89bff8", "#408ee0", "#1b75be", "#0051a2",
#                                "#ffe5ab", "#ffb87f", "#e88c56", "#b8622f", "#883906")) +
#   guides(fill = guide_legend(ncol=3, title = "Algorithm and Pseudo-Absence Technique")) +
#   scale_x_continuous(limits = c(0,1.01), breaks = seq(0, 1, 0.2)) +
#   ylab("Count") + xlab("Continuous Boyce Index Score")
# p3


#plot both together
rm(list = setdiff(ls(), c("p1", "p2")))
#rm(list=setdiff(ls(), c("p1", "p2", "p3")))

# get legend
legend <- get_legend(p2 + theme(legend.position = "right"))
plot(legend)

#remove legends
p1 <- p1 + guides(fill="none")
p2 <- p2 + guides(fill="none")

#add titles
p1 <- p1 + ggtitle("Temporal Transfer")
p2 <- p2 + ggtitle("Spatial Transfer")

#customise font sizes
p1 <- p1 + theme_classic(base_size = 10)
p2 <- p2 + theme_classic(base_size = 10)

#plot without legend
plots <- plot_grid(p1, p2 + theme(legend.position="none"), ncol=2, align="v")
plots + ggview::canvas(width = 300, height = 150, units = "mm")

#plots and legend together
p3 <- plot_grid(legend, plots, ncol=1, rel_heights = c(2,6))
p3 + ggview::canvas(width = 300, height = 200, units = "mm")

#export
ggsave("text/figures/ecography/s04_spat_temp_transfer_spatiotemp.png", plots,
       width = 300, height = 150, units = "mm", dpi = 300)
