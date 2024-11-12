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
boyce <- readRDS("output/leave-year-out/boyce_filtered.RDS")

#average rank (9 Highest 1 Lowest)
avg_rank <- boyce %>% arrange(score) %>% group_by(species, site, stage, season) %>%
  mutate(rank=rank(score)) %>% group_by(algopseudo) %>% summarise(average_rank = mean(rank)) %>%
  arrange(average_rank)

#individual ranks for each row
ranked <- boyce %>% arrange(score) %>% group_by(species, site, stage, season) %>%
  mutate(rank=round(rank(score))) 
ggplot(ranked, aes(x=rank, fill=algopseudo)) + geom_bar(position="dodge", stat="count") +
  scale_fill_viridis_d() 

# % of tests where model is top ranking
ranked$top <- if_else(ranked$rank == 9, 1, 0)
top <- ranked %>% group_by(algopseudo) %>% summarise(prop = sum(top)/84) 


# % where models are within 0.05 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, season) %>% 
  nest()

within_5 <- function(data) {
  data <- data %>% mutate(lim = max(score)-0.05)
  data <- data %>% mutate(within = if_else(score > lim, 1, 0))
  print(data)
}

test <- comp[1,] %>% unnest(cols = c(data))
test <- within_5(data = test)

all <- NULL

for(i in 1:84){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_5(data = test)
  all <- rbind(all, test)
}

within_0.05 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/84) 


# % where models are within 0.1 of top (i.e. comparable)
comp <- ranked %>% group_by(species, site, stage, season) %>% 
  nest()

within_10 <- function(data) {
  data <- data %>% mutate(lim = max(score)-0.1)
  data <- data %>% mutate(within = if_else(score > lim, 1, 0))
  print(data)
}

test <- comp[1,] %>% unnest(cols = c(data))
test <- within_10(data = test)

all <- NULL

for(i in 1:84){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_10(data = test)
  all <- rbind(all, test)
}

within_0.1 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/84) 

#check examples where RF background isn't within 0.1
rfback <- all %>% filter(algopseudo == "RF back", within == 0)
topscores <- all %>% filter(top == 1) %>% 
  select(species, stage, site, season, algopseudo, score) %>%
  rename(topscore = score, topconf = algopseudo)
rfback <- rfback %>% left_join(topscores)

gamback <- all %>% filter(algopseudo == "GAM back") %>% 
  select(species, site, stage, season, score) %>%
  rename(gamback = score)
rfback <- rfback %>% left_join(gamback)
fail <- rfback %>% filter(score < 0.5 & gamback < 0.5)

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
                              "GAM Background", "GAM Buffer", "GAM CRW",
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
                 `RF Background`, `GAM Background`, `BRT Background`,
                 `RF CRW`, `GAM CRW`, `BRT CRW`,
                 `RF Buffer`, `GAM Buffer`, `BRT Buffer`)
tab <- props %>% gt() %>% cols_align(align = "center") %>%
  data_color(direction = "row", palette = "Purples", method = "numeric") 
tab <- tab %>%
  tab_spanner(label = "Background", columns = 2:4) %>%
  tab_spanner(label = "CRW", columns = 5:7) %>%
  tab_spanner(label = "Buffer", columns = 8:10) %>%
  tab_spanner(label = "Pseudo-Absence", columns = 1) %>%
  cols_label(algopseudo = "Algorithm",
    `RF Background` = "RF", `GAM Background` = "GAM", `BRT Background` = "BRT",
    `RF CRW` = "RF", `GAM CRW` = "GAM", `BRT CRW` = "BRT",
    `RF Buffer` = "RF", `GAM Buffer` = "GAM", `BRT Buffer` = "BRT") %>%
  cols_width(algopseudo ~ px(140), everything() ~ px(70)) 
tab

gtsave(tab, filename = "text/figures/temp table.png")


#top ranked model scores
best <- filter(ranked, top == 1) 
min(best$score, na.rm=T)

levels(best$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW", 
                             "GAM Background", "GAM Buffer", "GAM CRW",
                             "RF Background", "RF Buffer", "RF CRW")

p1 <- ggplot(best, aes(x=score)) + 
  geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) + 
  theme_classic() +
  scale_fill_manual(values = c("red1", "red3", "red4",
                               "steelblue1", "steelblue", "steelblue4",
                               "lightsalmon", "darkorange2", "darkorange3"),
                    breaks = c("BRT Background", "GAM Background", "RF Background",
                               "BRT Buffer", "GAM Buffer", "RF Buffer",
                               "BRT CRW", "GAM CRW", "RF CRW")) +
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
boyce <- readRDS("output/spatial/boyce_final.RDS")

#average rank (9 Highest 1 Lowest)
avg_rank <- boyce %>% arrange(score) %>% group_by(species, site, stage, train_site) %>%
  mutate(rank=rank(score)) %>% group_by(algopseudo) %>% summarise(average_rank = mean(rank)) %>%
  arrange(average_rank)

#individual ranks for each row
ranked <- boyce %>% arrange(score) %>% group_by(species, site, stage, train_site) %>%
  mutate(rank=round(rank(score))) 
ggplot(ranked, aes(x=rank, fill=algopseudo)) + geom_bar(position="dodge", stat="count") +
  scale_fill_viridis_d() 

# % of tests where model is top ranking
ranked$top <- if_else(ranked$rank == 9, 1, 0)
top <- ranked %>% group_by(algopseudo) %>% summarise(prop = sum(top)/42) 


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

for(i in 1:42){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_5(data = test)
  all <- rbind(all, test)
}

within_0.05 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/42) 


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

for(i in 1:42){
  test <- comp[i,] %>% unnest(cols = c(data))
  test <- within_10(data = test)
  all <- rbind(all, test)
}

within_0.1 <- all %>% group_by(algopseudo) %>% summarise(prop = sum(within, na.rm=T)/42) 


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
                              "GAM Background", "GAM Buffer", "GAM CRW",
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
                          `RF Background`, `GAM Background`, `BRT Background`,
                          `RF CRW`, `GAM CRW`, `BRT CRW`,
                          `RF Buffer`, `GAM Buffer`, `BRT Buffer`)
tab <- props %>% gt() %>% cols_align(align = "center") %>%
  data_color(direction = "row", palette = "Purples", method = "numeric") 
tab <- tab %>%
  tab_spanner(label = "Background", columns = 2:4) %>%
  tab_spanner(label = "CRW", columns = 5:7) %>%
  tab_spanner(label = "Buffer", columns = 8:10) %>%
  tab_spanner(label = "Pseudo-Absence", columns = 1) %>%
  cols_label(algopseudo = "Algorithm",
             `RF Background` = "RF", `GAM Background` = "GAM", `BRT Background` = "BRT",
             `RF CRW` = "RF", `GAM CRW` = "GAM", `BRT CRW` = "BRT",
             `RF Buffer` = "RF", `GAM Buffer` = "GAM", `BRT Buffer` = "BRT") %>%
  cols_width(algopseudo ~ px(140), everything() ~ px(70)) 
tab

gtsave(tab, filename = "text/figures/spat table.png")


#top ranked model scores
best <- filter(ranked, top == 1) 
min(best$score, na.rm=T)

levels(best$algopseudo) <- c("BRT Background", "BRT Buffer", "BRT CRW", 
                             "GAM Background", "GAM Buffer", "GAM CRW",
                             "RF Background", "RF Buffer", "RF CRW")

p2 <- ggplot(best, aes(x=score)) + 
  geom_histogram(binwidth = 0.05, boundary = 1, aes(fill=algopseudo)) + 
  theme_classic() +
  scale_fill_manual(values = c("red1", "red3", "red4",
                               "steelblue1", "steelblue", "steelblue4",
                               "lightsalmon", "darkorange2", "darkorange3"),
                    breaks = c("BRT Background", "GAM Background", "RF Background",
                               "BRT Buffer", "GAM Buffer", "RF Buffer",
                               "BRT CRW", "GAM CRW", "RF CRW")) +
  guides(fill = guide_legend(ncol=3, title = "Algorithm and Pseudo-Absence Technique")) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) + 
  ylab("Count") + xlab("Continuous Boyce Index Score")
p2 

#plot both together
rm(list=setdiff(ls(), c("p1", "p2")))

#remove legends
p1 <- p1 + guides(fill="none")

#get legend
legend <- get_legend(p2)
legend

#add titles
p1 <- p1 + ggtitle("Temporal Transfer")
p2 <- p2 + ggtitle("Spatial Transfer")

plots <- plot_grid(p1, p2 + theme(legend.position="none", ncol=2, align="v"))
plots

#plots and legend together
plot_grid(legend, plots, ncol=1, rel_heights = c(1,4))
