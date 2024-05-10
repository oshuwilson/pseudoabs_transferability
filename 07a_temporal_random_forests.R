#automated script to run Random Forests for leave-year-out validation
#Change cv_scheme to 10?
#cross validation over individuals? how for pseudo-abs?
rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

{
  library(dplyr)
  library(lubridate)
  library(ranger)
  library(caret)
  library(miceRanger)
  library(enmSdmX)
  library(flexsdm)
}

meta <- read.csv("data/species_site_stage_metadata.csv")
predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope")

for(z in 1:19){
try({
  
#define parameters in loop
rm(list=setdiff(ls(), c("meta", "predictors", "z")))
this.species <- meta[z, 1]
this.site <- meta[z, 2]
this.stage <- meta[z, 3]
season <- meta[z, 4]

# 1. Formatting
#read in presences and pseudo-absences
tracks <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/presences.csv"))
buff <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/buffers.csv"))
back <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/background.csv"))
crw <- read.csv(paste0("output/extraction/", this.species, "/", this.site, "/", this.stage, "/CRWs.csv"))

#create presence/absence column
tracks$pa <- 1
buff$pa <- 0
back$pa <- 0
crw$pa <- 0

#remove extra columns to allow rbind
columns <- c("date", "depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "chl", "wind", "slope", "x", "y", "pa")
tracks <- tracks %>% select(all_of(columns))
buff <- buff %>% select(all_of(columns))
back <- back %>% select(all_of(columns))
crw <- crw %>% select(all_of(columns))

#keep background points for testing
back_test <- back

#rbind for models
buff <- rbind(tracks, buff)
back <- rbind(tracks, back)
crw <- rbind(tracks, crw)

#setup date column
buff$date <- as_date(buff$date)
back$date <- as_date(back$date)
crw$date <- as_date(crw$date)
tracks$date <- as_date(tracks$date)


# 2. Isolate training and testing data

#if season = FALSE, separate by year
if(season == FALSE){
  buff$season <- year(buff$date)
  back$season <- year(back$date)
  crw$season <- year(crw$date)
  tracks$season <- year(tracks$date)
}

#if season = TRUE, separate by season
if(season == TRUE){
  buff$season <- year(round_date(buff$date, unit="year"))
  back$season <- year(round_date(back$date, unit="year"))
  crw$season <- year(round_date(crw$date, unit="year"))
  tracks$season <- year(round_date(tracks$date, unit="year"))
}

#extract seasons for loop
seasons <- levels(as.factor(tracks$season))

#null tables for output
rf_mtry_meta <- NULL
rf_boyce_final <- NULL
counts <- NULL

#loop starts here
for(i in seasons){
this.test <- i

#extract training data
buff_train <- buff %>% filter(season != this.test)
back_train <- back %>% filter(season != this.test)
crw_train <- crw %>% filter(season != this.test)
tracks_train <- tracks %>% filter(season != this.test)

#extract testing data
back_test <- back %>% filter(season == this.test, pa == 0)
tracks_test <- tracks %>% filter(season == this.test)

#if test data is missing over 10% of a predictor, remove predictor from models
pred_check <- back %>% filter(season == this.test) %>% select(all_of(predictors))

if(sum(is.na(pred_check)) > 0.1*nrow(pred_check)){
  pred_check <- pred_check[colSums(is.na(pred_check)) < 0.1*nrow(pred_check)]
}

predictors <- names(pred_check)


# 3. RF Predictions - Tune mtry

#create parameter grid to vary mtry between 2, 3, and 4
param_grid <- expand.grid(mtry=2:4, splitrule = "gini", min.node.size=1)

#setup 5-fold cross-validation
cv_scheme <- trainControl(method = "cv", number = 10, verboseIter = FALSE,
                          summaryFunction = twoClassSummary, classProbs = TRUE)

#remove NAs from test data
back_test <- back_test %>% select(all_of(predictors)) %>% na.omit()
tracks_test <- tracks_test %>% select(all_of(predictors)) %>% na.omit()


#BUFFERS
#remove non-predictor columns
buff_sel <- buff_train %>% select(all_of(predictors), pa)

#make presence-absence a character name
buff_sel$pa <- if_else(buff_sel$pa == 1, "presence", "absence")

#check for NA - less than 10% of training data okay for imputing
if(sum(is.na(buff_sel)) < 0.1*nrow(buff_sel) & sum(is.na(buff_sel)) > 0){
  buff_mice <- miceRanger(buff_sel, m=1)
  buff_sel <- completeData(buff_mice)[[1]]
}

#remove columns where missing data is over 10% of rows then impute
if(sum(is.na(buff_sel)) > 0.1*nrow(buff_sel)){
  buff_sel <- buff_sel[colSums(is.na(buff_sel)) < 0.1*nrow(buff_sel)]
  buff_mice <- miceRanger(buff_sel, m=1)
  buff_sel <- completeData(buff_mice)[[1]]
}

#perform tuning search
X <- buff_sel %>% select(-pa)
Y <- as.factor(buff_sel$pa)
buff_rf <- train(x = X, y = Y, method = "ranger", metric = "ROC", trControl = cv_scheme, 
                 tuneGrid = param_grid, num.trees = 1000, importance = "impurity")

#save mtry results
buff_mtry <- buff_rf$results[,c(1,4)]
buff_mtry$pseudo <- "buffer"
buff_mtry$season <- i

#predict and evaluate
p1 <- predict(buff_rf, tracks_test, type = "prob")[,2]
p2 <- predict(buff_rf, back_test, type = "prob")[,2]
buff_rf_boyce <- evalContBoyce(p1, p2)

#save model
saveRDS(buff_rf, 
        file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/buff_rf_", this.test, ".RDS"))


#remove unnecessary parameters to continue
rm(buff_rf, buff_mice, buff_sel, X, Y, p1, p2)


#BACKGROUND
#remove non-predictor columns
back_sel <- back_train %>% select(all_of(predictors), pa)

#make presence-absence a character name
back_sel$pa <- if_else(back_sel$pa == 1, "presence", "absence")

#check for NA - less than 10% of training data okay for imputing
if(sum(is.na(back_sel)) < 0.1*nrow(back_sel) & sum(is.na(back_sel)) > 0){
  back_mice <- miceRanger(back_sel, m=1)
  back_sel <- completeData(back_mice)[[1]]
}

#remove columns where missing data is over 10% of rows then impute
if(sum(is.na(back_sel)) > 0.1*nrow(back_sel)){
  back_sel <- back_sel[colSums(is.na(back_sel)) < 0.1*nrow(back_sel)]
  back_mice <- miceRanger(back_sel, m=1)
  back_sel <- completeData(back_mice)[[1]]
}

#perform tuning search
X <- back_sel %>% select(-pa)
Y <- as.factor(back_sel$pa)
back_rf <- train(x = X, y = Y, method = "ranger", metric = "ROC", trControl = cv_scheme, 
                 tuneGrid = param_grid, num.trees = 1000, importance = "impurity")

#save mtry results
back_mtry <- back_rf$results[,c(1,4)]
back_mtry$pseudo <- "background"
back_mtry$season <- i

#predict and evaluate
p1 <- predict(back_rf, tracks_test, type = "prob")[,2]
p2 <- predict(back_rf, back_test, type = "prob")[,2]
back_rf_boyce <- evalContBoyce(p1, p2)

#save model
saveRDS(back_rf, 
        file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/back_rf_", this.test, ".RDS"))

#remove unnecessary parameters to continue
rm(back_rf, back_mice, back_sel, X, Y, p1, p2)


#CRWs
#remove non-predictor columns
crw_sel <- crw_train %>% select(all_of(predictors), pa)

#make presence-absence a character name
crw_sel$pa <- if_else(crw_sel$pa == 1, "presence", "absence")

#check for NA - less than 10% of training data okay for imputing
if(sum(is.na(crw_sel)) < 0.1*nrow(crw_sel) & sum(is.na(crw_sel)) > 0){
  crw_mice <- miceRanger(crw_sel, m=1)
  crw_sel <- completeData(crw_mice)[[1]]
}

#remove columns where missing data is over 10% of rows then impute
if(sum(is.na(crw_sel)) > 0.1*nrow(crw_sel)){
  crw_sel <- crw_sel[colSums(is.na(crw_sel)) < 0.1*nrow(crw_sel)]
  crw_mice <- miceRanger(crw_sel, m=1)
  crw_sel <- completeData(crw_mice)[[1]]
}

#perform tuning search
X <- crw_sel %>% select(-pa)
Y <- as.factor(crw_sel$pa)
crw_rf <- train(x = X, y = Y, method = "ranger", metric = "ROC", trControl = cv_scheme, 
                 tuneGrid = param_grid, num.trees = 1000, importance = "impurity")

#save mtry results
crw_mtry <- crw_rf$results[,c(1,4)]
crw_mtry$pseudo <- "crw"
crw_mtry$season <- i

#predict and evaluate
p1 <- predict(crw_rf, tracks_test, type = "prob")[,2]
p2 <- predict(crw_rf, back_test, type = "prob")[,2]
crw_rf_boyce <- evalContBoyce(p1, p2)

#save model to scratch
saveRDS(crw_rf, 
        file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/crw_rf_", this.test, ".RDS"))

#remove unnecessary parameters to continue
rm(crw_rf, crw_mice, crw_sel, X, Y, p1, p2)


#FINAL DATA
rf_boyce <- expand.grid(buff = buff_rf_boyce, back = back_rf_boyce, crw = crw_rf_boyce)
rf_boyce$season <- i
rf_boyce_final <- rbind(rf_boyce_final, rf_boyce)

mtry_values <- rbind(buff_mtry, back_mtry, crw_mtry)
rf_mtry_meta <- rbind(rf_mtry_meta, mtry_values)

n <- length(tracks_train$depth)
ndays <- n_distinct(tracks_train$date)
nyears <- n_distinct(tracks_train$season)

season_count <- expand.grid(n = n, ndays = ndays, nyears = nyears)
season_count$season <- i
counts <- rbind(counts, season_count)
}


# 4. Export Boyce, Mtry, and Metadata
saveRDS(rf_boyce_final, 
        file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/boyce_scores_rf.RDS"))
saveRDS(rf_mtry_meta, 
        file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/mtry_values.RDS"))
saveRDS(counts,
        file = paste0("output/leave-year-out/", this.species, "/", this.site, "/", this.stage, "/sample_sizes.RDS"))

print(paste0(this.species, " ", this.site, " ", this.stage, " success"))

}) 

}
