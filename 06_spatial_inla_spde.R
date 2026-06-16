#-------------------------------------------------------------------------------
# Spatial INLA-SPDE model
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Chapter 01")

library(INLA)
library(raster)
library(terra)
library(sp)
library(tidyverse)
library(tidysdm)
library(terra)
library(tidyterra)
library(fmesher)
library(miceRanger)

# read in table with info for each species, site and stage
meta <- read.csv("data/species_site_stage_metadata.csv")
meta2 <- read.csv("data/spatial_site_metadata.csv")

#remove ANPE, EMPE, SUFS, and MAPE incubation/post-breeding (no spatial transfer)
meta <- meta %>% filter(Species != "ANPE" & Species != "EMPE" & Species != "SUFS")
meta <- meta %>% filter(Species != "MAPE" | 
                          (Species == "MAPE" & Stage != "incubation" & Stage != "post-breeding"))

# run over each row of metadata
for(z in 14){
  print(z)
  
  #define parameters in loop
  rm(list=setdiff(ls(), c("meta", "meta2", "predictors", "z", "cores", "missing")))
  this.species <- meta[z, 1]
  this.site <- meta[z, 2]
  this.stage <- meta[z, 3]
  
  # read in presences with environmental data
  pres <- read.csv(paste0("output_old/extraction/", this.species, "/", this.site, "/", this.stage, "/presences.csv"))
  
  # read in pseudo-absences with environmental data
  pseudo <- readRDS(paste0("output/extraction/", this.species, "_", this.site, "_", this.stage, "_extracted.rds")) %>%
    filter(pb != "presence")
  
  # filter spatial meta to test sites
  spatial_tests <- meta2 %>%
    filter(Species == this.species, Stage == this.stage)
  
  # define missing var options
  missing_options <- spatial_tests %>%
    pull(Missing) %>%
    unique()
  
  # for each missing option
  for(missing in missing_options){
    
    #define initial predictors
    predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "slope")
    
    #remove predictors if missing
    if(missing == "chl"){
      predictors <- c("depth", "dshelf", "sst", "mld", "sal", "ssh", "sic", "curr", "eke", "slope")
    }
    
    # combine presences and pesudoabsences
    data <- bind_rows(
      pseudo %>% select(pb, all_of(predictors), individual_id, date, x, y, test_year),
      pres %>% mutate(pb = "presence", date = as_datetime(date)) %>% 
        select(pb, all_of(predictors), individual_id, date, x, y)
    )
    
    # isolate self-testing dataset from within training data
    self_test_inds <- data %>%
      filter(pb == "presence") %>%
      pull(individual_id) %>%
      unique() %>%
      sample(size = length(.)*0.2)
    self_test <- data %>% 
      filter(pb == "background" & is.na(test_year) | pb == "presence") %>%
      select(all_of(predictors), pb, individual_id, x, y) %>%
      mutate(pb = ordered(pb, levels = c("presence", "background"))) %>%
      filter(individual_id %in% self_test_inds)
    
    # # impute if any missing values
    # if(sum(is.na(self_test)) > 0){
    #   self_mice <- miceRanger(self_test, m=1)
    #   self_test <- completeData(self_mice)[[1]]
    # }
    
    # remove self-test data from training data
    train <- data %>%
      filter(!individual_id %in% self_test_inds)
    
    
    # get minimum number of psuedo-absences
    min_n <- train %>% filter(pb != "presence") %>%
      group_by(pb) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      pull(n) %>%
      min()
    
    # subsample all pseudoabsences
    pres <- train %>% filter(pb == "presence")
    train <- train %>%  filter(pb != "presence") %>%
      group_by(pb) %>%
      sample_n(min_n) %>%
      ungroup() %>%
      bind_rows(pres)
    
    #filter spatial metadata to this species and stage
    meta3 <- meta2 %>% filter(Species == this.species & Stage == this.stage & Missing == missing)
    
    #extract list of sites for this species and stage
    sites <- unique(meta3$Site)
    
    #null table for output
    sdpe_boyce_final <- NULL
    
    # test for each site
    for(test.site in sites){
      
      #load in test data
      test <- readRDS(paste0("output/extraction/spatial/", this.species, "_", test.site, "_", this.stage, "_extracted.rds"))
      
      #only select predictors for testing
      test <- test %>%
        select(all_of(predictors), pb, x, y) %>%
        mutate(individual_id = "test")
      
      #remove NA values from test data
      test <- test %>% drop_na()
      
      #order PA
      test$pb <- ordered(as.factor(test$pb), levels = c("presence", "background"))
      
      # isolate background and presence
      bg <- train %>%
        filter(pb %in% c("background", "presence"))
      
      # list coordinates
      coords <- bg %>% select(x, y)
      
      # make boundary for mesh
      bg_sf <- bg %>% vect(geom = c("x", "y"), crs = "EPSG:4326") %>% sf::st_as_sf()
      boundary <- fm_nonconvex_hull(bg_sf, format = "fm")
      
      # create mesh
      mesh <- fm_mesh_2d_inla(boundary = boundary, offset = c(0.5, 1), max.edge = c(45, 150))
      
      # create A matrix
      A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coords));dim(A)
      
      # create spatial structure
      spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
      
      # create indexes for SPDE model
      iset <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)
      
      # convert pb to 1-0
      bg$use <- ifelse(bg$pb == "presence", 1, 0)
      
      # create stack for INLA
      stk <- inla.stack(data = list(y = bg$use),
                        A = list(A, 1),
                        effects = list(c(list(Intercept = 1),
                                         iset),
                                       list(depth = bg$depth,
                                            dshelf = bg$dshelf,
                                            sst = bg$sst,
                                            mld = bg$mld,
                                            sal = bg$sal,
                                            ssh = bg$ssh,
                                            sic = bg$sic,
                                            curr = bg$curr,
                                            eke = bg$eke,
                                            slope = bg$slope)))
      
      # model formula
      form <- y ~ -1 + Intercept + depth + dshelf + sst + mld + sal + ssh + sic + curr + eke + slope +
        f(spatial.field, model = spde)
      
      # fit model
      m0 <- inla(form,
                 data = inla.stack.data(stk),
                 family = "binomial",
                 control.predictor = list(A = inla.stack.A(stk)),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                 verbose = F)
      
      # to predict we also need the A matrix for the prediction coordinates
      pred_coords <- test %>% select(x, y)
      Apred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pred_coords))
      
      # create prediction stack
      stk.pred <- inla.stack(data = list(y = NA),
                             A = list(Apred, 1),
                             effects = list(c(list(Intercept = 1),
                                              iset),
                                            list(depth = test$depth,
                                                 dshelf = test$dshelf,
                                                 sst = test$sst,
                                                 mld = test$mld,
                                                 sal = test$sal,
                                                 ssh = test$ssh,
                                                 sic = test$sic,
                                                 curr = test$curr,
                                                 eke = test$eke,
                                                 slope = test$slope)),
                             tag = "pred")
      
      # join prediction stack with full stack
      stk.full <- inla.stack(stk, stk.pred)
      
      # prediction inla
      p.res.pred <- inla(form, data = inla.stack.data(stk.full, spde = spde), 
                         family = "binomial", 
                         control.predictor = list(A = inla.stack.A(stk.full), compute = F),
                         control.compute = list(config = T),
                         control.inla(strategy = "simplified.laplace", huge = T),
                         verbose = F)
      
      # extract indices of prediction nodes and posterior mean of response
      index.pred <- inla.stack.index(stk.full, tag = "pred")$data
      post.mean.pred.logit <- p.res.pred$summary.fitted.values[index.pred, "mean"]
      p.pred <- exp(post.mean.pred.logit) / (1 + exp(post.mean.pred.logit))
      
      # assign to test data
      test$bg_pred <- p.pred
      
      # calculate Boyce index
      back_boyce <- boyce_cont(test, pb, bg_pred) %>%
        pull(.estimate)
      
      # repeat for self test data
      self_test_coords <- self_test %>% select(x, y)
      Aself_test <- inla.spde.make.A(mesh = mesh, loc = as.matrix(self_test_coords))
      stk.self_test <- inla.stack(data = list(y = NA),
                                  A = list(Aself_test, 1),
                                  effects = list(c(list(Intercept = 1),
                                                   iset),
                                                 list(depth = self_test$depth,
                                                      dshelf = self_test$dshelf,
                                                      sst = self_test$sst,
                                                      mld = self_test$mld,
                                                      sal = self_test$sal,
                                                      ssh = self_test$ssh,
                                                      sic = self_test$sic,
                                                      curr = self_test$curr,
                                                      eke = self_test$eke,
                                                      slope = self_test$slope)),
                                  tag = "self_test")
      stk.full_self_test <- inla.stack(stk, stk.self_test)
      p.res.self_test <- inla(form, data = inla.stack.data(stk.full_self_test, spde = spde), 
                              family = "binomial", 
                              control.predictor = list(A = inla.stack.A(stk.full_self_test), compute = F),
                              control.compute = list(config = T),
                              control.inla(strategy = "simplified.laplace", huge = T),
                              verbose = F)
      index.self_test <- inla.stack.index(stk.full_self_test, tag = "self_test")$data
      post.mean.self_test.logit <- p.res.self_test$summary.fitted.values[index.self_test, "mean"]
      self_test$bg_pred <- exp(post.mean.self_test.logit) / (1 + exp(post.mean.self_test.logit))
      self_back_boyce <- boyce_cont(self_test, pb, bg_pred) %>%
        pull(.estimate)
      
      
      # 2. buffer
      
      # isolate buffer and presence
      buff <- train %>%
        filter(pb %in% c("buffer", "presence"))
      
      # list coordinates
      coords <- buff %>% select(x, y)
      
      # make boundary for mesh
      buff_sf <- buff %>% vect(geom = c("x", "y"), crs = "EPSG:4326") %>% sf::st_as_sf()
      boundary <- fm_nonconvex_hull(buff_sf, format = "fm")
      
      # create mesh
      mesh <- fm_mesh_2d_inla(boundary = boundary, offset = c(0.5, 1), max.edge = c(45, 150))
      
      # create A matrix
      A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coords));dim(A)
      
      # create spatial structure
      spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
      
      # create indexes for SPDE model
      iset <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)
      
      # convert pb to 1-0
      buff$use <- ifelse(buff$pb == "presence", 1, 0)
      
      # create stack for INLA
      stk <- inla.stack(data = list(y = buff$use),
                        A = list(A, 1),
                        effects = list(c(list(Intercept = 1),
                                         iset),
                                       list(depth = buff$depth,
                                            dshelf = buff$dshelf,
                                            sst = buff$sst,
                                            mld = buff$mld,
                                            sal = buff$sal,
                                            ssh = buff$ssh,
                                            sic = buff$sic,
                                            curr = buff$curr,
                                            eke = buff$eke,
                                            slope = buff$slope)))
      
      # model formula
      form <- y ~ -1 + Intercept + depth + dshelf + sst + mld + sal + ssh + sic + curr + eke + slope +
        f(spatial.field, model = spde)
      
      # fit model
      m0 <- inla(form,
                 data = inla.stack.data(stk),
                 family = "binomial",
                 control.predictor = list(A = inla.stack.A(stk)),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                 verbose = F)
      
      # to predict we also need the A matrix for the prediction coordinates
      pred_coords <- test %>% select(x, y)
      Apred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pred_coords))
      
      # create prediction stack
      stk.pred <- inla.stack(data = list(y = NA),
                             A = list(Apred, 1),
                             effects = list(c(list(Intercept = 1),
                                              iset),
                                            list(depth = test$depth,
                                                 dshelf = test$dshelf,
                                                 sst = test$sst,
                                                 mld = test$mld,
                                                 sal = test$sal,
                                                 ssh = test$ssh,
                                                 sic = test$sic,
                                                 curr = test$curr,
                                                 eke = test$eke,
                                                 slope = test$slope)),
                             tag = "pred")
      
      # join prediction stack with full stack
      stk.full <- inla.stack(stk, stk.pred)
      
      # prediction inla
      p.res.pred <- inla(form, data = inla.stack.data(stk.full, spde = spde), 
                         family = "binomial", 
                         control.predictor = list(A = inla.stack.A(stk.full), compute = F),
                         control.compute = list(config = T),
                         control.inla(strategy = "simplified.laplace", huge = T),
                         verbose = F)
      
      # extract indices of prediction nodes and posterior mean of response
      index.pred <- inla.stack.index(stk.full, tag = "pred")$data
      post.mean.pred.logit <- p.res.pred$summary.fitted.values[index.pred, "mean"]
      p.pred <- exp(post.mean.pred.logit) / (1 + exp(post.mean.pred.logit))
      
      # assign to test data
      test$buff_pred <- p.pred
      
      # calculate Boyce index
      buff_boyce <- boyce_cont(test, pb, buff_pred) %>%
        pull(.estimate)
      
      # repeat for self test data
      self_test_coords <- self_test %>% select(x, y)
      Aself_test <- inla.spde.make.A(mesh = mesh, loc = as.matrix(self_test_coords))
      stk.self_test <- inla.stack(data = list(y = NA),
                                  A = list(Aself_test, 1),
                                  effects = list(c(list(Intercept = 1),
                                                   iset),
                                                 list(depth = self_test$depth,
                                                      dshelf = self_test$dshelf,
                                                      sst = self_test$sst,
                                                      mld = self_test$mld,
                                                      sal = self_test$sal,
                                                      ssh = self_test$ssh,
                                                      sic = self_test$sic,
                                                      curr = self_test$curr,
                                                      eke = self_test$eke,
                                                      slope = self_test$slope)),
                                  tag = "self_test")
      stk.full_self_test <- inla.stack(stk, stk.self_test)
      p.res.self_test <- inla(form, data = inla.stack.data(stk.full_self_test, spde = spde), 
                              family = "binomial", 
                              control.predictor = list(A = inla.stack.A(stk.full_self_test), compute = F),
                              control.compute = list(config = T),
                              control.inla(strategy = "simplified.laplace", huge = T),
                              verbose = F)
      index.self_test <- inla.stack.index(stk.full_self_test, tag = "self_test")$data
      post.mean.self_test.logit <- p.res.self_test$summary.fitted.values[index.self_test, "mean"]
      self_test$buff_pred <- exp(post.mean.self_test.logit) / (1 + exp(post.mean.self_test.logit))
      self_buff_boyce <- boyce_cont(self_test, pb, buff_pred) %>%
        pull(.estimate)
      
      
      
      # 3. crw
      
      # isolate crw and presence
      crw <- train %>%
        filter(pb %in% c("crw", "presence"))
      
      # list coordinates
      coords <- crw %>% select(x, y)
      
      # make boundary for mesh
      crw_sf <- crw %>% vect(geom = c("x", "y"), crs = "EPSG:4326") %>% sf::st_as_sf()
      boundary <- fm_nonconvex_hull(crw_sf, format = "fm")
      
      # create mesh
      mesh <- fm_mesh_2d_inla(boundary = boundary, offset = c(0.5, 1), max.edge = c(45, 150))
      
      # create A matrix
      A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coords));dim(A)
      
      # create spatial structure
      spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
      
      # create indexes for SPDE model
      iset <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)
      
      # convert pb to 1-0
      crw$use <- ifelse(crw$pb == "presence", 1, 0)
      
      # create stack for INLA
      stk <- inla.stack(data = list(y = crw$use),
                        A = list(A, 1),
                        effects = list(c(list(Intercept = 1),
                                         iset),
                                       list(depth = crw$depth,
                                            dshelf = crw$dshelf,
                                            sst = crw$sst,
                                            mld = crw$mld,
                                            sal = crw$sal,
                                            ssh = crw$ssh,
                                            sic = crw$sic,
                                            curr = crw$curr,
                                            eke = crw$eke,
                                            slope = crw$slope)))
      
      # model formula
      form <- y ~ -1 + Intercept + depth + dshelf + sst + mld + sal + ssh + sic + curr + eke + slope +
        f(spatial.field, model = spde)
      
      # fit model
      m0 <- inla(form,
                 data = inla.stack.data(stk),
                 family = "binomial",
                 control.predictor = list(A = inla.stack.A(stk)),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                 verbose = F)
      
      # to predict we also need the A matrix for the prediction coordinates
      pred_coords <- test %>% select(x, y)
      Apred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pred_coords))
      
      # create prediction stack
      stk.pred <- inla.stack(data = list(y = NA),
                             A = list(Apred, 1),
                             effects = list(c(list(Intercept = 1),
                                              iset),
                                            list(depth = test$depth,
                                                 dshelf = test$dshelf,
                                                 sst = test$sst,
                                                 mld = test$mld,
                                                 sal = test$sal,
                                                 ssh = test$ssh,
                                                 sic = test$sic,
                                                 curr = test$curr,
                                                 eke = test$eke,
                                                 slope = test$slope)),
                             tag = "pred")
      
      # join prediction stack with full stack
      stk.full <- inla.stack(stk, stk.pred)
      
      # prediction inla
      p.res.pred <- inla(form, data = inla.stack.data(stk.full, spde = spde), 
                         family = "binomial", 
                         control.predictor = list(A = inla.stack.A(stk.full), compute = F),
                         control.compute = list(config = T),
                         control.inla(strategy = "simplified.laplace", huge = T),
                         verbose = F)
      
      # extract indices of prediction nodes and posterior mean of response
      index.pred <- inla.stack.index(stk.full, tag = "pred")$data
      post.mean.pred.logit <- p.res.pred$summary.fitted.values[index.pred, "mean"]
      p.pred <- exp(post.mean.pred.logit) / (1 + exp(post.mean.pred.logit))
      
      # assign to test data
      test$crw_pred <- p.pred
      
      # calculate Boyce index
      crw_boyce <- boyce_cont(test, pb, crw_pred) %>%
        pull(.estimate)
      
      # repeat for self test data
      self_test_coords <- self_test %>% select(x, y)
      Aself_test <- inla.spde.make.A(mesh = mesh, loc = as.matrix(self_test_coords))
      stk.self_test <- inla.stack(data = list(y = NA),
                                  A = list(Aself_test, 1),
                                  effects = list(c(list(Intercept = 1),
                                                   iset),
                                                 list(depth = self_test$depth,
                                                      dshelf = self_test$dshelf,
                                                      sst = self_test$sst,
                                                      mld = self_test$mld,
                                                      sal = self_test$sal,
                                                      ssh = self_test$ssh,
                                                      sic = self_test$sic,
                                                      curr = self_test$curr,
                                                      eke = self_test$eke,
                                                      slope = self_test$slope)),
                                  tag = "self_test")
      stk.full_self_test <- inla.stack(stk, stk.self_test)
      p.res.self_test <- inla(form, data = inla.stack.data(stk.full_self_test, spde = spde), 
                              family = "binomial", 
                              control.predictor = list(A = inla.stack.A(stk.full_self_test), compute = F),
                              control.compute = list(config = T),
                              control.inla(strategy = "simplified.laplace", huge = T),
                              verbose = F)
      index.self_test <- inla.stack.index(stk.full_self_test, tag = "self_test")$data
      post.mean.self_test.logit <- p.res.self_test$summary.fitted.values[index.self_test, "mean"]
      self_test$crw_pred <- exp(post.mean.self_test.logit) / (1 + exp(post.mean.self_test.logit))
      self_crw_boyce <- boyce_cont(self_test, pb, crw_pred) %>%
        pull(.estimate)
      
      #combine boyce scores
      sdpe_boyce <- expand.grid(buff = buff_boyce, back = back_boyce, crw = crw_boyce)
      sdpe_boyce$test_site <- test.site
      sdpe_boyce_final <- rbind(sdpe_boyce_final, sdpe_boyce)
      
    }
    
    # collate self test scores
    self <- expand.grid(buff = self_buff_boyce, back = self_back_boyce, crw = self_crw_boyce) %>%
      mutate(missing_covar = missing)
    
    # join to other self test scores
    if(missing == missing_options[1]){
      all_self <- self
    } else {
      all_self <- rbind(all_self, self)
    }
    
    # join together boyce scores
    if(missing == missing_options[1]){
      all_boyce <- sdpe_boyce_final
    } else {
      all_boyce <- rbind(all_boyce, sdpe_boyce_final)
    }
  }
  
  # export boyce scores
  saveRDS(all_boyce, 
          file = paste0("output/spatial/sdpe/", this.species, "_", this.site, "_", this.stage, "_boyce_scores_sdpe.RDS"))
  
  # export self test 
  saveRDS(all_self,
          file = paste0("output/spatial/sdpe/", this.species, "_", this.site, "_", this.stage, "_self_test_sdpe.RDS"))
  
  #show species has finished
  print(paste0(this.species, " ", this.site, " ", this.stage, " success"))
}
