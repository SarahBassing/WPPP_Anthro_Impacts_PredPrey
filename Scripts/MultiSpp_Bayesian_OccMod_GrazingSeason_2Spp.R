  #'  ==========================================================================
  #'  Bayesian Multispecies Occupancy Models with 2 species: GRAZING SEASON DATA 
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  July 2022
  #'  ===========================================================================
  #'  Script to source data formatting scripts & run multi-species, single-season
  #'  occupancy models within a Bayesian framework. Code adapted from Kery & Royle
  #'  AHM2 book (Ch 8.2.3). Code runs 3-spp models using detection/non-detection
  #'  data for deer, elk, moose, black bears, cougars, wolves, coyotes, & bobcats
  #'  during the grazing season 2018-2020 (July - Sept). Grazing season co-occurrence 
  #'  models include 13 7-day sampling occasions comprising the peak of livestock
  #'  activity detected on camera. Models test whether predator-prey co-occurrence 
  #'  is not independent and whether their occurrence, co-occurrence, and detection 
  #'  are influenced by livestock activity at camera sites.
  #'  
  #'  Encounter histories are generated with CameratTrap_DetectionHistories_for_unmarked.R
  #'  and Cattle_Hunter_Activity.R scripts. Covariate data formatted with the
  #'  Data_Formatting_3SppX_OccMod.R script.
  #'  ==========================================================================
  
  #'  Clean workspace & load libraries
  rm(list = ls())
  
  library(jagsUI)
  library(abind)
  library(mcmcplots)
  #library(loo)
  library(tidyverse)
  
  
  ####  Source detection history & covariate scripts  ####
  #'  ================================================
  #'  Detection histories come trimmed to desired season length based on unique
  #'  detection events requiring 30 min interval elapse between detections of 
  #'  same species at a given camera
  source("./Scripts/Detection_Histories_for_unmarked.R")
  
  #'  Detection histories representing grazing and hunting activity at each camera 
  #'  during each sampling occasion --> based on summing number of unique detection 
  #'  events, requiring 5 min intervals between detections of cattle or humans, 
  #'  per sampling occasion
  source("./Scripts/Cattle_Hunter_Activity.R")
  
  #'  Format covariate data and detection histories for multi-species occupancy 
  #'  models with THREE-SPECIES interactions
  source("./Scripts/Data_Formatting_2SppX_OccMods.R")
  
  
  ####  Bundle detection histories  ####
  #'  ==============================
  #'  Bundle species detections in an array (site x survey x species)
  str(coug_md_graze_DH)
  coug_md_graze_dat <- abind(coug_md_graze_DH, along = 3)
  coug_elk_graze_dat <- abind(coug_elk_graze_DH, along = 3)
  coug_wtd_graze_dat <- abind(coug_wtd_graze_DH, along = 3) 
  coug_moose_graze_dat <- abind(coug_moose_graze_DH, along = 3)
  
  wolf_md_graze_dat <- abind(wolf_md_graze_DH, along = 3)
  wolf_elk_graze_dat <- abind(wolf_elk_graze_DH, along = 3)
  wolf_wtd_graze_dat <- abind(wolf_wtd_graze_DH, along = 3)
  wolf_moose_graze_dat <- abind(wolf_moose_graze_DH, along = 3) 
  
  bear_md_graze_dat <- abind(bear_md_graze_DH, along = 3) 
  bear_elk_graze_dat <- abind(bear_elk_graze_DH, along = 3) 
  bear_wtd_graze_dat <- abind(bear_wtd_graze_DH, along = 3) 
  bear_moose_graze_dat <- abind(bear_moose_graze_DH, along = 3) 
  
  bob_md_graze_dat <- abind(bob_md_graze_DH, along = 3)
  bob_wtd_graze_dat <- abind(bob_wtd_graze_DH, along = 3)
  
  coy_md_graze_dat <- abind(coy_md_graze_DH, along = 3)
  coy_wtd_graze_dat <- abind(coy_wtd_graze_DH, along = 3)
  
  #'  List lists of detection histories for faster formatting in functions below
  megaList_graze <- list(coug_md_graze_dat, coug_elk_graze_dat, coug_wtd_graze_dat, coug_moose_graze_dat,
                         wolf_md_graze_dat, wolf_elk_graze_dat, wolf_wtd_graze_dat, wolf_moose_graze_dat,
                         bear_md_graze_dat, bear_elk_graze_dat, bear_wtd_graze_dat, bear_moose_graze_dat,
                         bob_md_graze_dat, bob_wtd_graze_dat, coy_md_graze_dat, coy_wtd_graze_dat)
  
  
  ####  Identify survey dimensions  ####
  #'  ==============================
  nsites_graze <- dim(coug_md_graze_dat)[1]
  nsurveys_graze <- dim(coug_md_graze_dat)[2]
  nspp_graze <- dim(coug_md_graze_dat)[3]
  #'  Number of possible community states (species interactions)
  ncat_graze <- 2^nspp_graze
  
  
  ####  Format detection histories for JAGS  ####
  #'  =======================================
  #'  Function to convert species-specific detection histories to multi-species 
  #'  detection array to be site x survey
  detection_array <- function(DH_list) {
    ycat <- apply(DH_list, c(1,2), paste, collapse = "")
    ycat[ycat == "00"] <- 1        # Unoccupied (U)
    ycat[ycat == "10"] <- 2        # Only predator (spp1) detected
    ycat[ycat == "01"] <- 3        # Only prey (spp2) detected
    ycat[ycat == "11"] <- 4        # both detected (all)
    ycat[ycat == "NANA"] <- NA     # Not sampled, no data
    #'  Make each column numeric so JAGS knows to handle it as a response variable
    ycat <- apply(ycat, 2, as.numeric)
    #' #'  Drop all row and column names generated by camTrapR
    #' attr(ycat, "dimnames") <- NULL
    return(ycat)
  } 
  #'  Run list of detection histories through function to create detection arrays
  #'  for each model
  det_array_graze <- lapply(megaList_graze, detection_array)
  names(det_array_graze) <- c("coug_md", "coug_elk", "coug_wtd", "coug_moose",
                              "wolf_md", "wolf_elk", "wolf_wtd", "wolf_moose",
                              "bear_md", "bear_elk", "bear_wtd", "bear_moose",
                              "bob_md", "bob_wtd", "coy_md", "coy_wtd")
  
  
  ####  Format covariates for JAGS  ####
  #'  ==============================
  #'  Prepare site-level covariates for detection sub-models
  #'  study area covariate: 0 = NE, 1 = OK
  #'  trail covariate: 0 = trail, 1 = dirt road
  covs_graze <- stations_graze %>%
    mutate(Study_Area = ifelse(Study_Area == "NE", 0, 1),
           Trail = ifelse(Trail == "Trail", 0, 1), 
           Public = as.integer(as.character(Public)))
  table(gTrail <- covs_graze[,"Trail"])
  
  #'  Janky patch to fill NAs for weekly grazing/hunter activity --> MUST find a better solution for this
  cattle_week_cov <- cattle_week_scaled
  cattle_week_cov[is.na(cattle_week_cov)] <- -0.1546436 # fill NAs with scaled value that equals 0
  
  #'  Matrix for first order occupancy (psi): main effects
  psi_cov_graze <- matrix(NA, ncol = 6, nrow = nsites_graze)
  psi_cov_graze[,1] <- 1
  psi_cov_graze[,2] <- covs_graze$Study_Area
  psi_cov_graze[,3] <- covs_graze$Elev
  psi_cov_graze[,4] <- (covs_graze$Elev)^2
  psi_cov_graze[,5] <- covs_graze$PercForest
  psi_cov_graze[,6] <- covs_graze$GrazingActivity
  head(psi_cov_graze)
  
  #'  Matrix for second order occupancy (psi): 2-way interactions
  psi_inxs_cov_graze <- matrix(NA, ncol = 2, nrow = nsites_graze)
  psi_inxs_cov_graze[,1] <- 1
  psi_inxs_cov_graze[,2] <- covs_graze$GrazingActivity
  head(psi_inxs_cov_graze)
  
  #'  Matrix for first order detection (rho): main effects
  rho_cov_graze <- array(NA, dim = c(nsites_graze, nsurveys_graze, 3)) # last digit = number of covariates
  rho_cov_graze[,,1] <- 1
  rho_cov_graze[,1:nsurveys_graze,2] <- gTrail   # repeat site-level covariate for each survey
  rho_cov_graze[,,3] <- cattle_week_cov
  head(rho_cov_graze)
  
  #'  Matrix for second order detection (rho): 2-way interactions
  rho_inxs_cov_graze <- rep(1, nsites_graze)
  
  
  ####  Bundle all data for JAGS  ####
  #'  ============================
  #'  Function to bundle detection and covariate data
  bundle_data <- function(det_array, psi_covs, psi_inxs, rho_covs, rho_inxs, 
                          sites, surveys, psi_1order, psi_2order, rho_1order, 
                          rho_2order, ncats) {
    #'  list all pieces of data together
    bundled <- list(y = det_array, psi_cov = psi_covs, psi_inxs_cov = psi_inxs,
                    rho_cov = rho_covs, rho_inxs_cov = rho_inxs, nsites = sites,
                    nsurveys = surveys, nfirst_order_psi = ncol(psi_1order), 
                    nsecond_order_psi = ncol(psi_2order), 
                    nfirst_order_rho = dim(rho_1order)[3], 
                    nsecond_order_rho = rho_2order, ncat = ncats)
    #'  Summarize to make sure it looks right
    str(bundled)
    return(bundled)
  }
  #'  Run giant list of detection arrays through function to bundle with covariate
  #'  data for each model
  bundled_graze_list <- lapply(det_array_graze, bundle_data, psi_covs = psi_cov_graze, 
                               psi_inxs = psi_inxs_cov_graze, rho_covs = rho_cov_graze, 
                               rho_inxs = rho_inxs_cov_graze, sites = nsites_graze, 
                               surveys = nsurveys_graze, psi_1order = psi_cov_graze, 
                               psi_2order = psi_inxs_cov_graze, rho_1order = rho_cov_graze, 
                               rho_2order = 1, ncats = ncat_graze)
 
  
  ####  Initial values, Monitor params, MCMC settings  ####
  #'  =================================================
  #'  Provide initial values for model
  #'  Get maximum possible state across all potential surveys at a site (site x spp matrix)
  initial_z <- function(bundled_dat) {
    zinit <- apply(bundled_dat, c(1,3), sum, na.rm = TRUE)
    zinit[zinit > 1] <- 1
    return(zinit)
  }
  #'  Run arrays of species-specific detection data through function
  zinit_graze <- lapply(megaList_graze, initial_z)
  
  #'  Convert initial values to a category
  initial_z_categories <- function(zlist) {
    #'  Collapse species-specific detection histories into 8 categories
    zcat <- apply(zlist, 1, paste, collapse = "")
    zcat[zcat == "00"] <- 1
    zcat[zcat == "10"] <- 2
    zcat[zcat == "01"] <- 3
    zcat[zcat == "11"] <- 4
    #'  Make z numeric again
    zcat <- as.numeric(zcat)
    return(zcat)
  }
  #'  Run list of detection histories through function to create detection arrays
  #'  for each model
  zcat_graze <- lapply(zinit_graze, initial_z_categories)
  names(zcat_graze) <- c("coug_md", "coug_elk", "coug_wtd", "coug_moose",
                         "wolf_md", "wolf_elk", "wolf_wtd", "wolf_moose",
                         "bear_md", "bear_elk", "bear_wtd", "bear_moose",
                         "bob_md", "bob_wtd", "coy_md", "coy_wtd")
  
  #'  Parameters monitored
  params <- c("betaSpp1", "betaSpp2", "betaSpp3", "betaSpp12", "betaSpp23", "betaSpp13",
              "alphaSpp1", "alphaSpp2", "alphaSpp3", "mean.psiSpp1", "mean.psiSpp2", "mean.psiSpp3",
              "mean.pSpp1", "mean.pSpp2", "mean.pSpp3", "z")
  
  #'  MCMC settings
  #nc <- 3; ni <- 500; nb <- 300; nt <- 1; na <- 100
  #nc <- 3; ni <- 5000; nb <- 3000; nt <- 2; na <- 1000
  nc <- 3; ni <- 50000; nb <- 30000; nt <- 4; na <- 10000
  
  
  #'  ===================
  ####  RUN JAGS MODELS  ####
  #'  ===================
  #'  For each species combination: 
  #'    1. set initial values with correct detection data, 
  #'    2. source and run each model in JAGS 
  #'    3. visually inspect traceplots
  #'    4. review model summary and any parameters that did not converge well
  #'    5. save results
  
  ####  Cougar-Mule Deer-Cattle  ####
  #'  ---------------------------
  inits <- function(){list(z = zcat_graze[[1]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  coug.md.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.md.mod0$samples)
  which(coug.md.mod0$summary[,"Rhat"] > 1.1) 
  print(coug.md.mod0$summary[1:36, -c(4:6)], 3)
  save(coug.md.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_coug.md.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  coug.md.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.md.mod1$samples)
  which(coug.md.mod1$summary[,"Rhat"] > 1.1) 
  print(coug.md.mod1$summary[1:36, -c(4:6)], 3)
  save(coug.md.mod1, file="./Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_coug.md.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.md.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.md.mod2$samples)
  which(coug.md.mod2$summary[,"Rhat"] > 1.1) 
  print(coug.md.mod2$summary[1:36, -c(4:6)], 3)
  save(coug.md.mod2, file="./Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.)_coug.md.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  coug.md.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.md.mod3$samples)
  which(coug.md.mod3$summary[,"Rhat"] > 1.1) 
  print(coug.md.mod3$summary[1:36, -c(4:6)], 3)
  save(coug.md.mod3, file="./Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_coug.md.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.md.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.md.mod4$samples)
  which(coug.md.mod4$summary[,"Rhat"] > 1.1) 
  print(coug.md.mod4$summary[1:36, -c(4:6)], 3)
  save(coug.md.mod4, file="./Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_coug.md.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  coug.md.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.md.mod5$samples)
  which(coug.md.mod5$summary[,"Rhat"] > 1.1) 
  print(coug.md.mod5$summary[1:36, -c(4:6)], 3)
  save(coug.md.mod5, file="./Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_coug.md.Rdata")
  
  
  ####  Cougar-Elk-Cattle  ####
  #'  ---------------------
  inits <- function(){list(z = zcat_graze[[2]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  coug.elk.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.elk.mod0$samples)
  which(coug.elk.mod0$summary[,"Rhat"] > 1.1) 
  print(coug.elk.mod0$summary[1:36, -c(4:6)], 3)
  save(coug.elk.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_coug.elk.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  coug.elk.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.elk.mod1$samples)
  which(coug.elk.mod1$summary[,"Rhat"] > 1.1) 
  print(coug.elk.mod1$summary[1:36, -c(4:6)], 3)
  save(coug.elk.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_coug.elk.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.elk.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.elk.mod2$samples)
  which(coug.elk.mod2$summary[,"Rhat"] > 1.1) 
  print(coug.elk.mod2$summary[1:36, -c(4:6)], 3)
  save(coug.elk.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_coug.elk.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  coug.elk.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.elk.mod3$samples)
  which(coug.elk.mod3$summary[,"Rhat"] > 1.1) 
  print(coug.elk.mod3$summary[1:36, -c(4:6)], 3)
  save(coug.elk.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_coug.elk.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.elk.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.elk.mod4$samples)
  which(coug.elk.mod4$summary[,"Rhat"] > 1.1) 
  print(coug.elk.mod4$summary[1:36, -c(4:6)], 3)
  save(coug.elk.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_coug.elk.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  coug.elk.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.elk.mod5$samples)
  which(coug.elk.mod5$summary[,"Rhat"] > 1.1) 
  print(coug.elk.mod5$summary[1:36, -c(4:6)], 3)
  save(coug.elk.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_coug.elk.Rdata")  
  
  
  ####  Cougar-White-tailed Deer-Cattle  ####
  #'  -----------------------------------
  inits <- function(){list(z = zcat_graze[[3]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  coug.wtd.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.wtd.mod0$samples)
  which(coug.wtd.mod0$summary[,"Rhat"] > 1.1) 
  print(coug.wtd.mod0$summary[1:36, -c(4:6)], 3)
  save(coug.wtd.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_coug.wtd.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  coug.wtd.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.wtd.mod1$samples)
  which(coug.wtd.mod1$summary[,"Rhat"] > 1.1) 
  print(coug.wtd.mod1$summary[1:36, -c(4:6)], 3)
  save(coug.wtd.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_coug.wtd.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.wtd.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.wtd.mod2$samples)
  which(coug.wtd.mod2$summary[,"Rhat"] > 1.1) 
  print(coug.wtd.mod2$summary[1:36, -c(4:6)], 3)
  save(coug.wtd.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_coug.wtd.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  coug.wtd.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.wtd.mod3$samples)
  which(coug.wtd.mod3$summary[,"Rhat"] > 1.1) 
  print(coug.wtd.mod3$summary[1:36, -c(4:6)], 3)
  save(coug.wtd.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_coug.wtd.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.wtd.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.wtd.mod4$samples)
  which(coug.wtd.mod4$summary[,"Rhat"] > 1.1) 
  print(coug.wtd.mod4$summary[1:36, -c(4:6)], 3)
  save(coug.wtd.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_coug.wtd.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  coug.wtd.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.wtd.mod5$samples)
  which(coug.wtd.mod5$summary[,"Rhat"] > 1.1) 
  print(coug.wtd.mod5$summary[1:36, -c(4:6)], 3)
  save(coug.wtd.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_coug.wtd.Rdata")
  
  
  ####  Cougar-Moose-Cattle  ####
  #'  -----------------------------------
  inits <- function(){list(z = zcat_graze[[3]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  coug.moose.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.moose.mod0$samples)
  which(coug.moose.mod0$summary[,"Rhat"] > 1.1) 
  print(coug.moose.mod0$summary[1:36, -c(4:6)], 3)
  save(coug.moose.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_coug.moose.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  coug.moose.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.moose.mod1$samples)
  which(coug.moose.mod1$summary[,"Rhat"] > 1.1) 
  print(coug.moose.mod1$summary[1:36, -c(4:6)], 3)
  save(coug.moose.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_coug.moose.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.moose.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.moose.mod2$samples)
  which(coug.moose.mod2$summary[,"Rhat"] > 1.1) 
  print(coug.moose.mod2$summary[1:36, -c(4:6)], 3)
  save(coug.moose.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_coug.moose.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  coug.moose.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.moose.mod3$samples)
  which(coug.moose.mod3$summary[,"Rhat"] > 1.1) 
  print(coug.moose.mod3$summary[1:36, -c(4:6)], 3)
  save(coug.moose.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_coug.moose.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  coug.moose.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.moose.mod4$samples)
  which(coug.moose.mod4$summary[,"Rhat"] > 1.1) 
  print(coug.moose.mod4$summary[1:36, -c(4:6)], 3)
  save(coug.moose.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_coug.moose.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  coug.moose.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coug.moose.mod5$samples)
  which(coug.moose.mod5$summary[,"Rhat"] > 1.1) 
  print(coug.moose.mod5$summary[1:36, -c(4:6)], 3)
  save(coug.moose.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_coug.moose.Rdata")
  
  
  ####  Wolf-Mule Deer-Cattle  ####
  #'  -------------------------
  inits <- function(){list(z = zcat_graze[[5]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  wolf.md.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.md.mod0$samples)
  which(wolf.md.mod0$summary[,"Rhat"] > 1.1) 
  print(wolf.md.mod0$summary[1:36, -c(4:6)], 3)
  save(wolf.md.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_wolf.md.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  wolf.md.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.md.mod1$samples)
  which(wolf.md.mod1$summary[,"Rhat"] > 1.1) 
  print(wolf.md.mod1$summary[1:36, -c(4:6)], 3)
  save(wolf.md.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_wolf.md.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.md.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.md.mod2$samples)
  which(wolf.md.mod2$summary[,"Rhat"] > 1.1) 
  print(wolf.md.mod2$summary[1:36, -c(4:6)], 3)
  save(wolf.md.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_wolf.md.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  wolf.md.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.md.mod3$samples)
  which(wolf.md.mod3$summary[,"Rhat"] > 1.1) 
  print(wolf.md.mod3$summary[1:36, -c(4:6)], 3)
  save(wolf.md.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_wolf.md.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.md.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.md.mod4$samples)
  which(wolf.md.mod4$summary[,"Rhat"] > 1.1) 
  print(wolf.md.mod4$summary[1:36, -c(4:6)], 3)
  save(wolf.md.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_wolf.md.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  wolf.md.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.md.mod5$samples)
  which(wolf.md.mod5$summary[,"Rhat"] > 1.1) 
  print(wolf.md.mod5$summary[1:36, -c(4:6)], 3)
  save(wolf.md.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_wolf.md.Rdata")
  
  
  ####  Wolf-Elk-Cattle  ####
  #'  -------------------
  inits <- function(){list(z = zcat_graze[[6]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  wolf.elk.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.elk.mod0$samples)
  which(wolf.elk.mod0$summary[,"Rhat"] > 1.1) 
  print(wolf.elk.mod0$summary[1:36, -c(4:6)], 3)
  save(wolf.elk.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_wolf.elk.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  wolf.elk.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.elk.mod1$samples)
  which(wolf.elk.mod1$summary[,"Rhat"] > 1.1) 
  print(wolf.elk.mod1$summary[1:36, -c(4:6)], 3)
  save(wolf.elk.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_wolf.elk.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.elk.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.elk.mod2$samples)
  which(wolf.elk.mod2$summary[,"Rhat"] > 1.1) 
  print(wolf.elk.mod2$summary[1:36, -c(4:6)], 3)
  save(wolf.elk.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_wolf.elk.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  wolf.elk.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.elk.mod3$samples)
  which(wolf.elk.mod3$summary[,"Rhat"] > 1.1) 
  print(wolf.elk.mod3$summary[1:36, -c(4:6)], 3)
  save(wolf.elk.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_wolf.elk.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.elk.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.elk.mod4$samples)
  which(wolf.elk.mod4$summary[,"Rhat"] > 1.1) 
  print(wolf.elk.mod4$summary[1:36, -c(4:6)], 3)
  save(wolf.elk.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_wolf.elk.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  wolf.elk.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.elk.mod5$samples)
  which(wolf.elk.mod5$summary[,"Rhat"] > 1.1) 
  print(wolf.elk.mod5$summary[1:36, -c(4:6)], 3)
  save(wolf.elk.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_wolf.elk.Rdata")  
  
  
  
  ####  Wolf-White-tailed Deer-Cattle  ####
  #'  ---------------------------------
  inits <- function(){list(z = zcat_graze[[7]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  wolf.wtd.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.wtd.mod0$samples)
  which(wolf.wtd.mod0$summary[,"Rhat"] > 1.1) 
  print(wolf.wtd.mod0$summary[1:36, -c(4:6)], 3)
  save(wolf.wtd.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_wolf.wtd.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  wolf.wtd.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.wtd.mod1$samples)
  which(wolf.wtd.mod1$summary[,"Rhat"] > 1.1) 
  print(wolf.wtd.mod1$summary[1:36, -c(4:6)], 3)
  save(wolf.wtd.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_wolf.wtd.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.wtd.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.wtd.mod2$samples)
  which(wolf.wtd.mod2$summary[,"Rhat"] > 1.1) 
  print(wolf.wtd.mod2$summary[1:36, -c(4:6)], 3)
  save(wolf.wtd.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_wolf.wtd.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  wolf.wtd.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.wtd.mod3$samples)
  which(wolf.wtd.mod3$summary[,"Rhat"] > 1.1) 
  print(wolf.wtd.mod3$summary[1:36, -c(4:6)], 3)
  save(wolf.wtd.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_wolf.wtd.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.wtd.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.wtd.mod4$samples)
  which(wolf.wtd.mod4$summary[,"Rhat"] > 1.1) 
  print(wolf.wtd.mod4$summary[1:36, -c(4:6)], 3)
  save(wolf.wtd.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_wolf.wtd.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  wolf.wtd.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.wtd.mod5$samples)
  which(wolf.wtd.mod5$summary[,"Rhat"] > 1.1) 
  print(wolf.wtd.mod5$summary[1:36, -c(4:6)], 3)
  save(wolf.wtd.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_wolf.wtd.Rdata")
  
  
  
  
  ####  Wolf-Moose-Cattle  ####
  #'  ---------------------
  inits <- function(){list(z = zcat_graze[[8]])}
  
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  wolf.moose.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.moose.mod0$samples)
  which(wolf.moose.mod0$summary[,"Rhat"] > 1.1) 
  print(wolf.moose.mod0$summary[1:36, -c(4:6)], 3)
  save(wolf.moose.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_wolf.moose.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  wolf.moose.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.moose.mod1$samples)
  which(wolf.moose.mod1$summary[,"Rhat"] > 1.1) 
  print(wolf.moose.mod1$summary[1:36, -c(4:6)], 3)
  save(wolf.moose.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_wolf.moose.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.moose.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.moose.mod2$samples)
  which(wolf.moose.mod2$summary[,"Rhat"] > 1.1) 
  print(wolf.moose.mod2$summary[1:36, -c(4:6)], 3)
  save(wolf.moose.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_wolf.moose.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  wolf.moose.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.moose.mod3$samples)
  which(wolf.moose.mod3$summary[,"Rhat"] > 1.1) 
  print(wolf.moose.mod3$summary[1:36, -c(4:6)], 3)
  save(wolf.moose.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_wolf.moose.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  wolf.moose.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.moose.mod4$samples)
  which(wolf.moose.mod4$summary[,"Rhat"] > 1.1) 
  print(wolf.moose.mod4$summary[1:36, -c(4:6)], 3)
  save(wolf.moose.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_wolf.moose.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  wolf.moose.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(wolf.moose.mod5$samples)
  which(wolf.moose.mod5$summary[,"Rhat"] > 1.1) 
  print(wolf.moose.mod5$summary[1:36, -c(4:6)], 3)
  save(wolf.moose.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_wolf.moose.Rdata")
  
  
  
  
  ####  Bear-Mule Deer-Cattle  ####
  #'  -------------------------
  inits <- function(){list(z = zcat_graze[[9]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  bear.md.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.md.mod0$samples)
  which(bear.md.mod0$summary[,"Rhat"] > 1.1) 
  print(bear.md.mod0$summary[1:36, -c(4:6)], 3)
  save(bear.md.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_bear.md.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  bear.md.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.md.mod1$samples)
  which(bear.md.mod1$summary[,"Rhat"] > 1.1) 
  print(bear.md.mod1$summary[1:36, -c(4:6)], 3)
  save(bear.md.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_bear.md.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.md.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.md.mod2$samples)
  which(bear.md.mod2$summary[,"Rhat"] > 1.1) 
  print(bear.md.mod2$summary[1:36, -c(4:6)], 3)
  save(bear.md.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_bear.md.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  bear.md.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.md.mod3$samples)
  which(bear.md.mod3$summary[,"Rhat"] > 1.1) 
  print(bear.md.mod3$summary[1:36, -c(4:6)], 3)
  save(bear.md.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_bear.md.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.md.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.md.mod4$samples)
  which(bear.md.mod4$summary[,"Rhat"] > 1.1) 
  print(bear.md.mod4$summary[1:36, -c(4:6)], 3)
  save(bear.md.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_bear.md.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  bear.md.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.md.mod5$samples)
  which(bear.md.mod5$summary[,"Rhat"] > 1.1) 
  print(bear.md.mod5$summary[1:36, -c(4:6)], 3)
  save(bear.md.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_bear.md.Rdata")
  
  
  
  ####  Bear-Elk-Cattle  ####
  #'  -------------------
  inits <- function(){list(z = zcat_graze[[10]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  bear.elk.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.elk.mod0$samples)
  which(bear.elk.mod0$summary[,"Rhat"] > 1.1) 
  print(bear.elk.mod0$summary[1:36, -c(4:6)], 3)
  save(bear.elk.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_bear.elk.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  bear.elk.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.elk.mod1$samples)
  which(bear.elk.mod1$summary[,"Rhat"] > 1.1) 
  print(bear.elk.mod1$summary[1:36, -c(4:6)], 3)
  save(bear.elk.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_bear.elk.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.elk.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.elk.mod2$samples)
  which(bear.elk.mod2$summary[,"Rhat"] > 1.1) 
  print(bear.elk.mod2$summary[1:36, -c(4:6)], 3)
  save(bear.elk.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_bear.elk.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  bear.elk.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.elk.mod3$samples)
  which(bear.elk.mod3$summary[,"Rhat"] > 1.1) 
  print(bear.elk.mod3$summary[1:36, -c(4:6)], 3)
  save(bear.elk.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_bear.elk.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.elk.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.elk.mod4$samples)
  which(bear.elk.mod4$summary[,"Rhat"] > 1.1) 
  print(bear.elk.mod4$summary[1:36, -c(4:6)], 3)
  save(bear.elk.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_bear.elk.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  bear.elk.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.elk.mod5$samples)
  which(bear.elk.mod5$summary[,"Rhat"] > 1.1) 
  print(bear.elk.mod5$summary[1:36, -c(4:6)], 3)
  save(bear.elk.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_bear.elk.Rdata")  
  
  
  
  
  ####  Bear-White-tailed Deer-Cattle  ####
  #'  ---------------------------------
  inits <- function(){list(z = zcat_graze[[11]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  bear.wtd.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.wtd.mod0$samples)
  which(bear.wtd.mod0$summary[,"Rhat"] > 1.1) 
  print(bear.wtd.mod0$summary[1:36, -c(4:6)], 3)
  save(bear.wtd.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_bear.wtd.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  bear.wtd.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.wtd.mod1$samples)
  which(bear.wtd.mod1$summary[,"Rhat"] > 1.1) 
  print(bear.wtd.mod1$summary[1:36, -c(4:6)], 3)
  save(bear.wtd.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_bear.wtd.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.wtd.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.wtd.mod2$samples)
  which(bear.wtd.mod2$summary[,"Rhat"] > 1.1) 
  print(bear.wtd.mod2$summary[1:36, -c(4:6)], 3)
  save(bear.wtd.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_bear.wtd.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  bear.wtd.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.wtd.mod3$samples)
  which(bear.wtd.mod3$summary[,"Rhat"] > 1.1) 
  print(bear.wtd.mod3$summary[1:36, -c(4:6)], 3)
  save(bear.wtd.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_bear.wtd.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.wtd.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.wtd.mod4$samples)
  which(bear.wtd.mod4$summary[,"Rhat"] > 1.1) 
  print(bear.wtd.mod4$summary[1:36, -c(4:6)], 3)
  save(bear.wtd.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_bear.wtd.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  bear.wtd.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.wtd.mod5$samples)
  which(bear.wtd.mod5$summary[,"Rhat"] > 1.1) 
  print(bear.wtd.mod5$summary[1:36, -c(4:6)], 3)
  save(bear.wtd.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_bear.wtd.Rdata")
  
  
  
  
  ####  Bear-Moose-Cattle  ####
  #'  ---------------------
  inits <- function(){list(z = zcat_graze[[12]])}
  
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  bear.moose.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.moose.mod0$samples)
  which(bear.moose.mod0$summary[,"Rhat"] > 1.1) 
  print(bear.moose.mod0$summary[1:36, -c(4:6)], 3)
  save(bear.moose.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_bear.moose.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  bear.moose.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.moose.mod1$samples)
  which(bear.moose.mod1$summary[,"Rhat"] > 1.1) 
  print(bear.moose.mod1$summary[1:36, -c(4:6)], 3)
  save(bear.moose.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_bear.moose.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.moose.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.moose.mod2$samples)
  which(bear.moose.mod2$summary[,"Rhat"] > 1.1) 
  print(bear.moose.mod2$summary[1:36, -c(4:6)], 3)
  save(bear.moose.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_bear.moose.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  bear.moose.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.moose.mod3$samples)
  which(bear.moose.mod3$summary[,"Rhat"] > 1.1) 
  print(bear.moose.mod3$summary[1:36, -c(4:6)], 3)
  save(bear.moose.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_bear.moose.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  bear.moose.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.moose.mod4$samples)
  which(bear.moose.mod4$summary[,"Rhat"] > 1.1) 
  print(bear.moose.mod4$summary[1:36, -c(4:6)], 3)
  save(bear.moose.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_bear.moose.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  bear.moose.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                          n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bear.moose.mod5$samples)
  which(bear.moose.mod5$summary[,"Rhat"] > 1.1) 
  print(bear.moose.mod5$summary[1:36, -c(4:6)], 3)
  save(bear.moose.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_bear.moose.Rdata")
  
  
  
  
  ####  Bobcat-Mule Deer-Cattle  ####
  #'  ---------------------------
  inits <- function(){list(z = zcat_graze[[13]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  bob.md.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.md.mod0$samples)
  which(bob.md.mod0$summary[,"Rhat"] > 1.1) 
  print(bob.md.mod0$summary[1:36, -c(4:6)], 3)
  save(bob.md.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_bob.md.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  bob.md.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.md.mod1$samples)
  which(bob.md.mod1$summary[,"Rhat"] > 1.1) 
  print(bob.md.mod1$summary[1:36, -c(4:6)], 3)
  save(bob.md.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_bob.md.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  bob.md.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.md.mod2$samples)
  which(bob.md.mod2$summary[,"Rhat"] > 1.1) 
  print(bob.md.mod2$summary[1:36, -c(4:6)], 3)
  save(bob.md.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_bob.md.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  bob.md.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.md.mod3$samples)
  which(bob.md.mod3$summary[,"Rhat"] > 1.1) 
  print(bob.md.mod3$summary[1:36, -c(4:6)], 3)
  save(bob.md.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_bob.md.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  bob.md.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.md.mod4$samples)
  which(bob.md.mod4$summary[,"Rhat"] > 1.1) 
  print(bob.md.mod4$summary[1:36, -c(4:6)], 3)
  save(bob.md.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_bob.md.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  bob.md.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.md.mod5$samples)
  which(bob.md.mod5$summary[,"Rhat"] > 1.1) 
  print(bob.md.mod5$summary[1:36, -c(4:6)], 3)
  save(bob.md.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_bob.md.Rdata")
  
  
  
  ####  Bobcat-White-tailed Deer-Cattle  ####
  #'  -----------------------------------
  inits <- function(){list(z = zcat_graze[[14]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  bob.wtd.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.wtd.mod0$samples)
  which(bob.wtd.mod0$summary[,"Rhat"] > 1.1) 
  print(bob.wtd.mod0$summary[1:36, -c(4:6)], 3)
  save(bob.wtd.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_bob.wtd.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  bob.wtd.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.wtd.mod1$samples)
  which(bob.wtd.mod1$summary[,"Rhat"] > 1.1) 
  print(bob.wtd.mod1$summary[1:36, -c(4:6)], 3)
  save(bob.wtd.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_bob.wtd.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  bob.wtd.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.wtd.mod2$samples)
  which(bob.wtd.mod2$summary[,"Rhat"] > 1.1) 
  print(bob.wtd.mod2$summary[1:36, -c(4:6)], 3)
  save(bob.wtd.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_bob.wtd.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  bob.wtd.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.wtd.mod3$samples)
  which(bob.wtd.mod3$summary[,"Rhat"] > 1.1) 
  print(bob.wtd.mod3$summary[1:36, -c(4:6)], 3)
  save(bob.wtd.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_bob.wtd.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  bob.wtd.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.wtd.mod4$samples)
  which(bob.wtd.mod4$summary[,"Rhat"] > 1.1) 
  print(bob.wtd.mod4$summary[1:36, -c(4:6)], 3)
  save(bob.wtd.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_bob.wtd.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  bob.wtd.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(bob.wtd.mod5$samples)
  which(bob.wtd.mod5$summary[,"Rhat"] > 1.1) 
  print(bob.wtd.mod5$summary[1:36, -c(4:6)], 3)
  save(bob.wtd.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_bob.wtd.Rdata")
  
  
  
  
  ####  Coyote-Mule Deer-Cattle  ####
  #'  ---------------------------
  inits <- function(){list(z = zcat_graze[[15]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  coy.md.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.md.mod0$samples)
  which(coy.md.mod0$summary[,"Rhat"] > 1.1) 
  print(coy.md.mod0$summary[1:36, -c(4:6)], 3)
  save(coy.md.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_coy.md.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  coy.md.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.md.mod1$samples)
  which(coy.md.mod1$summary[,"Rhat"] > 1.1) 
  print(coy.md.mod1$summary[1:36, -c(4:6)], 3)
  save(coy.md.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_coy.md.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  coy.md.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.md.mod2$samples)
  which(coy.md.mod2$summary[,"Rhat"] > 1.1) 
  print(coy.md.mod2$summary[1:36, -c(4:6)], 3)
  save(coy.md.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_coy.md.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  coy.md.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.md.mod3$samples)
  which(coy.md.mod3$summary[,"Rhat"] > 1.1) 
  print(coy.md.mod3$summary[1:36, -c(4:6)], 3)
  save(coy.md.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_coy.md.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  coy.md.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.md.mod4$samples)
  which(coy.md.mod4$summary[,"Rhat"] > 1.1) 
  print(coy.md.mod4$summary[1:36, -c(4:6)], 3)
  save(coy.md.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_coy.md.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  coy.md.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                       n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.md.mod5$samples)
  which(coy.md.mod5$summary[,"Rhat"] > 1.1) 
  print(coy.md.mod5$summary[1:36, -c(4:6)], 3)
  save(coy.md.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_coy.md.Rdata")
  
  
  
  ####  Coyote-White-tailed Deer-Cattle  ####
  #'  -----------------------------------
  inits <- function(){list(z = zcat_graze[[16]])}
  
  #'  Null model
  source("./Scripts/JAGS_models_2spp/Grazing_p(.)_psi(.).R")
  start.time <- Sys.time()
  coy.wtd.mod0 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.wtd.mod0$samples)
  which(coy.wtd.mod0$summary[,"Rhat"] > 1.1) 
  print(coy.wtd.mod0$summary[1:36, -c(4:6)], 3)
  save(coy.wtd.mod0, file="./Outputs/JAGS_models_2spp/Grazing_p(.)_psi(.)_coy.wtd.Rdata")
  
  #'  Conditional occupancy model, no inxs - baseline habitat model
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).R")
  start.time <- Sys.time()
  coy.wtd.mod1 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.wtd.mod1$samples)
  which(coy.wtd.mod1$summary[,"Rhat"] > 1.1) 
  print(coy.wtd.mod1$summary[1:36, -c(4:6)], 3)
  save(coy.wtd.mod1, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_coy.wtd.Rdata")
  
  #'  Pairwise interaction with habitat on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).R")
  start.time <- Sys.time()
  coy.wtd.mod2 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail)_psi(habitat)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.wtd.mod2$samples)
  which(coy.wtd.mod2$summary[,"Rhat"] > 1.1) 
  print(coy.wtd.mod2$summary[1:36, -c(4:6)], 3)
  save(coy.wtd.mod2, file="./Outputs/JAGS_models/Grazing_p(trail)_psi(habitat)_inxpsi(.)_coy.wtd.Rdata")
  
  #'  Conditional occupancy model with cattle activity, no inxs
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).R")
  start.time <- Sys.time()
  coy.wtd.mod3 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.wtd.mod3$samples)
  which(coy.wtd.mod3$summary[,"Rhat"] > 1.1) 
  print(coy.wtd.mod3$summary[1:36, -c(4:6)], 3)
  save(coy.wtd.mod3, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_coy.wtd.Rdata")
  
  #'  Pairwise interaction with habitat and cattle activity on conditional occupancy
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).R")
  start.time <- Sys.time()
  coy.wtd.mod4 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.wtd.mod4$samples)
  which(coy.wtd.mod4$summary[,"Rhat"] > 1.1) 
  print(coy.wtd.mod4$summary[1:36, -c(4:6)], 3)
  save(coy.wtd.mod4, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(.)_coy.wtd.Rdata")
  
  #'  Pairwise interaction with cattle activity
  source("./Scripts/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).R")
  start.time <- Sys.time()
  coy.wtd.mod5 <- jags(bundled_graze_list[[1]], inits = inits, params, './Outputs/JAGS_models_2spp/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle).txt',
                        n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, n.adapt = na, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  mcmcplot(coy.wtd.mod5$samples)
  which(coy.wtd.mod5$summary[,"Rhat"] > 1.1) 
  print(coy.wtd.mod5$summary[1:36, -c(4:6)], 3)
  save(coy.wtd.mod5, file="./Outputs/JAGS_models/Grazing_p(trail_cattle)_psi(habitat_cattle)_inxpsi(cattle)_coy.wtd.Rdata")
  
  
  
  

  
  # QUESTIONS/TO-DO List
  # 1. did I correctly set the 3-way interaction parameter to 0 I don't estimate it?
  # 2. what to do about these missing survey-level covs

  
  
  
  
