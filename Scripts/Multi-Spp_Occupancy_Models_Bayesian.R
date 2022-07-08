  #'  =====================================================
  #'  Bayesian Multi-Species Co-Occurrence Occupancy Models 
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  July 2022
  #'  =====================================================
  #'  Script to source data formatting scripts & run multi-species, single-season
  #'  occupancy models for deer, elk, moose, black bears, cougars, wolves, coyotes, 
  #'  and bobcats during the grazing season 2018-2020 (July - Sept) and hunting
  #'  season 2018-2020 (Oct - Nov), respectively. Grazing season co-occurrence 
  #'  models include 13 7-day sampling occasions comprising the peak of livestock
  #'  activity detected on camera. Hunting season co-occurrence models include
  #'  8 7-day sampling occasions comprising the two general rifle hunting seasons
  #'  in eastern Washington. Co-occurrence models test whether predator-prey 
  #'  co-occurrence is not independent and whether their occurrence, co-occurrence,
  #'  and detection are influenced by livestock and hunting activity at camera sites.
  #'  
  #'  Encounter histories are generated with the CameratTrap_DetectionHistories_for_unmarked.R
  #'  and Cattle_Hunter_Activity.R scripts. Covariate data formatted with the
  #'  Data_Formatting_3SppX_OccMod.R script.
  #'  ============================================
  
  #'  Clean workspace & load libraries
  rm(list = ls())
  
  # library(unmarked)
  library(jagsUI)
  library(abind)
  # library(MuMIn)
  # library(condformat)
  library(tidyverse)
  
  #'  Source scripts that generates detection histories
  #'  Detection histories come trimmed to desired season length based on unique
  #'  detection events requiring 30 min interval elapse between detections of 
  #'  same species at a given camera
  source("./Scripts/Detection_Histories_for_unmarked.R")
  
  #'  Detection histories representing grazing and hunting activity at each camera 
  #'  during each sampling occasion --> based on summing number of unique detection 
  #'  events, requiring 5 min intervals between detections of cattle or humans, 
  #'  per sampling occasion
  source("./Scripts/Cattle_Hunter_Activity.R")
  
  #'  Formats covariate data and detection histories for multi-species occupancy 
  #'  models in unmarked with THREE-SPECIES interactions
  source("./Scripts/Data_Formatting_3SppX_OccMods.R")
  
  #'  Bundle species detections in an array (site x survey x species)
  #'  Grazing season data
  str(coug_md_cattle_graze_DH)
  coug_md_cattle_graze_dat <- abind(coug_md_cattle_graze_DH, along = 3)
  coug_elk_cattle_graze_dat <- abind(coug_elk_cattle_graze_DH, along = 3)
  coug_wtd_cattle_graze_dat <- abind(coug_wtd_cattle_graze_DH, along = 3) 
  coug_moose_cattle_graze_dat <- abind(coug_moose_cattle_graze_DH, along = 3)
  
  wolf_md_cattle_graze_dat <- abind(wolf_md_cattle_graze_DH, along = 3)
  wolf_elk_cattle_graze_dat <- abind(wolf_elk_cattle_graze_DH, along = 3)
  wolf_wtd_cattle_graze_dat <- abind(wolf_wtd_cattle_graze_DH, along = 3)
  wolf_moose_cattle_graze_dat <- abind(wolf_moose_cattle_graze_DH, along = 3) 
  
  bear_md_cattle_graze_dat <- abind(bear_md_cattle_graze_DH, along = 3) 
  bear_elk_cattle_graze_dat <- abind(bear_elk_cattle_graze_DH, along = 3) 
  bear_wtd_cattle_graze_dat <- abind(bear_wtd_cattle_graze_DH, along = 3) 
  bear_moose_cattle_graze_dat <- abind(bear_moose_cattle_graze_DH, along = 3) 
  
  bob_md_cattle_graze_dat <- abind(bob_md_cattle_graze_DH, along = 3)
  bob_wtd_cattle_graze_dat <- abind(bob_wtd_cattle_graze_DH, along = 3)
  
  coy_md_cattle_graze_dat <- abind(coy_md_cattle_graze_DH, along = 3)
  coy_wtd_cattle_graze_dat <- abind(coy_wtd_cattle_graze_DH, along = 3)
  
  #'  List lists of detection histories for faster formatting in functions below
  megaList_graze <- list(coug_md_cattle_graze_dat, coug_elk_cattle_graze_dat, coug_wtd_cattle_graze_dat, coug_moose_cattle_graze_dat,
                         wolf_md_cattle_graze_dat, wolf_elk_cattle_graze_dat, wolf_wtd_cattle_graze_dat, wolf_moose_cattle_graze_dat,
                         bear_md_cattle_graze_dat, bear_elk_cattle_graze_dat, bear_wtd_cattle_graze_dat, bear_moose_cattle_graze_dat,
                         bob_md_cattle_graze_dat, bob_wtd_cattle_graze_dat, coy_md_cattle_graze_dat, coy_wtd_cattle_graze_dat)
  
  #'  Hunting season data
  str(coug_md_hunter_hunt_DH)
  coug_md_hunter_hunt_dat <- abind(coug_md_hunter_hunt_DH, along = 3)
  coug_elk_hunter_hunt_dat <- abind(coug_elk_hunter_hunt_DH, along = 3)
  coug_wtd_hunter_hunt_dat <- abind(coug_wtd_hunter_hunt_DH, along = 3)
  coug_moose_hunter_hunt_dat <- abind(coug_moose_hunter_hunt_DH, along = 3)

  wolf_md_hunter_hunt_dat <- abind(wolf_md_hunter_hunt_DH, along = 3)
  wolf_elk_hunter_hunt_dat <- abind(wolf_elk_hunter_hunt_DH, along = 3)
  wolf_wtd_hunter_hunt_dat <- abind(wolf_wtd_hunter_hunt_DH, along = 3) 
  wolf_moose_hunter_hunt_dat <- abind(wolf_moose_hunter_hunt_DH, along = 3)

  bear_md_hunter_hunt_dat <- abind(bear_md_hunter_hunt_DH, along = 3) 
  bear_elk_hunter_hunt_dat <- abind(bear_elk_hunter_hunt_DH, along = 3)
  bear_wtd_hunter_hunt_dat <- abind(bear_wtd_hunter_hunt_DH, along = 3) 
  bear_moose_hunter_hunt_dat <- abind(bear_moose_hunter_hunt_DH, along = 3) 
  
  bob_md_hunter_hunt_dat <- abind(bob_md_hunter_hunt_DH, along = 3)
  bob_wtd_hunter_hunt_dat <- abind(bob_wtd_hunter_hunt_DH, along = 3)
  
  coy_md_hunter_hunt_dat <- abind(coy_md_hunter_hunt_DH, along = 3)
  coy_wtd_hunter_hunt_dat <- abind(coy_wtd_hunter_hunt_DH, along = 3)
  
  #'  List lists of detection histories for faster formatting in functions below
  megaList_hunt <- list(coug_md_hunter_hunt_dat, coug_elk_hunter_hunt_dat, coug_wtd_hunter_hunt_dat, coug_moose_hunter_hunt_dat,
                        wolf_md_hunter_hunt_dat, wolf_elk_hunter_hunt_dat, wolf_wtd_hunter_hunt_dat, wolf_moose_hunter_hunt_dat,
                        bear_md_hunter_hunt_dat, bear_elk_hunter_hunt_dat, bear_wtd_hunter_hunt_dat, bear_moose_hunter_hunt_dat,
                        bob_md_hunter_hunt_dat, bob_wtd_hunter_hunt_dat, coy_md_hunter_hunt_dat, coy_wtd_hunter_hunt_dat)
  
  #'  Identify survey dimensions
  nsites_graze <- dim(coug_md_cattle_graze_dat)[1]
  nsurvey_graze <- dim(coug_md_cattle_graze_dat)[2]
  nspp_graze <- dim(coug_md_cattle_graze_dat)[3]
  #'  Number of possible community states (species interactions)
  ncat_graze <- 2^nspp_graze
  
  nsites_hunt <- dim(coug_md_hunter_hunt_dat)[1]
  nsurvey_hunt <- dim(coug_md_hunter_hunt_dat)[2]
  nspp_hunt <- dim(coug_md_hunter_hunt_dat)[3]
  #'  Number of possible community states (species interactions)
  ncat_hunt <- 2^nspp_hunt
  
  #'  Prepare site-level covariates for detection sub-model
  #'  trail covariate: 1 = trail, 2 = dirt road
  #'  public covariate: 1 = private, 2 = public
  table(gTrail <- stations_graze[,"Trail"])
  table(hTrail <- stations_hunt[,"Trail"])
  table(Public <- stations_hunt[,"Public"])
  
  #'  Function to convert species-specific detection histories to multi-species 
  #'  detection array to be site x survey
  detection_array <- function(DH_list) {
    ycat <- apply(DH_list, c(1,2), paste, collapse = "")
    ycat[ycat == "000"] <- 1        # Unoccupied (U)
    ycat[ycat == "100"] <- 2        # Only predator (spp1) detected
    ycat[ycat == "010"] <- 3        # Only prey (spp2) detected
    ycat[ycat == "001"] <- 4        # Only anthropogenic acticity (spp3) detected
    ycat[ycat == "110"] <- 5        # spp1 & spp2 detected
    ycat[ycat == "011"] <- 6        # spp2 & spp3 detected
    ycat[ycat == "101"] <- 7        # spp1 & spp3 detected
    ycat[ycat == "111"] <- 8        # All three detected (all)
    ycat[ycat == "NANANA"] <- NA    # Not sampled, no data
    #'  Make each column numeric so JAGS knows to handle it as a response variable
    ycat <- apply(ycat, 2, as.numeric)
    return(ycat)
  } 
  #'  Run list of detection histories through function to create detection arrays
  #'  for each model
  det_array_graze <- lapply(megaList_graze, detection_array)
  names(det_array_graze) <- c("coug_md_cattle", "coug_elk_cattle", "coug_wtd_cattle", "coug_moose_cattle",
                              "wolf_md_cattle", "wolf_elk_cattle", "wolf_wtd_cattle", "wolf_moose_cattle",
                              "bear_md_cattle", "bear_elk_cattle", "bear_wtd_cattle", "bear_moose_cattle",
                              "bob_md_cattle", "bob_wtd_cattle", "coy_md_cattle", "coy_wtd_cattle")
  det_array_hunt <- lapply(megaList_hunt, detection_array)
  names(det_array_hunt) <- c("coug_md_hunter", "coug_elk_hunter", "coug_wtd_hunter", "coug_moose_hunter",
                              "wolf_md_hunter", "wolf_elk_hunter", "wolf_wtd_hunter", "wolf_moose_hunter",
                              "bear_md_hunter", "bear_elk_hunter", "bear_wtd_hunter", "bear_moose_hunter",
                              "bob_md_hunter", "bob_wtd_hunter", "coy_md_hunter", "coy_wtd_hunter")
  
  ####  Covariate matrices for different parts of model  ####
  #'  ===================================================
  #'  Matrix for first order occupancy (psi): main effects
  psi_cov_graze <- matrix(NA, ncol = 5, nrow = nsites_graze)
  psi_cov_graze[,1] <- 1
  psi_cov_graze[,2] <- stations_graze$Study_Area
  psi_cov_graze[,3] <- stations_graze$Elev
  psi_cov_graze[,4] <- (stations_graze$Elev)^2
  psi_cov_graze[,5] <- stations_graze$PercForest
  
  psi_cov_hunt <- matrix(NA, ncol = 5, nrow = nsites_hunt)
  psi_cov_hunt[,1] <- 1
  psi_cov_hunt[,2] <- stations_hunt$Study_Area
  psi_cov_hunt[,3] <- stations_hunt$Elev
  psi_cov_hunt[,4] <- (stations_hunt$Elev)^2
  psi_cov_hunt[,5] <- stations_hunt$PercForest
  
  #'  Matrix for second order occupancy (psi): 2-way interactions
  psi_inxs_covs_graze <- matrix(NA, ncol = 2, nrow = nsites_graze)
  psi_inxs_covs_graze[,1] <- 1
  psi_inxs_covs_graze[,2] <- stations_graze$GrazingActivity
  
  psi_inxs_covs_hunt <- matrix(NA, ncol = 3, nrow = nsites_hunt)
  psi_inxs_covs_hunt[,1] <- 1
  psi_inxs_covs_hunt[,2] <- stations_hunt$GrazingActivity
  psi_inxs_covs_hunt[,3] <- stations_hunt$Public
  
  #'  Matrix for first order detection (rho): main effects
  rho_cov_graze <- array(NA, dim = c(nsites_graze, nsurvey_graze, 3)) # last digit = number of covariates
  rho_cov_graze[,,1] <- 1
  rho_cov_graze[,1:nsurvey_graze,2] <- gTrail   # repeat site-level covariate for each survey
  rho_cov_graze[,,3] <- cattle_week_scaled
  head(rho_cov_graze)
  
  rho_cov_hunt <- array(NA, dim = c(nsites_hunt, nsurvey_hunt, 4))
  rho_cov_hunt[,,1] <- 1
  rho_cov_hunt[,1:nsurvey_hunt,2] <- hTrail    # repeat site-level covariate for each survey
  rho_cov_hunt[,1:nsurvey_hunt,3] <- Public    # repeat site-level covariate for each survey
  rho_cov_hunt[,,4] <- hunter_week_scaled
  head(rho_cov_hunt)
  
  #'  Matrix for second order detection (rho): 2-way interactions
  rho_inxs_cov_graze <- rep(1, nsites_graze)
  rho_inxs_cov_hunt <- rep(1, nsites_hunt)
  
  #'  Function to bundle detection and covariate data
  bundle_data <- function(det_array, psi_covs, psi_inxs, rho_covs, rho_inxs, 
                          sites, surveys, psi_1order, psi_2order, rho_1order, 
                          rho_2order, ncats) {
    #'  list all pieces of data together
    bundled <- list(y = det_array, psi_cov = psi_covs, psi_inxs_covs = psi_inxs,
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
                               psi_inxs = psi_inxs_covs_graze, rho_covs = rho_cov_graze, 
                               rho_inxs = rho_inxs_cov_graze, sites = nsites_graze, 
                               surveys = nsurvey_graze, psi_1order = psi_cov_graze, 
                               psi_2order = psi_inxs_covs_graze, rho_1order = rho_cov_graze, 
                               rho_2order = 1, ncats = ncat_graze)
  bundled_hunt_list <- lapply(det_array_hunt, bundle_data, psi_covs = psi_cov_hunt, 
                              psi_inxs = psi_inxs_covs_hunt, rho_covs = rho_cov_hunt, 
                              rho_inxs = rho_inxs_cov_hunt, sites = nsites_hunt, 
                              surveys = nsurvey_hunt, psi_1order = psi_cov_hunt, 
                              psi_2order = psi_inxs_covs_hunt, rho_1order = rho_cov_hunt, 
                              rho_2order = 1, ncats = ncat_hunt)
  
  
  ####  Specify model  ####
  #'  =================
  #'  1. Define basic hierarchical model that separates latent and observed processes.
  #'  Each process uses a Categorical distribution since this is essentially a 
  #'  multi-state occupancy model.
  
  #'  Specify model in JAGS language
  cat(file = "multi-spp_OccMod.txt", "
      model{
      
      ##### Define Priors  ####
      #'  First order psi priors
      betaSpp1[1] <- logit(mean.psiSpp1)        # fo occupancy intercepts
      betaSpp2[1] <- logit(mean.psiSpp2)
      betaSpp3[1] <- logit(mean.psiSpp3)
      mean.psiSpp1[1] ~ dunif(0, 1)
      mean.psiSpp2[1] ~ dunif(0, 1)
      mean.psiSpp3[1] ~ dunif(0, 1)
      for(fo_psi in 2:nsecond_order_psi) {      # fo occupancy slopes
        betaSpp1[fo_psi] ~ dnorm(0, 0.1)    ####----- Sarah- make sure this precision makes sense
        betaSpp2[fo_psi] ~ dnorm(0, 0.1)
        betaSpp3[fo_psi] ~ dnorm(0, 0.1)
      }
      #'  Second order psi priors
      for[so_psi in 1:nsecond_order_psi] {
        betaSpp12[so_psi] ~ dnorm(0, 0.1)
        betaSpp23[so_psi] ~ dnorm(0, 0.1)
        betaSpp13[so_psi] ~ dnorm(0, 0.1)
      }
      #'  First order detection priors (rho)
      alphaSpp1[1] <- logit(mean.pSpp1)        # fo detection intercepts
      alphaSpp2[1] <- logit(mean.pSpp2)
      alphaSpp3[1] <- logit(mean.pSpp3)
      mean.pSpp1 ~ dunif(0, 1)
      mean.pSpp2 ~ dunif(0, 1)
      mean.pSpp3 ~ dunif(0, 1)
      for(fo_rho in 2:nfirst_order_rho) {      # fo detection slopes
        alphaSpp1[fo_rho] ~ dnorm(0, 0.01)
        alphaSpp2[fo_rho] ~ dnorm(0, 0.01)
        alphaSpp3[fo_rho] ~ dnorm(0, 0.01)
      }
      #'  Second order detection priors (rho)
      #'  none for now
      
      ####  Define Likelihood  ####
      #'  (1) Set up basic hierarchical model
      #'  Latent state model
      #'  For each site, true occupancy (z) is drawn from a categorical distribution
      #'  with 8 mututally exclusive occupancy probabilities
      for(i in 1:nsites) {
        z[i] ~ dcat(lsv[i, (1:ncat)])
      }
      #'  Detection model
      #'  For each site and survey occasion, the deteciton data are drawn from a
      #'  categorical distribution with 8 mututally exclusive latent states (z)
      for(i in 1:nsites) {
        for(j in 1:nsurveys) {
          y[i,j] ~ dcat(rdm[i, j, (1:ncat), z])
        }
      }
      
      #'  (2) Define latent state vector and observation matrices
      
      
      
      
      
      
      #'  2. Define arrays containing cell probabilities for categorical distributions.
      #'  3. Define linear models for each fundamental parameter that governs the cell probs.

      }")
  
    

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  