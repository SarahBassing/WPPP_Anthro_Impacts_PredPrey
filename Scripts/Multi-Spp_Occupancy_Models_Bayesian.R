  #'  =====================================================
  #'  Bayesian Multi-Species Co-Occurrence Occupancy Models 
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  July 2022
  #'  =====================================================
  #'  Script to source data formatting scripts & run multi-species, single-season
  #'  occupancy models within a Bayesian framework. Code adapted from Kery & Royle
  #'  AHM2 book (Ch 8.2.3). Code runs 3-spp models using detection/non-detection
  #'  data for deer, elk, moose, black bears, cougars, wolves, coyotes, & bobcats
  #'  during the grazing season 2018-2020 (July - Sept) and hunting season 
  #'  2018-2020 (Oct - Nov), respectively. Grazing season co-occurrence 
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
  
  #library(jagsUI)
  library(R2jags)
  library(abind)
  library(mcmcplots)
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
  source("./Scripts/Data_Formatting_3SppX_OccMods.R")
  
  
  ####  Bundle detection histories  ####
  #'  ==============================
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
  
  
  ####  Identify survey dimensions  ####
  #'  ==============================
  nsites_graze <- dim(coug_md_cattle_graze_dat)[1]
  nsurveys_graze <- dim(coug_md_cattle_graze_dat)[2]
  nspp_graze <- dim(coug_md_cattle_graze_dat)[3]
  #'  Number of possible community states (species interactions)
  ncat_graze <- 2^nspp_graze
  
  nsites_hunt <- dim(coug_md_hunter_hunt_dat)[1]
  nsurveys_hunt <- dim(coug_md_hunter_hunt_dat)[2]
  nspp_hunt <- dim(coug_md_hunter_hunt_dat)[3]
  #'  Number of possible community states (species interactions)
  ncat_hunt <- 2^nspp_hunt
  
  
  ####  Format detection histories for JAGS  ####
  #'  =======================================
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
    #' #'  Drop all row and column names generated by camTrapR
    #' attr(ycat, "dimnames") <- NULL
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
  
  
  ####  Format covariates for JAGS  ####
  #'  ==============================
  #'  Prepare site-level covariates for detection sub-models
  #'  study area covariate: 0 = NE, 1 = OK
  #'  trail covariate: 0 = trail, 1 = dirt road
  #'  public covariate: 0 = private, 1 = public
  covs_graze <- stations_graze %>%
    mutate(Study_Area = ifelse(Study_Area == "NE", 0, 1),
           Trail = ifelse(Trail == "Trail", 0, 1), 
           Public = as.integer(as.character(Public)))
  covs_hunt <- stations_hunt %>%
    mutate(Study_Area = ifelse(Study_Area == "NE", 0, 1),
           Trail = ifelse(Trail == "Trail", 0, 1), 
           Public = as.integer(as.character(Public)))
  table(gTrail <- covs_graze[,"Trail"])
  table(hTrail <- covs_hunt[,"Trail"])
  table(Public <- covs_hunt[,"Public"])
  
  #'  Janky patch to fill NAs for weekly grazing/hunter activity --> MUST find a better solution for this
  cattle_week_cov <- cattle_week_scaled
  cattle_week_cov[is.na(cattle_week_cov)] <- -0.1546436 # fill NAs with scaled value that equals 0
  hunter_week_cov <- hunter_week_scaled
  hunter_week_cov[is.na(hunter_week_cov)] <- -0.1565421 # fill NAs with scaled value that equals 0
  
  #'  Matrix for first order occupancy (psi): main effects
  psi_cov_graze <- matrix(NA, ncol = 5, nrow = nsites_graze)
  psi_cov_graze[,1] <- 1
  psi_cov_graze[,2] <- covs_graze$Study_Area
  psi_cov_graze[,3] <- covs_graze$Elev
  psi_cov_graze[,4] <- (covs_graze$Elev)^2
  psi_cov_graze[,5] <- covs_graze$PercForest
  head(psi_cov_graze)
  
  psi_cov_hunt <- matrix(NA, ncol = 5, nrow = nsites_hunt)
  psi_cov_hunt[,1] <- 1
  psi_cov_hunt[,2] <- covs_hunt$Study_Area
  psi_cov_hunt[,3] <- covs_hunt$Elev
  psi_cov_hunt[,4] <- (covs_hunt$Elev)^2
  psi_cov_hunt[,5] <- covs_hunt$PercForest
  head(psi_cov_hunt)
  
  #'  Matrix for second order occupancy (psi): 2-way interactions
  psi_inxs_cov_graze <- matrix(NA, ncol = 2, nrow = nsites_graze)
  psi_inxs_cov_graze[,1] <- 1
  psi_inxs_cov_graze[,2] <- covs_graze$GrazingActivity
  head(psi_inxs_cov_graze)
  
  psi_inxs_cov_hunt <- matrix(NA, ncol = 3, nrow = nsites_hunt) 
  psi_inxs_cov_hunt[,1] <- 1
  psi_inxs_cov_hunt[,2] <- covs_hunt$HuntingActivity
  psi_inxs_cov_hunt[,3] <- covs_hunt$Public
  head(psi_inxs_cov_hunt)
  
  #'  Matrix for first order detection (rho): main effects
  rho_cov_graze <- array(NA, dim = c(nsites_graze, nsurveys_graze, 3)) # last digit = number of covariates
  rho_cov_graze[,,1] <- 1
  rho_cov_graze[,1:nsurveys_graze,2] <- gTrail   # repeat site-level covariate for each survey
  rho_cov_graze[,,3] <- cattle_week_cov
  head(rho_cov_graze)
  
  rho_cov_hunt <- array(NA, dim = c(nsites_hunt, nsurveys_hunt, 4)) # last digit = number of covariates
  rho_cov_hunt[,,1] <- 1
  rho_cov_hunt[,1:nsurveys_hunt,2] <- hTrail    # repeat site-level covariate for each survey
  rho_cov_hunt[,1:nsurveys_hunt,3] <- Public    # repeat site-level covariate for each survey
  rho_cov_hunt[,,4] <- hunter_week_cov
  head(rho_cov_hunt)
  
  #'  Matrix for second order detection (rho): 2-way interactions
  rho_inxs_cov_graze <- rep(1, nsites_graze)
  rho_inxs_cov_hunt <- rep(1, nsites_hunt)
  
  
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
  bundled_hunt_list <- lapply(det_array_hunt, bundle_data, psi_covs = psi_cov_hunt, 
                              psi_inxs = psi_inxs_cov_hunt, rho_covs = rho_cov_hunt, 
                              rho_inxs = rho_inxs_cov_hunt, sites = nsites_hunt, 
                              surveys = nsurveys_hunt, psi_1order = psi_cov_hunt, 
                              psi_2order = psi_inxs_cov_hunt, rho_1order = rho_cov_hunt, 
                              rho_2order = 1, ncats = ncat_hunt)
  
  
  ####  Specify model in JAGS language  ####
  #'  ==================================
  #sink("./Outputs/JAGS_models/multi-spp_OccMod.txt")
  cat(file = 'multi-spp_OccMod.txt', "
  model{
  
    ##### Define Priors  ####
    #'  =================
    #'  Priors for parameters of interest 
    #'  Intercepts and slopes for linear models associated with each natural parameter
  
    betaSpp1[1] <- logit(mean.psiSpp1)          # fo occupancy intercepts (logit scale)
    betaSpp2[1] <- logit(mean.psiSpp2)
    betaSpp3[1] <- logit(mean.psiSpp3)
    mean.psiSpp1[1] ~ dunif(0, 1)               # on the probability scale
    mean.psiSpp2[1] ~ dunif(0, 1)
    mean.psiSpp3[1] ~ dunif(0, 1)
    
    for(fo_psi in 2:nfirst_order_psi){          # fo occupancy slopes (logit scale)
      betaSpp1[fo_psi] ~ dnorm(0, 0.1)
      betaSpp2[fo_psi] ~ dnorm(0, 0.1)
      betaSpp3[fo_psi] ~ dnorm(0, 0.1)
    }    
    
    #'  Second order psi priors
    for(so_psi in 1:nsecond_order_psi){
      betaSpp12[so_psi] ~ dnorm(0, 0.1)
      betaSpp13[so_psi] ~ dnorm(0, 0.1)
      betaSpp23[so_psi] ~ dnorm(0, 0.1)
    }

    #'  First order detection priors (rho)
    alphaSpp1[1] <- logit(mean.pSpp1)           # fo detection intercepts (logit scale)
    alphaSpp2[1] <- logit(mean.pSpp2)
    alphaSpp3[1] <- logit(mean.pSpp3)
    mean.pSpp1 ~ dunif(0, 1)                    # on the probability scale
    mean.pSpp2 ~ dunif(0, 1)
    mean.pSpp3 ~ dunif(0, 1)
    
    for(fo_rho in 2:nfirst_order_rho){          # fo detection slopes (logit scale)
      alphaSpp1[fo_rho] ~ dnorm(0, 0.1)
      alphaSpp2[fo_rho] ~ dnorm(0, 0.1)
      alphaSpp3[fo_rho] ~ dnorm(0, 0.1)
    }
  
    #'  Second order detection priors (rho)
    #'  none for now
        
    
    ####  Define Likelihood  ####
    #'  =====================

    #'  1. Set up basic hierarchical model
    #'  Latent state and observation processes uses a Categorical distribution 
    #'  since this is essentially a multi-state occupancy model.

    #'  Latent state model
    #'  ------------------
    #'  For each site, true occupancy (z) is drawn from a categorical distribution
    #'  with 8 mututally exclusive occupancy probabilities
    
    for(i in 1:nsites) {
      z[i] ~ dcat(lsv[i, (1:ncat)])
    }
  
    #'  Observation model
    #'  -----------------
    #'  For each site and survey occasion, the deteciton data are drawn from a
    #'  categorical distribution with 8 latent states (z)
  
    for(i in 1:nsites) {
      for(j in 1:nsurveys) {
        y[i,j] ~ dcat(rdm[i, j, (1:ncat), z[i]])
      }
    }
  
    #'  2. Define arrays containing cell probabilities for categorical distributions
      
    for(i in 1:nsites) {
      #'  Probabilities for each latent state (z), held in latent state vector (lsv)
      lsv[i, 1] <- 1                   # Unoccupied
      lsv[i, 2] <- exp(psiSpp1[i])     # Pr(Spp1 present)
      lsv[i, 3] <- exp(psiSpp2[i])     # Pr(Spp2 present)
      lsv[i, 4] <- exp(psiSpp3[i])     # Pr(Spp3 present)
      lsv[i, 5] <- exp(psiSpp12[i])    # Pr(Spp1 & Spp2 present)
      lsv[i, 6] <- exp(psiSpp23[i])    # Pr(Spp2 & Spp3 present)
      lsv[i, 7] <- exp(psiSpp13[i])    # Pr(Spp1 & Spp3 present)
      lsv[i, 8] <- exp(psiSpp123[i])   # Pr(all spp present)
 
      for(j in 1:nsurveys) {
        #'  Probabilities for each detection array, held in rho detection matrix (rdm) 
        #'  where OS = observed state, TS = true state and each row sums to 1. 
        #'  Exponentiating log odds so rdm holds estimates on probability scale.???
        #'  Reminder - this model assumes NO false positives in the data so
        #'  probability is 0 when OS x TS combinations are not possible.
        #'  Mmmk don't freak out over this section!
        #'  Example 1: when only Spp1 is observed and only Spp1 is truly 
        #'  present, the detection probability is rhoSpp1.
        #'  Example 2: when only Spp1 is observed by in reality Spp1 & Spp2
        #'  are truly present, the detection probability is 
        #'  True state = unoccupied (z = 1 --> 000)
        rdm[i, j, 1, 1] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 1] <- 0 # ------------------------------------ OS = Spp1 present
        rdm[i, j, 3, 1] <- 0 # ------------------------------------ OS = Spp2 present
        rdm[i, j, 4, 1] <- 0 # ------------------------------------ OS = Spp3 present
        rdm[i, j, 5, 1] <- 0 # ------------------------------------ OS = Spp12 present
        rdm[i, j, 6, 1] <- 0 # ------------------------------------ OS = Spp23 present
        rdm[i, j, 7, 1] <- 0 # ------------------------------------ OS = Spp13 present
        rdm[i, j, 8, 1] <- 0 # ------------------------------------ OS = all species present
        #'  True state = Spp1 present (z = 2 --> 100)
        rdm[i, j, 1, 2] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 2] <- exp(rhoSpp1[i, j]) # ------------------- OS = Spp1 present
        rdm[i, j, 3, 2] <- 0 # ------------------------------------ OS = Spp2 present
        rdm[i, j, 4, 2] <- 0 # ------------------------------------ OS = Spp3 present
        rdm[i, j, 5, 2] <- 0 # ------------------------------------ OS = Spp12 present
        rdm[i, j, 6, 2] <- 0 # ------------------------------------ OS = Spp23 present
        rdm[i, j, 7, 2] <- 0 # ------------------------------------ OS = Spp13 present
        rdm[i, j, 8, 2] <- 0 # ------------------------------------ OS = all species present
        #'  True state = Spp2 present (z = 3 --> 010 )
        rdm[i, j, 1, 3] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 3] <- 0 # ------------------------------------ OS = Spp1 present
        rdm[i, j, 3, 3] <- exp(rhoSpp2[i, j]) # ------------------- OS = Spp2 present
        rdm[i, j, 4, 3] <- 0 # ------------------------------------ OS = Spp3 present
        rdm[i, j, 5, 3] <- 0 # ------------------------------------ OS = Spp12 present
        rdm[i, j, 6, 3] <- 0 # ------------------------------------ OS = Spp23 present
        rdm[i, j, 7, 3] <- 0 # ------------------------------------ OS = Spp13 present
        rdm[i, j, 8, 3] <- 0 # ------------------------------------ OS = all species present
        #'  True state = Spp3 present (z = 4 --> 001)
        rdm[i, j, 1, 4] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 4] <- 0 # ------------------------------------ OS = Spp1 present
        rdm[i, j, 3, 4] <- 0 # ------------------------------------ OS = Spp2 present
        rdm[i, j, 4, 4] <- exp(rhoSpp3[i, j]) # ------------------- OS = Spp3 present
        rdm[i, j, 5, 4] <- 0 # ------------------------------------ OS = Spp12 present
        rdm[i, j, 6, 4] <- 0 # ------------------------------------ OS = Spp23 present
        rdm[i, j, 7, 4] <- 0 # ------------------------------------ OS = Spp13 present
        rdm[i, j, 8, 4] <- 0 # ------------------------------------ OS = all species present
        #'  True state = Spp1 & Spp2 present (z = 5 --> 110)
        rdm[i, j, 1, 5] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 5] <- exp(rhoSpp12[i, j]) # ------------------ OS = Spp1 present
        rdm[i, j, 3, 5] <- exp(rhoSpp21[i, j]) # ------------------ OS = Spp2 present
        rdm[i, j, 4, 5] <- 0 # ------------------------------------ OS = Spp3 present
        rdm[i, j, 5, 5] <- exp(rhoSpp12[i, j] + rhoSpp21[i, j]) # - OS = Spp12 present
        rdm[i, j, 6, 5] <- 0 # ------------------------------------ OS = Spp23 present
        rdm[i, j, 7, 5] <- 0 # ------------------------------------ OS = Spp13 present
        rdm[i, j, 8, 5] <- 0 # ------------------------------------ OS = all species present
        #'  True state = Spp2 & Spp3 present (z = 6 --> 011)
        rdm[i, j, 1, 6] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 6] <- 0 # ------------------------------------ OS = Spp1 present
        rdm[i, j, 3, 6] <- exp(rhoSpp23[i, j]) # ------------------ OS = Spp2 present
        rdm[i, j, 4, 6] <- exp(rhoSpp32[i, j]) # ------------------ OS = Spp3 present
        rdm[i, j, 5, 6] <- 0 # ------------------------------------ OS = Spp12 present
        rdm[i, j, 6, 6] <- exp(rhoSpp23[i, j] + rhoSpp32[i, j]) # - OS = Spp23 present
        rdm[i, j, 7, 6] <- 0 # ------------------------------------ OS = Spp13 present
        rdm[i, j, 8, 6] <- 0 # ------------------------------------ OS = all species present
        #'  True state = Spp1 & Spp3 present (z = 7 --> 101)
        rdm[i, j, 1, 7] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 7] <- exp(rhoSpp13[i, j]) # ------------------ OS = Spp1 present
        rdm[i, j, 3, 7] <- 0 # ------------------------------------ OS = Spp2 present
        rdm[i, j, 4, 7] <- exp(rhoSpp31[i, j]) # ------------------ OS = Spp3 present
        rdm[i, j, 5, 7] <- 0 # ------------------------------------ OS = Spp12 present
        rdm[i, j, 6, 7] <- 0 # ------------------------------------ OS = Spp23 present
        rdm[i, j, 7, 7] <- exp(rhoSpp13[i, j] + rhoSpp31[i, j]) # - OS = Spp13 present
        rdm[i, j, 8, 7] <- 0 # ------------------------------------ OS = all species present
        #'  True state = All three species present (z = 8 --> 111)
        rdm[i, j, 1, 8] <- 1 # ------------------------------------ OS = unoccupied
        rdm[i, j, 2, 8] <- exp(rhoSpp123[i, j]) # ----------------- OS = Spp1 present
        rdm[i, j, 3, 8] <- exp(rhoSpp213[i, j]) # ----------------- OS = Spp2 present
        rdm[i, j, 4, 8] <- exp(rhoSpp312[i, j]) # ----------------- OS = Spp3 present
        rdm[i, j, 5, 8] <- exp(rhoSpp123[i, j] + rhoSpp213[i, j]) # OS = Spp12 present
        rdm[i, j, 6, 8] <- exp(rhoSpp213[i, j] + rhoSpp312[i, j]) # OS = Spp23 present
        rdm[i, j, 7, 8] <- exp(rhoSpp123[i, j] + rhoSpp312[i, j]) # OS = Spp13 present
        rdm[i, j, 8, 8] <- exp(rhoSpp123[i, j] + rhoSpp213[i, j] + rhoSpp312[i, j]) # OS = all species present
      }
      
      #'  3. Define linear models for each fundamental parameter that governs the cell probs
      #'  These are my natural parameters (f1, f2, f3, f12, f13, f23, f123)!
      #'  Linear models for the occupancy parameters on the logit scale
      
      #'  ...for states Spp1, Spp2, Spp3
      # psiSpp1[i] <- inprod(betaSpp1, psi_cov[i, ])
      # psiSpp2[i] <- inprod(betaSpp2, psi_cov[i, ])
      # psiSpp3[i] <- inprod(betaSpp3, psi_cov[i, ])
    
      #'  Covariate order: Intercept + Study_Area + Elevation + Elevation^2 + Forest
      psiSpp1[i] <- betaSpp1[1]*psi_cov[i,1] + betaSpp1[2]*psi_cov[i,2] + betaSpp1[3]*psi_cov[i,3] + betaSpp1[4]*psi_cov[i,4] + betaSpp1[5]*psi_cov[i,5]
      psiSpp2[i] <- betaSpp2[1]*psi_cov[i,1] + betaSpp2[2]*psi_cov[i,2] + betaSpp2[3]*psi_cov[i,3] + betaSpp2[4]*psi_cov[i,4] + betaSpp2[5]*psi_cov[i,5]
      psiSpp3[i] <- betaSpp3[1]*psi_cov[i,1] + betaSpp3[2]*psi_cov[i,2] + betaSpp3[3]*psi_cov[i,3] + betaSpp3[4]*psi_cov[i,4] + betaSpp3[5]*psi_cov[i,5]
  
      #'  ...for states Spp12, Spp23, Spp13 (in this specific order!)
      # psiSpp12[i] <- psiSpp1[i] + psiSpp2[i] + inprod(betaSpp12, psi_inxs_cov[i, ])
      # psiSpp23[i] <- psiSpp2[i] + psiSpp3[i] + inprod(betaSpp23, psi_inxs_cov[i, ])
      # psiSpp13[i] <- psiSpp1[i] + psiSpp3[i] + inprod(betaSpp13, psi_inxs_cov[i, ])
  
      #'  Covariate order: Intercept + Trail + GrazingActivity
      psiSpp12[i] <- psiSpp1[i] + psiSpp2[i] + betaSpp12[1]*psi_inxs_cov[i,1] + betaSpp12[2]*psi_inxs_cov[i,2]
      psiSpp23[i] <- psiSpp2[i] + psiSpp3[i] + betaSpp23[1]*psi_inxs_cov[i,1] # beta2 on this interaction
      psiSpp13[i] <- psiSpp1[i] + psiSpp3[i] + betaSpp13[1]*psi_inxs_cov[i,1] # beta2 on this interaction
  
      #'  ...for state Spp123
      # psiSpp123[i] <- psiSpp1[i] + psiSpp2[i] + psiSpp3[i] + inprod(betaSpp12, psi_inxs_cov[i, ]) +
      #                 inprod(betaSpp23, psi_inxs_cov[i, ]) + inprod(betaSpp13, psi_inxs_cov[i, ])
      psiSpp123[i] <- 0  # my attempt to set this natural parameter to 0
  
      #'  Linear models for the detection parameters on the logit scale
      #'  (can specific interactions on detection here as well)
      for(j in 1:nsurveys) {
        #'  Baseline detection linear predictors
        # rhoSpp1[i, j] <- inprod(alphaSpp1, rho_cov[i, j, ])
        # rhoSpp2[i, j] <- inprod(alphaSpp2, rho_cov[i, j, ])
        # rhoSpp3[i, j] <- inprod(alphaSpp3, rho_cov[i, j, ])
  
        rhoSpp1[i, j] <- alphaSpp1[1] + alphaSpp1[2]*rho_cov[i,j,2] + alphaSpp1[3]*0 #rho_cov[i,j,3]
        rhoSpp2[i, j] <- alphaSpp2[1] + alphaSpp2[2]*rho_cov[i,j,2] + alphaSpp2[3]*0 #rho_cov[i,j,3]
        rhoSpp3[i, j] <- alphaSpp3[1] + alphaSpp3[2]*rho_cov[i,j,2] # no anthro activity on this spp
  
        #'  Asymetric interactirons between all 3 species
        #'  Currently forcing interactions to equal baseline detection probs above
        rhoSpp12[i, j] <- rhoSpp1[i, j]
        rhoSpp13[i, j] <- rhoSpp1[i, j]
        rhoSpp21[i, j] <- rhoSpp2[i, j]
        rhoSpp23[i, j] <- rhoSpp2[i, j]
        rhoSpp31[i, j] <- rhoSpp3[i, j]
        rhoSpp32[i, j] <- rhoSpp3[i, j]
  
        #'  Asymetric interactions when all three species are present
        #'  Currently forcing interactions to equal baseline detection probs above
        rhoSpp123[i, j] <- rhoSpp1[i, j]
        rhoSpp213[i, j] <- rhoSpp2[i, j]
        rhoSpp312[i, j] <- rhoSpp3[i, j]
      }
    }
  }
  ")
  
  #'  Provide initial values for model
  #'  Get maximum possible state across all potential surveys at a site (site x spp matrix)
  initial_z <- function(bundled_dat) {
    zinit <- apply(bundled_dat, c(1,3), sum, na.rm = TRUE)
    zinit[zinit > 1] <- 1
    return(zinit)
  }
  #'  Run arrays of species-specific detection data through function
  zinit_graze <- lapply(megaList_graze, initial_z)
  zinit_hunt <- lapply(megaList_hunt, initial_z)
  
  #'  Convert initial values to a category
  initial_z_categories <- function(zlist) {
    #'  Collapse species-specific detection histories into 8 categories
    zcat <- apply(zlist, 1, paste, collapse = "")
    zcat[zcat == "000"] <- 1
    zcat[zcat == "100"] <- 2
    zcat[zcat == "010"] <- 3
    zcat[zcat == "001"] <- 4
    zcat[zcat == "110"] <- 5
    zcat[zcat == "011"] <- 6
    zcat[zcat == "101"] <- 7
    zcat[zcat == "111"] <- 8
    #'  Make z numeric again
    zcat <- as.numeric(zcat)
    return(zcat)
  }
  #'  Run list of detection histories through function to create detection arrays
  #'  for each model
  zcat_graze <- lapply(zinit_graze, initial_z_categories)
  names(zcat_graze) <- c("coug_md_cattle", "coug_elk_cattle", "coug_wtd_cattle", "coug_moose_cattle",
                         "wolf_md_cattle", "wolf_elk_cattle", "wolf_wtd_cattle", "wolf_moose_cattle",
                         "bear_md_cattle", "bear_elk_cattle", "bear_wtd_cattle", "bear_moose_cattle",
                         "bob_md_cattle", "bob_wtd_cattle", "coy_md_cattle", "coy_wtd_cattle")
  zcat_hunt <- lapply(zinit_hunt, initial_z_categories)
  names(zcat_hunt) <- c("coug_md_hunter", "coug_elk_hunter", "coug_wtd_hunter", "coug_moose_hunter",
                        "wolf_md_hunter", "wolf_elk_hunter", "wolf_wtd_hunter", "wolf_moose_hunter",
                        "bear_md_hunter", "bear_elk_hunter", "bear_wtd_hunter", "bear_moose_hunter",
                        "bob_md_hunter", "bob_wtd_hunter", "coy_md_hunter", "coy_wtd_hunter")
  
  #'  Inits function
  inits <- function() list(z = zcat_graze[[1]])
  
  #'  Parameters monitored
  params <- c("betaSpp1", "betaSpp2", "betaSpp3", "betaSpp12", "betaSpp23", "betaSpp13",
              "alphaSpp1", "alphaSpp2", "alphaSpp3", "mean.psiSpp1", "mean.psiSpp2", "mean.psiSpp3",
              "mean.pSpp1", "mean.pSpp2", "mean.pSpp3", "z")

  #'  MCMC settings
  na <- 10000; nc <- 3; ni <- 50000; nb <- 30000; nt <- 20
  
  start.time <- Sys.time()
  #'  Call JAGS, check convergence and summarize posteriors
  out1 <- jags.parallel(bundled_graze_list[[1]], inits, params, 'multi-spp_OccMod.txt', n.chains = nc,
               n.burnin = nb, n.iter = ni, n.thin = nt, jags.module = c("glm","dic"), n.cluster= n.chains) #n.adapt = na, parallel = TRUE
  end.time <- Sys.time()
  (run.time <- end.time - start.time)
  
  out2 <- out1
  mcmcplot(out2)
  jag.sum <- out1$BUGSoutput$summary
  print(jag.sum[1:37,-c(4:6)], 3)
  which(jag.sum[,"Rhat"] > 1.1)
  write.table(x=jag.sum, file="Test_MultiSpp_Model.RData", sep="\t")
  
  # traceplot(jag.sum)
  # which(out1$summary[,"Rhat"] > 1.1) #summary[,8]
  # print(out1$summary[1:36, -c(4:6)], 3)
  # save(list=ls(), file="Test_MultiSpp_Model.Rdata")
  
  
  
  # QUESTIONS/TO-DO List
  # 1. did I correctly set the 3-way interaction parameter to 0 I don't estimate it?
  # 2. how to functionalize this so I can run models with different covariate per species 
  # and for multiple species without rewriting whole thing over and over
  # 3. what to do about these missing survey-level covs
  
  
  
  
  