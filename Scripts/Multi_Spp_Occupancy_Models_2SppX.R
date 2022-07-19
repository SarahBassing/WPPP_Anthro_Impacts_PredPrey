  #'  ============================================
  #'  Multi-Species Co-Occurance Occupancy Models
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2022
  #'  ============================================
  #'  Script to create unmarked data frames and run multi-species, single-season
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
  #'  Data_Formatting_2SppX_OccMod.R script.
  #'  ============================================

  #'  Clean workspace & load libraries
  rm(list = ls())
  
  library(unmarked)
  library(MuMIn)
  library(condformat)
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
  #'  models in unmarked with TWO-SPECIES interactions
  source("./Scripts/Data_Formatting_2SppX_OccMods.R")
  
  
  ####  Multi-Species Occupancy models  ####
  #'  ==================================
  #'  Multi-species occupancy model --> occuMulti (pg 83) in unmarked manual
  #'  Occupancy formulas should match number/order of columns in fDesign matrix
  #'  i.e., one formula for each species/interaction of interest
  #'  Detection formulas should match number/order of species in list of DH
  #'  Use ~1 to estimate intercept without covariates
  #'  Use ~0 to fix a natural parameter to 0
  #'  E.g., occFormulas <- c("~1", "~1", "~1", "~0", "~1", "~1", "~0") estimates 
  #'  intercept for 1st order natural parameters (3 spp), fixes first 2nd order  
  #'  natural parameter to 0 but estimates intercept for other 2nd order parameters, 
  #'  and fixes 3rd order natural parameter to 0 (only applies if 3+ spp)
  #'  Covariates: Can use different covariates on different natural parameters, 
  #'  E.g., covs on 1st order parameters to explain single-spp occurrence 
  #'  regardless of other spp, covs on 2nd order parameters to explain co-occ
  #'  
  #'  Testing hypothesis that co-occurrence is non-independent and that cattle/
  #'  hunter activity impacts occurrence and/or co-occurrence patterns between
  #'  predators and prey.
  #'  
  #'  Include a consistent set of additional covariates to account for habitat
  #'  variation and other factors we know influence occurrence and deteciton.
  #'  Use AIC for model selection
  #'  =============================
  
  ####  GRAZING SEASON MODELS  ####
  #'  Detection and occupancy formulas
  #'  Question 1: Does grazing activity affect detection?
  detFormulas_null <- c("~1", "~1")
  detFormulas_trail <- c("~Trail", "~Trail")
  detFormulas_graze <- c("~Trail + WeeklyGrazing", "~Trail + WeeklyGrazing")
  
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null <- c("~1", "~1", "~0")
  #occFormulas_null1 <- c("~1", "~1", "~1")
  occFormulas_hab0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~0")
  occFormulas_hab1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~1")
  #'  Question 2: Does grazing activity affect occurrence and/or co-occurrence patterns
  occFormulas_graze0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity",
                         "~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity",
                         "~0")
  occFormulas_graze1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~1")
  occFormulas_graze2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + GrazingActivity", 
                          "~GrazingActivity")

  
  ####  Cougar-Mule Deer Grazing Season  ####
  (gs_cougmd_null <- occuMulti(detFormulas_null, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougmd_fl <- fitList(gs_cougmd_null, gs_cougmd_hab0, gs_cougmd_hab1, gs_cougmd_graze0, gs_cougmd_graze1, gs_cougmd_graze2)
  #' Model selection
  modSel(gs_cougmd_fl)
  summary(gs_cougmd_hab0)
  summary(gs_cougmd_hab1)
  
  #' (gs_cougmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  #' (gs_cougmd_graze <- occuMulti(detFormulas_graze, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougmd_fld <- fitList(gs_cougmd_trail, gs_cougmd_graze)
  #' #' Model selection
  #' modSel(gs_cougmd_fld)
  #' 
  #' (gs_cougmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  #' (gs_cougmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_md_grazing_UMF, silent = TRUE))
  #' (gs_cougmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_md_grazing_UMF, silent = TRUE))
  #' (gs_cougmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_md_grazing_UMF, silent = TRUE))
  #' (gs_cougmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, coug_md_grazing_UMF, silent = TRUE))
  #' (gs_cougmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, coug_md_grazing_UMF, silent = TRUE))
  #' (gs_cougmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, coug_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougmd_fl <- fitList(gs_cougmd_null0, gs_cougmd_null1, gs_cougmd_hab0, gs_cougmd_hab1, gs_cougmd_graze0, gs_cougmd_graze1, gs_cougmd_graze2)
  #' #' Model selection
  #' modSel(gs_cougmd_fl)
  #' summary(gs_cougmd_hab0)
  #' summary(gs_cougmd_hab1) # f12 not significant
  
  ####  Cougar-ELK Grazing Season  ####
  (gs_cougelk_null <- occuMulti(detFormulas_null, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_elk_grazing_UMF, silent = TRUE)) #fail
  (gs_cougelk_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougelk_fl <- fitList(gs_cougelk_null, gs_cougelk_hab0, gs_cougelk_hab1, gs_cougelk_graze0, gs_cougelk_graze2)
  #' Model selection
  modSel(gs_cougelk_fl)
  summary(gs_cougelk_hab0)
  summary(gs_cougelk_hab1)
  
  #' (gs_cougelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  #' (gs_cougelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougelk_fld <- fitList(gs_cougelk_trail, gs_cougelk_dgraze)
  #' #' Model selection
  #' modSel(gs_cougelk_fld)
  #' 
  #' (gs_cougelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  #' (gs_cougelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_elk_grazing_UMF, silent = TRUE))
  #' (gs_cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_grazing_UMF, silent = TRUE))
  #' (gs_cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_grazing_UMF, silent = TRUE))
  #' (gs_cougelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, coug_elk_grazing_UMF, silent = TRUE))
  #' (gs_cougelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, coug_elk_grazing_UMF, silent = TRUE)) # fails
  #' (gs_cougelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, coug_elk_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougelk_fl <- fitList(gs_cougelk_null0, gs_cougelk_null1, gs_cougelk_hab0, gs_cougelk_hab1, gs_cougelk_graze0, gs_cougelk_graze2)
  #' #' Model selection
  #' modSel(gs_cougelk_fl)
  #' summary(gs_cougelk_hab0)
  #' # summary(gs_cougelk_hab1) # f12 not significant
  
  ####  Cougar-White-tailed Deer Grazing Season  ####
  (gs_cougwtd_null <- occuMulti(detFormulas_null, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougwtd_fl <- fitList(gs_cougwtd_null, gs_cougwtd_hab0, gs_cougwtd_hab1, gs_cougwtd_graze0, gs_cougwtd_graze1, gs_cougwtd_graze2)
  #' Model selection
  modSel(gs_cougwtd_fl)
  summary(gs_cougwtd_graze2)
  summary(gs_cougwtd_graze0)
  
  #' (gs_cougwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  #' (gs_cougwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougwtd_fld <- fitList(gs_cougwtd_trail, gs_cougwtd_dgraze)
  #' #' Model selection
  #' modSel(gs_cougwtd_fld)
  #' 
  #' (gs_cougwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  #' (gs_cougwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_wtd_grazing_UMF, silent = TRUE))
  #' (gs_cougwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_wtd_grazing_UMF, silent = TRUE))
  #' (gs_cougwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_wtd_grazing_UMF, silent = TRUE))
  #' (gs_cougwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_wtd_grazing_UMF, silent = TRUE))
  #' (gs_cougwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_wtd_grazing_UMF, silent = TRUE))
  #' (gs_cougwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougwtd_fl <- fitList(gs_cougwtd_null0, gs_cougwtd_null1, gs_cougwtd_hab0, gs_cougwtd_hab1, gs_cougwtd_graze0, gs_cougwtd_graze1, gs_cougwtd_graze2)
  #' #' Model selection
  #' modSel(gs_cougwtd_fl)
  #' #summary(gs_cougwtd_graze2) # f12 not significant
  #' summary(gs_cougwtd_graze0)
  #' # summary(gs_cougwtd_graze1) # f12 not significant
  
  ####  Cougar-Moose Grazing Season  ####
  (gs_cougmoose_null <- occuMulti(detFormulas_null, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougmoose_fl <- fitList(gs_cougmoose_null, gs_cougmoose_hab0, gs_cougmoose_hab1, gs_cougmoose_graze0, gs_cougmoose_graze1, gs_cougmoose_graze2)
  #' Model selection
  modSel(gs_cougmoose_fl)
  summary(gs_cougmoose_graze0)
  summary(gs_cougmoose_hab0)
  
  #' (gs_cougmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  #' (gs_cougmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougmoose_fld <- fitList(gs_cougmoose_trail, gs_cougmoose_dgraze)
  #' #' Model selection
  #' modSel(gs_cougmoose_fld)
  #' 
  #' (gs_cougmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  #' (gs_cougmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_moose_grazing_UMF, silent = TRUE))
  #' (gs_cougmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_moose_grazing_UMF, silent = TRUE))
  #' (gs_cougmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_moose_grazing_UMF, silent = TRUE))
  #' (gs_cougmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_moose_grazing_UMF, silent = TRUE))
  #' (gs_cougmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_moose_grazing_UMF, silent = TRUE))
  #' (gs_cougmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_moose_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_cougmoose_fl <- fitList(gs_cougmoose_null0, gs_cougmoose_null1, gs_cougmoose_hab0, gs_cougmoose_hab1, gs_cougmoose_graze0, gs_cougmoose_graze1, gs_cougmoose_graze2)
  #' #' Model selection
  #' modSel(gs_cougmoose_fl)
  #' summary(gs_cougmoose_hab0)
  #' # summary(gs_cougmoose_hab1) # f12 not significant
  #' summary(gs_cougmoose_graze0)
  
  
  ####  Wolf-Mule Deer Grazing Season  ####
  (gs_wolfmd_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfmd_fl <- fitList(gs_wolfmd_null, gs_wolfmd_hab0, gs_wolfmd_hab1, gs_wolfmd_graze0, gs_wolfmd_graze1, gs_wolfmd_graze2)
  #' Model selection
  modSel(gs_wolfmd_fl)
  summary(gs_wolfmd_hab0)
  summary(gs_wolfmd_hab1)
  
  #' (gs_wolfmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  #' (gs_wolfmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_wolfmd_fld <- fitList(gs_wolfmd_trail, gs_wolfmd_dgraze)
  #' #' Model selection
  #' modSel(gs_wolfmd_fld)
  #' 
  #' (gs_wolfmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  #' (gs_wolfmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_md_grazing_UMF, silent = TRUE))
  #' (gs_wolfmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_md_grazing_UMF, silent = TRUE))
  #' (gs_wolfmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_md_grazing_UMF, silent = TRUE))
  #' (gs_wolfmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, wolf_md_grazing_UMF, silent = TRUE))
  #' (gs_wolfmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, wolf_md_grazing_UMF, silent = TRUE))
  #' (gs_wolfmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, wolf_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_wolfmd_fl <- fitList(gs_wolfmd_null0, gs_wolfmd_null1, gs_wolfmd_hab0, gs_wolfmd_hab1, gs_wolfmd_graze0, gs_wolfmd_graze1, gs_wolfmd_graze2)
  #' #' Model selection
  #' modSel(gs_wolfmd_fl)
  #' summary(gs_wolfmd_hab0)
  #' # summary(gs_wolfmd_hab1) # f12 not significant
  
  
  ####  Wolf-Elk Grazing Season  ####
  (gs_wolfelk_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_elk_grazing_UMF, silent = TRUE)) 
  (gs_wolfelk_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_elk_grazing_UMF, silent = TRUE)) # FAILS
  #' List of fitted models
  gs_wolfelk_fl <- fitList(gs_wolfelk_null, gs_wolfelk_hab0, gs_wolfelk_hab1, gs_wolfelk_graze0, gs_wolfelk_graze1)
  #' Model selection
  modSel(gs_wolfelk_fl)
  summary(gs_wolfelk_hab0)
  summary(gs_wolfelk_hab1)
  
  #' (gs_wolfelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  #' (gs_wolfelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_wolfelk_fld <- fitList(gs_wolfelk_trail, gs_wolfelk_dgraze)
  #' #' Model selection
  #' modSel(gs_wolfelk_fld)
  #' 
  #' (gs_wolfelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  #' (gs_wolfelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_elk_grazing_UMF, silent = TRUE))
  #' (gs_wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_grazing_UMF, silent = TRUE))
  #' (gs_wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_grazing_UMF, silent = TRUE))
  #' (gs_wolfelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, wolf_elk_grazing_UMF, silent = TRUE))
  #' (gs_wolfelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, wolf_elk_grazing_UMF, silent = TRUE)) 
  #' (gs_wolfelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, wolf_elk_grazing_UMF, silent = TRUE)) # FAILS
  #' #' List of fitted models
  #' gs_wolfelk_fl <- fitList(gs_wolfelk_null0, gs_wolfelk_null1, gs_wolfelk_hab0, gs_wolfelk_hab1, gs_wolfelk_graze0, gs_wolfelk_graze1)
  #' #' Model selection
  #' modSel(gs_wolfelk_fl)
  #' summary(gs_wolfelk_hab0)
  #' summary(gs_wolfelk_graze0) # not adding any significant info
  #' summary(gs_wolfelk_hab1) # f1 not significant
  
  ####  Wolf-White-tailed Deer Grazing Season  ####
  (gs_wolfwtd_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfwtd_fl <- fitList(gs_wolfwtd_null, gs_wolfwtd_hab0, gs_wolfwtd_hab1, gs_wolfwtd_graze0, gs_wolfwtd_graze1, gs_wolfwtd_graze2)
  #' Model selection
  modSel(gs_wolfwtd_fl)
  summary(gs_wolfwtd_graze0)
  summary(gs_wolfelk_graze1)
  
  #' (gs_wolfwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  #' (gs_wolfwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_wolfwtd_fld <- fitList(gs_wolfwtd_trail, gs_wolfwtd_dgraze)
  #' #' Model selection
  #' modSel(gs_wolfwtd_fld)
  #' 
  #' (gs_wolfwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  #' (gs_wolfwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_wtd_grazing_UMF, silent = TRUE))
  #' (gs_wolfwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_wtd_grazing_UMF, silent = TRUE))
  #' (gs_wolfwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_wtd_grazing_UMF, silent = TRUE))
  #' (gs_wolfwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_wtd_grazing_UMF, silent = TRUE))
  #' (gs_wolfwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_wtd_grazing_UMF, silent = TRUE))
  #' (gs_wolfwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_wolfwtd_fl <- fitList(gs_wolfwtd_null0, gs_wolfwtd_null1, gs_wolfwtd_hab0, gs_wolfwtd_hab1, gs_wolfwtd_graze0, gs_wolfwtd_graze1, gs_wolfwtd_graze2)
  #' #' Model selection
  #' modSel(gs_wolfwtd_fl)
  #' summary(gs_wolfwtd_graze0)
  #' # summary(gs_wolfwtd_graze1) # f12 not significant
  
  ####  Wolf-Moose Grazing Season  ####
  (gs_wolfmoose_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfmoose_fl <- fitList(gs_wolfmoose_null, gs_wolfmoose_hab0, gs_wolfmoose_hab1, gs_wolfmoose_graze0, gs_wolfmoose_graze1, gs_wolfmoose_graze2)
  #' Model selection
  modSel(gs_wolfmoose_fl)
  summary(gs_wolfmoose_graze0)
  summary(gs_wolfmoose_graze1)
  
  #' (gs_wolfmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  #' (gs_wolfmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_wolfmoose_fld <- fitList(gs_wolfmoose_trail, gs_wolfmoose_dgraze)
  #' #' Model selection
  #' modSel(gs_wolfmoose_fld)
  #' 
  #' (gs_wolfmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  #' (gs_wolfmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_moose_grazing_UMF, silent = TRUE))
  #' (gs_wolfmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_moose_grazing_UMF, silent = TRUE))
  #' (gs_wolfmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_moose_grazing_UMF, silent = TRUE))
  #' (gs_wolfmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_moose_grazing_UMF, silent = TRUE))
  #' (gs_wolfmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_moose_grazing_UMF, silent = TRUE))
  #' (gs_wolfmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_moose_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_wolfmoose_fl <- fitList(gs_wolfmoose_null0, gs_wolfmoose_null1, gs_wolfmoose_hab0, gs_wolfmoose_hab1, gs_wolfmoose_graze0, gs_wolfmoose_graze1, gs_wolfmoose_graze2)
  #' #' Model selection
  #' modSel(gs_wolfmoose_fl)
  #' summary(gs_wolfmoose_graze0)
  #' summary(gs_wolfmoose_hab0)
  #' # summary(gs_wolfmoose_graze1) # f12 not significant
  
  
  
  ####  Bear-Mule Deer Grazing Season  ####
  (gs_bearmd_null <- occuMulti(detFormulas_null, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearmd_fl <- fitList(gs_bearmd_null, gs_bearmd_hab0, gs_bearmd_hab1, gs_bearmd_graze0, gs_bearmd_graze1, gs_bearmd_graze2)
  #' Model selection
  modSel(gs_bearmd_fl)
  summary(gs_bearmd_hab1) # f12 not significant
  summary(gs_bearmd_graze1) # f12 not significant
  summary(gs_bearmd_hab0)
  
  #' (gs_bearmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  #' (gs_bearmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bearmd_fld <- fitList(gs_bearmd_trail, gs_bearmd_dgraze)
  #' #' Model selection
  #' modSel(gs_bearmd_fld)
  #' 
  #' (gs_bearmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  #' (gs_bearmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_md_grazing_UMF, silent = TRUE))
  #' (gs_bearmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_md_grazing_UMF, silent = TRUE))
  #' (gs_bearmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_md_grazing_UMF, silent = TRUE))
  #' (gs_bearmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bear_md_grazing_UMF, silent = TRUE))
  #' (gs_bearmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bear_md_grazing_UMF, silent = TRUE))
  #' (gs_bearmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bear_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bearmd_fl <- fitList(gs_bearmd_null0, gs_bearmd_null1, gs_bearmd_hab0, gs_bearmd_hab1, gs_bearmd_graze0, gs_bearmd_graze1, gs_bearmd_graze2)
  #' #' Model selection
  #' modSel(gs_bearmd_fl)
  #' # summary(gs_bearmd_graze1) # f12 not significant
  #' # summary(gs_bearmd_graze0) # grazing not actually significant
  #' # summary(gs_bearmd_graze2) # f12 and grazing not significant
  #' # summary(gs_bearmd_hab1) # f12 not significant
  #' summary(gs_bearmd_hab0)
  
  ####  Bear-Elk Grazing Season  ####
  (gs_bearelk_null <- occuMulti(detFormulas_null, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_elk_grazing_UMF, silent = TRUE)) 
  (gs_bearelk_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_elk_grazing_UMF, silent = TRUE)) # not converging well
  #' List of fitted models
  gs_bearelk_fl <- fitList(gs_bearelk_null, gs_bearelk_hab0, gs_bearelk_hab1, gs_bearelk_graze0, gs_bearelk_graze1)
  #' Model selection
  modSel(gs_bearelk_fl)
  summary(gs_bearelk_hab1) # f12 not significant
  summary(gs_bearelk_graze1) # f12 not significant
  summary(gs_bearelk_hab0)
  summary(gs_bearelk_graze0)
  
  #' (gs_bearelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  #' (gs_bearelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bearelk_fld <- fitList(gs_bearelk_trail, gs_bearelk_dgraze)
  #' #' Model selection
  #' modSel(gs_bearelk_fld)
  #' 
  #' (gs_bearelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  #' (gs_bearelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_elk_grazing_UMF, silent = TRUE))
  #' (gs_bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_grazing_UMF, silent = TRUE))
  #' (gs_bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_grazing_UMF, silent = TRUE))
  #' (gs_bearelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bear_elk_grazing_UMF, silent = TRUE))
  #' (gs_bearelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bear_elk_grazing_UMF, silent = TRUE)) 
  #' (gs_bearelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bear_elk_grazing_UMF, silent = TRUE)) # FAILS
  #' #' List of fitted models
  #' gs_bearelk_fl <- fitList(gs_bearelk_null0, gs_bearelk_null1, gs_bearelk_hab0, gs_bearelk_hab1, gs_bearelk_graze0, gs_bearelk_graze1)
  #' #' Model selection
  #' modSel(gs_bearelk_fl)
  #' # summary(gs_bearelk_graze1) # f12 not significant
  #' # summary(gs_bearelk_graze0) # grazing not significant
  #' # summary(gs_bearelk_hab1) # f12 not signigicant
  #' summary(gs_bearelk_hab0)
  
  ####  Bear-White-tailed Deer Grazing Season  ####
  (gs_bearwtd_null <- occuMulti(detFormulas_null, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearwtd_fl <- fitList(gs_bearwtd_null, gs_bearwtd_hab0, gs_bearwtd_hab1, gs_bearwtd_graze0, gs_bearwtd_graze1, gs_bearwtd_graze2)
  #' Model selection
  modSel(gs_bearwtd_fl)
  summary(gs_bearwtd_graze0)
  summary(gs_bearwtd_graze1)
  
  #' (gs_bearwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bearwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bearwtd_fld <- fitList(gs_bearwtd_trail, gs_bearwtd_dgraze)
  #' #' Model selection
  #' modSel(gs_bearwtd_fld)
  #' 
  #' (gs_bearwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bearwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bearwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bearwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bearwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bearwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bearwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bearwtd_fl <- fitList(gs_bearwtd_null0, gs_bearwtd_null1, gs_bearwtd_hab0, gs_bearwtd_hab1, gs_bearwtd_graze0, gs_bearwtd_graze1, gs_bearwtd_graze2)
  #' #' Model selection
  #' modSel(gs_bearwtd_fl)
  #' summary(gs_bearwtd_graze0)
  #' # summary(gs_bearwtd_graze1) # f12 not significant
  
  
  ####  Bear-Moose Grazing Season  ####
  (gs_bearmoose_null <- occuMulti(detFormulas_null, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_moose_grazing_UMF, silent = TRUE)) # not converging well
  #' List of fitted models
  gs_bearmoose_fl <- fitList(gs_bearmoose_null, gs_bearmoose_hab0, gs_bearmoose_hab1, gs_bearmoose_graze0, gs_bearmoose_graze1)
  #' Model selection
  modSel(gs_bearmoose_fl)
  summary(gs_bearmoose_graze0)
  summary(gs_bearmoose_graze1) 
  
  #' (gs_bearmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  #' (gs_bearmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bearmoose_fld <- fitList(gs_bearmoose_trail, gs_bearmoose_dgraze)
  #' #' Model selection
  #' modSel(gs_bearmoose_fld)
  #' 
  #' (gs_bearmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  #' (gs_bearmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_moose_grazing_UMF, silent = TRUE))
  #' (gs_bearmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_moose_grazing_UMF, silent = TRUE))
  #' (gs_bearmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_moose_grazing_UMF, silent = TRUE))
  #' (gs_bearmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_moose_grazing_UMF, silent = TRUE))
  #' (gs_bearmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_moose_grazing_UMF, silent = TRUE))
  #' (gs_bearmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_moose_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bearmoose_fl <- fitList(gs_bearmoose_null0, gs_bearmoose_null1, gs_bearmoose_hab0, gs_bearmoose_hab1, gs_bearmoose_graze0, gs_bearmoose_graze1, gs_bearmoose_graze2)
  #' #' Model selection
  #' modSel(gs_bearmoose_fl)
  #' summary(gs_bearmoose_graze0)
  #' # summary(gs_bearmoose_graze2) # f12 not significant
  #' # summary(gs_bearmoose_graze1) # f12 not significant
  
  
  
  ####  Coyote-Mule Deer Grazing Season  ####
  (gs_coymd_null <- occuMulti(detFormulas_null, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_coymd_fl <- fitList(gs_coymd_null, gs_coymd_hab0, gs_coymd_hab1, gs_coymd_graze0, gs_coymd_graze1, gs_coymd_graze2)
  #' Model selection
  modSel(gs_coymd_fl)
  summary(gs_coymd_graze2)
  
  #' (gs_coymd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  #' (gs_coymd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_coymd_fld <- fitList(gs_coymd_trail, gs_coymd_dgraze)
  #' #' Model selection
  #' modSel(gs_coymd_fld)
  #' 
  #' (gs_coymd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  #' (gs_coymd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_md_grazing_UMF, silent = TRUE))
  #' (gs_coymd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_md_grazing_UMF, silent = TRUE))
  #' (gs_coymd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_md_grazing_UMF, silent = TRUE))
  #' (gs_coymd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_md_grazing_UMF, silent = TRUE))
  #' (gs_coymd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_md_grazing_UMF, silent = TRUE))
  #' (gs_coymd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_coymd_fl <- fitList(gs_coymd_null0, gs_coymd_null1, gs_coymd_hab0, gs_coymd_hab1, gs_coymd_graze0, gs_coymd_graze1, gs_coymd_graze2)
  #' #' Model selection
  #' modSel(gs_coymd_fl)
  #' summary(gs_coymd_graze2)
  
  ####  Coyote-White-tailed Deer Grazing Season  ####
  (gs_coywtd_null <- occuMulti(detFormulas_null, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_coywtd_fl <- fitList(gs_coywtd_null, gs_coywtd_hab0, gs_coywtd_hab1, gs_coywtd_graze0, gs_coywtd_graze1, gs_coywtd_graze2)
  #' Model selection
  modSel(gs_coywtd_fl)
  summary(gs_coywtd_graze2) 
  
  #' (gs_coywtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  #' (gs_coywtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_coywtd_fld <- fitList(gs_coywtd_trail, gs_coywtd_dgraze)
  #' #' Model selection
  #' modSel(gs_coywtd_fld)
  #' 
  #' (gs_coywtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  #' (gs_coywtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_wtd_grazing_UMF, silent = TRUE))
  #' (gs_coywtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_wtd_grazing_UMF, silent = TRUE))
  #' (gs_coywtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_wtd_grazing_UMF, silent = TRUE))
  #' (gs_coywtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_wtd_grazing_UMF, silent = TRUE))
  #' (gs_coywtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_wtd_grazing_UMF, silent = TRUE))
  #' (gs_coywtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_coywtd_fl <- fitList(gs_coywtd_null0, gs_coywtd_null1, gs_coywtd_hab0, gs_coywtd_hab1, gs_coywtd_graze0, gs_coywtd_graze1, gs_coywtd_graze2)
  #' #' Model selection
  #' modSel(gs_coywtd_fl)
  #' summary(gs_coywtd_graze2) 
  #' summary(gs_coywtd_graze1)
  
  
  ####  Bobcat-Mule Deer Grazing Season  ####
  (gs_bobmd_null <- occuMulti(detFormulas_null, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bob_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bobmd_fl <- fitList(gs_bobmd_null, gs_bobmd_hab0, gs_bobmd_hab1, gs_bobmd_graze0, gs_bobmd_graze1, gs_bobmd_graze2)
  #' Model selection
  modSel(gs_bobmd_fl)
  summary(gs_bobmd_graze0)
  summary(gs_bobmd_graze2)
  summary(gs_bobmd_graze1)
  
  #' (gs_bobmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  #' (gs_bobmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bobmd_fld <- fitList(gs_bobmd_trail, gs_bobmd_dgraze)
  #' #' Model selection
  #' modSel(gs_bobmd_fld)
  #' 
  #' (gs_bobmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  #' (gs_bobmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_md_grazing_UMF, silent = TRUE))
  #' (gs_bobmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_md_grazing_UMF, silent = TRUE))
  #' (gs_bobmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_md_grazing_UMF, silent = TRUE))
  #' (gs_bobmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bob_md_grazing_UMF, silent = TRUE))
  #' (gs_bobmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bob_md_grazing_UMF, silent = TRUE))
  #' (gs_bobmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bob_md_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bobmd_fl <- fitList(gs_bobmd_null0, gs_bobmd_null1, gs_bobmd_hab0, gs_bobmd_hab1, gs_bobmd_graze0, gs_bobmd_graze1, gs_bobmd_graze2)
  #' #' Model selection
  #' modSel(gs_bobmd_fl)
  #' summary(gs_bobmd_graze0)
  #' # summary(gs_bobmd_graze2) # f12 not significant
  #' # summary(gs_bobmd_graze1) # f12 not significant
  
  ####  Bobcat-White-tailed Deer Grazing Season  ####
  (gs_bobwtd_null <- occuMulti(detFormulas_null, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bob_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bobwtd_fl <- fitList(gs_bobwtd_null, gs_bobwtd_hab0, gs_bobwtd_hab1, gs_bobwtd_graze0, gs_bobwtd_graze1, gs_bobwtd_graze2)
  #' Model selection
  modSel(gs_bobwtd_fl)
  summary(gs_bobwtd_graze2) # f12 not significant
  summary(gs_bobwtd_graze0)
  summary(gs_bobwtd_graze1) # f12 not significant
  
  #' (gs_bobwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bobwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bobwtd_fld <- fitList(gs_bobwtd_trail, gs_bobwtd_dgraze)
  #' #' Model selection
  #' modSel(gs_bobwtd_fld)
  #' 
  #' (gs_bobwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bobwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bob_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bobwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bob_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bobwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bob_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bobwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bob_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bobwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bob_wtd_grazing_UMF, silent = TRUE))
  #' (gs_bobwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bob_wtd_grazing_UMF, silent = TRUE))
  #' #' List of fitted models
  #' gs_bobwtd_fl <- fitList(gs_bobwtd_null0, gs_bobwtd_null1, gs_bobwtd_hab0, gs_bobwtd_hab1, gs_bobwtd_graze0, gs_bobwtd_graze1, gs_bobwtd_graze2)
  #' #' Model selection
  #' modSel(gs_bobwtd_fl)
  #' # summary(gs_bobwtd_graze2) # f12 and grazing not actually significant
  #' summary(gs_bobwtd_graze0)
  #' # summary(gs_bobwtd_graze1) # f12 not significant

  
  
  ####  HUNTING SEASON MODELS  ####
  #'  Detection and occupancy formulas
  #'  Question 1: Does hunting activity affect detection?
  detFormulas_null <- c("~1", "~1")
  detFormulas_trail <- c("~Trail", "~Trail")
  detFormulas_pub <- c("~Trail + Public", "~Trail + Public")
  detFormulas_hunt <- c("~Trail + WeeklyHunting", "~Trail + WeeklyHunting") #' not significant for any species
  #detFormulas_pubhunt <- c("~Trail + Public + WeeklyHunting", "~Trail + Public + WeeklyHunting")
  #'  Remove Public vs Private covariate from predator sub-model
  detFormulas_pubish <- c("~Trail", "~Trail + Public")
  #detFormulas_pubhuntish <- c("~Trail + WeeklyHunting", "~Trail + Public + WeeklyHunting")
  
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null <- c("~1", "~1", "~0")
  #occFormulas_null1 <- c("~1", "~1", "~1")
  occFormulas_hab0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~0")
  occFormulas_hab1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~1")
  #'  Question 2: Does hunting activity affect occurrence and/or co-occurrence patterns
  occFormulas_pub0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~0")
  occFormulas_pub1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~1")
  occFormulas_pub2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Public")
  occFormulas_hunt0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~0")
  occFormulas_hunt1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~1")
  occFormulas_hunt2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~HuntingActivity")
  
  #'  Remove Public vs Private covariate from predator sub-model
  occFormulas_pubish0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~0")
  occFormulas_pubish1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~1")
  occFormulas_pubish2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                           "~Public")
  
  
  ####  Cougar-Mule Deer Hunting Season  ####
  (hs_cougmd_null <- occuMulti(detFormulas_null, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_md_hunting_UMF, silent = TRUE)) 
  (hs_cougmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_md_hunting_UMF, silent = TRUE)) 
  (hs_cougmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_md_hunting_UMF, silent = TRUE)) 
  (hs_cougmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_md_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_cougmd_fl <- fitList(hs_cougmd_null, hs_cougmd_hab0, hs_cougmd_hab1, hs_cougmd_pub0, hs_cougmd_pub1, hs_cougmd_pub2, hs_cougmd_hunt0, hs_cougmd_hunt1) 
  #' Model selection
  modSel(hs_cougmd_fl)
  summary(hs_cougmd_pub2)
  summary(hs_cougmd_hunt0)
  summary(hs_cougmd_hab1)
  
  
  #' #'  Note: using occupancy model formulas that drop Public covariate from predator sub-model
  #' (hs_cougmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_cougmd_fld <- fitList(hs_cougmd_trail, hs_cougmd_hunt, hs_cougmd_pub, hs_cougmd_pubhunt)
  #' #' Model selection
  #' modSel(hs_cougmd_fld)
  #' 
  #' (hs_cougmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_md_hunting_UMF, silent = TRUE)) 
  #' (hs_cougmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_md_hunting_UMF, silent = TRUE)) 
  #' (hs_cougmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_md_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_md_hunting_UMF, silent = TRUE)) 
  #' (hs_cougmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_md_hunting_UMF, silent = TRUE))
  #' (hs_cougmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_md_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_cougmd_fl <- fitList(hs_cougmd_null0, hs_cougmd_null1, hs_cougmd_hab0, hs_cougmd_hab1, hs_cougmd_hunt0, hs_cougmd_hunt1, hs_cougmd_pub0, hs_cougmd_pub1, hs_cougmd_pub2) 
  #' #' Model selection
  #' modSel(hs_cougmd_fl)
  #' summary(hs_cougmd_pub2) 
  #' summary(hs_cougmd_hunt0)
  #' summary(hs_cougmd_hunt1)

  ####  Cougar-Elk Hunting Season  ####
  (hs_cougelk_null <- occuMulti(detFormulas_null, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_hunting_UMF, silent = TRUE)) 
  (hs_cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_elk_hunting_UMF, silent = TRUE)) # not converging well 
  (hs_cougelk_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_elk_hunting_UMF, silent = TRUE)) # not converging well
  (hs_cougelk_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_cougelk_fl <- fitList(hs_cougelk_null, hs_cougelk_hab0, hs_cougelk_hunt0, hs_cougelk_hunt1) 
  #' Model selection
  modSel(hs_cougelk_fl)
  summary(hs_cougelk_hunt1)
  
  #' (hs_cougelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  #' (hs_cougelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  #' (hs_cougelk_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  #' (hs_cougelk_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_cougelk_fld <- fitList(hs_cougelk_trail, hs_cougelk_hunt, hs_cougelk_pub, hs_cougelk_pubhunt)
  #' #' Model selection
  #' modSel(hs_cougelk_fld)
  #' 
  #' (hs_cougelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  #' (hs_cougelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_hunting_UMF, silent = TRUE)) 
  #' (hs_cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, coug_elk_hunting_UMF, silent = TRUE))
  #' (hs_cougelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pubish0, coug_elk_hunting_UMF, silent = TRUE))
  #' (hs_cougelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pubish1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pubish2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' #' List of fitted models
  #' hs_cougelk_fl <- fitList(hs_cougelk_null0, hs_cougelk_hab0, hs_cougelk_hunt0, hs_cougelk_pub0) 
  #' #' Model selection
  #' modSel(hs_cougelk_fl)
  #' summary(hs_cougelk_pub0)
  
  ####  Cougar-White-tailed Deer Hunting Season  ####
  (hs_cougwtd_null <- occuMulti(detFormulas_null, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_wtd_hunting_UMF, silent = TRUE)) 
  (hs_cougwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_wtd_hunting_UMF, silent = TRUE)) 
  (hs_cougwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_wtd_hunting_UMF, silent = TRUE)) 
  (hs_cougwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_wtd_hunting_UMF, silent = TRUE)) # not converging well
  (hs_cougwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougwtd_fl <- fitList(hs_cougwtd_null, hs_cougwtd_hab0, hs_cougwtd_hab1, hs_cougwtd_pub0, hs_cougwtd_pub1, hs_cougwtd_pub2, hs_cougwtd_hunt0, hs_cougwtd_hunt1, hs_cougwtd_hunt2) 
  #' Model selection
  modSel(hs_cougwtd_fl)
  summary(hs_cougwtd_pub2)
  
  #' (hs_cougwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_cougwtd_fld <- fitList(hs_cougwtd_trail, hs_cougwtd_hunt, hs_cougwtd_pub, hs_cougwtd_pubhunt)
  #' #' Model selection
  #' modSel(hs_cougwtd_fld)
  #' 
  #' (hs_cougwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_cougwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_cougwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_wtd_hunting_UMF, silent = TRUE))
  #' (hs_cougwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_cougwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_wtd_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_cougwtd_fl <- fitList(hs_cougwtd_null0, hs_cougwtd_null1, hs_cougwtd_hab0, hs_cougwtd_hab1, hs_cougwtd_hunt0, hs_cougwtd_hunt1, hs_cougwtd_hunt2, hs_cougwtd_pub0, hs_cougwtd_pub1, hs_cougwtd_pub2) 
  #' #' Model selection
  #' modSel(hs_cougwtd_fl)
  #' summary(hs_cougwtd_pub2)
  #' summary(hs_cougwtd_pub0)
  
  ####  Cougar-Moose Hunting Season  ####
  (hs_cougmoose_null <- occuMulti(detFormulas_null, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_moose_hunting_UMF, silent = TRUE)) 
  (hs_cougmoose_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_moose_hunting_UMF, silent = TRUE)) # not converging well
  (hs_cougmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_moose_hunting_UMF, silent = TRUE)) # not converging well
  (hs_cougmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_cougmoose_fl <- fitList(hs_cougmoose_null, hs_cougmoose_hab0, hs_cougmoose_hab1, hs_cougmoose_pub0, hs_cougmoose_pub1, hs_cougmoose_hunt0) 
  #' Model selection
  modSel(hs_cougmoose_fl)
  summary(hs_cougmoose_hunt0)
  summary(hs_cougmoose_pub1) # f12 not significant
  
  #' (hs_cougmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_cougmoose_fld <- fitList(hs_cougmoose_trail, hs_cougmoose_hunt, hs_cougmoose_pub, hs_cougmoose_pubhunt)
  #' #' Model selection
  #' modSel(hs_cougmoose_fld)
  #' 
  #' (hs_cougmoose_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_moose_hunting_UMF, silent = TRUE)) 
  #' (hs_cougmoose_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_cougmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_moose_hunting_UMF, silent = TRUE))
  #' (hs_cougmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_moose_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_cougmoose_fl <- fitList(hs_cougmoose_null0, hs_cougmoose_null1, hs_cougmoose_hab0, hs_cougmoose_hunt0, hs_cougmoose_pub0, hs_cougmoose_pub1, hs_cougmoose_pub2) 
  #' #' Model selection
  #' modSel(hs_cougmoose_fl)
  #' summary(hs_cougmoose_hunt0)
  

  ####  Wolf-Mule Deer Hunting Season  ####
  (hs_wolfmd_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_md_hunting_UMF, silent = TRUE)) 
  (hs_wolfmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_md_hunting_UMF, silent = TRUE)) 
  (hs_wolfmd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_md_hunting_UMF, silent = TRUE)) # FAIL
  (hs_wolfmd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_wolfmd_fl <- fitList(hs_wolfmd_null, hs_wolfmd_hab0, hs_wolfmd_hab1, hs_wolfmd_pub0, hs_wolfmd_pub1, hs_wolfmd_hunt0, hs_wolfmd_hunt1, hs_wolfmd_hunt2) 
  #' Model selection
  modSel(hs_wolfmd_fl)
  summary(hs_wolfmd_hunt1)
  summary(hs_wolfmd_hunt2)
  
  #' #'  Note: using detection & occupancy model formulas that drop Public covariate from predator sub-models
  #' (hs_wolfmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_wolfmd_fld <- fitList(hs_wolfmd_trail, hs_wolfmd_hunt, hs_wolfmd_pub, hs_wolfmd_pubhunt) 
  #' #' Model selection
  #' modSel(hs_wolfmd_fld)
  #' 
  #' (hs_wolfmd_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_md_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfmd_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_md_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfmd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_md_hunting_UMF, silent = TRUE))
  #' (hs_wolfmd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_md_hunting_UMF, silent = TRUE)) # FAIL
  #' #' List of fitted models
  #' hs_wolfmd_fl <- fitList(hs_wolfmd_null0, hs_wolfmd_null1, hs_wolfmd_hab0, hs_wolfmd_hab1, hs_wolfmd_hunt0, hs_wolfmd_hunt1, hs_wolfmd_hunt2, hs_wolfmd_pub0, hs_wolfmd_pub1) 
  #' #' Model selection
  #' modSel(hs_wolfmd_fl)
  #' summary(hs_wolfmd_hunt1)
  #' summary(hs_wolfmd_hunt2)
  
  
  ####  Wolf-Elk Hunting Season  ####
  (hs_wolfelk_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_hunting_UMF, silent = TRUE)) 
  (hs_wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_hunting_UMF, silent = TRUE)) 
  (hs_wolfelk_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_wolfelk_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_wolfelk_fl <- fitList(hs_wolfelk_null, hs_wolfelk_hab0, hs_wolfelk_hab1, hs_wolfelk_pub0, hs_wolfelk_pub1, hs_wolfelk_hunt0, hs_wolfelk_hunt1, hs_wolfelk_hunt2) 
  #' Model selection
  modSel(hs_wolfelk_fl)
  summary(hs_wolfelk_hunt0)
  
  #' (hs_wolfelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_wolfelk_fld <- fitList(hs_wolfelk_trail, hs_wolfelk_hunt, hs_wolfelk_pub, hs_wolfelk_pubhunt) 
  #' #' Model selection
  #' modSel(hs_wolfelk_fld)
  #' 
  #' (hs_wolfelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pubish0, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pubish1, wolf_elk_hunting_UMF, silent = TRUE))
  #' (hs_wolfelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pubish2, wolf_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' #' List of fitted models
  #' hs_wolfelk_fl <- fitList(hs_wolfelk_null0, hs_wolfelk_null1, hs_wolfelk_hab0, hs_wolfelk_hab1, hs_wolfelk_hunt0, hs_wolfelk_hunt1, hs_wolfelk_hunt2, hs_wolfelk_pub0, hs_wolfelk_pub1) 
  #' #' Model selection
  #' modSel(hs_wolfelk_fl)
  #' summary(hs_wolfelk_hunt0)
  
  ####  Wolf-White-tailed Deer Hunting Season  ####
  (hs_wolfwtd_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_wtd_hunting_UMF, silent = TRUE)) 
  (hs_wolfwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_wtd_hunting_UMF, silent = TRUE)) 
  (hs_wolfwtd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  (hs_wolfwtd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  (hs_wolfwtd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_wtd_hunting_UMF, silent = TRUE))  
  #' List of fitted models
  hs_wolfwtd_fl <- fitList(hs_wolfwtd_null, hs_wolfwtd_hab0, hs_wolfwtd_hab1, hs_wolfwtd_pub0, hs_wolfwtd_pub1, hs_wolfwtd_hunt0, hs_wolfwtd_hunt2) 
  #' Model selection
  modSel(hs_wolfwtd_fl)
  summary(hs_wolfwtd_hunt0)
  
  #' (hs_wolfwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_wolfwtd_fld <- fitList(hs_wolfwtd_trail, hs_wolfwtd_hunt, hs_wolfwtd_pub, hs_wolfwtd_pubhunt)
  #' #' Model selection
  #' modSel(hs_wolfwtd_fld)
  #' 
  #' (hs_wolfwtd_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfwtd_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfwtd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_wolfwtd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_wtd_hunting_UMF, silent = TRUE))  
  #' (hs_wolfwtd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_wtd_hunting_UMF, silent = TRUE))
  #' (hs_wolfwtd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  #' #' List of fitted models
  #' hs_wolfwtd_fl <- fitList(hs_wolfwtd_null0, hs_wolfwtd_null1, hs_wolfwtd_hab0, hs_wolfwtd_hab1, hs_wolfwtd_hunt0, hs_wolfwtd_hunt2, hs_wolfwtd_pub0, hs_wolfwtd_pub1) 
  #' #' Model selection
  #' modSel(hs_wolfwtd_fl)
  #' summary(hs_wolfwtd_hunt0)
  
  ####  Wolf-Moose Hunting Season  ####
  (hs_wolfmoose_null <- occuMulti(detFormulas_null, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_moose_hunting_UMF, silent = TRUE)) 
  (hs_wolfmoose_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_moose_hunting_UMF, silent = TRUE)) 
  (hs_wolfmoose_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_moose_hunting_UMF, silent = TRUE)) # FAIL
  (hs_wolfmoose_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_wolfmoose_fl <- fitList(hs_wolfmoose_null, hs_wolfmoose_hab0, hs_wolfmoose_hab1, hs_wolfmoose_pub0, hs_wolfmoose_pub1, hs_wolfmoose_hunt0, hs_wolfmoose_hunt1, hs_wolfmoose_hunt2)
  #' Model selection
  modSel(hs_wolfmoose_fl)
  summary(hs_wolfmoose_hunt2)
  
  #' (hs_wolfmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_wolfmoose_fld <- fitList(hs_wolfmoose_trail, hs_wolfmoose_hunt, hs_wolfmoose_pub, hs_wolfmoose_pubhunt)
  #' #' Model selection
  #' modSel(hs_wolfmoose_fld)
  #' 
  #' (hs_wolfmoose_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_moose_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfmoose_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_moose_hunting_UMF, silent = TRUE)) 
  #' (hs_wolfmoose_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_moose_hunting_UMF, silent = TRUE))
  #' (hs_wolfmoose_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_moose_hunting_UMF, silent = TRUE)) # FAIL
  #' #' List of fitted models
  #' hs_wolfmoose_fl <- fitList(hs_wolfmoose_null0, hs_wolfmoose_null1, hs_wolfmoose_hab0, hs_wolfmoose_hab1, hs_wolfmoose_hunt0, hs_wolfmoose_hunt1, hs_wolfmoose_hunt2, hs_wolfmoose_pub0, hs_wolfmoose_pub1)
  #' #' Model selection
  #' modSel(hs_wolfmoose_fl)
  #' summary(hs_wolfmoose_hunt2)
  
  
  ####  Bear-Mule Deer Hunting Season  ####
  (hs_bearmd_null <- occuMulti(detFormulas_null, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_md_hunting_UMF, silent = TRUE)) 
  (hs_bearmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_md_hunting_UMF, silent = TRUE)) 
  (hs_bearmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_md_hunting_UMF, silent = TRUE)) 
  (hs_bearmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearmd_fl <- fitList(hs_bearmd_null, hs_bearmd_hab0, hs_bearmd_hab1, hs_bearmd_pub0, hs_bearmd_pub1, hs_bearmd_pub2, hs_bearmd_hunt0, hs_bearmd_hunt1, hs_bearmd_hunt2) 
  #' Model selection
  modSel(hs_bearmd_fl)
  summary(hs_bearmd_hunt0)
  summary(hs_bearmd_hunt1)
  summary(hs_bearmd_hunt2)
  summary(hs_bearmd_pub0)
  summary(hs_bearmd_pub1)
  
  #' (hs_bearmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bearmd_fld <- fitList(hs_bearmd_trail, hs_bearmd_hunt, hs_bearmd_pub, hs_bearmd_pubhunt)
  #' #' Model selection
  #' modSel(hs_bearmd_fld)
  #' 
  #' (hs_bearmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bearmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bearmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bearmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_md_hunting_UMF, silent = TRUE))
  #' (hs_bearmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_md_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bearmd_fl <- fitList(hs_bearmd_null0, hs_bearmd_null1, hs_bearmd_hab0, hs_bearmd_hab1, hs_bearmd_hunt0, hs_bearmd_hunt1, hs_bearmd_hunt2, hs_bearmd_pub0, hs_bearmd_pub1, hs_bearmd_pub2) 
  #' #' Model selection
  #' modSel(hs_bearmd_fl)
  #' summary(hs_bearmd_hab0)
  #' summary(hs_bearmd_hab1)
  
  ####  Bear-Elk Hunting Season  ####
  (hs_bearelk_null <- occuMulti(detFormulas_null, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_hunting_UMF, silent = TRUE)) 
  (hs_bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_hunting_UMF, silent = TRUE)) 
  (hs_bearelk_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_elk_hunting_UMF, silent = TRUE)) #FAIL
  (hs_bearelk_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearelk_fl <- fitList(hs_bearelk_null, hs_bearelk_hab0, hs_bearelk_hab1, hs_bearelk_pub0, hs_bearelk_pub1, hs_bearelk_hunt0, hs_bearelk_hunt1, hs_bearelk_hunt2)
  #' Model selection
  modSel(hs_bearelk_fl)
  summary(hs_bearelk_hab1) 
  
  #' (hs_bearelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bearelk_fld <- fitList(hs_bearelk_trail, hs_bearelk_hunt, hs_bearelk_pub, hs_bearelk_pubhunt)
  #' #' Model selection
  #' modSel(hs_bearelk_fld)
  #' 
  #' (hs_bearelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_hunting_UMF, silent = TRUE)) 
  #' (hs_bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_hunting_UMF, silent = TRUE)) 
  #' (hs_bearelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pub0, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pub1, bear_elk_hunting_UMF, silent = TRUE))
  #' (hs_bearelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pub2, bear_elk_hunting_UMF, silent = TRUE)) #FAIL
  #' #' List of fitted models
  #' hs_bearelk_fl <- fitList(hs_bearelk_null0, hs_bearelk_null1, hs_bearelk_hab0, hs_bearelk_hab1, hs_bearelk_hunt0, hs_bearelk_hunt1, hs_bearelk_hunt2, hs_bearelk_pub0, hs_bearelk_pub1)
  #' #' Model selection
  #' modSel(hs_bearelk_fl)
  #' summary(hs_bearelk_pub1) # f12 not significant
  #' summary(hs_bearelk_hab1) # f12 not significant
  
  
  ####  Bear-White-tailed Deer Hunting Season  ####
  (hs_bearwtd_null <- occuMulti(detFormulas_null, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearwtd_fl <- fitList(hs_bearwtd_null, hs_bearwtd_hab0, hs_bearwtd_hab1, hs_bearwtd_pub0, hs_bearwtd_pub1, hs_bearwtd_pub2, hs_bearwtd_hunt0, hs_bearwtd_hunt1, hs_bearwtd_hunt2) 
  #' Model selection
  modSel(hs_bearwtd_fl)
  summary(hs_bearwtd_pub0)
  summary(hs_bearwtd_pub1)
  
  #' (hs_bearwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bearwtd_fld <- fitList(hs_bearwtd_trail, hs_bearwtd_hunt, hs_bearwtd_pub, hs_bearwtd_pubhunt)
  #' #' Model selection
  #' modSel(hs_bearwtd_fld)
  #' 
  #' (hs_bearwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bearwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bearwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bearwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bearwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bearwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_wtd_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bearwtd_fl <- fitList(hs_bearwtd_null0, hs_bearwtd_null1, hs_bearwtd_hab0, hs_bearwtd_hab1, hs_bearwtd_hunt0, hs_bearwtd_hunt1, hs_bearwtd_hunt2, hs_bearwtd_pub0, hs_bearwtd_pub1, hs_bearwtd_pub2) 
  #' #' Model selection
  #' modSel(hs_bearwtd_fl)
  #' summary(hs_bearwtd_pub0)
  #' summary(hs_bearwtd_hab0)
  
  ####  Bear-Moose Hunting Season  ####
  (hs_bearmoose_null <- occuMulti(detFormulas_null, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_moose_hunting_UMF, silent = TRUE)) 
  (hs_bearmoose_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_moose_hunting_UMF, silent = TRUE)) 
  (hs_bearmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  (hs_bearmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_moose_hunting_UMF, silent = TRUE)) 
  (hs_bearmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  hs_bearmoose_fl <- fitList(hs_bearmoose_null, hs_bearmoose_hab0, hs_bearmoose_hab1, hs_bearmoose_pub0, hs_bearmoose_pub1, hs_bearmoose_hunt0, hs_bearmoose_hunt2) 
  #' Model selection
  modSel(hs_bearmoose_fl)
  summary(hs_bearmoose_hunt0)
  
  #' (hs_bearmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bearmoose_fld <- fitList(hs_bearmoose_trail, hs_bearmoose_hunt, hs_bearmoose_pub, hs_bearmoose_pubhunt)
  #' #' Model selection
  #' modSel(hs_bearmoose_fld)
  #' 
  #' (hs_bearmoose_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_moose_hunting_UMF, silent = TRUE)) 
  #' (hs_bearmoose_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_moose_hunting_UMF, silent = TRUE)) 
  #' (hs_bearmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  #' (hs_bearmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_moose_hunting_UMF, silent = TRUE)) 
  #' (hs_bearmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_moose_hunting_UMF, silent = TRUE))
  #' (hs_bearmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  #' #' List of fitted models
  #' hs_bearmoose_fl <- fitList(hs_bearmoose_null0, hs_bearmoose_null1, hs_bearmoose_hab0, hs_bearmoose_hab1, hs_bearmoose_hunt0, hs_bearmoose_hunt2, hs_bearmoose_pub0, hs_bearmoose_pub1) 
  #' #' Model selection
  #' modSel(hs_bearmoose_fl)
  #' summary(hs_bearmoose_hunt0)
  #' summary(hs_bearmoose_hab0)

  
  ####  Coyote-Mule Deer Hunting Season  ####
  (hs_coymd_null <- occuMulti(detFormulas_null, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_md_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_coymd_fl <- fitList(hs_coymd_null, hs_coymd_hab0, hs_coymd_hab1, hs_coymd_pub0, hs_coymd_pub1, hs_coymd_pub2, hs_coymd_hunt0, hs_coymd_hunt1)
  #' Model selection
  modSel(hs_coymd_fl)
  summary(hs_coymd_hunt0)
  summary(hs_coymd_hunt1)
  
  #' (hs_coymd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_hunting_UMF, silent = TRUE))
  #' (hs_coymd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_md_hunting_UMF, silent = TRUE))
  #' (hs_coymd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  #' (hs_coymd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  #' #' List of fitted models
  #' hs_coymd_fld <- fitList(hs_coymd_trail, hs_coymd_hunt, hs_coymd_pub, hs_coymd_pubhunt)
  #' #' Model selection
  #' modSel(hs_coymd_fld)
  #' 
  #' (hs_coymd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  #' (hs_coymd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_md_hunting_UMF, silent = TRUE))
  #' (hs_coymd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coy_md_hunting_UMF, silent = TRUE)) 
  #' (hs_coymd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coy_md_hunting_UMF, silent = TRUE))
  #' (hs_coymd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_md_hunting_UMF, silent = TRUE)) 
  #' (hs_coymd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_md_hunting_UMF, silent = TRUE))
  #' (hs_coymd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_md_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_coymd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_md_hunting_UMF, silent = TRUE)) 
  #' (hs_coymd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_md_hunting_UMF, silent = TRUE))
  #' (hs_coymd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_md_hunting_UMF, silent = TRUE))
  #' # (coymd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, coy_md_hunting_UMF, silent = TRUE)) 
  #' # (coymd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coy_md_hunting_UMF, silent = TRUE))
  #' # (coymd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coy_md_hunting_UMF, silent = TRUE)) #FAIL
  #' #' List of fitted models
  #' hs_coymd_fl <- fitList(hs_coymd_null0, hs_coymd_null1, hs_coymd_hab0, hs_coymd_hab1, hs_coymd_hunt0, hs_coymd_hunt1, hs_coymd_pub0, hs_coymd_pub1, hs_coymd_pub2)
  #' #' Model selection
  #' modSel(hs_coymd_fl)
  #' summary(hs_coymd_hunt0)
  #' summary(hs_coymd_hunt1)
  
  ####  Coyote-White-tailed Deer Hunting Season  ####
  (hs_coywtd_null <- occuMulti(detFormulas_null, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_wtd_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_coywtd_fl <- fitList(hs_coywtd_null, hs_coywtd_hab0, hs_coywtd_hab1, hs_coywtd_pub0, hs_coywtd_pub1, hs_coywtd_pub2, hs_coywtd_hunt0, hs_coywtd_hunt1)
  #' Model selection
  modSel(hs_coywtd_fl)
  summary(hs_coywtd_hunt1)
  
  #' (hs_coywtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE))
  #' (hs_coywtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE))
  #' (hs_coywtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_coywtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' #' List of fitted models
  #' hs_coywtd_fld <- fitList(hs_coywtd_trail, hs_coywtd_hunt, hs_coywtd_pub, hs_coywtd_pubhunt)
  #' #' Model selection
  #' modSel(hs_coywtd_fld)
  #' 
  #' (hs_coywtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_coywtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_wtd_hunting_UMF, silent = TRUE))
  #' (hs_coywtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_coywtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coy_wtd_hunting_UMF, silent = TRUE))
  #' (hs_coywtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_coywtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_wtd_hunting_UMF, silent = TRUE))
  #' (hs_coywtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_wtd_hunting_UMF, silent = TRUE)) # FAIL
  #' (hs_coywtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_coywtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_wtd_hunting_UMF, silent = TRUE))
  #' (hs_coywtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_wtd_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_coywtd_fl <- fitList(hs_coywtd_null0, hs_coywtd_null1, hs_coywtd_hab0, hs_coywtd_hab1, hs_coywtd_hunt0, hs_coywtd_hunt1, hs_coywtd_pub0, hs_coywtd_pub1, hs_coywtd_pub2)
  #' #' Model selection
  #' modSel(hs_coywtd_fl)
  #' summary(hs_coywtd_hunt1)
  
  
  ####  Bobcat-Mule Deer Hunting Season  ####
  (hs_bobmd_null <- occuMulti(detFormulas_null, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_md_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  hs_bobmd_fl <- fitList(hs_bobmd_null, hs_bobmd_hab0, hs_bobmd_hab1, hs_bobmd_pub0, hs_bobmd_pub1, hs_bobmd_pub2, hs_bobmd_hunt0, hs_bobmd_hunt1, hs_bobmd_hunt2)
  #' Model selection
  modSel(hs_bobmd_fl)
  summary(hs_bobmd_hunt0)
  summary(hs_bobmd_hunt1)
  summary(hs_bobmd_pub0)
  
  #' (hs_bobmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_hunting_UMF, silent = TRUE))
  #' (hs_bobmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_md_hunting_UMF, silent = TRUE))
  #' (hs_bobmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bobmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  #' #' List of fitted models
  #' hs_bobmd_fld <- fitList(hs_bobmd_trail, hs_bobmd_hunt, hs_bobmd_pub, hs_bobmd_pubhunt)
  #' #' Model selection
  #' modSel(hs_bobmd_fld)
  #' 
  #' (hs_bobmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bobmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_md_hunting_UMF, silent = TRUE))
  #' (hs_bobmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bob_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bobmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bob_md_hunting_UMF, silent = TRUE))
  #' (hs_bobmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bobmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bobmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bobmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_md_hunting_UMF, silent = TRUE)) 
  #' (hs_bobmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_md_hunting_UMF, silent = TRUE))
  #' (hs_bobmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_md_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bobmd_fl <- fitList(hs_bobmd_null0, hs_bobmd_null1, hs_bobmd_hab0, hs_bobmd_hab1, hs_bobmd_hunt0, hs_bobmd_hunt1, hs_bobmd_hunt2, hs_bobmd_pub0, hs_bobmd_pub1, hs_bobmd_pub2)
  #' #' Model selection
  #' modSel(hs_bobmd_fl)
  #' summary(hs_bobmd_hunt0) # ok but the hunting activity effect on bobcat sub-model has unusually large effect size...
  #' summary(hs_bobmd_hunt1)
  #' summary(hs_bobmd_pub0)
  
  
  ####  Bobcat-White-tailed Deer Hunting Season  ####
  (hs_bobwtd_null <- occuMulti(detFormulas_null, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_wtd_hunting_UMF, silent = TRUE)) # not converging well
  (hs_bobwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bobwtd_fl <- fitList(hs_bobwtd_null, hs_bobwtd_hab0, hs_bobwtd_hab1, hs_bobwtd_pub0, hs_bobwtd_pub1, hs_bobwtd_pub2, hs_bobwtd_hunt0, hs_bobwtd_hunt1)
  #' Model selection
  modSel(hs_bobwtd_fl)
  summary(hs_bobwtd_pub1) # f12 not significant
  summary(hs_bobwtd_pub0)
  summary(hs_bobwtd_pub2)
  summary(hs_bobwtd_hunt1)
  
  
  #' (hs_bobwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bobwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bobwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bobwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' #' List of fitted models
  #' hs_bobwtd_fld <- fitList(hs_bobwtd_trail, hs_bobwtd_hunt, hs_bobwtd_pub, hs_bobwtd_pubhunt)
  #' #' Model selection
  #' modSel(hs_bobwtd_fld)
  #' 
  #' (hs_bobwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bobwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bobwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bobwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bob_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bobwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bobwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bobwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_wtd_hunting_UMF, silent = TRUE)) # FAIL massive effect sizes and sign switches on intercept
  #' (hs_bobwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' (hs_bobwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_wtd_hunting_UMF, silent = TRUE))
  #' (hs_bobwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_wtd_hunting_UMF, silent = TRUE))
  #' #' List of fitted models
  #' hs_bobwtd_fl <- fitList(hs_bobwtd_null0, hs_bobwtd_null1, hs_bobwtd_hab0, hs_bobwtd_hab1, hs_bobwtd_hunt0, hs_bobwtd_hunt1, hs_bobwtd_pub0, hs_bobwtd_pub1, hs_bobwtd_pub2)
  #' #' Model selection
  #' modSel(hs_bobwtd_fl)
  #' summary(hs_bobwtd_pub1) # f12 not significant
  #' summary(hs_bobwtd_pub0)
  
  
  #' Save model outputs in one giant R image
  save.image(file = "./Outputs/MultiSpp_CoOcc_Models_2spp.RData")
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info

  #'  Function to save occupancy results
  rounddig <- 2
  occ_out <- function(mod, spp1, spp2, season) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Season = rep(season, nrow(.)),
        Estimate = round(Estimate, rounddig),
        SE = round(SE, rounddig),
        z = round(z, rounddig),
        Pval = round(`P(>|z|)`, rounddig)
      ) %>%
      dplyr::select(-`P(>|z|)`) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species1, .before = Parameter) %>%
      relocate(Species2, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
      
    return(out)
  }
  
  #'  Run each model through function
  #'  Grazing season models
  occ_cougmd_grazing <- occ_out(gs_cougmd_hab0, "Cougar", "Mule Deer", "Grazing")
  occ_cougelk_grazing <- occ_out(gs_cougelk_hab0, "Cougar", "Elk", "Grazing")
  occ_cougwtd_grazing <- occ_out(gs_cougwtd_graze2, "Cougar", "White-tailed Deer", "Grazing")
  occ_cougmoose_grazing <- occ_out(gs_cougmoose_graze0, "Cougar", "Moose", "Grazing")
  occ_wolfmd_grazing <- occ_out(gs_wolfmd_hab0, "Wolf", "Mule Deer", "Grazing")
  occ_wolfelk_grazing <- occ_out(gs_wolfelk_hab0, "Wolf", "Elk", "Grazing")
  occ_wolfwtd_grazing <- occ_out(gs_wolfwtd_graze0, "Wolf", "White-tailed Deer", "Grazing")
  occ_wolfmoose_grazing <- occ_out(gs_wolfmoose_graze0, "Wolf", "Moose", "Grazing")
  occ_bearmd_grazing <- occ_out(gs_bearmd_hab1, "Black Bear", "Mule Deer", "Grazing")
  occ_bearelk_grazing <- occ_out(gs_bearelk_hab1, "Black Bear", "Elk", "Grazing")
  occ_bearwtd_grazing <- occ_out(gs_bearwtd_graze0, "Black Bear", "White-tailed Deer", "Grazing")
  occ_bearmoose_grazing <- occ_out(gs_bearmoose_graze0, "Black Bear", "Moose", "Grazing")
  occ_coymd_grazing <- occ_out(gs_coymd_graze2, "Coyote", "Mule Deer", "Grazing")
  occ_coyelk_grazing <- occ_out(gs_coywtd_graze2, "Coyote", "White-tailed Deer", "Grazing")
  occ_bobwtd_grazing <- occ_out(gs_bobmd_graze0, "Bobcat", "Mule Deer Deer", "Grazing")
  occ_bobmoose_grazing <- occ_out(gs_bobwtd_graze2, "Bobcat", "White-tailed Deer", "Grazing")
  #'  Hunting season models
  occ_cougmd_hunting <- occ_out(hs_cougmd_pub2, "Cougar", "Mule Deer", "Hunting")
  occ_cougelk_hunting <- occ_out(hs_cougelk_hunt1, "Cougar", "Elk", "Hunting")
  occ_cougwtd_hunting <- occ_out(hs_cougwtd_pub2, "Cougar", "White-tailed Deer", "Hunting")
  occ_cougmoose_hunting <- occ_out(hs_cougmoose_hunt0, "Cougar", "Moose", "Hunting")
  occ_wolfmd_hunting <- occ_out(hs_wolfmd_hunt1, "Wolf", "Mule Deer", "Hunting")
  occ_wolfelk_hunting <- occ_out(hs_wolfelk_hunt0, "Wolf", "Elk", "Hunting")
  occ_wolfwtd_hunting <- occ_out(hs_wolfwtd_hunt0, "Wolf", "White-tailed Deer", "Hunting")
  occ_wolfmoose_hunting <- occ_out(hs_wolfmoose_hunt2, "Wolf", "Moose", "Hunting")
  occ_bearmd_hunting <- occ_out(hs_bearmd_hunt0, "Black Bear", "Mule Deer", "Hunting")
  occ_bearelk_hunting <- occ_out(hs_bearelk_hab1, "Black Bear", "Elk", "Hunting")
  occ_bearwtd_hunting <- occ_out(hs_bearwtd_pub0, "Black Bear", "White-tailed Deer", "Hunting")
  occ_bearmoose_hunting <- occ_out(hs_bearmoose_hunt0, "Black Bear", "Moose", "Hunting")
  occ_coymd_hunting <- occ_out(hs_coymd_hunt0, "Coyote", "Mule Deer", "Hunting")
  occ_coyelk_hunting <- occ_out(hs_coywtd_hunt1, "Coyote", "White-tailed Deer", "Hunting")
  occ_bobwtd_hunting <- occ_out(hs_bobmd_hunt0, "Bobcat", "Mule Deer Deer", "Hunting")
  occ_bobmoose_hunting <- occ_out(hs_bobwtd_pub1, "Bobcat", "White-tailed Deer", "Hunting")
  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  graze_results <- rbind(occ_cougmd_grazing, occ_cougelk_grazing, occ_cougwtd_grazing, occ_cougmoose_grazing,
                     occ_wolfmd_grazing, occ_wolfelk_grazing, occ_wolfwtd_grazing, occ_wolfmoose_grazing,
                     occ_bearmd_grazing, occ_bearelk_grazing, occ_bearwtd_grazing, occ_bearmoose_grazing,
                     occ_coymd_grazing, occ_coyelk_grazing, occ_bobwtd_grazing, occ_bobmoose_grazing)
  hunt_results <- rbind(occ_cougmd_hunting, occ_cougelk_hunting, occ_cougwtd_hunting, occ_cougmoose_hunting,
                    occ_wolfmd_hunting, occ_wolfelk_hunting, occ_wolfwtd_hunting, occ_wolfmoose_hunting,
                    occ_bearmd_hunting, occ_bearelk_hunting, occ_bearwtd_hunting, occ_bearmoose_hunting,
                    occ_coymd_hunting, occ_coyelk_hunting, occ_bobwtd_hunting, occ_bobmoose_hunting)
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  #'  Grazing results tables
  results_graze <- graze_results %>%  
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) 
  results_graze_wide <- results_graze %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    #'  Change species names to general classes
    mutate(
      Parameter = gsub("cougar", "Species 1", Parameter),
      Parameter = gsub("wolf", "Species 1", Parameter),
      Parameter = gsub("blackbear", "Species 1", Parameter),
      Parameter = gsub("coyote", "Species 1", Parameter),
      Parameter = gsub("bobcat", "Species 1", Parameter),
      Parameter = gsub("muledeer", "Species 2", Parameter),
      Parameter = gsub("mule_deer", "Species 2", Parameter),
      Parameter = gsub("elk", "Species 2", Parameter),
      Parameter = gsub("wtd", "Species 2", Parameter),
      Parameter = gsub("moose", "Species 2", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Study_AreaOK") %>%
    relocate("[Species 1:Species 2] GrazingActivity", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1] GrazingActivity", c("[Species 1] GrazingActivity (SE)", "[Species 1] GrazingActivity Pval"), sep = "_") %>%
    separate("[Species 2] GrazingActivity", c("[Species 2] GrazingActivity (SE)", "[Species 2] GrazingActivity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] GrazingActivity", c("[Species 1:Species 2] GrazingActivity (SE)", "[Species 1:Species 2] GrazingActivity Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
    
  #'  Hunting results tables
  results_hunt <- hunt_results %>%  
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) 
  results_hunt_wide <- results_hunt %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    #'  Change species names to general classes
    mutate(
      Parameter = gsub("cougar", "Species 1", Parameter),
      Parameter = gsub("wolf", "Species 1", Parameter),
      Parameter = gsub("blackbear", "Species 1", Parameter),
      Parameter = gsub("coyote", "Species 1", Parameter),
      Parameter = gsub("bobcat", "Species 1", Parameter),
      Parameter = gsub("muledeer", "Species 2", Parameter),
      Parameter = gsub("mule_deer", "Species 2", Parameter),
      Parameter = gsub("elk", "Species 2", Parameter),
      Parameter = gsub("wtd", "Species 2", Parameter),
      Parameter = gsub("moose", "Species 2", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 2] Study_AreaOK") %>%
    relocate("[Species 1:Species 2] HuntingActivity", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] Public1", .after = "[Species 1:Species 2] HuntingActivity") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1] HuntingActivity", c("[Species 1] HuntingActivity (SE)", "[Species 1] HuntingActivity Pval"), sep = "_") %>%
    separate("[Species 2] HuntingActivity", c("[Species 2] HuntingActivity (SE)", "[Species 2] HuntingActivity Pval"), sep = "_") %>%
    separate("[Species 1] Public1", c("[Species 1] Public1 (SE)", "[Species 1] Public1 Pval"), sep = "_") %>%
    separate("[Species 2] Public1", c("[Species 2] Public1 (SE)", "[Species 2] Public1 Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] HuntingActivity", c("[Species 1:Species 2] HuntingActivity (SE)", "[Species 1:Species 2] HuntingActivity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Public1", c("[Species 1:Species 2] Public1 (SE)", "[Species 1:Species 2] Public1 Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  #'  Save!
  write.csv(results_graze, paste0("./Outputs/Tables/CoOcc_OccProb_GrazingResults_", Sys.Date(), ".csv"))  
  write.csv(results_graze_wide, paste0("./Outputs/Tables/CoOcc_OccProb_GrazingResults_wide_", Sys.Date(), ".csv"))
  write.csv(results_hunt, paste0("./Outputs/Tables/CoOcc_OccProb_HuntingResults_", Sys.Date(), ".csv"))  
  write.csv(results_hunt_wide, paste0("./Outputs/Tables/CoOcc_OccProb_HuntingResults_wide_", Sys.Date(), ".csv"))
  
 
  #'  Function to save detection results
  det_out <- function(mod, spp1, spp2, season) {
    out <- summary(mod@estimates)$det %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$det),
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Season = rep(season, nrow(.)),
        Estimate = round(Estimate, 2),
        SE = round(SE, 2),
        z = round(z, 2),
        Pval = round(`P(>|z|)`, 2)
      ) %>%
      dplyr::select(-`P(>|z|)`) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species1, .before = Parameter) %>%
      relocate(Species2, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
    return(out)
  }
  
  
  #'  Run each model through detection function
  #'  Grazing season models
  det_cougmd_grazing <- det_out(gs_cougmd_hab0, "Cougar", "Mule Deer", "Grazing")
  det_cougelk_grazing <- det_out(gs_cougelk_hab0, "Cougar", "Elk", "Grazing")
  det_cougwtd_grazing <- det_out(gs_cougwtd_graze2, "Cougar", "White-tailed Deer", "Grazing")
  det_cougmoose_grazing <- det_out(gs_cougmoose_graze0, "Cougar", "Moose", "Grazing")
  det_wolfmd_grazing <- det_out(gs_wolfmd_hab0, "Wolf", "Mule Deer", "Grazing")
  det_wolfelk_grazing <- det_out(gs_wolfelk_hab0, "Wolf", "Elk", "Grazing")
  det_wolfwtd_grazing <- det_out(gs_wolfwtd_graze0, "Wolf", "White-tailed Deer", "Grazing")
  det_wolfmoose_grazing <- det_out(gs_wolfmoose_graze0, "Wolf", "Moose", "Grazing")
  det_bearmd_grazing <- det_out(gs_bearmd_hab1, "Black Bear", "Mule Deer", "Grazing")
  det_bearelk_grazing <- det_out(gs_bearelk_hab1, "Black Bear", "Elk", "Grazing")
  det_bearwtd_grazing <- det_out(gs_bearwtd_graze0, "Black Bear", "White-tailed Deer", "Grazing")
  det_bearmoose_grazing <- det_out(gs_bearmoose_graze0, "Black Bear", "Moose", "Grazing")
  det_coymd_grazing <- det_out(gs_coymd_graze2, "Coyote", "Mule Deer", "Grazing")
  det_coyelk_grazing <- det_out(gs_coywtd_graze2, "Coyote", "White-tailed Deer", "Grazing")
  det_bobwtd_grazing <- det_out(gs_bobmd_graze0, "Bobcat", "Mule Deer Deer", "Grazing")
  det_bobmoose_grazing <- det_out(gs_bobwtd_graze2, "Bobcat", "White-tailed Deer", "Grazing")
  #'  Hunting season models
  det_cougmd_hunting <- det_out(hs_cougmd_pub2, "Cougar", "Mule Deer", "Hunting")
  det_cougelk_hunting <- det_out(hs_cougelk_hunt1, "Cougar", "Elk", "Hunting")
  det_cougwtd_hunting <- det_out(hs_cougwtd_pub2, "Cougar", "White-tailed Deer", "Hunting")
  det_cougmoose_hunting <- det_out(hs_cougmoose_hunt0, "Cougar", "Moose", "Hunting")
  det_wolfmd_hunting <- det_out(hs_wolfmd_hunt1, "Wolf", "Mule Deer", "Hunting")
  det_wolfelk_hunting <- det_out(hs_wolfelk_hunt0, "Wolf", "Elk", "Hunting")
  det_wolfwtd_hunting <- det_out(hs_wolfwtd_hunt0, "Wolf", "White-tailed Deer", "Hunting")
  det_wolfmoose_hunting <- det_out(hs_wolfmoose_hunt2, "Wolf", "Moose", "Hunting")
  det_bearmd_hunting <- det_out(hs_bearmd_hunt0, "Black Bear", "Mule Deer", "Hunting")
  det_bearelk_hunting <- det_out(hs_bearelk_hab1, "Black Bear", "Elk", "Hunting")
  det_bearwtd_hunting <- det_out(hs_bearwtd_pub0, "Black Bear", "White-tailed Deer", "Hunting")
  det_bearmoose_hunting <- det_out(hs_bearmoose_hunt0, "Black Bear", "Moose", "Hunting")
  det_coymd_hunting <- det_out(hs_coymd_hunt0, "Coyote", "Mule Deer", "Hunting")
  det_coyelk_hunting <- det_out(hs_coywtd_hunt1, "Coyote", "White-tailed Deer", "Hunting")
  det_bobwtd_hunting <- det_out(hs_bobmd_hunt0, "Bobcat", "Mule Deer Deer", "Hunting")
  det_bobmoose_hunting <- det_out(hs_bobwtd_pub1, "Bobcat", "White-tailed Deer", "Hunting")

  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  graze_det_results <- rbind(det_cougmd_grazing, det_cougelk_grazing, det_cougwtd_grazing, det_cougmoose_grazing,
                         det_wolfmd_grazing, det_wolfelk_grazing, det_wolfwtd_grazing, det_wolfmoose_grazing,
                         det_bearmd_grazing, det_bearelk_grazing, det_bearwtd_grazing, det_bearmoose_grazing,
                         det_coymd_grazing, det_coyelk_grazing, det_bobwtd_grazing, det_bobmoose_grazing)
  hunt_det_results <- rbind(det_cougmd_hunting, det_cougelk_hunting, det_cougwtd_hunting, det_cougmoose_hunting,
                        det_wolfmd_hunting, det_wolfelk_hunting, det_wolfwtd_hunting, det_wolfmoose_hunting,
                        det_bearmd_hunting, det_bearelk_hunting, det_bearwtd_hunting, det_bearmoose_hunting,
                        det_coymd_hunting, det_coyelk_hunting, det_bobwtd_hunting, det_bobmoose_hunting)

  #'  Round so numbers are a little easier to interpret
  results_det_graze <- graze_det_results %>%
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    )
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_det_graze_wide <- results_det_graze %>% 
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    #'  Change species names to general classes
    mutate(
      Parameter = gsub("cougar", "Species 1", Parameter),
      Parameter = gsub("wolf", "Species 1", Parameter),
      Parameter = gsub("blackbear", "Species 1", Parameter),
      Parameter = gsub("coyote", "Species 1", Parameter),
      Parameter = gsub("bobcat", "Species 1", Parameter),
      Parameter = gsub("muledeer", "Species 2", Parameter),
      Parameter = gsub("mule_deer", "Species 2", Parameter),
      Parameter = gsub("elk", "Species 2", Parameter),
      Parameter = gsub("wtd", "Species 2", Parameter),
      Parameter = gsub("moose", "Species 2", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 1] WeeklyGrazing", c("[Species 1] WeeklyGrazing (SE)", "[Species 1] WeeklyGrazing Pval"), sep = "_") %>%
    separate("[Species 2] WeeklyGrazing", c("[Species 2] WeeklyGrazing (SE)", "[Species 2] WeeklyGrazing Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 

  #'  Hunting season detection results
  results_det_hunt <- hunt_det_results %>%
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    )
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_det_hunt_wide <- results_det_hunt %>% 
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    #'  Change species names to general classes
    mutate(
      Parameter = gsub("cougar", "Species 1", Parameter),
      Parameter = gsub("wolf", "Species 1", Parameter),
      Parameter = gsub("blackbear", "Species 1", Parameter),
      Parameter = gsub("coyote", "Species 1", Parameter),
      Parameter = gsub("bobcat", "Species 1", Parameter),
      Parameter = gsub("muledeer", "Species 2", Parameter),
      Parameter = gsub("mule_deer", "Species 2", Parameter),
      Parameter = gsub("elk", "Species 2", Parameter),
      Parameter = gsub("wtd", "Species 2", Parameter),
      Parameter = gsub("moose", "Species 2", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 1] Public1", c("[Species 1] Public1 (SE)", "[Species 1] Public1 Pval"), sep = "_") %>%
    separate("[Species 2] Public1", c("[Species 2] Public1 (SE)", "[Species 2] Public1 Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  #'  Save!
  write.csv(results_det_graze, paste0("./Outputs/Tables/CoOcc_DetProb_GrazingResults_", Sys.Date(), ".csv"))  #'  KEEP TRACK of whether human mod was excluded from models
  write.csv(results_det_graze_wide, paste0("./Outputs/Tables/CoOcc_DetProb_GrazingResults_wide", Sys.Date(), ".csv"))
  write.csv(results_det_hunt, paste0("./Outputs/Tables/CoOcc_DetProb_HuntingResults_", Sys.Date(), ".csv"))  #'  KEEP TRACK of whether human mod was excluded from models
  write.csv(results_det_hunt_wide, paste0("./Outputs/Tables/CoOcc_DetProb_HuntingResults_wide", Sys.Date(), ".csv"))
  

  #' #'  Predict probability of occupancy across sites
  #' mu_occ <- function(mod, species, season) {
  #'   #'  Predict occupancy probability for all camera sties
  #'   occu_mean <- predict(object = mod, type = "state") %>%  # do I provide newdata = sitecovariates if it's the same values as what went into the model?
  #'     #'  Average occupancy probabilities across sites for mean psi
  #'     summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE)
  #'   #'  Predict occupancy probability for all camera sties
  #'   det_mean <- predict(object = mod, type = "det") %>%
  #'     #'  Average occupancy probabilities across sites for mean psi
  #'     summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE)
  #'   predicted <- as.data.frame(rbind(occu_mean, det_mean))
  #'   colnames(predicted) <- c("Mean", "SE")
  #'   Parameter <- c("Occupancy", "Detection")
  #'   Species <- species
  #'   Season <- season
  #'   predicted <- cbind(predicted, Parameter)
  #'   predicted <- cbind(predicted, Species)
  #'   predicted <- cbind(predicted, Season)
  #'   return(predicted)
  #' }
  #' #'  Estimate mean probability of occupancy and detection per species and season
  #' md_predict_smr <- mu_occ(md_s1819_global, "Mule Deer", "Summer")
  #' md_predict_wtr <- mu_occ(md_w1820_global, "Mule Deer", "Winter")
  #' elk_predict_smr <- mu_occ(elk_s1819_global2, "Elk", "Summer")
  #' elk_predict_wtr <- mu_occ(elk_w1820_global2, "Elk", "Winter")
  #' wtd_predict_smr <- mu_occ(wtd_s1819_global2, "White-tailed Deer", "Summer")
  #' wtd_predict_wtr <- mu_occ(wtd_w1820_global2, "White-tailed Deer", "Winter")
  #' coug_predict_smr <- mu_occ(coug_s1819_global, "Cougar", "Summer")
  #' coug_predict_wtr <- mu_occ(coug_w1820_global, "Cougar", "Winter")
  #' wolf_predict_smr <- mu_occ(wolf_s1819_global2, "Wolf", "Summer")
  #' wolf_predict_wtr <- mu_occ(wolf_w1820_global2, "Wolf", "Winter")
  #' bob_predict_smr <- mu_occ(bob_s1819_global, "Bobcat", "Summer")
  #' bob_predict_wtr <- mu_occ(bob_w1820_global, "Bobcat", "Winter")
  #' coy_predict_smr <- mu_occ(coy_s1819_global, "Coyote", "Summer")
  #' coy_predict_wtr <- mu_occ(coy_w1820_global, "Coyote", "Winter")
  #' 
  #' #'  Make a pretty table
  #' Mean_tbl <- bind_rows(md_predict_smr, md_predict_wtr, elk_predict_smr, elk_predict_wtr, wtd_predict_smr,
  #'             wtd_predict_wtr, coug_predict_smr, coug_predict_wtr, wolf_predict_smr,
  #'             wolf_predict_wtr, bob_predict_smr, bob_predict_wtr, coy_predict_smr,
  #'             coy_predict_wtr) %>%
  #'   relocate(Species, .before = Mean) %>%
  #'   relocate(Season, .after = Species) %>%
  #'   relocate(Parameter, .after = Season) %>%
  #'   arrange(Parameter, Mean, Species)
  #' 
  #' #'  Save
  #' write.csv(Mean_tbl, paste0("./Outputs/Tables/OccMod_Mean_Estimates_NoHM_", Sys.Date(), ".csv"))  #"  Keep track of whether human mod was excluded from analyses

  