  #'  ============================================
  #'  Multi-Species Co-Occurrence Occupancy Models 
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
  #'  and Cattle_Hunter_Activity.R scripts. Covariate data included in occupancy 
  #'  models were collected at each camera site or extracted from remotely sensed data.
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
  
  #'  Read in covariate data
  #'  Land mgnt, habitat, roads, sum anthro activity, etc. => site-level occupancy covs
  #'  Dist. to focal pt, height, monitoring, weekly anthro activity = > survey/site-level detection covs
  stations <- read.csv("./Outputs/CameraLocation_Covariates18-21_2022-05-02.csv") %>% #2022-04-29
    #'  Get rid of mysterious space after one of the NEs (ugh)
    mutate(
      Study_Area = ifelse(Study_Area == "NE ", "NE", as.character(Study_Area)),
    ) %>%
    dplyr::select(-X) %>%
    #'  Consolidate trail types
    mutate(
      Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Decommissioned road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Game trail", "Trail", as.character(Monitoring)),
    ) %>%
    #'  Consolidate landcover into fewer categories
    mutate(
      landcover_type = ifelse(landcover_type == "Water", "Other", landcover_type),
      landcover_type = ifelse(landcover_type == "Glacier", "Other", landcover_type),
      landcover_type = ifelse(landcover_type == "Barren", "Other", landcover_type),
      landcover_type = ifelse(landcover_type == "Wetland", "Other", landcover_type),
      landcover_type = ifelse(landcover_type == "Wetland", "Other", landcover_type),
      landcover_type = ifelse(landcover_type == "Woody Wetland", "Open Grass", landcover_type),
      landcover_type = ifelse(landcover_type == "Mesic Grass", "Open Grass", landcover_type),
      landcover_type = ifelse(landcover_type == "Xeric Grass", "Open Grass", landcover_type),
      landcover_type = ifelse(landcover_type == "Mesic Shrub", "Shrub Mix", landcover_type),
      landcover_type = ifelse(landcover_type == "Xeric Shrub", "Shrub Mix", landcover_type),
      landcover_type = ifelse(landcover_type == "Forest", "Forest", landcover_type),
      landcover_type = ifelse(landcover_type == "Agriculture", "Other", landcover_type),
      landcover_type = ifelse(landcover_type == "Commercial", "Other", landcover_type),
      landcover_type = ifelse(landcover_type == "Developed", "Other", landcover_type)
      # landcover_type = ifelse(landcover_type == "310", "Other", landcover_type)
    ) %>%
    #'  Consolidate Land_Owner to binary private (0) vs public (1)
    #'  Note this allows private timberlands to be considered "public" for access
    #'  Use Land_Mgnt if want to keep private timberlands private
    mutate(
      Public = ifelse(Land_Owner == "Private", 0, 1)
    ) %>%
    arrange(Year, CameraLocation) #NECESSARY TO MATCH DH's CAMERALOCATION ORDER
  
  
  #'  Remove rows when camera was completely inoperable so missing detection data
  #'  for ALL sampling occasions at that site
  #'  Grazing season missing all data at sites: 10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336
  #'  Hunting season missing all data at sites: 16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283
  stations_graze <- stations[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  stations_hunt <- stations[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  
  
  #'  Split site-level covaraites up by study area
  stations_grazeNE <- filter(stations_graze, Study_Area == "NE")
  stations_grazeOK <- filter(stations_graze, Study_Area == "OK")
  stations_huntNE <- filter(stations_hunt, Study_Area == "NE")
  stations_huntOK <- filter(stations_hunt, Study_Area == "OK")
  
  #'  Scale and format site-level covariates
  format_covs <- function(cams) {
    #'  Rename, format, and scale as needed
    formatted <- transmute(cams,
                           CameraLocation = as.factor(CameraLocation),
                           Year = as.factor(Year),
                           Study_Area = as.factor(Study_Area),
                           Distance = scale(Distance_Focal_Point),
                           Height = scale(Height_frm_grnd),
                           Trail = as.factor(Monitoring),
                           Land_Mgnt = as.factor(Land_Mgnt),
                           Land_Owner = as.factor(Land_Owner),
                           Public = as.factor(Public),
                           GrazingActivity = scale(GrazingActivity), # Summed over 13 week period
                           HuntingActivity = scale(HuntingActivity), # Summed over 8 week period
                           VehicleActivity = scale(VehicleActivity), # Summed over 8 week period
                           Landcover = as.factor(landcover_type),
                           TreeCover = scale(TreeCover),   
                           PercForest = scale(PercForest),    
                           PercXericShrub = scale(PercXericShrub),
                           PercXericGrass = scale(PercXericGrass),
                           Dist2Edge = scale(Dist2Edge),
                           Elev = scale(Elev), 
                           Slope = scale(Slope),
                           TRI = scale(TRI),
                           TPI = scale(TPI),
                           Dist2Water = scale(Dist2Water),
                           NearestRd = scale(NearestRd),
                           RoadDensity = scale(RoadDen), 
                           HumanModified = scale(HumanMod),
                           HumanDensity = scale(HumanDen)
    ) %>%
      arrange(Year, CameraLocation) #NECESSARY TO MATCH DH's CAMERALOCATION ORDER
    
    #'  Adjust reference category for Trail factors
    order_trail <- c("Trail", "Dirt road", "Decommissioned road")
    formatted <- formatted %>%
      mutate(
        Trail = fct_relevel(Trail, order_trail)
      )
    
    #'  Identify which sites have missing data for distance & height covariates
    noDist <- which(is.na(formatted$Distance)); print(noDist)
    noHgt <- which(is.na(formatted$Height)); print(noHgt)
    
    #'  Replace missing values with mean once covariates are z-transformed
    formatted$Distance[is.na(formatted$Distance),] <- 0
    formatted$Height[is.na(formatted$Height),] <- 0
    
    return(formatted)
  }
  stations_graze <- format_covs(stations_graze)
  stations_hunt <- format_covs(stations_hunt)
  stations_grazeNE <- format_covs(stations_grazeNE)
  stations_grazeOK <- format_covs(stations_grazeOK)
  stations_huntNE <- format_covs(stations_huntNE)
  stations_huntOK <- format_covs(stations_huntOK)
  
  #'  Make sure order of rows for DHs match order of rows for covariates
  head(DH_wtd_graze1820)
  head(DH_cattle_graze1820)
  head(stations_graze)
  
  DH_bear_hunt1820[115:120,] # should be end of Year1, start of Year2 cams
  DH_all_hunt1820[115:120,]
  stations_hunt[115:120,]
  
  head(DH_md_graze1820_OK)
  head(DH_cattle_graze1820_OK)
  head(stations_grazeOK)
  
  DH_moose_hunt1820_NE[52:57,] # should be end of Year1, start of Year2 NE cams
  DH_all_hunt1820_NE[52:57,]
  stations_huntNE[52:57,]
  
  
  #'  Check for correlation among covariates
  #'  Watch out for NAs (use = "complete.obs")
  site_covs <- stations_graze[ , c("Elev", "Slope", "TRI", "TPI", "TreeCover", "PercForest", 
                                   "PercXericGrass", "PercXericShrub", "Dist2Water", 
                                   "NearestRd", "RoadDensity", "HumanModified", "HumanDensity", 
                                   "GrazingActivity", "HuntingActivity", "VehicleActivity",
                                   "Distance", "Height")]
  (cov_matrix <- cor(site_covs, use = "complete.obs")) 
  # Elev & HM -0.629; Slope & TRI 0.951; TRI & TPI 0.984; Forest & Tree 0.718
  site_covs_NE <- stations_huntNE[ , c("Elev", "Slope", "TRI", "TPI", "TreeCover", "PercForest",
                                       "PercXericGrass", "PercXericShrub", "Dist2Water",
                                       "NearestRd", "RoadDensity", "HumanModified", "HumanDensity",
                                       "GrazingActivity", "HuntingActivity", "VehicleActivity",
                                       "Distance", "Height")]
  (cov_matrix_NE <- cor(site_covs_NE, use = "complete.obs"))
  site_covs_OK <- stations_huntOK[ , c("Elev", "Slope", "TRI", "TPI", "TreeCover", "PercForest",
                                       "PercXericGrass", "PercXericShrub", "Dist2Water",
                                       "NearestRd", "RoadDensity", "HumanModified", "HumanDensity",
                                       "GrazingActivity", "HuntingActivity", "VehicleActivity",
                                       "Distance", "Height")]
  (cov_matrix_OK <- cor(site_covs_OK, use = "complete.obs"))
  
  #' #'  Survey covariates
  #' #'  Daily max precip and mean temp; extracted from NARR by O.Sanderfoot
  #' #'  Will need to summarize/average based on desired date range & sampling occasions
  #' #'  Note I've added "b" to 2020-2021 cameras to mirror camera station/detection data
  #' narr <- read.csv("./Data/NARR_data_Sanderfoot.csv") 
  #' temp <- narr %>%
  #'   dplyr::select(-Daily_Accumulated_Precipitation.mm.) %>%
  #'   # group_by(Year) %>%
  #'   # pivot_wider(names_from = numbers, values_from = Daily_Mean_Air_Temperature.K.)
  #'   # ungroup() %>%
  #'   arrange(Camera_Location)
  #' temp <- as.data.frame(temp)
  #' 
  #' #'  Break up temp data in summer and winter data frames
  #' temp_smr <- filter(temp, Season == "Summer18" | Season == "Summer19") 
  #' #'  Create rows for cameras that are missing from seasonal temp data
  #' #'  Happens when a camera is deployed after specific data range
  #' missing_smr <- subset(stations, !(CameraLocation %in% temp_smr$CameraLocation)) %>%
  #'   dplyr::select(CameraLocation, Study_Area, Year) %>%
  #'   mutate(
  #'     Season = ifelse(Year == "Year1", "Summer18", "Summer19"),
  #'   )
  #' missing_occ <- matrix(NA, nrow = nrow(missing_smr), ncol = 13)
  #' missing_smr <- cbind(missing_smr, missing_occ)
  #' temp_smr <- rbind(temp_smr, missing_smr) %>%
  #'   arrange(Year, CameraLocation)
  #' #'  Format for UMFs and standardize each column
  #' Temp_smr <- temp_smr %>%
  #'   dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
  #'   scale(.)
  #'   
  #' temp_wtr <- filter(temp, Season == "Winter1819" | Season == "Winter1920") 
  #' #'  Create rows for cameras that are missing from seasonal temp data
  #' #'  Happens when a camera was moved part way through the year
  #' missing_wtr <- subset(stations, !(CameraLocation %in% temp_wtr$CameraLocation)) %>%
  #'   dplyr::select(CameraLocation, Study_Area, Year) %>%
  #'   mutate(
  #'     Season = ifelse(Year == "Year1", "Winter1819", "Winter1920"),
  #'   )
  #' missing_occ <- matrix(NA, nrow = nrow(missing_wtr), ncol = 13)
  #' missing_wtr <- cbind(missing_wtr, missing_occ)
  #' temp_wtr <- rbind(temp_wtr, missing_wtr) %>%
  #'   arrange(Year, CameraLocation)
  #' #'  Format for UMFs and standardize each column
  #' Temp_wtr <- temp_wtr %>%
  #'   dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
  #'   scale(.)
  
  
  #'  Scale survey-level covariates
  scale_srvy_cov <- function(DH) {
    #'  Find mean & standard deviation of covariates across all sites & occasions
    mu <- mean(as.matrix(DH), na.rm = TRUE)
    sd <- sd(as.matrix(DH), na.rm = TRUE)
    
    #'  Z-transform (center observations around mean & scale by 1 SD)
    scaled <- ((DH - mu) / sd)
    
    return(scaled)
  }
  #'  Scale each survey-level covariate
  cattle_week_scaled <- scale_srvy_cov(DH_cattle_graze1820)
  cattle_week_scaled_NE <- scale_srvy_cov(DH_cattle_graze1820_NE)
  cattle_week_scaled_OK <- scale_srvy_cov(DH_cattle_graze1820_OK)
  
  hunter_week_scaled <- scale_srvy_cov(DH_all_hunt1820)
  hunter_week_scaled_NE <- scale_srvy_cov(DH_all_hunt1820_NE)
  hunter_week_scaled_OK <- scale_srvy_cov(DH_all_hunt1820_OK)
  
  vehicle_week_scaled <- scale_srvy_cov(DH_vehicle1820)
  vehicle_week_scaled_NE <- scale_srvy_cov(DH_vehicle1820_NE)
  vehicle_week_scaled_OK <- scale_srvy_cov(DH_vehicle1820_OK)
  
  Effort_graze_scaled <- scale_srvy_cov(Effort_graze1820)
  Effort_graze_scaled_NE <- scale_srvy_cov(Effort_graze1820_NE)
  Effort_graze_scaled_OK <- scale_srvy_cov(Effort_graze1820_OK)
  
  Effort_hunt_scaled <- scale_srvy_cov(Effort_hunt1820)
  Effort_hunt_scaled_NE <- scale_srvy_cov(Effort_hunt1820_NE)
  Effort_hunt_scaled_OK <- scale_srvy_cov(Effort_hunt1820_OK)
  
  #'  Create list of survey-level covariates for unmarked
  srvy_grazing_covs <- list(WeeklyGrazing = cattle_week_scaled,
                            SamplingEffort = Effort_graze_scaled)
  srvy_grazing_covs_NE <- list(WeeklyGrazing = cattle_week_scaled_NE,
                               SamplingEffort = Effort_graze_scaled_NE)
  srvy_grazing_covs_OK <- list(WeeklyGrazing = cattle_week_scaled_OK,
                               SamplingEffort = Effort_graze_scaled_OK)
  
  srvy_hunting_covs <- list(WeeklyHunting = hunter_week_scaled,
                            WeeklyVehicles = vehicle_week_scaled,
                            SamplingEffort = Effort_hunt_scaled)
  srvy_hunting_covs_NE <- list(WeeklyHunting = hunter_week_scaled_NE,
                               WeeklyVehicles = vehicle_week_scaled_NE,
                               SamplingEffort = Effort_hunt_scaled_NE)
  srvy_hunting_covs_OK <- list(WeeklyHunting = hunter_week_scaled_OK,
                               WeeklyVehicles = vehicle_week_scaled_OK,
                               SamplingEffort = Effort_hunt_scaled_OK)
  
  #'  Double check each list looks correct
  head(srvy_grazing_covs[[1]])
  head(srvy_grazing_covs[[2]])
  head(srvy_hunting_covs[[1]])
  head(srvy_hunting_covs[[2]])
  head(srvy_hunting_covs[[3]])
  
  #' #'  Survey-level covariates by study area  ----- NEED TO SCALE BY STUDY AREA IF DOING THIS
  #' Temp_smr_NE <- temp_smr[temp_smr$Study_Area == "NE",] %>%
  #'   dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
  #'   scale(.)
  #' Temp_smr_OK <- temp_smr[temp_smr$Study_Area == "OK",] %>%
  #'   dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
  #'   scale(.)
  #' Temp_wtr_NE <- temp_wtr[temp_wtr$Study_Area == "NE",] %>%
  #'   dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
  #'   scale(.)
  #' Temp_wtr_OK <- temp_wtr[temp_wtr$Study_Area == "OK",] %>%
  #'   dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
  #'   scale(.)
  #' srvy_covs_NE <- list(
  #' )
  #' srvy_covs_OK <- list(
  #' )
  #' 
  #' #'  Double check these
  #' head(srvy_covs_NE[[1]])
  #' tail(srvy_covs_OK[[1]])
  
  
  
  
  
  ####  DON'T FORGET TO DROP GEBBERS CAMERAS FROM GRAZING ANALYSES!!!  #### 
  
  
  
  #'  Create unmarked dataframes for multi-species occupancy models
  #'  FYI: unmarkedFrameOccu is NOT smart enough to match up DH and covs by CameraLocation
  #'  so BE SURE DH AND STATIONS CAMERA SITES ARE ORDERED IN THE SAME WAY!!!!
  
  
  ####  Multi-species unmarkedDF --> unmarkedFrameOccuMulti (pg 151 of unmarked manual)
  ####  Multi-species occupancy model --> occuMulti (pg 83)
  
  
  ####  Setup data for unmarked  ####
  #'  ---------------------------
  #'  Multi-species unmarkedDF --> unmarkedFrameOccuMulti (pg 151 of unmarked manual)
  #'  List relevant detection histories and create unmarked data frame for each 
  #'  multi-species occupancy model. Currently running 2-species occupancy models.
  #'  Define maximum interaction order with maxOrder. Defaults to all possible
  #'  interactions if not defined. 
  
  ####  COUGAR-MULE DEER-CATTLE UMFs  ####
  #'  Grazing Season
  coug_md_cattle_graze_DH <- list(cougar = DH_coug_graze1820, muledeer = DH_md_graze1820, cattle = DH_livestock_graze1820)
  coug_md_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_md_cattle_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_md_cattle_grazing_UMF)
  #'  Review covariates
  summary(coug_md_cattle_grazing_UMF)
  #'  Look at natural parameters f design matrix
  coug_md_cattle_grazing_UMF@fDesign
  
  #'  Hunting Season
  coug_md_hunter_hunt_DH <- list(cougar = DH_coug_hunt1820, muledeer = DH_md_hunt1820, hunter = DH_hunters_hunt1820)
  coug_md_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_md_hunter_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_md_hunter_hunting_UMF)
  #'  Review covariates
  summary(coug_md_hunter_hunting_UMF)
  #'  Look at natural parameters f design matrix
  coug_md_hunter_hunting_UMF@fDesign
  
  ####  COUGAR-ELK UMFs  ####
  #'  Grazing Season
  coug_elk_cattle_graze_DH <- list(cougar = DH_coug_graze1820, elk = DH_elk_graze1820, cattle = DH_livestock_graze1820)
  coug_elk_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_elk_cattle_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_elk_cattle_grazing_UMF)
  
  #'  Hunting Season
  coug_elk_hunter_hunt_DH <- list(cougar = DH_coug_hunt1820, elk = DH_elk_hunt1820, hunter = DH_hunters_hunt1820)
  coug_elk_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_elk_hunter_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_elk_hunter_hunting_UMF)
  
  ####  COUGAR-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  coug_wtd_cattle_graze_DH <- list(cougar = DH_coug_graze1820, wtd = DH_wtd_graze1820, cattle = DH_livestock_graze1820)
  coug_wtd_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_wtd_cattle_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_wtd_cattle_grazing_UMF)
  
  #'  Hunting Season
  coug_wtd_hunter_hunt_DH <- list(cougar = DH_coug_hunt1820, wtd = DH_wtd_hunt1820, hunter = DH_hunters_hunt1820)
  coug_wtd_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_wtd_hunter_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_wtd_hunter_hunting_UMF)
  
  ####  COUGAR-MOOSE UMFs  ####
  #'  Grazing Season
  coug_moose_cattle_graze_DH <- list(cougar = DH_coug_graze1820, moose = DH_moose_graze1820, cattle = DH_livestock_graze1820)
  coug_moose_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_moose_cattle_graze_DH,
                                                   siteCovs = stations_graze,
                                                   obsCovs = srvy_grazing_covs,
                                                   maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_moose_cattle_grazing_UMF)
  
  #'  Hunting Season
  coug_moose_hunter_hunt_DH <- list(cougar = DH_coug_hunt1820, moose = DH_moose_hunt1820, hunter = DH_hunters_hunt1820)
  coug_moose_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_moose_hunter_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coug_moose_hunter_hunting_UMF)
  
  
  ####  WOLF-MULE DEER UMFs  ####
  #'  Grazing Season
  wolf_md_cattle_graze_DH <- list(wolf = DH_wolf_graze1820, muledeer = DH_md_graze1820, cattle = DH_livestock_graze1820)
  wolf_md_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_md_cattle_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_md_cattle_grazing_UMF)
  
  #'  Hunting Season
  wolf_md_hunter_hunt_DH <- list(wolf = DH_wolf_hunt1820, muledeer = DH_md_hunt1820, hunter = DH_hunters_hunt1820)
  wolf_md_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_md_hunter_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_md_hunter_hunting_UMF)
  
  ####  WOLF-ELK UMFs  ####
  #'  Grazing Season
  wolf_elk_cattle_graze_DH <- list(wolf = DH_wolf_graze1820, elk = DH_elk_graze1820, cattle = DH_livestock_graze1820)
  wolf_elk_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_elk_cattle_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_elk_cattle_grazing_UMF)
  
  #'  Hunting Season
  wolf_elk_hunter_hunt_DH <- list(wolf = DH_wolf_hunt1820, elk = DH_elk_hunt1820, hunter = DH_hunters_hunt1820)
  wolf_elk_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_elk_hunter_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_elk_hunter_hunting_UMF)
  
  ####  WOLF-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  wolf_wtd_cattle_graze_DH <- list(wolf = DH_wolf_graze1820, wtd = DH_wtd_graze1820, cattle = DH_livestock_graze1820)
  wolf_wtd_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_wtd_cattle_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_wtd_cattle_grazing_UMF)
  
  #'  Hunting Season
  wolf_wtd_hunter_hunt_DH <- list(wolf = DH_wolf_hunt1820, wtd = DH_wtd_hunt1820, hunter = DH_hunters_hunt1820)
  wolf_wtd_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_wtd_hunter_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_wtd_hunter_hunting_UMF)
  
  ####  WOLF-MOOSE UMFs  ####
  #'  Grazing Season
  wolf_moose_cattle_graze_DH <- list(wolf = DH_wolf_graze1820, moose = DH_moose_graze1820, cattle = DH_livestock_graze1820)
  wolf_moose_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_moose_cattle_graze_DH,
                                                   siteCovs = stations_graze,
                                                   obsCovs = srvy_grazing_covs,
                                                   maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_moose_cattle_grazing_UMF)
  
  #'  Hunting Season
  wolf_moose_hunter_hunt_DH <- list(wolf = DH_wolf_hunt1820, moose = DH_moose_hunt1820, hunter = DH_hunters_hunt1820)
  wolf_moose_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_moose_hunter_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(wolf_moose_hunter_hunting_UMF)
  
  
  
  ####  BLACK Bear-MULE DEER UMFs  ####
  #'  Grazing Season
  bear_md_cattle_graze_DH <- list(blackbear = DH_bear_graze1820, muledeer = DH_md_graze1820, cattle = DH_livestock_graze1820)
  bear_md_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_md_cattle_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_md_cattle_grazing_UMF)
  
  #'  Hunting Season
  bear_md_hunter_hunt_DH <- list(blackbear = DH_bear_hunt1820, muledeer = DH_md_hunt1820, hunter = DH_hunters_hunt1820)
  bear_md_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_md_hunter_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_md_hunter_hunting_UMF)
  
  ####  BLACK BEAR-ELK UMFs  ####
  #'  Grazing Season
  bear_elk_cattle_graze_DH <- list(blackbear = DH_bear_graze1820, elk = DH_elk_graze1820, cattle = DH_livestock_graze1820)
  bear_elk_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_elk_cattle_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_elk_cattle_grazing_UMF)
  
  #'  Hunting Season
  bear_elk_hunter_hunt_DH <- list(blackbear = DH_bear_hunt1820, elk = DH_elk_hunt1820, hunter = DH_hunters_hunt1820)
  bear_elk_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_elk_hunter_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_elk_hunter_hunting_UMF)
  
  ####  BLACK BEAR-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  bear_wtd_cattle_graze_DH <- list(blackbear = DH_bear_graze1820, wtd = DH_wtd_graze1820, cattle = DH_livestock_graze1820)
  bear_wtd_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_wtd_cattle_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_wtd_cattle_grazing_UMF)
  
  #'  Hunting Season
  bear_wtd_hunter_hunt_DH <- list(blackbear = DH_bear_hunt1820, wtd = DH_wtd_hunt1820, hunter = DH_hunters_hunt1820)
  bear_wtd_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_wtd_hunter_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_wtd_hunter_hunting_UMF)
  
  ####  BLACK BEAR-MOOSE UMFs  ####
  #'  Grazing Season
  bear_moose_cattle_graze_DH <- list(blackbear = DH_bear_graze1820, moose = DH_moose_graze1820, cattle = DH_livestock_graze1820)
  bear_moose_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_moose_cattle_graze_DH,
                                                   siteCovs = stations_graze,
                                                   obsCovs = srvy_grazing_covs,
                                                   maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_moose_cattle_grazing_UMF)
  
  #'  Hunting Season
  bear_moose_hunter_hunt_DH <- list(blackbear = DH_bear_hunt1820, moose = DH_moose_hunt1820, hunter = DH_hunters_hunt1820)
  bear_moose_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_moose_hunter_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bear_moose_hunter_hunting_UMF)
  
  
  
  ####  COYOTE-MULE DEER UMFs  ####
  #'  Grazing Season
  coy_md_cattle_graze_DH <- list(coyote = DH_coy_graze1820, mule_deer = DH_md_graze1820, cattle = DH_livestock_graze1820)
  coy_md_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = coy_md_cattle_graze_DH,
                                               siteCovs = stations_graze,
                                               obsCovs = srvy_grazing_covs,
                                               maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coy_md_cattle_grazing_UMF)
  
  #'  Hunting Season
  coy_md_hunter_hunt_DH <- list(coyote = DH_coy_hunt1820, mule_deer = DH_md_hunt1820, hunter = DH_hunters_hunt1820)
  coy_md_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = coy_md_hunter_hunt_DH,
                                               siteCovs = stations_hunt,
                                               obsCovs = srvy_hunting_covs,
                                               maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coy_md_hunter_hunting_UMF)
  
  ####  COYOTE-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  coy_wtd_cattle_graze_DH <- list(coyote = DH_coy_graze1820, wtd = DH_wtd_graze1820, cattle = DH_livestock_graze1820)
  coy_wtd_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = coy_wtd_cattle_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coy_wtd_cattle_grazing_UMF)
  
  #'  Hunting Season
  coy_wtd_hunter_hunt_DH <- list(coyote = DH_coy_hunt1820, wtd = DH_wtd_hunt1820, hunter = DH_hunters_hunt1820)
  coy_wtd_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = coy_wtd_hunter_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(coy_wtd_hunter_hunting_UMF)
  
  
  
  ####  BOBCAT-MULE DEER UMFs  ####
  #'  Grazing Season
  bob_md_cattle_graze_DH <- list(bobcat = DH_bob_graze1820, mule_deer = DH_md_graze1820, cattle = DH_livestock_graze1820)
  bob_md_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = bob_md_cattle_graze_DH,
                                               siteCovs = stations_graze,
                                               obsCovs = srvy_grazing_covs,
                                               maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bob_md_cattle_grazing_UMF)
  
  #'  Hunting Season
  bob_md_hunter_hunt_DH <- list(bobcat = DH_bob_hunt1820, mule_deer = DH_md_hunt1820, hunter = DH_hunters_hunt1820)
  bob_md_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = bob_md_hunter_hunt_DH,
                                               siteCovs = stations_hunt,
                                               obsCovs = srvy_hunting_covs,
                                               maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bob_md_hunter_hunting_UMF)
  
  ####  BOBCAT-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  bob_wtd_cattle_graze_DH <- list(bobcat = DH_bob_graze1820, wtd = DH_wtd_graze1820, cattle = DH_livestock_graze1820)
  bob_wtd_cattle_grazing_UMF <- unmarkedFrameOccuMulti(y = bob_wtd_cattle_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bob_wtd_cattle_grazing_UMF)
  
  #'  Hunting Season
  bob_wtd_hunter_hunt_DH <- list(bobcat = DH_bob_hunt1820, wtd = DH_wtd_hunt1820, hunter = DH_hunters_hunt1820)
  bob_wtd_hunter_hunting_UMF <- unmarkedFrameOccuMulti(y = bob_wtd_hunter_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 3)
  #'  Visualize detection/non-detection data
  plot(bob_wtd_hunter_hunting_UMF)
  
  
  
  ####  Occupancy models  ####
  #'  ====================
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
  #'  Hypothesis testing with 4 models: 
  #'    1) anthro activity on detection, assume independent co-occurrence/detection
  #'    2) anthro activity on detection, allow non-independent co-occurrence/detection
  #'    3) anthro activity on occupancy, assume independent co-occurrence/detection
  #'    4) anthro acticity on occupancy, allow non-independent co-occurrence/detection
  #'  Include a consistent set of additional covariates to account for habitat
  #'  variation and other factors we know influence occurrence and deteciton.
  #'  Use AIC for model selection
  #'  =============================
  
  ####  GRAZING SEASON MODELS  ####
  #'  Detection and occupancy formulas
  #'  Question 1: Does grazing activity affect detection?
  detFormulas_trail <- c("~Trail", "~Trail", "~Trail")
  detFormulas_graze <- c("~Trail + WeeklyGrazing", "~Trail + WeeklyGrazing", "~Trail") # no grazing on cattle detection
  
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null <- c("~1", "~1", "~1", "~0", "~0", "~0", "~0") # no interactions
  occFormulas_null1 <- c("~1", "~1", "~1", "~1", "~1", "~1", "~0") # 2-spp interactions
  occFormulas_null2 <- c("~1", "~1", "~1", "~1", "~1", "~1", "~1") # 3-spp interaction
  
  #'  Question 2: Does grazing activity and/or presence of cattle influence 
  #'  co-occurrence or conditional occupancy of predators and prey?
  occFormulas_hab0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~0", "~0", "~0", "~0")
  occFormulas_hab1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~1", "~1", "~1",
                        "~0")
  occFormulas_grz <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~GrazingActivity", "~1", "~1",
                        "~0")
  occFormulas_pub <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~Public", "~1", "~1",
                        "~0")
  occFormulas_grzpub <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~GrazingActivity + Public", "~1", "~1",
                        "~0")

  
  ####  Cougar-Mule Deer-Cattle Grazing Season  ####
  (coug.md.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, coug_md_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coug.md.cow_fld <- fitList(coug.md.cow_trail, coug.md.cow_graze)
  #' Model selection
  modSel(coug.md.cow_fld)
  
  (coug.md.cow_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_grz <- occuMulti(detFormulas_trail, occFormulas_grz, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_pub <- occuMulti(detFormulas_trail, occFormulas_pub, coug_md_cattle_grazing_UMF, silent = TRUE))
  (coug.md.cow_grzpub <- occuMulti(detFormulas_trail, occFormulas_grzpub, coug_md_cattle_grazing_UMF, silent = TRUE))
  coug.md.cow_fl <- fitList(coug.md.cow_hab0, coug.md.cow_hab1, coug.md.cow_grz, coug.md.cow_pub, coug.md.cow_grzpub)
  #' Model selection
  modSel(coug.md.cow_fl)
  summary(coug.md.cow_hab0)
  
  ####  Cougar-Elk-Cattle Grazing Season  ####
  (coug.elk.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, coug_elk_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coug.elk.cow_fld <- fitList(coug.elk.cow_trail, coug.elk.cow_graze)
  #' Model selection
  modSel(coug.elk.cow_fld)
  
  (coug.elk.cow_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_grz <- occuMulti(detFormulas_trail, occFormulas_grz, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_pub <- occuMulti(detFormulas_trail, occFormulas_pub, coug_elk_cattle_grazing_UMF, silent = TRUE))
  (coug.elk.cow_grzpub <- occuMulti(detFormulas_trail, occFormulas_grzpub, coug_elk_cattle_grazing_UMF, silent = TRUE))
  coug.elk.cow_fl <- fitList(coug.elk.cow_hab0, coug.elk.cow_hab1, coug.elk.cow_grz, coug.elk.cow_pub, coug.elk.cow_grzpub)
  #' Model selection
  modSel(coug.elk.cow_fl)
  summary(coug.elk.cow_hab0)
  
  ####  Cougar-White-tailed Deer-Cattle Grazing Season  ####
  (coug.wtd.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_cattle_grazing_UMF, silent = TRUE))
  (coug.wtd.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coug.wtd.cow_fld <- fitList(coug.wtd.cow_trail, coug.wtd.cow_graze)
  #' Model selection
  modSel(coug.wtd.cow_fld)
  
  (coug.wtd.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_cattle_grazing_UMF, silent = TRUE))
  (coug.wtd.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_wtd_cattle_grazing_UMF, silent = TRUE))
  (coug.wtd.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, coug_wtd_cattle_grazing_UMF, silent = TRUE))
  (coug.wtd.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_wtd_cattle_grazing_UMF, silent = TRUE))
  (coug.wtd.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_wtd_cattle_grazing_UMF, silent = TRUE))
  (coug.wtd.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, coug_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (coug.wtd.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, coug_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (coug.wtd.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, coug_wtd_cattle_grazing_UMF, silent = TRUE)) 
  coug.wtd.cow_fl <- fitList(coug.wtd.cow_hab0, coug.wtd.cow_hab1, coug.wtd.cow_grz, coug.wtd.cow_pub, coug.wtd.cow_grzpub)
  #' Model selection
  modSel(coug.wtd.cow_fl)
  summary(coug.wtd.cow_grz)
  summary(coug.wtd.cow_grzpub)
  
  ####  Cougar-Moose-Cattle Grazing Season  ####
  (coug.moose.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_cattle_grazing_UMF, silent = TRUE))
  (coug.moose.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coug.moose.cow_fld <- fitList(coug.moose.cow_trail, coug.moose.cow_graze)
  #' Model selection
  modSel(coug.moose.cow_fld)
  
  (coug.moose.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_cattle_grazing_UMF, silent = TRUE))
  (coug.moose.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_moose_cattle_grazing_UMF, silent = TRUE))
  (coug.moose.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, coug_moose_cattle_grazing_UMF, silent = TRUE))
  (coug.moose.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_moose_cattle_grazing_UMF, silent = TRUE))
  (coug.moose.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_moose_cattle_grazing_UMF, silent = TRUE))
  (coug.moose.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, coug_moose_cattle_grazing_UMF, silent = TRUE)) 
  (coug.moose.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, coug_moose_cattle_grazing_UMF, silent = TRUE)) 
  (coug.moose.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, coug_moose_cattle_grazing_UMF, silent = TRUE)) 
  coug.moose.cow_fl <- fitList(coug.moose.cow_hab0, coug.moose.cow_hab1, coug.moose.cow_grz, coug.moose.cow_pub, coug.moose.cow_grzpub)
  #' Model selection
  modSel(coug.moose.cow_fl)
  summary(coug.moose.cow_grz)
  summary(coug.moose.cow_hab0)
  summary(coug.moose.cow_grzpub)
  
  
  ####  Wolf-Mule Deer-Cattle Grazing Season  ####
  (wolf.md.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_cattle_grazing_UMF, silent = TRUE))
  (wolf.md.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_md_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolf.md.cow_fld <- fitList(wolf.md.cow_trail, wolf.md.cow_graze)
  #' Model selection
  modSel(wolf.md.cow_fld)
  
  (wolf.md.cow_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_cattle_grazing_UMF, silent = TRUE))
  (wolf.md.cow_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_md_cattle_grazing_UMF, silent = TRUE))
  (wolf.md.cow_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_md_cattle_grazing_UMF, silent = TRUE))
  (wolf.md.cow_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_md_cattle_grazing_UMF, silent = TRUE))
  (wolf.md.cow_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_md_cattle_grazing_UMF, silent = TRUE))
  (wolf.md.cow_grz <- occuMulti(detFormulas_trail, occFormulas_grz, wolf_md_cattle_grazing_UMF, silent = TRUE))
  (wolf.md.cow_pub <- occuMulti(detFormulas_trail, occFormulas_pub, wolf_md_cattle_grazing_UMF, silent = TRUE)) # FAIL
  (wolf.md.cow_grzpub <- occuMulti(detFormulas_trail, occFormulas_grzpub, wolf_md_cattle_grazing_UMF, silent = TRUE)) 
  wolf.md.cow_fl <- fitList(wolf.md.cow_hab0, wolf.md.cow_hab1, wolf.md.cow_grz, wolf.md.cow_grzpub)
  #' Model selection
  modSel(wolf.md.cow_fl)
  summary(wolf.md.cow_grz)
  summary(wolf.md.cow_grzpub)
  
  
  ####  Wolf-Elk-Cattle Grazing Season  ####
  (wolf.elk.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_cattle_grazing_UMF, silent = TRUE))
  (wolf.elk.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_elk_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolf.elk.cow_fld <- fitList(wolf.elk.cow_trail, wolf.elk.cow_graze)
  #' Model selection
  modSel(wolf.elk.cow_fld)
  
  (wolf.elk.cow_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_cattle_grazing_UMF, silent = TRUE))
  (wolf.elk.cow_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_elk_cattle_grazing_UMF, silent = TRUE))
  (wolf.elk.cow_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, wolf_elk_cattle_grazing_UMF, silent = TRUE)) # FAIL
  (wolf.elk.cow_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_cattle_grazing_UMF, silent = TRUE))
  (wolf.elk.cow_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_cattle_grazing_UMF, silent = TRUE))
  (wolf.elk.cow_grz <- occuMulti(detFormulas_trail, occFormulas_grz, wolf_elk_cattle_grazing_UMF, silent = TRUE)) # FAIL
  (wolf.elk.cow_pub <- occuMulti(detFormulas_trail, occFormulas_pub, wolf_elk_cattle_grazing_UMF, silent = TRUE)) # FAIL
  (wolf.elk.cow_grzpub <- occuMulti(detFormulas_trail, occFormulas_grzpub, wolf_elk_cattle_grazing_UMF, silent = TRUE)) # FAIL 
  wolf.elk.cow_fl <- fitList(wolf.elk.cow_hab0, wolf.elk.cow_hab1)
  #' Model selection
  modSel(wolf.elk.cow_fl)
  summary(wolf.elk.cow_hab0)
  
  ####  Wolf-White-tailed Deer-Cattle Grazing Season  ####
  (wolf.wtd.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_cattle_grazing_UMF, silent = TRUE))
  (wolf.wtd.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolf.wtd.cow_fld <- fitList(wolf.wtd.cow_trail, wolf.wtd.cow_graze)
  #' Model selection
  modSel(wolf.wtd.cow_fld)
  
  (wolf.wtd.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_cattle_grazing_UMF, silent = TRUE))
  (wolf.wtd.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_wtd_cattle_grazing_UMF, silent = TRUE))
  (wolf.wtd.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, wolf_wtd_cattle_grazing_UMF, silent = TRUE))
  (wolf.wtd.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_wtd_cattle_grazing_UMF, silent = TRUE))
  (wolf.wtd.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_wtd_cattle_grazing_UMF, silent = TRUE))
  (wolf.wtd.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, wolf_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (wolf.wtd.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, wolf_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (wolf.wtd.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, wolf_wtd_cattle_grazing_UMF, silent = TRUE)) 
  wolf.wtd.cow_fl <- fitList(wolf.wtd.cow_hab0, wolf.wtd.cow_hab1, wolf.wtd.cow_grz, wolf.wtd.cow_pub, wolf.wtd.cow_grzpub)
  #' Model selection
  modSel(wolf.wtd.cow_fl)
  summary(wolf.wtd.cow_grz)
  summary(wolf.wtd.cow_grzpub)
  summary(wolf.wtd.cow_hab1)
  
  ####  Wolf-Moose-Cattle Grazing Season  ####
  (wolf.moose.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_cattle_grazing_UMF, silent = TRUE))
  (wolf.moose.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolf.moose.cow_fld <- fitList(wolf.moose.cow_trail, wolf.moose.cow_graze)
  #' Model selection
  modSel(wolf.moose.cow_fld)
  
  (wolf.moose.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_cattle_grazing_UMF, silent = TRUE))
  (wolf.moose.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_moose_cattle_grazing_UMF, silent = TRUE))
  (wolf.moose.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, wolf_moose_cattle_grazing_UMF, silent = TRUE))
  (wolf.moose.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_moose_cattle_grazing_UMF, silent = TRUE))
  (wolf.moose.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_moose_cattle_grazing_UMF, silent = TRUE))
  (wolf.moose.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, wolf_moose_cattle_grazing_UMF, silent = TRUE)) 
  (wolf.moose.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, wolf_moose_cattle_grazing_UMF, silent = TRUE)) # FAIL
  (wolf.moose.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, wolf_moose_cattle_grazing_UMF, silent = TRUE)) # FAIL
  wolf.moose.cow_fl <- fitList(wolf.moose.cow_hab0, wolf.moose.cow_hab1, wolf.moose.cow_grz)
  #' Model selection
  modSel(wolf.moose.cow_fl)
  summary(wolf.moose.cow_hab0)
  summary(wolf.moose.cow_grz)
  
  
  ####  Bear-Mule Deer-Cattle Grazing Season  ####
  (bear.md.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_cattle_grazing_UMF, silent = TRUE))
  (bear.md.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, bear_md_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bear.md.cow_fld <- fitList(bear.md.cow_trail, bear.md.cow_graze)
  #' Model selection
  modSel(bear.md.cow_fld)
  
  (bear.md.cow_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_cattle_grazing_UMF, silent = TRUE))
  (bear.md.cow_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_md_cattle_grazing_UMF, silent = TRUE))
  (bear.md.cow_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bear_md_cattle_grazing_UMF, silent = TRUE))
  (bear.md.cow_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_md_cattle_grazing_UMF, silent = TRUE))
  (bear.md.cow_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_md_cattle_grazing_UMF, silent = TRUE))
  (bear.md.cow_grz <- occuMulti(detFormulas_trail, occFormulas_grz, bear_md_cattle_grazing_UMF, silent = TRUE))
  (bear.md.cow_pub <- occuMulti(detFormulas_trail, occFormulas_pub, bear_md_cattle_grazing_UMF, silent = TRUE)) 
  (bear.md.cow_grzpub <- occuMulti(detFormulas_trail, occFormulas_grzpub, bear_md_cattle_grazing_UMF, silent = TRUE)) 
  bear.md.cow_fl <- fitList(bear.md.cow_hab0, bear.md.cow_hab1, bear.md.cow_grz, bear.md.cow_pub, bear.md.cow_grzpub)
  #' Model selection
  modSel(bear.md.cow_fl)
  # summary(bear.md.cow_grzpub) # but basically none of the interactions are significant
  # summary(bear.md.cow_grz) # but basically none of the interactions are significant
  # summary(bear.md.cow_pub) # but basically none of the interactions are significant
  # summary(bear.md.cow_hab1) # but basically none of the interactions are significant
  summary(bear.md.cow_hab0)
  
  ####  Bear-Elk-Cattle Grazing Season  ####
  (bear.elk.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_cattle_grazing_UMF, silent = TRUE))
  (bear.elk.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, bear_elk_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bear.elk.cow_fld <- fitList(bear.elk.cow_trail, bear.elk.cow_graze)
  #' Model selection
  modSel(bear.elk.cow_fld)
  
  (bear.elk.cow_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_cattle_grazing_UMF, silent = TRUE))
  (bear.elk.cow_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_elk_cattle_grazing_UMF, silent = TRUE))
  (bear.elk.cow_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bear_elk_cattle_grazing_UMF, silent = TRUE)) 
  (bear.elk.cow_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_cattle_grazing_UMF, silent = TRUE))
  (bear.elk.cow_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_cattle_grazing_UMF, silent = TRUE))
  (bear.elk.cow_grz <- occuMulti(detFormulas_trail, occFormulas_grz, bear_elk_cattle_grazing_UMF, silent = TRUE))
  (bear.elk.cow_pub <- occuMulti(detFormulas_trail, occFormulas_pub, bear_elk_cattle_grazing_UMF, silent = TRUE)) 
  (bear.elk.cow_grzpub <- occuMulti(detFormulas_trail, occFormulas_grzpub, bear_elk_cattle_grazing_UMF, silent = TRUE)) 
  bear.elk.cow_fl <- fitList(bear.elk.cow_hab0, bear.elk.cow_hab1, bear.elk.cow_grz, bear.elk.cow_pub, bear.elk.cow_grzpub)
  #' Model selection
  modSel(bear.elk.cow_fl)
  summary(bear.elk.cow_pub) # but none of the interactions are significant
  summary(bear.elk.cow_hab1) # but none of the interactions are significant
  summary(bear.elk.cow_hab0)
  summary(bear.elk.cow_grzpub) # but none of the interactions are significant

  
  ####  Bear-White-tailed Deer-Cattle Grazing Season  ####
  (bear.wtd.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_cattle_grazing_UMF, silent = TRUE))
  (bear.wtd.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bear.wtd.cow_fld <- fitList(bear.wtd.cow_trail, bear.wtd.cow_graze)
  #' Model selection
  modSel(bear.wtd.cow_fld)
  
  (bear.wtd.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_cattle_grazing_UMF, silent = TRUE))
  (bear.wtd.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_wtd_cattle_grazing_UMF, silent = TRUE))
  (bear.wtd.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, bear_wtd_cattle_grazing_UMF, silent = TRUE))
  (bear.wtd.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_wtd_cattle_grazing_UMF, silent = TRUE))
  (bear.wtd.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_wtd_cattle_grazing_UMF, silent = TRUE))
  (bear.wtd.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, bear_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (bear.wtd.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, bear_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (bear.wtd.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, bear_wtd_cattle_grazing_UMF, silent = TRUE)) 
  bear.wtd.cow_fl <- fitList(bear.wtd.cow_hab0, bear.wtd.cow_hab1, bear.wtd.cow_grz, bear.wtd.cow_pub, bear.wtd.cow_grzpub)
  #' Model selection
  modSel(bear.wtd.cow_fl)
  summary(bear.wtd.cow_grz)
  summary(bear.wtd.cow_grzpub)
  
  
  ####  Bear-Moose-Cattle Grazing Season  ####
  (bear.moose.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_cattle_grazing_UMF, silent = TRUE))
  (bear.moose.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bear.moose.cow_fld <- fitList(bear.moose.cow_trail, bear.moose.cow_graze)
  #' Model selection
  modSel(bear.moose.cow_fld)
  
  (bear.moose.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_cattle_grazing_UMF, silent = TRUE))
  (bear.moose.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_moose_cattle_grazing_UMF, silent = TRUE))
  (bear.moose.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, bear_moose_cattle_grazing_UMF, silent = TRUE))
  (bear.moose.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_moose_cattle_grazing_UMF, silent = TRUE))
  (bear.moose.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_moose_cattle_grazing_UMF, silent = TRUE))
  (bear.moose.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, bear_moose_cattle_grazing_UMF, silent = TRUE)) 
  (bear.moose.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, bear_moose_cattle_grazing_UMF, silent = TRUE))
  (bear.moose.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, bear_moose_cattle_grazing_UMF, silent = TRUE)) 
  bear.moose.cow_fl <- fitList(bear.moose.cow_hab0, bear.moose.cow_hab1, bear.moose.cow_grz, bear.moose.cow_pub, bear.moose.cow_grzpub)
  #' Model selection
  modSel(bear.moose.cow_fl)
  summary(bear.moose.cow_grzpub)
  summary(bear.moose.cow_grz)
  
  
  
  ####  Coyote-Mule Deer-Cattle Grazing Season  ####
  (coy.md.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_cattle_grazing_UMF, silent = TRUE))
  (coy.md.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coy.md.cow_fld <- fitList(coy.md.cow_trail, coy.md.cow_graze)
  #' Model selection
  modSel(coy.md.cow_fld)
  
  (coy.md.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_cattle_grazing_UMF, silent = TRUE))
  (coy.md.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_md_cattle_grazing_UMF, silent = TRUE))
  (coy.md.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, coy_md_cattle_grazing_UMF, silent = TRUE))
  (coy.md.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_md_cattle_grazing_UMF, silent = TRUE))
  (coy.md.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_md_cattle_grazing_UMF, silent = TRUE))
  (coy.md.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, coy_md_cattle_grazing_UMF, silent = TRUE))
  (coy.md.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, coy_md_cattle_grazing_UMF, silent = TRUE)) 
  (coy.md.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, coy_md_cattle_grazing_UMF, silent = TRUE)) 
  coy.md.cow_fl <- fitList(coy.md.cow_hab0, coy.md.cow_hab1, coy.md.cow_grz, coy.md.cow_pub, coy.md.cow_grzpub)
  #' Model selection
  modSel(coy.md.cow_fl)
  summary(coy.md.cow_grz)
  summary(coy.md.cow_grzpub)
  
  
  ####  Coyote-White-tailed Deer-Cattle Grazing Season  ####
  (coy.wtd.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_cattle_grazing_UMF, silent = TRUE))
  (coy.wtd.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coy.wtd.cow_fld <- fitList(coy.wtd.cow_trail, coy.wtd.cow_graze)
  #' Model selection
  modSel(coy.wtd.cow_fld)
  
  (coy.wtd.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_cattle_grazing_UMF, silent = TRUE))
  (coy.wtd.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_wtd_cattle_grazing_UMF, silent = TRUE))
  (coy.wtd.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, coy_wtd_cattle_grazing_UMF, silent = TRUE))
  (coy.wtd.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_wtd_cattle_grazing_UMF, silent = TRUE))
  (coy.wtd.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_wtd_cattle_grazing_UMF, silent = TRUE))
  (coy.wtd.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, coy_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (coy.wtd.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, coy_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (coy.wtd.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, coy_wtd_cattle_grazing_UMF, silent = TRUE)) 
  coy.wtd.cow_fl <- fitList(coy.wtd.cow_hab0, coy.wtd.cow_hab1, coy.wtd.cow_grz, coy.wtd.cow_pub, coy.wtd.cow_grzpub)
  #' Model selection
  modSel(coy.wtd.cow_fl)
  summary(coy.wtd.cow_grzpub)
  
  
  
  ####  Bobcat-Mule Deer-Cattle Grazing Season  ####
  (bob.md.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_cattle_grazing_UMF, silent = TRUE))
  (bob.md.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, bob_md_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bob.md.cow_fld <- fitList(bob.md.cow_trail, bob.md.cow_graze)
  #' Model selection
  modSel(bob.md.cow_fld)
  
  (bob.md.cow_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_cattle_grazing_UMF, silent = TRUE))
  (bob.md.cow_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_md_cattle_grazing_UMF, silent = TRUE))
  (bob.md.cow_null2 <- occuMulti(detFormulas_trail, occFormulas_null2, bob_md_cattle_grazing_UMF, silent = TRUE))
  (bob.md.cow_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_md_cattle_grazing_UMF, silent = TRUE))
  (bob.md.cow_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_md_cattle_grazing_UMF, silent = TRUE))
  (bob.md.cow_grz <- occuMulti(detFormulas_trail, occFormulas_grz, bob_md_cattle_grazing_UMF, silent = TRUE))
  (bob.md.cow_pub <- occuMulti(detFormulas_trail, occFormulas_pub, bob_md_cattle_grazing_UMF, silent = TRUE)) 
  (bob.md.cow_grzpub <- occuMulti(detFormulas_trail, occFormulas_grzpub, bob_md_cattle_grazing_UMF, silent = TRUE)) 
  bob.md.cow_fl <- fitList(bob.md.cow_hab0, bob.md.cow_hab1, bob.md.cow_grz, bob.md.cow_pub, bob.md.cow_grzpub)
  #' Model selection
  modSel(bob.md.cow_fl)
  summary(bob.md.cow_grzpub)
  summary(bob.md.cow_grz)
  
  
  ####  Bobcat-White-tailed Deer-Cattle Grazing Season  ####
  (bob.wtd.cow_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_cattle_grazing_UMF, silent = TRUE))
  (bob.wtd.cow_graze <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_cattle_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bob.wtd.cow_fld <- fitList(bob.wtd.cow_trail, bob.wtd.cow_graze)
  #' Model selection
  modSel(bob.wtd.cow_fld)
  
  (bob.wtd.cow_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_cattle_grazing_UMF, silent = TRUE))
  (bob.wtd.cow_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bob_wtd_cattle_grazing_UMF, silent = TRUE))
  (bob.wtd.cow_null2 <- occuMulti(detFormulas_graze, occFormulas_null2, bob_wtd_cattle_grazing_UMF, silent = TRUE))
  (bob.wtd.cow_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bob_wtd_cattle_grazing_UMF, silent = TRUE))
  (bob.wtd.cow_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bob_wtd_cattle_grazing_UMF, silent = TRUE))
  (bob.wtd.cow_grz <- occuMulti(detFormulas_graze, occFormulas_grz, bob_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (bob.wtd.cow_pub <- occuMulti(detFormulas_graze, occFormulas_pub, bob_wtd_cattle_grazing_UMF, silent = TRUE)) 
  (bob.wtd.cow_grzpub <- occuMulti(detFormulas_graze, occFormulas_grzpub, bob_wtd_cattle_grazing_UMF, silent = TRUE)) 
  bob.wtd.cow_fl <- fitList(bob.wtd.cow_hab0, bob.wtd.cow_hab1, bob.wtd.cow_grz, bob.wtd.cow_pub, bob.wtd.cow_grzpub)
  #' Model selection
  modSel(bob.wtd.cow_fl)
  summary(bob.wtd.cow_grz)
  summary(bob.wtd.cow_grzpub)
  
  
  
  ####  HUNTING SEASON MODELS  ####
  #'  Detection and occupancy formulas
  #'  Question 1: Does hunting activity affect detection?
  detFormulas_trail <- c("~Trail", "~Trail", "~Trail")
  detFormulas_hunt <- c("~Trail + WeeklyHunting", "~Trail + WeeklyHunting", "~Trail")
  detFormulas_pub <- c("~Trail + Public", "~Trail + Public", "~Trail + Public")
  #'  Remove Public vs Private covariate from predator sub-model
  detFormulas_pubish <- c("~Trail", "~Trail + Public", "~Trail + Public")
  
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null <- c("~1", "~1", "~1", "~0", "~0", "~0", "~0") # no interactions
  occFormulas_null1 <- c("~1", "~1", "~1", "~1", "~1", "~1", "~0") # 2-spp interactions
  occFormulas_null2 <- c("~1", "~1", "~1", "~1", "~1", "~1", "~1") # 3-spp interaction
  
  occFormulas_habpub0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                           "~0", "~0", "~0", "~0")
  occFormulas_habpub1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                           "~1", "~1", "~1",
                           "~0")
  occFormulas_habpub2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                           "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                           "~Public", "~1", "~1",
                           "~0")
  occFormulas_huntpub1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                            "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                            "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                            "~HuntingActivity",  "~1", "~1",
                            "~0")
  occFormulas_huntpub2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                            "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                            "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                            "~HuntingActivity + Public",  "~1", "~1",
                            "~0")

  #'  Remove Public vs Private covariate from predator sub-model and interactions
  #'  Necessary for cougar & wolf models- this binary covariate blows things up
  occFormulas_habpubish0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                              "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                              "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                              "~0", "~0", "~0", "~0")
  occFormulas_habpubish1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                              "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                              "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                              "~1", "~1", "~1",
                              "~0")
  occFormulas_huntpubish1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest",
                               "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                               "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                               "~HuntingActivity",  "~1", "~1",
                               "~0")
  occFormulas_huntpubish2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest",
                               "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                               "~Study_Area + Elev + I(Elev^2) + PercForest + Public",
                               "~HuntingActivity + Public",  "~1", "~1",
                               "~0")
  
  ####  Cougar-Mule Deer-Hunter Hunting Season  ####
  (coug.md.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_hunter_hunting_UMF, silent = TRUE))
  (coug.md.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_md_hunter_hunting_UMF, silent = TRUE))
  (coug.md.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coug.md.hunt_fld <- fitList(coug.md.hunt_trail, coug.md.hunt_hunt, coug.md.hunt_pub)
  #' Model selection
  modSel(coug.md.hunt_fld)
  
  (coug.md.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunter_hunting_UMF, silent = TRUE))
  (coug.md.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_md_hunter_hunting_UMF, silent = TRUE))
  (coug.md.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, coug_md_hunter_hunting_UMF, silent = TRUE)) # FAIL 
  (coug.md.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpubish0, coug_md_hunter_hunting_UMF, silent = TRUE)) 
  (coug.md.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpubish1, coug_md_hunter_hunting_UMF, silent = TRUE)) 
  (coug.md.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coug_md_hunter_hunting_UMF, silent = TRUE)) 
  (coug.md.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coug_md_hunter_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  coug.md.hunt_fl <- fitList(coug.md.hunt_habpub0, coug.md.hunt_habpub1, coug.md.hunt_huntpub1, coug.md.hunt_huntpub2) 
  #' Model selection
  modSel(coug.md.hunt_fl)
  summary(coug.md.hunt_huntpub1)
  summary(coug.md.hunt_huntpub2)
  
  ####  Cougar-Elk-Hunter Hunting Season  ####
  (coug.elk.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_hunter_hunting_UMF, silent = TRUE))
  (coug.elk.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_elk_hunter_hunting_UMF, silent = TRUE))
  (coug.elk.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_elk_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coug.elk.hunt_fld <- fitList(coug.elk.hunt_trail, coug.elk.hunt_hunt, coug.elk.hunt_pub)
  #' Model selection
  modSel(coug.elk.hunt_fld)
  
  (coug.elk.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_elk_hunter_hunting_UMF, silent = TRUE))
  (coug.elk.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_elk_hunter_hunting_UMF, silent = TRUE))
  (coug.elk.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, coug_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL 
  (coug.elk.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpubish0, coug_elk_hunter_hunting_UMF, silent = TRUE)) 
  (coug.elk.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpubish1, coug_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (coug.elk.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coug_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL Error in optim(starts, nll, method = method, hessian = se, ...) : non-finite finite-difference value [19]
  (coug.elk.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coug_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  #' Including null models in model selection process b/c only one model with 
  #' site-level covariates did not fail to converge
  coug.elk.hunt_fl <- fitList(coug.elk.hunt_null0, coug.elk.hunt_null1, coug.elk.hunt_habpub0) 
  #' Model selection
  modSel(coug.elk.hunt_fl)
  summary(coug.elk.hunt_habpub0)
  
  ####  Cougar-White-tailed Deer-Hunter Hunting Season  ####
  (coug.wtd.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_hunter_hunting_UMF, silent = TRUE))
  (coug.wtd.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_wtd_hunter_hunting_UMF, silent = TRUE))
  (coug.wtd.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coug.wtd.hunt_fld <- fitList(coug.wtd.hunt_trail, coug.wtd.hunt_hunt, coug.wtd.hunt_pub)
  #' Model selection
  modSel(coug.wtd.hunt_fld)
  
  (coug.wtd.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunter_hunting_UMF, silent = TRUE))
  (coug.wtd.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_wtd_hunter_hunting_UMF, silent = TRUE))
  (coug.wtd.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, coug_wtd_hunter_hunting_UMF, silent = TRUE))  
  (coug.wtd.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpubish0, coug_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (coug.wtd.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpubish1, coug_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (coug.wtd.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coug_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (coug.wtd.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coug_wtd_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coug.wtd.hunt_fl <- fitList(coug.wtd.hunt_habpub0, coug.wtd.hunt_habpub1, coug.wtd.hunt_huntpub1, coug.wtd.hunt_huntpub2) 
  #' Model selection
  modSel(coug.elk.hunt_fl)
  summary(coug.wtd.hunt_habpub0)
  
  ####  Cougar-Moose-Hunter Hunting Season  ####
  (coug.moose.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_hunter_hunting_UMF, silent = TRUE))
  (coug.moose.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_moose_hunter_hunting_UMF, silent = TRUE))
  (coug.moose.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coug.moose.hunt_fld <- fitList(coug.moose.hunt_trail, coug.moose.hunt_hunt, coug.moose.hunt_pub)
  #' Model selection
  modSel(coug.moose.hunt_fld)
  
  (coug.moose.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunter_hunting_UMF, silent = TRUE))
  (coug.moose.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_moose_hunter_hunting_UMF, silent = TRUE))
  (coug.moose.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, coug_moose_hunter_hunting_UMF, silent = TRUE))  
  (coug.moose.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpubish0, coug_moose_hunter_hunting_UMF, silent = TRUE)) 
  (coug.moose.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpubish1, coug_moose_hunter_hunting_UMF, silent = TRUE)) 
  (coug.moose.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coug_moose_hunter_hunting_UMF, silent = TRUE))
  (coug.moose.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coug_moose_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coug.moose.hunt_fl <- fitList(coug.moose.hunt_habpub0, coug.moose.hunt_habpub1, coug.moose.hunt_huntpub1, coug.moose.hunt_huntpub2) 
  #' Model selection
  modSel(coug.moose.hunt_fl)
  summary(coug.moose.hunt_huntpub1)
  summary(coug.moose.hunt_huntpub2)

  
  
  ####  Wolf-Mule Deer-Hunter Hunting Season  ####
  (wolf.md.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_hunter_hunting_UMF, silent = TRUE))
  (wolf.md.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_md_hunter_hunting_UMF, silent = TRUE))
  (wolf.md.hunt_pubish <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolf.md.hunt_fld <- fitList(wolf.md.hunt_trail, wolf.md.hunt_hunt, wolf.md.hunt_pubish)
  #' Model selection
  modSel(wolf.md.hunt_fld)
  
  (wolf.md.hunt_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunter_hunting_UMF, silent = TRUE))
  (wolf.md.hunt_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_md_hunter_hunting_UMF, silent = TRUE))
  (wolf.md.hunt_null2 <- occuMulti(detFormulas_pubish, occFormulas_null2, wolf_md_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.md.hunt_habpub0 <- occuMulti(detFormulas_pubish, occFormulas_habpubish0, wolf_md_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.md.hunt_habpub1 <- occuMulti(detFormulas_pubish, occFormulas_habpubish1, wolf_md_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.md.hunt_huntpub1 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish1, wolf_md_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (wolf.md.hunt_huntpub2 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish2, wolf_md_hunter_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  wolf.md.hunt_fl <- fitList(wolf.md.hunt_habpub0, wolf.md.hunt_habpub1) 
  #' Model selection
  modSel(wolf.md.hunt_fl)
  summary(wolf.md.hunt_habpub1)
  summary(wolf.md.hunt_habpub2)
  
  
  ####  Wolf-Elk-Hunter Hunting Season  ####
  (wolf.elk.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_hunter_hunting_UMF, silent = TRUE))
  (wolf.elk.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_elk_hunter_hunting_UMF, silent = TRUE))
  (wolf.elk.hunt_pubish <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_elk_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolf.elk.hunt_fld <- fitList(wolf.elk.hunt_trail, wolf.elk.hunt_hunt, wolf.elk.hunt_pubish)
  #' Model selection
  modSel(wolf.elk.hunt_fld)
  
  (wolf.elk.hunt_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_elk_hunter_hunting_UMF, silent = TRUE))
  (wolf.elk.hunt_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_elk_hunter_hunting_UMF, silent = TRUE))
  (wolf.elk.hunt_null2 <- occuMulti(detFormulas_pubish, occFormulas_null2, wolf_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL 
  (wolf.elk.hunt_habpub0 <- occuMulti(detFormulas_pubish, occFormulas_habpubish0, wolf_elk_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.elk.hunt_habpub1 <- occuMulti(detFormulas_pubish, occFormulas_habpubish1, wolf_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (wolf.elk.hunt_huntpub1 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish1, wolf_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL Error in optim(starts, nll, method = method, hessian = se, ...) : non-finite finite-difference value [19]
  (wolf.elk.hunt_huntpub2 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish2, wolf_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  #' Including null models in model selection process b/c only one model with 
  #' site-level covariates did not fail to converge
  wolf.elk.hunt_fl <- fitList(wolf.elk.hunt_null0, wolf.elk.hunt_null1, wolf.elk.hunt_habpub0) 
  #' Model selection
  modSel(wolf.elk.hunt_fl)
  summary(wolf.elk.hunt_habpub0)
  
  ####  Wolf-White-tailed Deer-Hunter Hunting Season  ####
  (wolf.wtd.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_hunter_hunting_UMF, silent = TRUE))
  (wolf.wtd.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_wtd_hunter_hunting_UMF, silent = TRUE))
  (wolf.wtd.hunt_pubish <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolf.wtd.hunt_fld <- fitList(wolf.elk.hunt_trail, wolf.elk.hunt_hunt, wolf.elk.hunt_pubish)
  #' Model selection
  modSel(wolf.wtd.hunt_fld)
  
  (wolf.wtd.hunt_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunter_hunting_UMF, silent = TRUE))
  (wolf.wtd.hunt_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_wtd_hunter_hunting_UMF, silent = TRUE))
  (wolf.wtd.hunt_null2 <- occuMulti(detFormulas_pubish, occFormulas_null2, wolf_wtd_hunter_hunting_UMF, silent = TRUE))  
  (wolf.wtd.hunt_habpub0 <- occuMulti(detFormulas_pubish, occFormulas_habpubish0, wolf_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.wtd.hunt_habpub1 <- occuMulti(detFormulas_pubish, occFormulas_habpubish1, wolf_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.wtd.hunt_huntpub1 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish1, wolf_wtd_hunter_hunting_UMF, silent = TRUE))
  (wolf.wtd.hunt_huntpub2 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish2, wolf_wtd_hunter_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  wolf.wtd.hunt_fl <- fitList(wolf.wtd.hunt_habpub0, wolf.wtd.hunt_habpub1, wolf.wtd.hunt_huntpub1) 
  #' Model selection
  modSel(wolf.wtd.hunt_fl)
  summary(wolf.wtd.hunt_huntpub1)
  summary(wolf.wtd.hunt_habpub1)
  
  ####  Wolf-Moose-Hunter Hunting Season  ####
  (wolf.moose.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_hunter_hunting_UMF, silent = TRUE))
  (wolf.moose.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_moose_hunter_hunting_UMF, silent = TRUE))
  (wolf.moose.hunt_pubish <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolf.moose.hunt_fld <- fitList(wolf.moose.hunt_trail, wolf.moose.hunt_hunt, wolf.moose.hunt_pubish)
  #' Model selection
  modSel(wolf.moose.hunt_fld)
  
  (wolf.moose.hunt_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunter_hunting_UMF, silent = TRUE))
  (wolf.moose.hunt_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_moose_hunter_hunting_UMF, silent = TRUE))
  (wolf.moose.hunt_null2 <- occuMulti(detFormulas_pubish, occFormulas_null2, wolf_moose_hunter_hunting_UMF, silent = TRUE))  
  (wolf.moose.hunt_habpub0 <- occuMulti(detFormulas_pubish, occFormulas_habpubish0, wolf_moose_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.moose.hunt_habpub1 <- occuMulti(detFormulas_pubish, occFormulas_habpubish1, wolf_moose_hunter_hunting_UMF, silent = TRUE)) 
  (wolf.moose.hunt_huntpub1 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish1, wolf_moose_hunter_hunting_UMF, silent = TRUE)) # FAIL Error in optim(starts, nll, method = method, hessian = se, ...) : non-finite finite-difference value [19]
  (wolf.moose.hunt_huntpub2 <- occuMulti(detFormulas_pubish, occFormulas_huntpubish2, wolf_moose_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolf.moose.hunt_fl <- fitList(wolf.moose.hunt_habpub0, wolf.moose.hunt_habpub1) 
  #' Model selection
  modSel(wolf.moose.hunt_fl)
  summary(wolf.moose.hunt_habpub1)
  
  
  ####  Bear-Mule Deer-Hunter Hunting Season  ####
  (bear.md.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_hunter_hunting_UMF, silent = TRUE))
  (bear.md.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_md_hunter_hunting_UMF, silent = TRUE))
  (bear.md.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bear.md.hunt_fld <- fitList(bear.md.hunt_trail, bear.md.hunt_hunt, bear.md.hunt_pub)
  #' Model selection
  modSel(bear.md.hunt_fld)
  
  (bear.md.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunter_hunting_UMF, silent = TRUE))
  (bear.md.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_md_hunter_hunting_UMF, silent = TRUE))
  (bear.md.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, bear_md_hunter_hunting_UMF, silent = TRUE)) 
  (bear.md.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, bear_md_hunter_hunting_UMF, silent = TRUE)) 
  (bear.md.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, bear_md_hunter_hunting_UMF, silent = TRUE)) 
  (bear.md.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, bear_md_hunter_hunting_UMF, silent = TRUE)) 
  (bear.md.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, bear_md_hunter_hunting_UMF, silent = TRUE))
  (bear.md.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, bear_md_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bear.md.hunt_fl <- fitList(bear.md.hunt_habpub0, bear.md.hunt_habpub1, bear.md.hunt_habpub2, bear.md.hunt_huntpub1, bear.md.hunt_huntpub2) 
  #' Model selection
  modSel(bear.md.hunt_fl)
  summary(bear.md.hunt_habpub0)
  
  ####  Bear-Elk-Hunter Hunting Season  ####
  (bear.elk.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_hunter_hunting_UMF, silent = TRUE))
  (bear.elk.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_elk_hunter_hunting_UMF, silent = TRUE))
  (bear.elk.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_elk_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bear.elk.hunt_fld <- fitList(bear.elk.hunt_trail, bear.elk.hunt_hunt, bear.elk.hunt_pub)
  #' Model selection
  modSel(bear.elk.hunt_fld)
  
  (bear.elk.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_elk_hunter_hunting_UMF, silent = TRUE))
  (bear.elk.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_elk_hunter_hunting_UMF, silent = TRUE))
  (bear.elk.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, bear_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (bear.elk.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, bear_elk_hunter_hunting_UMF, silent = TRUE)) 
  (bear.elk.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, bear_elk_hunter_hunting_UMF, silent = TRUE)) 
  (bear.elk.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, bear_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (bear.elk.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, bear_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (bear.elk.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, bear_elk_hunter_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  bear.elk.hunt_fl <- fitList(bear.elk.hunt_habpub0, bear.elk.hunt_habpub1) 
  #' Model selection
  modSel(bear.elk.hunt_fl)
  summary(bear.elk.hunt_habpub1) # but f12s aren't significant at 0.05
  
  
  ####  Bear-White-tailed Deer-Hunter Hunting Season  ####
  (bear.wtd.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_hunter_hunting_UMF, silent = TRUE))
  (bear.wtd.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_wtd_hunter_hunting_UMF, silent = TRUE))
  (bear.wtd.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bear.wtd.hunt_fld <- fitList(bear.wtd.hunt_trail, bear.wtd.hunt_hunt, bear.wtd.hunt_pub)
  #' Model selection
  modSel(bear.wtd.hunt_fld)
  
  (bear.wtd.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunter_hunting_UMF, silent = TRUE))
  (bear.wtd.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_wtd_hunter_hunting_UMF, silent = TRUE))
  (bear.wtd.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, bear_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bear.wtd.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, bear_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bear.wtd.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, bear_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bear.wtd.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, bear_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bear.wtd.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, bear_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bear.wtd.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, bear_wtd_hunter_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  bear.wtd.hunt_fl <- fitList(bear.wtd.hunt_habpub0, bear.wtd.hunt_habpub1, bear.wtd.hunt_huntpub1, bear.wtd.hunt_huntpub2) 
  #' Model selection
  modSel(bear.wtd.hunt_fl)
  summary(bear.wtd.hunt_habpub0)
  summary(bear.wtd.hunt_habpub1)
  
  ####  Bear-Moose-Hunter Hunting Season  ####
  (bear.moose.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_hunter_hunting_UMF, silent = TRUE))
  (bear.moose.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_moose_hunter_hunting_UMF, silent = TRUE))
  (bear.moose.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bear.moose.hunt_fld <- fitList(bear.moose.hunt_trail, bear.moose.hunt_hunt, bear.moose.hunt_pub)
  #' Model selection
  modSel(bear.moose.hunt_fld)
  
  (bear.moose.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunter_hunting_UMF, silent = TRUE))
  (bear.moose.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_moose_hunter_hunting_UMF, silent = TRUE))
  (bear.moose.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, bear_moose_hunter_hunting_UMF, silent = TRUE)) 
  (bear.moose.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, bear_moose_hunter_hunting_UMF, silent = TRUE)) 
  (bear.moose.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, bear_moose_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (bear.moose.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, bear_moose_hunter_hunting_UMF, silent = TRUE)) # FAIL
  (bear.moose.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, bear_moose_hunter_hunting_UMF, silent = TRUE)) # FAIL 
  (bear.moose.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, bear_moose_hunter_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  #' Including null models in model selection process b/c only one model with 
  #' site-level covariates did not fail to converge
  bear.moose.hunt_fl <- fitList(bear.moose.hunt_null0, bear.moose.hunt_null1, bear.moose.hunt_null2, bear.moose.hunt_habpub0) 
  #' Model selection
  modSel(bear.moose.hunt_fl)
  summary(bear.moose.hunt_habpub0)
  
  
  ####  Coyote-Mule Deer-Hunters Hunting Season  ####
  (coy.md.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_hunter_hunting_UMF, silent = TRUE))
  (coy.md.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_md_hunter_hunting_UMF, silent = TRUE))
  (coy.md.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coy.md.hunt_fld <- fitList(coy.md.hunt_trail, coy.md.hunt_hunt, coy.md.hunt_pub)
  #' Model selection
  modSel(coy.md.hunt_fld)
  
  (coy.md.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunter_hunting_UMF, silent = TRUE))
  (coy.md.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_md_hunter_hunting_UMF, silent = TRUE))
  (coy.md.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, coy_md_hunter_hunting_UMF, silent = TRUE))  
  (coy.md.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, coy_md_hunter_hunting_UMF, silent = TRUE)) 
  (coy.md.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, coy_md_hunter_hunting_UMF, silent = TRUE)) 
  (coy.md.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, coy_md_hunter_hunting_UMF, silent = TRUE)) 
  (coy.md.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, coy_md_hunter_hunting_UMF, silent = TRUE)) 
  (coy.md.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, coy_md_hunter_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  coy.md.hunt_fl <- fitList(coy.md.hunt_habpub0, coy.md.hunt_habpub1, coy.md.hunt_habpub2, coy.md.hunt_huntpub1, coy.md.hunt_huntpub2) 
  #' Model selection
  modSel(coy.md.hunt_fl)
  summary(coy.md.hunt_huntpub1)
  summary(coy.md.hunt_huntpub2)
  summary(coy.md.hunt_habpub1)
  
  ####  Coyote-White-tailed Deer-Hunters Hunting Season  ####
  (coy.wtd.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_hunter_hunting_UMF, silent = TRUE))
  (coy.wtd.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_wtd_hunter_hunting_UMF, silent = TRUE))
  (coy.wtd.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  coy.wtd.hunt_fld <- fitList(coy.wtd.hunt_trail, coy.wtd.hunt_hunt, coy.wtd.hunt_pub)
  #' Model selection
  modSel(coy.wtd.hunt_fld)
  
  (coy.wtd.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunter_hunting_UMF, silent = TRUE))
  (coy.wtd.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_wtd_hunter_hunting_UMF, silent = TRUE))
  (coy.wtd.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, coy_wtd_hunter_hunting_UMF, silent = TRUE))  
  (coy.wtd.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, coy_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (coy.wtd.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, coy_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (coy.wtd.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, coy_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (coy.wtd.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, coy_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (coy.wtd.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, coy_wtd_hunter_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  coy.wtd.hunt_fl <- fitList(coy.wtd.hunt_habpub0, coy.wtd.hunt_habpub1, coy.wtd.hunt_habpub2, coy.wtd.hunt_huntpub1, coy.wtd.hunt_huntpub2) 
  #' Model selection
  modSel(coy.wtd.hunt_fl)
  summary(coy.wtd.hunt_habpub1)
  summary(coy.wtd.hunt_huntpub1)
  summary(coy.wtd.hunt_habpub2)
  
  
  ####  Bobcat-Mule Deer-Hunters Hunting Season  ####
  (bob.md.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_hunter_hunting_UMF, silent = TRUE))
  (bob.md.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_md_hunter_hunting_UMF, silent = TRUE))
  (bob.md.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bob.md.hunt_fld <- fitList(bob.md.hunt_trail, bob.md.hunt_hunt, bob.md.hunt_pub)
  #' Model selection
  modSel(bob.md.hunt_fld)
  
  (bob.md.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunter_hunting_UMF, silent = TRUE))
  (bob.md.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_md_hunter_hunting_UMF, silent = TRUE))
  (bob.md.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, bob_md_hunter_hunting_UMF, silent = TRUE))  
  (bob.md.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, bob_md_hunter_hunting_UMF, silent = TRUE)) 
  (bob.md.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, bob_md_hunter_hunting_UMF, silent = TRUE)) 
  (bob.md.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, bob_md_hunter_hunting_UMF, silent = TRUE)) 
  (bob.md.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, bob_md_hunter_hunting_UMF, silent = TRUE)) 
  (bob.md.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, bob_md_hunter_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  bob.md.hunt_fl <- fitList(bob.md.hunt_habpub0, bob.md.hunt_habpub1, bob.md.hunt_habpub2, bob.md.hunt_huntpub1, bob.md.hunt_huntpub2) 
  #' Model selection
  modSel(bob.md.hunt_fl)
  summary(bob.md.hunt_huntpub1)
  summary(bob.md.hunt_habpub0)
  summary(bob.md.hunt_habpub1)
  summary(bob.md.hunt_huntpub2)
  
  
  
  ####  Bobcat-White-tailed Deer-Hunters Hunting Season  ####
  (bob.wtd.hunt_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_hunter_hunting_UMF, silent = TRUE))
  (bob.wtd.hunt_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_wtd_hunter_hunting_UMF, silent = TRUE))
  (bob.wtd.hunt_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunter_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bob.wtd.hunt_fld <- fitList(bob.wtd.hunt_trail, bob.wtd.hunt_hunt, bob.wtd.hunt_pub)
  #' Model selection
  modSel(bob.wtd.hunt_fld)
  
  (bob.wtd.hunt_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunter_hunting_UMF, silent = TRUE))
  (bob.wtd.hunt_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_wtd_hunter_hunting_UMF, silent = TRUE))
  (bob.wtd.hunt_null2 <- occuMulti(detFormulas_pub, occFormulas_null2, bob_wtd_hunter_hunting_UMF, silent = TRUE))  
  (bob.wtd.hunt_habpub0 <- occuMulti(detFormulas_pub, occFormulas_habpub0, bob_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bob.wtd.hunt_habpub1 <- occuMulti(detFormulas_pub, occFormulas_habpub1, bob_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bob.wtd.hunt_habpub2 <- occuMulti(detFormulas_pub, occFormulas_habpub2, bob_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bob.wtd.hunt_huntpub1 <- occuMulti(detFormulas_pub, occFormulas_huntpub1, bob_wtd_hunter_hunting_UMF, silent = TRUE)) 
  (bob.wtd.hunt_huntpub2 <- occuMulti(detFormulas_pub, occFormulas_huntpub2, bob_wtd_hunter_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  bob.wtd.hunt_fl <- fitList(bob.wtd.hunt_habpub0, bob.wtd.hunt_habpub1, bob.wtd.hunt_habpub2, bob.wtd.hunt_huntpub1, bob.wtd.hunt_huntpub2) 
  #' Model selection
  modSel(bob.wtd.hunt_fl)
  summary(bob.wtd.hunt_habpub1)
  summary(bob.wtd.hunt_habpub2)
  
  
  #' Save model outputs in one giant R image
  save.image(file = "./Outputs/MultiSpp_CoOcc_Models.RData")
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info
  
  #'  Function to save occupancy results
  rounddig <- 2
  occ_out <- function(mod, spp1, spp2, spp3, season) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Species3 = rep(spp3, nrow(.)),
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
      relocate(Species3, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
    
    return(out)
  }
  
  #'  Run each model through function
  #'  Grazing season models
  occ_coug.md.cow <- occ_out(coug.md.cow_hab0, "Cougar", "Mule Deer", "Cattle", "Grazing")
  occ_coug.elk.cow <- occ_out(coug.elk.cow_hab0, "Cougar", "Elk", "Cattle", "Grazing")
  occ_coug.wtd.cow <- occ_out(coug.wtd.cow_grz, "Cougar", "White-tailed Deer", "Cattle", "Grazing")
  occ_coug.moose.cow <- occ_out(coug.moose.cow_grz, "Cougar", "Moose", "Cattle", "Grazing")
  occ_wolf.md.cow <- occ_out(wolf.md.cow_grz, "Wolf", "Mule Deer", "Cattle", "Grazing")
  occ_wolf.elk.cow <- occ_out(wolf.elk.cow_hab0, "Wolf", "Elk", "Cattle", "Grazing")
  occ_wolf.wtd.cow <- occ_out(wolf.wtd.cow_grz, "Wolf", "White-tailed Deer", "Cattle", "Grazing")
  occ_wolf.moose.cow <- occ_out(wolf.moose.cow_hab0, "Wolf", "Moose", "Cattle", "Grazing")
  occ_bear.md.cow <- occ_out(bear.md.cow_hab0, "Black Bear", "Mule Deer", "Cattle", "Grazing")
  occ_bear.elk.cow <- occ_out(bear.elk.cow_hab0, "Black Bear", "Elk", "Cattle", "Grazing")
  occ_bear.wtd.cow <- occ_out(bear.wtd.cow_grz, "Black Bear", "White-tailed Deer", "Cattle", "Grazing")
  occ_bear.moose.cow <- occ_out(bear.moose.cow_grzpub, "Black Bear", "Moose", "Cattle", "Grazing")
  occ_coy.md.cow <- occ_out(coy.md.cow_grz, "Coyote", "Mule Deer", "Cattle", "Grazing")
  occ_coy.wtd.cow <- occ_out(coy.wtd.cow_grzpub, "Coyote", "White-tailed Deer", "Cattle", "Grazing")
  occ_bob.md.cow <- occ_out(bob.md.cow_grzpub, "Bobcat", "Mule Deer Deer", "Cattle", "Grazing")
  occ_bob.wtd.cow <- occ_out(bob.wtd.cow_grz, "Bobcat", "White-tailed Deer", "Cattle", "Grazing")
  #'  Hunting season models
  occ_coug.md.hunt <- occ_out(coug.md.hunt_huntpub1, "Cougar", "Mule Deer", "Hunter", "Hunting")
  occ_coug.elk.hunt <- occ_out(coug.elk.hunt_habpub0, "Cougar", "Elk", "Hunter", "Hunting")
  occ_coug.wtd.hunt <- occ_out(coug.wtd.hunt_habpub0, "Cougar", "White-tailed Deer", "Hunter", "Hunting")
  occ_coug.moose.hunt <- occ_out(coug.moose.hunt_huntpub1, "Cougar", "Moose", "Hunter", "Hunting")
  occ_wolf.md.hunt <- occ_out(wolf.md.hunt_habpub1, "Wolf", "Mule Deer", "Hunter", "Hunting")
  occ_wolf.elk.hunt <- occ_out(wolf.elk.hunt_habpub0, "Wolf", "Elk", "Hunter", "Hunting")
  occ_wolf.wtd.hunt <- occ_out(wolf.wtd.hunt_huntpub1, "Wolf", "White-tailed Deer", "Hunter", "Hunting")
  occ_wolf.moose.hunt <- occ_out(wolf.moose.hunt_habpub1, "Wolf", "Moose", "Hunter", "Hunting")
  occ_bear.md.hunt <- occ_out(bear.md.hunt_habpub0, "Black Bear", "Mule Deer", "Hunter", "Hunting")
  occ_bear.elk.hunt <- occ_out(bear.elk.hunt_habpub1, "Black Bear", "Elk", "Hunter", "Hunting")
  occ_bear.wtd.hunt <- occ_out(bear.wtd.hunt_habpub0, "Black Bear", "White-tailed Deer", "Hunter", "Hunting")
  occ_bear.moose.hunt <- occ_out(bear.moose.hunt_habpub0, "Black Bear", "Moose", "Hunter", "Hunting")
  occ_coy.md.hunt <- occ_out(coy.md.hunt_huntpub1, "Coyote", "Mule Deer", "Hunter", "Hunting")
  occ_coy.wtd.hunt <- occ_out(coy.wtd.hunt_habpub1, "Coyote", "White-tailed Deer", "Hunter", "Hunting")
  occ_bob.md.hunt <- occ_out(bob.md.hunt_huntpub1, "Bobcat", "Mule Deer Deer", "Hunter", "Hunting")
  occ_bob.wtd.hunt <- occ_out(bob.wtd.hunt_habpub1, "Bobcat", "White-tailed Deer", "Hunter", "Hunting")
  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  graze_results <- rbind(occ_coug.md.cow, occ_coug.elk.cow, occ_coug.wtd.cow, occ_coug.moose.cow,
                         occ_wolf.md.cow, occ_wolf.elk.cow, occ_wolf.wtd.cow, occ_wolf.moose.cow,
                         occ_bear.md.cow, occ_bear.elk.cow, occ_bear.wtd.cow, occ_bear.moose.cow,
                         occ_coy.md.cow, occ_coy.wtd.cow, occ_bob.md.cow, occ_bob.wtd.cow)
  hunt_results <- rbind(occ_coug.md.hunt, occ_coug.elk.hunt, occ_coug.wtd.hunt, occ_coug.moose.hunt,
                        occ_wolf.md.hunt, occ_wolf.elk.hunt, occ_wolf.wtd.hunt, occ_wolf.moose.hunt,
                        occ_bear.md.hunt, occ_bear.elk.hunt, occ_bear.wtd.hunt, occ_bear.moose.hunt,
                        occ_coy.md.hunt, occ_coy.wtd.hunt, occ_bob.md.hunt, occ_bob.wtd.hunt)
  
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
      Parameter = gsub("moose", "Species 2", Parameter),
      Parameter = gsub("cattle", "Species 3", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 3] Study_AreaOK") %>%
    relocate("[Species 1:Species 3] (Intercept)", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 2:Species 3] (Intercept)", .after = "[Species 1:Species 3] (Intercept)") %>%
    relocate("[Species 1:Species 2] GrazingActivity", .after = "[Species 1:Species 2] (Intercept)") %>%
    relocate("[Species 1:Species 2] Public1", .after = "[Species 1:Species 2] GrazingActivity") %>%
    # relocate("[Species 2:Species 3] GrazingActivity", .after = "[Species 2:Species 3] (Intercept)") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    relocate("[Species 3] I(Elev^2)", .after = "[Species 3] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 3] (Intercept)", c("[Species 3] Intercept (SE)", "[Species 3] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 3] Elev", c("[Species 3] Elevation (SE)", "[Species 3] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 3] I(Elev^2)", c("[Species 3] Elevation^2 (SE)", "[Species 3] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 3] PercForest", c("[Species 3] PercentForest (SE)", "[Species 3] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 3] Study_AreaOK", c("[Species 3] Study_AreaOK (SE)", "[Species 3] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 3] (Intercept)", c("[Species 1:Species 3] Intercept (SE)", "[Species 1:Species 3] Intercept Pval"), sep = "_") %>%
    separate("[Species 2:Species 3] (Intercept)", c("[Species 2:Species 3] Intercept (SE)", "[Species 2:Species 3] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] GrazingActivity", c("[Species 1:Species 2] GrazingActivity (SE)", "[Species 1:Species 2] GrazingActivity Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] Public1", c("[Species 1:Species 2] Public1 (SE)", "[Species 1:Species 2] Public1 Pval"), sep = "_") %>%
    # separate("[Species 2:Species 3] GrazingActivity", c("[Species 2:Species 3] GrazingActivity (SE)", "[Species 2:Species 3] GrazingActivity Pval"), sep = "_") %>%
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
      Parameter = gsub("moose", "Species 2", Parameter),
      Parameter = gsub("hunter", "Species 3", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    relocate("[Species 1:Species 2] (Intercept)", .after = "[Species 3] Study_AreaOK") %>%
    relocate("[Species 1:Species 2] HuntingActivity", .after = "[Species 1:Species 2] (Intercept)") %>%
    # relocate("[Species 1:Species 2] Public1", .after = "[Species 1:Species 2] HuntingActivity") %>%
    relocate("[Species 1:Species 3] (Intercept)", .after = "[Species 1:Species 2] HuntingActivity") %>%
    relocate("[Species 2:Species 3] (Intercept)", .after = "[Species 1:Species 2] HuntingActivity") %>%
    relocate("[Species 1] I(Elev^2)", .after = "[Species 1] Elev") %>%
    relocate("[Species 2] I(Elev^2)", .after = "[Species 2] Elev") %>%
    relocate("[Species 3] I(Elev^2)", .after = "[Species 3] Elev") %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 3] (Intercept)", c("[Species 3] Intercept (SE)", "[Species 3] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] Elev", c("[Species 1] Elevation (SE)", "[Species 1] Elevation Pval"), sep = "_") %>%
    separate("[Species 2] Elev", c("[Species 2] Elevation (SE)", "[Species 2] Elevation Pval"), sep = "_") %>%
    separate("[Species 3] Elev", c("[Species 3] Elevation (SE)", "[Species 3] Elevation Pval"), sep = "_") %>%
    separate("[Species 1] I(Elev^2)", c("[Species 1] Elevation^2 (SE)", "[Species 1] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 2] I(Elev^2)", c("[Species 2] Elevation^2 (SE)", "[Species 2] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 3] I(Elev^2)", c("[Species 3] Elevation^2 (SE)", "[Species 3] Elevation^2 Pval"), sep = "_") %>%
    separate("[Species 1] PercForest", c("[Species 1] PercentForest (SE)", "[Species 1] PercentForest Pval"), sep = "_") %>%
    separate("[Species 2] PercForest", c("[Species 2] PercentForest (SE)", "[Species 2] PercentForest Pval"), sep = "_") %>%
    separate("[Species 3] PercForest", c("[Species 3] PercentForest (SE)", "[Species 3] PercentForest Pval"), sep = "_") %>%
    separate("[Species 1] Study_AreaOK", c("[Species 1] Study_AreaOK (SE)", "[Species 1] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 2] Study_AreaOK", c("[Species 2] Study_AreaOK (SE)", "[Species 2] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 3] Study_AreaOK", c("[Species 3] Study_AreaOK (SE)", "[Species 3] Study_AreaOK Pval"), sep = "_") %>%
    separate("[Species 1] Public1", c("[Species 1] Public1 (SE)", "[Species 1] Public1 Pval"), sep = "_") %>%
    separate("[Species 2] Public1", c("[Species 2] Public1 (SE)", "[Species 2] Public1 Pval"), sep = "_") %>%
    separate("[Species 3] Public1", c("[Species 3] Public1 (SE)", "[Species 3] Public1 Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] (Intercept)", c("[Species 1:Species 2] Intercept (SE)", "[Species 1:Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 1:Species 2] HuntingActivity", c("[Species 1:Species 2] HuntingActivity (SE)", "[Species 1:Species 2] HuntingActivity Pval"), sep = "_") %>%
    # separate("[Species 1:Species 2] Public1", c("[Species 1:Species 2] Public1 (SE)", "[Species 1:Species 2] Public1 Pval"), sep = "_") %>%
    separate("[Species 1:Species 3] (Intercept)", c("[Species 1:Species 3] Intercept (SE)", "[Species 1:Species 3] Intercept Pval"), sep = "_") %>%
    separate("[Species 2:Species 3] (Intercept)", c("[Species 2:Species 3] Intercept (SE)", "[Species 2:Species 3] Intercept Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  #'  Save!
  write.csv(results_graze, paste0("./Outputs/Tables/CoOcc_OccProb_GrazingResults_", Sys.Date(), ".csv"))  
  write.csv(results_graze_wide, paste0("./Outputs/Tables/CoOcc_OccProb_GrazingResults_wide_", Sys.Date(), ".csv"))
  write.csv(results_hunt, paste0("./Outputs/Tables/CoOcc_OccProb_HuntingResults_", Sys.Date(), ".csv"))  
  write.csv(results_hunt_wide, paste0("./Outputs/Tables/CoOcc_OccProb_HuntingResults_wide_", Sys.Date(), ".csv"))
  
  
  #'  Function to save detection results
  det_out <- function(mod, spp1, spp2, spp3, season) {
    out <- summary(mod@estimates)$det %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$det),
        Species1 = rep(spp1, nrow(.)),
        Species2 = rep(spp2, nrow(.)),
        Species3 = rep(spp3, nrow(.)),
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
      relocate(Species3, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
    return(out)
  }
  
  
  #'  Run each model through detection function
  #'  Grazing season models
  det_coug.md.cow <- det_out(coug.md.cow_hab0, "Cougar", "Mule Deer", "Cattle", "Grazing")
  det_coug.elk.cow <- det_out(coug.elk.cow_hab0, "Cougar", "Elk", "Cattle", "Grazing")
  det_coug.wtd.cow <- det_out(coug.wtd.cow_grz, "Cougar", "White-tailed Deer", "Cattle", "Grazing")
  det_coug.moose.cow <- det_out(coug.moose.cow_grz, "Cougar", "Moose", "Cattle", "Grazing")
  det_wolf.md.cow <- det_out(wolf.md.cow_grz, "Wolf", "Mule Deer", "Cattle", "Grazing")
  det_wolf.elk.cow <- det_out(wolf.elk.cow_hab0, "Wolf", "Elk", "Cattle", "Grazing")
  det_wolf.wtd.cow <- det_out(wolf.wtd.cow_grz, "Wolf", "White-tailed Deer", "Cattle", "Grazing")
  det_wolf.moose.cow <- det_out(wolf.moose.cow_hab0, "Wolf", "Moose", "Cattle", "Grazing")
  det_bear.md.cow <- det_out(bear.md.cow_hab0, "Black Bear", "Mule Deer", "Cattle", "Grazing")
  det_bear.elk.cow <- det_out(bear.elk.cow_hab0, "Black Bear", "Elk", "Cattle", "Grazing")
  det_bear.wtd.cow <- det_out(bear.wtd.cow_grz, "Black Bear", "White-tailed Deer", "Cattle", "Grazing")
  det_bear.moose.cow <- det_out(bear.moose.cow_grzpub, "Black Bear", "Moose", "Cattle", "Grazing")
  det_coy.md.cow <- det_out(coy.md.cow_grz, "Coyote", "Mule Deer", "Cattle", "Grazing")
  det_coy.wtd.cow <- det_out(coy.wtd.cow_grzpub, "Coyote", "White-tailed Deer", "Cattle", "Grazing")
  det_bob.md.cow <- det_out(bob.md.cow_grzpub, "Bobcat", "Mule Deer Deer", "Cattle", "Grazing")
  det_bob.wtd.cow <- det_out(bob.wtd.cow_grz, "Bobcat", "White-tailed Deer", "Cattle", "Grazing")
  #'  Hunting season models
  det_coug.md.hunt <- det_out(coug.md.hunt_huntpub1, "Cougar", "Mule Deer", "Hunter", "Hunting")
  det_coug.elk.hunt <- det_out(coug.elk.hunt_habpub0, "Cougar", "Elk", "Hunter", "Hunting")
  det_coug.wtd.hunt <- det_out(coug.wtd.hunt_habpub0, "Cougar", "White-tailed Deer", "Hunter", "Hunting")
  det_coug.moose.hunt <- det_out(coug.moose.hunt_huntpub1, "Cougar", "Moose", "Hunter", "Hunting")
  det_wolf.md.hunt <- det_out(wolf.md.hunt_habpub1, "Wolf", "Mule Deer", "Hunter", "Hunting")
  det_wolf.elk.hunt <- det_out(wolf.elk.hunt_habpub0, "Wolf", "Elk", "Hunter", "Hunting")
  det_wolf.wtd.hunt <- det_out(wolf.wtd.hunt_huntpub1, "Wolf", "White-tailed Deer", "Hunter", "Hunting")
  det_wolf.moose.hunt <- det_out(wolf.moose.hunt_habpub1, "Wolf", "Moose", "Hunter", "Hunting")
  det_bear.md.hunt <- det_out(bear.md.hunt_habpub0, "Black Bear", "Mule Deer", "Hunter", "Hunting")
  det_bear.elk.hunt <- det_out(bear.elk.hunt_habpub1, "Black Bear", "Elk", "Hunter", "Hunting")
  det_bear.wtd.hunt <- det_out(bear.wtd.hunt_habpub0, "Black Bear", "White-tailed Deer", "Hunter", "Hunting")
  det_bear.moose.hunt <- det_out(bear.moose.hunt_habpub0, "Black Bear", "Moose", "Hunter", "Hunting")
  det_coy.md.hunt <- det_out(coy.md.hunt_huntpub1, "Coyote", "Mule Deer", "Hunter", "Hunting")
  det_coy.wtd.hunt <- det_out(coy.wtd.hunt_habpub1, "Coyote", "White-tailed Deer", "Hunter", "Hunting")
  det_bob.md.hunt <- det_out(bob.md.hunt_huntpub1, "Bobcat", "Mule Deer Deer", "Hunter", "Hunting")
  det_bob.wtd.hunt <- det_out(bob.wtd.hunt_habpub1, "Bobcat", "White-tailed Deer", "Hunter", "Hunting")
  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  graze_det_results <- rbind(det_coug.md.cow, det_coug.elk.cow, det_coug.wtd.cow, det_coug.moose.cow,
                         det_wolf.md.cow, det_wolf.elk.cow, det_wolf.wtd.cow, det_wolf.moose.cow,
                         det_bear.md.cow, det_bear.elk.cow, det_bear.wtd.cow, det_bear.moose.cow,
                         det_coy.md.cow, det_coy.wtd.cow, det_bob.md.cow, det_bob.wtd.cow)
  hunt_det_results <- rbind(det_coug.md.hunt, det_coug.elk.hunt, det_coug.wtd.hunt, det_coug.moose.hunt,
                        det_wolf.md.hunt, det_wolf.elk.hunt, det_wolf.wtd.hunt, det_wolf.moose.hunt,
                        det_bear.md.hunt, det_bear.elk.hunt, det_bear.wtd.hunt, det_bear.moose.hunt,
                        det_coy.md.hunt, det_coy.wtd.hunt, det_bob.md.hunt, det_bob.wtd.hunt)
 
  
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
      Parameter = gsub("moose", "Species 2", Parameter),
      Parameter = gsub("cattle", "Species 3", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 3] (Intercept)", c("[Species 3] Intercept (SE)", "[Species 3] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 3] TrailDirt road", c("[Species 3] Dirt road (SE)", "[Species 3] Dirt road Pval"), sep = "_") %>%
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
      Parameter = gsub("moose", "Species 2", Parameter),
      Parameter = gsub("hunter", "Species 3", Parameter)
    ) %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("[Species 1] (Intercept)", c("[Species 1] Intercept (SE)", "[Species 1] Intercept Pval"), sep = "_") %>%
    separate("[Species 2] (Intercept)", c("[Species 2] Intercept (SE)", "[Species 2] Intercept Pval"), sep = "_") %>%
    separate("[Species 3] (Intercept)", c("[Species 3] Intercept (SE)", "[Species 3] Intercept Pval"), sep = "_") %>%
    separate("[Species 1] TrailDirt road", c("[Species 1] Dirt road (SE)", "[Species 1] Dirt road Pval"), sep = "_") %>%
    separate("[Species 2] TrailDirt road", c("[Species 2] Dirt road (SE)", "[Species 2] Dirt road Pval"), sep = "_") %>%
    separate("[Species 3] TrailDirt road", c("[Species 3] Dirt road (SE)", "[Species 3] Dirt road Pval"), sep = "_") %>%
    separate("[Species 1] Public1", c("[Species 1] Public1 (SE)", "[Species 1] Public1 Pval"), sep = "_") %>%
    separate("[Species 2] Public1", c("[Species 2] Public1 (SE)", "[Species 2] Public1 Pval"), sep = "_") %>%
    separate("[Species 3] Public1", c("[Species 3] Public1 (SE)", "[Species 3] Public1 Pval"), sep = "_") %>%
    arrange(match(Species1, c("Black Bear", "Bobcat", "Cougar", "Coyote", "Wolf"))) 
  
  #'  Save!
  write.csv(results_det_graze, paste0("./Outputs/Tables/CoOcc_DetProb_GrazingResults_", Sys.Date(), ".csv"))  
  write.csv(results_det_graze_wide, paste0("./Outputs/Tables/CoOcc_DetProb_GrazingResults_wide", Sys.Date(), ".csv"))
  write.csv(results_det_hunt, paste0("./Outputs/Tables/CoOcc_DetProb_HuntingResults_", Sys.Date(), ".csv"))  
  write.csv(results_det_hunt_wide, paste0("./Outputs/Tables/CoOcc_DetProb_HuntingResults_wide", Sys.Date(), ".csv"))
  
  
  
