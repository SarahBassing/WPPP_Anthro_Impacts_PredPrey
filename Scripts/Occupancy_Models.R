  #'  ============================================
  #'  Occupancy Models (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  Februray 2021
  #'  ============================================
  #'  Script to create unmarked data frames and run single-species, single-season
  #'  occupancy models for deer, elk, cougars, wolves, coyotes, and bobcats for
  #'  summer 2018 (7/1/18 - 9/29/18) and winter 2018-2019 (12/1/18 - 3/1/19),
  #'  respectively. Each single-season occupancy model includes 13 7-day sampling 
  #'  occasions comprising the warmest months in summer and coldest months in 
  #'  winter with the most consistent snow.
  #'  
  #'  Encounter histories are generated with the CameratTrap_DetectionHistories.R
  #'  script. Covariate data included in occupancy models were collected at each
  #'  camera site or extracted from remotely sensed data.
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
  
  ####  COUGAR-MULE DEER UMFs  ####
  #'  Grazing Season
  coug_md_graze_DH <- list(cougar = DH_coug_graze1820, muledeer = DH_md_graze1820)
  coug_md_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_md_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_md_grazing_UMF)
  #'  Review covariates
  summary(coug_md_grazing_UMF)
  #'  Look at natural parameters f design matrix
  coug_md_grazing_UMF@fDesign
  
  #'  Hunting Season
  coug_md_hunt_DH <- list(cougar = DH_coug_hunt1820, muledeer = DH_md_hunt1820)
  coug_md_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_md_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_md_hunting_UMF)
  #'  Review covariates
  summary(coug_md_hunting_UMF)
  #'  Look at natural parameters f design matrix
  coug_md_hunting_UMF@fDesign
  
  ####  COUGAR-ELK UMFs  ####
  #'  Grazing Season
  coug_elk_graze_DH <- list(cougar = DH_coug_graze1820, elk = DH_elk_graze1820)
  coug_elk_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_elk_graze_DH,
                                                     siteCovs = stations_graze,
                                                     obsCovs = srvy_grazing_covs,
                                                     maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_elk_grazing_UMF)
  
  #'  Hunting Season
  coug_elk_hunt_DH <- list(cougar = DH_coug_hunt1820, elk = DH_elk_hunt1820)
  coug_elk_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_elk_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_elk_hunting_UMF)
  
  ####  COUGAR-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  coug_wtd_graze_DH <- list(cougar = DH_coug_graze1820, wtd = DH_wtd_graze1820)
  coug_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_wtd_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_wtd_grazing_UMF)
  
  #'  Hunting Season
  coug_wtd_hunt_DH <- list(cougar = DH_coug_hunt1820, wtd = DH_wtd_hunt1820)
  coug_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_wtd_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_wtd_hunting_UMF)
  
  ####  COUGAR-MOOSE UMFs  ####
  #'  Grazing Season
  coug_moose_graze_DH <- list(cougar = DH_coug_graze1820, moose = DH_moose_graze1820)
  coug_moose_grazing_UMF <- unmarkedFrameOccuMulti(y = coug_moose_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_moose_grazing_UMF)
  
  #'  Hunting Season
  coug_moose_hunt_DH <- list(cougar = DH_coug_hunt1820, moose = DH_moose_hunt1820)
  coug_moose_hunting_UMF <- unmarkedFrameOccuMulti(y = coug_moose_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coug_moose_hunting_UMF)
  
  
  ####  WOLF-MULE DEER UMFs  ####
  #'  Grazing Season
  wolf_md_graze_DH <- list(wolf = DH_wolf_graze1820, muledeer = DH_md_graze1820)
  wolf_md_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_md_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_md_grazing_UMF)
  
  #'  Hunting Season
  wolf_md_hunt_DH <- list(wolf = DH_wolf_hunt1820, muledeer = DH_md_hunt1820)
  wolf_md_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_md_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_md_hunting_UMF)
  
  ####  WOLF-ELK UMFs  ####
  #'  Grazing Season
  wolf_elk_graze_DH <- list(wolf = DH_wolf_graze1820, elk = DH_elk_graze1820)
  wolf_elk_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_elk_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_elk_grazing_UMF)
  
  #'  Hunting Season
  wolf_elk_hunt_DH <- list(wolf = DH_wolf_hunt1820, elk = DH_elk_hunt1820)
  wolf_elk_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_elk_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_elk_hunting_UMF)
  
  ####  WOLF-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  wolf_wtd_graze_DH <- list(wolf = DH_wolf_graze1820, wtd = DH_wtd_graze1820)
  wolf_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_wtd_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_wtd_grazing_UMF)
  
  #'  Hunting Season
  wolf_wtd_hunt_DH <- list(wolf = DH_wolf_hunt1820, wtd = DH_wtd_hunt1820)
  wolf_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_wtd_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_wtd_hunting_UMF)
  
  ####  WOLF-MOOSE UMFs  ####
  #'  Grazing Season
  wolf_moose_graze_DH <- list(wolf = DH_wolf_graze1820, moose = DH_moose_graze1820)
  wolf_moose_grazing_UMF <- unmarkedFrameOccuMulti(y = wolf_moose_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_moose_grazing_UMF)
  
  #'  Hunting Season
  wolf_moose_hunt_DH <- list(wolf = DH_wolf_hunt1820, moose = DH_moose_hunt1820)
  wolf_moose_hunting_UMF <- unmarkedFrameOccuMulti(y = wolf_moose_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(wolf_moose_hunting_UMF)
  
  
  
  ####  BLACK Bear-MULE DEER UMFs  ####
  #'  Grazing Season
  bear_md_graze_DH <- list(blackbear = DH_bear_graze1820, muledeer = DH_md_graze1820)
  bear_md_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_md_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_md_grazing_UMF)
  
  #'  Hunting Season
  bear_md_hunt_DH <- list(blackbear = DH_bear_hunt1820, muledeer = DH_md_hunt1820)
  bear_md_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_md_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_md_hunting_UMF)
  
  ####  BLACK BEAR-ELK UMFs  ####
  #'  Grazing Season
  bear_elk_graze_DH <- list(blackbear = DH_bear_graze1820, elk = DH_elk_graze1820)
  bear_elk_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_elk_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_elk_grazing_UMF)
  
  #'  Hunting Season
  bear_elk_hunt_DH <- list(blackbear = DH_bear_hunt1820, elk = DH_elk_hunt1820)
  bear_elk_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_elk_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_elk_hunting_UMF)
  
  ####  BLACK BEAR-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  bear_wtd_graze_DH <- list(blackbear = DH_bear_graze1820, wtd = DH_wtd_graze1820)
  bear_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_wtd_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_wtd_grazing_UMF)
  
  #'  Hunting Season
  bear_wtd_hunt_DH <- list(blackbear = DH_bear_hunt1820, wtd = DH_wtd_hunt1820)
  bear_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_wtd_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_wtd_hunting_UMF)
  
  ####  BLACK BEAR-MOOSE UMFs  ####
  #'  Grazing Season
  bear_moose_graze_DH <- list(blackbear = DH_bear_graze1820, moose = DH_moose_graze1820)
  bear_moose_grazing_UMF <- unmarkedFrameOccuMulti(y = bear_moose_graze_DH,
                                                   siteCovs = stations_graze,
                                                   obsCovs = srvy_grazing_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_moose_grazing_UMF)
  
  #'  Hunting Season
  bear_moose_hunt_DH <- list(blackbear = DH_bear_hunt1820, moose = DH_moose_hunt1820)
  bear_moose_hunting_UMF <- unmarkedFrameOccuMulti(y = bear_moose_hunt_DH,
                                                   siteCovs = stations_hunt,
                                                   obsCovs = srvy_hunting_covs,
                                                   maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bear_moose_hunting_UMF)
  
  

  ####  COYOTE-MULE DEER UMFs  ####
  #'  Grazing Season
  coy_md_graze_DH <- list(coyote = DH_coy_graze1820, mule_deer = DH_md_graze1820)
  coy_md_grazing_UMF <- unmarkedFrameOccuMulti(y = coy_md_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_md_grazing_UMF)
  
  #'  Hunting Season
  coy_md_hunt_DH <- list(coyote = DH_coy_hunt1820, mule_deer = DH_md_hunt1820)
  coy_md_hunting_UMF <- unmarkedFrameOccuMulti(y = coy_md_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_md_hunting_UMF)
  
  ####  COYOTE-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  coy_wtd_graze_DH <- list(coyote = DH_coy_graze1820, wtd = DH_wtd_graze1820)
  coy_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = coy_wtd_graze_DH,
                                                 siteCovs = stations_graze,
                                                 obsCovs = srvy_grazing_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_wtd_grazing_UMF)
  
  #'  Hunting Season
  coy_wtd_hunt_DH <- list(coyote = DH_coy_hunt1820, wtd = DH_wtd_hunt1820)
  coy_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = coy_wtd_hunt_DH,
                                                 siteCovs = stations_hunt,
                                                 obsCovs = srvy_hunting_covs,
                                                 maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(coy_wtd_hunting_UMF)
  
  
  
  ####  BOBCAT-MULE DEER UMFs  ####
  #'  Grazing Season
  bob_md_graze_DH <- list(bobcat = DH_bob_graze1820, mule_deer = DH_md_graze1820)
  bob_md_grazing_UMF <- unmarkedFrameOccuMulti(y = bob_md_graze_DH,
                                               siteCovs = stations_graze,
                                               obsCovs = srvy_grazing_covs,
                                               maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_md_grazing_UMF)
  
  #'  Hunting Season
  bob_md_hunt_DH <- list(bobcat = DH_bob_hunt1820, mule_deer = DH_md_hunt1820)
  bob_md_hunting_UMF <- unmarkedFrameOccuMulti(y = bob_md_hunt_DH,
                                               siteCovs = stations_hunt,
                                               obsCovs = srvy_hunting_covs,
                                               maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_md_hunting_UMF)
  
  ####  BOBCAT-WHITE-TAILED DEER UMFs  ####
  #'  Grazing Season
  bob_wtd_graze_DH <- list(bobcat = DH_bob_graze1820, wtd = DH_wtd_graze1820)
  bob_wtd_grazing_UMF <- unmarkedFrameOccuMulti(y = bob_wtd_graze_DH,
                                                siteCovs = stations_graze,
                                                obsCovs = srvy_grazing_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_wtd_grazing_UMF)
  
  #'  Hunting Season
  bob_wtd_hunt_DH <- list(bobcat = DH_bob_hunt1820, wtd = DH_wtd_hunt1820)
  bob_wtd_hunting_UMF <- unmarkedFrameOccuMulti(y = bob_wtd_hunt_DH,
                                                siteCovs = stations_hunt,
                                                obsCovs = srvy_hunting_covs,
                                                maxOrder = 2)
  #'  Visualize detection/non-detection data
  plot(bob_wtd_hunting_UMF)
  
  
  
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
  detFormulas_trail <- c("~Trail", "~Trail")
  detFormulas_graze <- c("~Trail + WeeklyGrazing", "~Trail + WeeklyGrazing")
  
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null <- c("~1", "~1", "~0")
  occFormulas_null1 <- c("~1", "~1", "~1")
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
  (cougmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  (cougmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougmd_fld <- fitList(cougmd_trail, cougmd_dgraze)
  #' Model selection
  modSel(cougmd_fld)
  
  (cougmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  (cougmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_md_grazing_UMF, silent = TRUE))
  (cougmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_md_grazing_UMF, silent = TRUE))
  (cougmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_md_grazing_UMF, silent = TRUE))
  (cougmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, coug_md_grazing_UMF, silent = TRUE))
  (cougmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, coug_md_grazing_UMF, silent = TRUE))
  (cougmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, coug_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougmd_fl <- fitList(cougmd_null0, cougmd_null1, cougmd_hab0, cougmd_hab1, cougmd_graze0, cougmd_graze1, cougmd_graze2)
  #' Model selection
  modSel(cougmd_fl)
  summary(cougmd_hab0)
  
  ####  Cougar-ELK Grazing Season  ####
  (cougelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  (cougelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougelk_fld <- fitList(cougelk_trail, cougelk_dgraze)
  #' Model selection
  modSel(cougelk_fld)
  
  (cougelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  (cougelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_elk_grazing_UMF, silent = TRUE))
  (cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_grazing_UMF, silent = TRUE))
  (cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_grazing_UMF, silent = TRUE))
  (cougelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, coug_elk_grazing_UMF, silent = TRUE))
  (cougelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, coug_elk_grazing_UMF, silent = TRUE)) # fails
  (cougelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, coug_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougelk_fl <- fitList(cougelk_null0, cougelk_null1, cougelk_hab0, cougelk_hab1, cougelk_graze0, cougelk_graze2)
  #' Model selection
  modSel(cougelk_fl)
  summary(cougelk_hab0)
  
  ####  Cougar-White-tailed Deer Grazing Season  ####
  (cougwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  (cougwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougwtd_fld <- fitList(cougwtd_trail, cougwtd_dgraze)
  #' Model selection
  modSel(cougwtd_fld)
  
  (cougwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  (cougwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_wtd_grazing_UMF, silent = TRUE))
  (cougwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_wtd_grazing_UMF, silent = TRUE))
  (cougwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_wtd_grazing_UMF, silent = TRUE))
  (cougwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_wtd_grazing_UMF, silent = TRUE))
  (cougwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_wtd_grazing_UMF, silent = TRUE))
  (cougwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougwtd_fl <- fitList(cougwtd_null0, cougwtd_null1, cougwtd_hab0, cougwtd_hab1, cougwtd_graze0, cougwtd_graze1, cougwtd_graze2)
  #' Model selection
  modSel(cougwtd_fl)
  summary(cougwtd_graze2)
  
  ####  Cougar-Moose Grazing Season  ####
  (cougmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  (cougmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougmoose_fld <- fitList(cougmoose_trail, cougmoose_dgraze)
  #' Model selection
  modSel(cougmoose_fld)
  
  (cougmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  (cougmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_moose_grazing_UMF, silent = TRUE))
  (cougmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_moose_grazing_UMF, silent = TRUE))
  (cougmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_moose_grazing_UMF, silent = TRUE))
  (cougmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_moose_grazing_UMF, silent = TRUE))
  (cougmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_moose_grazing_UMF, silent = TRUE))
  (cougmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  cougmoose_fl <- fitList(cougmoose_null0, cougmoose_null1, cougmoose_hab0, cougmoose_hab1, cougmoose_graze0, cougmoose_graze1, cougmoose_graze2)
  #' Model selection
  modSel(cougmoose_fl)
  summary(cougmoose_hab0)
  summary(cougmoose_graze0)
  summary(cougmoose_hab1)
  
  
  
  ####  Wolf-Mule Deer Grazing Season  ####
  (wolfmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  (wolfmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolfmd_fld <- fitList(wolfmd_trail, wolfmd_dgraze)
  #' Model selection
  modSel(wolfmd_fld)
  
  (wolfmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  (wolfmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_md_grazing_UMF, silent = TRUE))
  (wolfmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_md_grazing_UMF, silent = TRUE))
  (wolfmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_md_grazing_UMF, silent = TRUE))
  (wolfmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, wolf_md_grazing_UMF, silent = TRUE))
  (wolfmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, wolf_md_grazing_UMF, silent = TRUE))
  (wolfmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, wolf_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolfmd_fl <- fitList(wolfmd_null0, wolfmd_null1, wolfmd_hab0, wolfmd_hab1, wolfmd_graze0, wolfmd_graze1, wolfmd_graze2)
  #' Model selection
  modSel(wolfmd_fl)
  summary(wolfmd_hab0)
  
  ####  Wolf-Elk Grazing Season  ####
  (wolfelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  (wolfelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolfelk_fld <- fitList(wolfelk_trail, wolfelk_dgraze)
  #' Model selection
  modSel(wolfelk_fld)
  
  (wolfelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  (wolfelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_elk_grazing_UMF, silent = TRUE))
  (wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_grazing_UMF, silent = TRUE))
  (wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_grazing_UMF, silent = TRUE))
  (wolfelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, wolf_elk_grazing_UMF, silent = TRUE))
  (wolfelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, wolf_elk_grazing_UMF, silent = TRUE)) 
  (wolfelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, wolf_elk_grazing_UMF, silent = TRUE)) # FAILS
  #' List of fitted models
  wolfelk_fl <- fitList(wolfelk_null0, wolfelk_null1, wolfelk_hab0, wolfelk_hab1, wolfelk_graze0, wolfelk_graze1)
  #' Model selection
  modSel(wolfelk_fl)
  summary(wolfelk_hab0)
  
  ####  Wolf-White-tailed Deer Grazing Season  ####
  (wolfwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  (wolfwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolfwtd_fld <- fitList(wolfwtd_trail, wolfwtd_dgraze)
  #' Model selection
  modSel(wolfwtd_fld)
  
  (wolfwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  (wolfwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_wtd_grazing_UMF, silent = TRUE))
  (wolfwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_wtd_grazing_UMF, silent = TRUE))
  (wolfwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_wtd_grazing_UMF, silent = TRUE))
  (wolfwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_wtd_grazing_UMF, silent = TRUE))
  (wolfwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_wtd_grazing_UMF, silent = TRUE))
  (wolfwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolfwtd_fl <- fitList(wolfwtd_null0, wolfwtd_null1, wolfwtd_hab0, wolfwtd_hab1, wolfwtd_graze0, wolfwtd_graze1, wolfwtd_graze2)
  #' Model selection
  modSel(wolfwtd_fl)
  summary(wolfwtd_graze0)
  
  ####  Wolf-Moose Grazing Season  ####
  (wolfmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  (wolfmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolfmoose_fld <- fitList(wolfmoose_trail, wolfmoose_dgraze)
  #' Model selection
  modSel(wolfmoose_fld)
  
  (wolfmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  (wolfmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_moose_grazing_UMF, silent = TRUE))
  (wolfmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_moose_grazing_UMF, silent = TRUE))
  (wolfmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_moose_grazing_UMF, silent = TRUE))
  (wolfmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_moose_grazing_UMF, silent = TRUE))
  (wolfmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_moose_grazing_UMF, silent = TRUE))
  (wolfmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  wolfmoose_fl <- fitList(wolfmoose_null0, wolfmoose_null1, wolfmoose_hab0, wolfmoose_hab1, wolfmoose_graze0, wolfmoose_graze1, wolfmoose_graze2)
  #' Model selection
  modSel(wolfmoose_fl)
  summary(wolfmoose_graze0)
  
  
  
  ####  Bear-Mule Deer Grazing Season  ####
  (bearmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  (bearmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bearmd_fld <- fitList(bearmd_trail, bearmd_dgraze)
  #' Model selection
  modSel(bearmd_fld)
  
  (bearmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  (bearmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_md_grazing_UMF, silent = TRUE))
  (bearmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_md_grazing_UMF, silent = TRUE))
  (bearmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_md_grazing_UMF, silent = TRUE))
  (bearmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bear_md_grazing_UMF, silent = TRUE))
  (bearmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bear_md_grazing_UMF, silent = TRUE))
  (bearmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bear_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bearmd_fl <- fitList(bearmd_null0, bearmd_null1, bearmd_hab0, bearmd_hab1, bearmd_graze0, bearmd_graze1, bearmd_graze2)
  #' Model selection
  modSel(bearmd_fl)
  summary(bearmd_graze1) # this makes no sense- graze has 0 effect on either species and co-occ appears independent
  
  ####  Bear-Elk Grazing Season  ####
  (bearelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  (bearelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bearelk_fld <- fitList(bearelk_trail, bearelk_dgraze)
  #' Model selection
  modSel(bearelk_fld)
  
  (bearelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  (bearelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_elk_grazing_UMF, silent = TRUE))
  (bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_grazing_UMF, silent = TRUE))
  (bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_grazing_UMF, silent = TRUE))
  (bearelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bear_elk_grazing_UMF, silent = TRUE))
  (bearelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bear_elk_grazing_UMF, silent = TRUE)) 
  (bearelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bear_elk_grazing_UMF, silent = TRUE)) # FAILS
  #' List of fitted models
  bearelk_fl <- fitList(bearelk_null0, bearelk_null1, bearelk_hab0, bearelk_hab1, bearelk_graze0, bearelk_graze1)
  #' Model selection
  modSel(bearelk_fl)
  summary(bearelk_graze1) # but doesn't make sense b/c grazing has 0 effect sooo....
  
  ####  Bear-White-tailed Deer Grazing Season  ####
  (bearwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  (bearwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bearwtd_fld <- fitList(bearwtd_trail, bearwtd_dgraze)
  #' Model selection
  modSel(bearwtd_fld)
  
  (bearwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  (bearwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_wtd_grazing_UMF, silent = TRUE))
  (bearwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_wtd_grazing_UMF, silent = TRUE))
  (bearwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_wtd_grazing_UMF, silent = TRUE))
  (bearwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_wtd_grazing_UMF, silent = TRUE))
  (bearwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_wtd_grazing_UMF, silent = TRUE))
  (bearwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bearwtd_fl <- fitList(bearwtd_null0, bearwtd_null1, bearwtd_hab0, bearwtd_hab1, bearwtd_graze0, bearwtd_graze1, bearwtd_graze2)
  #' Model selection
  modSel(bearwtd_fl)
  summary(bearwtd_graze0)
  
  ####  Bear-Moose Grazing Season  ####
  (bearmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  (bearmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bearmoose_fld <- fitList(bearmoose_trail, bearmoose_dgraze)
  #' Model selection
  modSel(bearmoose_fld)
  
  (bearmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  (bearmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_moose_grazing_UMF, silent = TRUE))
  (bearmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_moose_grazing_UMF, silent = TRUE))
  (bearmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_moose_grazing_UMF, silent = TRUE))
  (bearmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_moose_grazing_UMF, silent = TRUE))
  (bearmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_moose_grazing_UMF, silent = TRUE))
  (bearmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bearmoose_fl <- fitList(bearmoose_null0, bearmoose_null1, bearmoose_hab0, bearmoose_hab1, bearmoose_graze0, bearmoose_graze1, bearmoose_graze2)
  #' Model selection
  modSel(bearmoose_fl)
  summary(bearmoose_graze0)
  
  
  
  ####  Coyote-Mule Deer Grazing Season  ####
  (coymd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  (coymd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coymd_fld <- fitList(coymd_trail, coymd_dgraze)
  #' Model selection
  modSel(coymd_fld)
  
  (coymd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  (coymd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_md_grazing_UMF, silent = TRUE))
  (coymd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_md_grazing_UMF, silent = TRUE))
  (coymd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_md_grazing_UMF, silent = TRUE))
  (coymd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_md_grazing_UMF, silent = TRUE))
  (coymd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_md_grazing_UMF, silent = TRUE))
  (coymd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coymd_fl <- fitList(coymd_null0, coymd_null1, coymd_hab0, coymd_hab1, coymd_graze0, coymd_graze1, coymd_graze2)
  #' Model selection
  modSel(coymd_fl)
  summary(coymd_graze2)
  
  ####  Coyote-White-tailed Deer Grazing Season  ####
  (coywtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  (coywtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coywtd_fld <- fitList(coywtd_trail, coywtd_dgraze)
  #' Model selection
  modSel(coywtd_fld)
  
  (coywtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  (coywtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_wtd_grazing_UMF, silent = TRUE))
  (coywtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_wtd_grazing_UMF, silent = TRUE))
  (coywtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_wtd_grazing_UMF, silent = TRUE))
  (coywtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_wtd_grazing_UMF, silent = TRUE))
  (coywtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_wtd_grazing_UMF, silent = TRUE))
  (coywtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  coywtd_fl <- fitList(coywtd_null0, coywtd_null1, coywtd_hab0, coywtd_hab1, coywtd_graze0, coywtd_graze1, coywtd_graze2)
  #' Model selection
  modSel(coywtd_fl)
  summary(coywtd_graze2)
  
  
  
  ####  Bobcat-Mule Deer Grazing Season  ####
  (bobmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  (bobmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bobmd_fld <- fitList(bobmd_trail, bobmd_dgraze)
  #' Model selection
  modSel(bobmd_fld)
  
  (bobmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  (bobmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_md_grazing_UMF, silent = TRUE))
  (bobmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_md_grazing_UMF, silent = TRUE))
  (bobmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_md_grazing_UMF, silent = TRUE))
  (bobmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bob_md_grazing_UMF, silent = TRUE))
  (bobmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bob_md_grazing_UMF, silent = TRUE))
  (bobmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bob_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bobmd_fl <- fitList(bobmd_null0, bobmd_null1, bobmd_hab0, bobmd_hab1, bobmd_graze0, bobmd_graze1, bobmd_graze2)
  #' Model selection
  modSel(bobmd_fl)
  summary(bobmd_graze0)
  
  ####  Bobcat-White-tailed Deer Grazing Season  ####
  (bobwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  (bobwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bobwtd_fld <- fitList(bobwtd_trail, bobwtd_dgraze)
  #' Model selection
  modSel(bobwtd_fld)
  
  (bobwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  (bobwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bob_wtd_grazing_UMF, silent = TRUE))
  (bobwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bob_wtd_grazing_UMF, silent = TRUE))
  (bobwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bob_wtd_grazing_UMF, silent = TRUE))
  (bobwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bob_wtd_grazing_UMF, silent = TRUE))
  (bobwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bob_wtd_grazing_UMF, silent = TRUE))
  (bobwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bob_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  bobwtd_fl <- fitList(bobwtd_null0, bobwtd_null1, bobwtd_hab0, bobwtd_hab1, bobwtd_graze0, bobwtd_graze1, bobwtd_graze2)
  #' Model selection
  modSel(bobwtd_fl)
  summary(bobwtd_graze2)

  
  
  ####  HUNTING SEASON MODELS  ####
  #'  Detection and occupancy formulas
  #'  Question 1: Does hunting activity affect detection?
  detFormulas_trail <- c("~Trail", "~Trail")
  detFormulas_hunt <- c("~Trail + WeeklyHunting", "~Trail + WeeklyHunting")
  # detFormulas_veh <- c("~Trail + WeeklyVehicles", "~Trail + WeeklyVehicles")
  detFormulas_pub <- c("~Trail + Public", "~Trail + Public")
  detFormulas_pubhunt <- c("~Trail + Public + WeeklyHunting", "~Trail + Public + WeeklyHunting")
  #'  Remove Public vs Private covariate from predator sub-model
  detFormulas_pubish <- c("~Trail", "~Trail + Public")
  detFormulas_pubhuntish <- c("~Trail + WeeklyHunting", "~Trail + Public + WeeklyHunting")
  
  #'  Question 1: Is co-occurrence dependent?
  occFormulas_null <- c("~1", "~1", "~0")
  occFormulas_null1 <- c("~1", "~1", "~1")
  occFormulas_hab0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~0")
  occFormulas_hab1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest",
                        "~1")
  #'  Question 2: Does hunting activity affect occurrence and/or co-occurrence patterns
  occFormulas_hunt0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~0")
  occFormulas_hunt1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                          "~1")
  occFormulas_hunt2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                         "~HuntingActivity")
  occFormulas_pub0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                         "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                         "~0")
  occFormulas_pub1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                         "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                         "~1")
  occFormulas_pub2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Study_Area + Elev + I(Elev^2) + PercForest + Public", 
                        "~Public")
  # occFormulas_huntpub0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                              "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                              "~0")
  # occFormulas_huntpub1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                              "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                              "~1")
  # occFormulas_huntpub2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                              "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                              "~HuntingActivity + Public")
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
  occFormulas_huntpubish0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
                           "~0")
  occFormulas_huntpubish1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
                           "~1")
  occFormulas_huntpubish2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
                           "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
                           "~HuntingActivity + Public")
  
  
  ####  Cougar-Mule Deer Hunting Season  ####
  #'  Note: using occupancy model formulas that drop Public covariate from predator sub-model
  (cougmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  cougmd_fld <- fitList(cougmd_trail, cougmd_hunt, cougmd_pub, cougmd_pubhunt)
  #' Model selection
  modSel(cougmd_fld)
  
  (cougmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_md_hunting_UMF, silent = TRUE)) 
  (cougmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_md_hunting_UMF, silent = TRUE)) 
  (cougmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_md_hunting_UMF, silent = TRUE)) # FAIL
  (cougmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_md_hunting_UMF, silent = TRUE)) 
  (cougmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_md_hunting_UMF, silent = TRUE))
  (cougmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  cougmd_fl <- fitList(cougmd_null0, cougmd_null1, cougmd_hab0, cougmd_hab1, cougmd_hunt0, cougmd_hunt1, cougmd_pub0, cougmd_pub1, cougmd_pub2) 
  #' Model selection
  modSel(cougmd_fl)
  summary(cougmd_pub2)

  ####  Cougar-Elk Hunting Season  ####
  (cougelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  cougelk_fld <- fitList(cougelk_trail, cougelk_hunt, cougelk_pub, cougelk_pubhunt)
  #' Model selection
  modSel(cougelk_fld)
  
  (cougelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_hunting_UMF, silent = TRUE)) 
  (cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_hunting_UMF, silent = TRUE)) 
  (cougelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (cougelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pubish0, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pubish1, coug_elk_hunting_UMF, silent = TRUE))
  (cougelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pubish2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  cougelk_fl <- fitList(cougelk_null0, cougelk_null1, cougelk_hab0, cougelk_hab1, cougelk_hunt0, cougelk_hunt1, cougelk_pub0, cougelk_pub1) 
  #' Model selection
  modSel(cougelk_fl)
  summary(cougelk_pub1)
  
  ####  Cougar-White-tailed Deer Hunting Season  ####
  (cougwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  cougwtd_fld <- fitList(cougwtd_trail, cougwtd_hunt, cougwtd_pub, cougwtd_pubhunt)
  #' Model selection
  modSel(cougwtd_fld)
  
  (cougwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_wtd_hunting_UMF, silent = TRUE)) 
  (cougwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_wtd_hunting_UMF, silent = TRUE)) 
  (cougwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_wtd_hunting_UMF, silent = TRUE))
  (cougwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_wtd_hunting_UMF, silent = TRUE)) 
  (cougwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  cougwtd_fl <- fitList(cougwtd_null0, cougwtd_null1, cougwtd_hab0, cougwtd_hab1, cougwtd_hunt0, cougwtd_hunt1, cougwtd_hunt2, cougwtd_pub0, cougwtd_pub1, cougwtd_pub2) 
  #' Model selection
  modSel(cougwtd_fl)
  summary(cougwtd_pub2)
  
  ####  Cougar-Moose Hunting Season  ####
  (cougmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  cougmoose_fld <- fitList(cougmoose_trail, cougmoose_hunt, cougmoose_pub, cougmoose_pubhunt)
  #' Model selection
  modSel(cougmoose_fld)
  
  (cougmoose_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_moose_hunting_UMF, silent = TRUE)) 
  (cougmoose_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (cougmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (cougmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (cougmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_moose_hunting_UMF, silent = TRUE))
  (cougmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  cougmoose_fl <- fitList(cougmoose_null0, cougmoose_null1, cougmoose_hab0, cougmoose_hunt0, cougmoose_pub0, cougmoose_pub1, cougmoose_pub2) 
  #' Model selection
  modSel(cougmoose_fl)
  summary(cougmoose_hunt0)
  
  
  ####  Wolf-Mule Deer Hunting Season  ####
  #'  Note: using detection & occupancy model formulas that drop Public covariate from predator sub-models
  (wolfmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolfmd_fld <- fitList(wolfmd_trail, wolfmd_hunt, wolfmd_pub, wolfmd_pubhunt) 
  #' Model selection
  modSel(wolfmd_fld)
  
  (wolfmd_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_md_hunting_UMF, silent = TRUE)) 
  (wolfmd_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_md_hunting_UMF, silent = TRUE)) 
  (wolfmd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_md_hunting_UMF, silent = TRUE))
  (wolfmd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_md_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  wolfmd_fl <- fitList(wolfmd_null0, wolfmd_null1, wolfmd_hab0, wolfmd_hab1, wolfmd_hunt0, wolfmd_hunt1, wolfmd_hunt2, wolfmd_pub0, wolfmd_pub1) 
  #' Model selection
  modSel(wolfmd_fl)
  summary(wolfmd_hunt1)
  
  
  ####  Wolf-Elk Hunting Season  ####
  (wolfelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolfelk_fld <- fitList(wolfelk_trail, wolfelk_hunt, wolfelk_pub, wolfelk_pubhunt) 
  #' Model selection
  modSel(wolfelk_fld)
  
  (wolfelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_hunting_UMF, silent = TRUE)) 
  (wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_hunting_UMF, silent = TRUE)) 
  (wolfelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pubish0, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pubish1, wolf_elk_hunting_UMF, silent = TRUE))
  (wolfelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pubish2, wolf_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  wolfelk_fl <- fitList(wolfelk_null0, wolfelk_null1, wolfelk_hab0, wolfelk_hab1, wolfelk_hunt0, wolfelk_hunt1, wolfelk_hunt2, wolfelk_pub0, wolfelk_pub1) 
  #' Model selection
  modSel(wolfelk_fl)
  summary(wolfelk_hunt0)
  
  ####  Wolf-White-tailed Deer Hunting Season  ####
  (wolfwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolfwtd_fld <- fitList(wolfwtd_trail, wolfwtd_hunt, wolfwtd_pub, wolfwtd_pubhunt)
  #' Model selection
  modSel(wolfwtd_fld)
  
  (wolfwtd_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_wtd_hunting_UMF, silent = TRUE)) 
  (wolfwtd_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_wtd_hunting_UMF, silent = TRUE)) 
  (wolfwtd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  (wolfwtd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_wtd_hunting_UMF, silent = TRUE))  
  (wolfwtd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_wtd_hunting_UMF, silent = TRUE))
  (wolfwtd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  wolfwtd_fl <- fitList(wolfwtd_null0, wolfwtd_null1, wolfwtd_hab0, wolfwtd_hab1, wolfwtd_hunt0, wolfwtd_hunt2, wolfwtd_pub0, wolfwtd_pub1) 
  #' Model selection
  modSel(wolfwtd_fl)
  summary(wolfwtd_hunt0)
  
  ####  Wolf-Moose Hunting Season  ####
  (wolfmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  wolfmoose_fld <- fitList(wolfmoose_trail, wolfmoose_hunt, wolfmoose_pub, wolfmoose_pubhunt)
  #' Model selection
  modSel(wolfmoose_fld)
  
  (wolfmoose_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_moose_hunting_UMF, silent = TRUE)) 
  (wolfmoose_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_moose_hunting_UMF, silent = TRUE)) 
  (wolfmoose_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_moose_hunting_UMF, silent = TRUE))
  (wolfmoose_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_moose_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  wolfmoose_fl <- fitList(wolfmoose_null0, wolfmoose_null1, wolfmoose_hab0, wolfmoose_hab1, wolfmoose_hunt0, wolfmoose_hunt1, wolfmoose_hunt2, wolfmoose_pub0, wolfmoose_pub1)
  #' Model selection
  modSel(wolfmoose_fl)
  summary(wolfmoose_hunt2)
  
  
  ####  Bear-Mule Deer Hunting Season  ####
  (bearmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bearmd_fld <- fitList(bearmd_trail, bearmd_hunt, bearmd_pub, bearmd_pubhunt)
  #' Model selection
  modSel(bearmd_fld)
  
  (bearmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_md_hunting_UMF, silent = TRUE)) 
  (bearmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_md_hunting_UMF, silent = TRUE)) 
  (bearmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_md_hunting_UMF, silent = TRUE)) 
  (bearmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_md_hunting_UMF, silent = TRUE))
  (bearmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bearmd_fl <- fitList(bearmd_null0, bearmd_null1, bearmd_hab0, bearmd_hab1, bearmd_hunt0, bearmd_hunt1, bearmd_hunt2, bearmd_pub0, bearmd_pub1, bearmd_pub2) 
  #' Model selection
  modSel(bearmd_fl)
  summary(bearmd_hab0)
  
  ####  Bear-Elk Hunting Season  ####
  (bearelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bearelk_fld <- fitList(bearelk_trail, bearelk_hunt, bearelk_pub, bearelk_pubhunt)
  #' Model selection
  modSel(bearelk_fld)
  
  (bearelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_hunting_UMF, silent = TRUE)) 
  (bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_hunting_UMF, silent = TRUE)) 
  (bearelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pub0, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pub1, bear_elk_hunting_UMF, silent = TRUE))
  (bearelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pub2, bear_elk_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  bearelk_fl <- fitList(bearelk_null0, bearelk_null1, bearelk_hab0, bearelk_hab1, bearelk_hunt0, bearelk_hunt1, bearelk_hunt2, bearelk_pub0, bearelk_pub1)
  #' Model selection
  modSel(bearelk_fl)
  summary(bearelk_pub1) # don't buy it though b/c public not imprtant to either spp
  summary(bearelk_hab1)
  
  ####  Bear-White-tailed Deer Hunting Season  ####
  (bearwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bearwtd_fld <- fitList(bearwtd_trail, bearwtd_hunt, bearwtd_pub, bearwtd_pubhunt)
  #' Model selection
  modSel(bearwtd_fld)
  
  (bearwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_wtd_hunting_UMF, silent = TRUE)) 
  (bearwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_wtd_hunting_UMF, silent = TRUE)) 
  (bearwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_wtd_hunting_UMF, silent = TRUE)) 
  (bearwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_wtd_hunting_UMF, silent = TRUE)) 
  (bearwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_wtd_hunting_UMF, silent = TRUE))
  (bearwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bearwtd_fl <- fitList(bearwtd_null0, bearwtd_null1, bearwtd_hab0, bearwtd_hab1, bearwtd_hunt0, bearwtd_hunt1, bearwtd_hunt2, bearwtd_pub0, bearwtd_pub1, bearwtd_pub2) 
  #' Model selection
  modSel(bearwtd_fl)
  summary(bearwtd_pub0)
  
  ####  Bear-Moose Hunting Season  ####
  (bearmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bearmoose_fld <- fitList(bearmoose_trail, bearmoose_hunt, bearmoose_pub, bearmoose_pubhunt)
  #' Model selection
  modSel(bearmoose_fld)

  (bearmoose_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_moose_hunting_UMF, silent = TRUE)) 
  (bearmoose_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_moose_hunting_UMF, silent = TRUE)) 
  (bearmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  (bearmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_moose_hunting_UMF, silent = TRUE)) 
  (bearmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_moose_hunting_UMF, silent = TRUE))
  (bearmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  bearmoose_fl <- fitList(bearmoose_null0, bearmoose_null1, bearmoose_hab0, bearmoose_hab1, bearmoose_hunt0, bearmoose_hunt2, bearmoose_pub0, bearmoose_pub1) 
  #' Model selection
  modSel(bearmoose_fl)
  summary(bearmoose_hunt0)
  
  
  ####  Coyote-Mule Deer Hunting Season  ####
  (coymd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_hunting_UMF, silent = TRUE))
  (coymd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_md_hunting_UMF, silent = TRUE))
  (coymd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  (coymd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  coymd_fld <- fitList(coymd_trail, coymd_hunt, coymd_pub, coymd_pubhunt)
  #' Model selection
  modSel(coymd_fld)
  
  (coymd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  (coymd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_md_hunting_UMF, silent = TRUE))
  (coymd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coy_md_hunting_UMF, silent = TRUE)) 
  (coymd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coy_md_hunting_UMF, silent = TRUE))
  (coymd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_md_hunting_UMF, silent = TRUE)) 
  (coymd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_md_hunting_UMF, silent = TRUE))
  (coymd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_md_hunting_UMF, silent = TRUE)) # FAIL
  (coymd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_md_hunting_UMF, silent = TRUE)) 
  (coymd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_md_hunting_UMF, silent = TRUE))
  (coymd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_md_hunting_UMF, silent = TRUE))
  (coymd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, coy_md_hunting_UMF, silent = TRUE)) 
  (coymd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coy_md_hunting_UMF, silent = TRUE))
  (coymd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coy_md_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  coymd_fl <- fitList(coymd_null0, coymd_null1, coymd_hab0, coymd_hab1, coymd_hunt0, coymd_hunt1, coymd_pub0, coymd_pub1, coymd_pub2, coymd_pubhunt0, coymd_pubhunt1)
  #' Model selection
  modSel(coymd_fl)
  summary(coymd_hunt0)
  
  ####  Coyote-White-tailed Deer Hunting Season  ####
  (coywtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  (coywtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  coywtd_fld <- fitList(coywtd_trail, coywtd_hunt, coywtd_pub, coywtd_pubhunt)
  #' Model selection
  modSel(coywtd_fld)
  
  (coywtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  (coywtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (coywtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (coywtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_wtd_hunting_UMF, silent = TRUE)) # FAIL
  (coywtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (coywtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (coywtd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coy_wtd_hunting_UMF, silent = TRUE))
  (coywtd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coy_wtd_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  coywtd_fl <- fitList(coywtd_null0, coywtd_null1, coywtd_hab0, coywtd_hab1, coywtd_hunt0, coywtd_hunt1, coywtd_pub0, coywtd_pub1, coywtd_pub2, coywtd_pubhunt0, coywtd_pubhunt1)
  #' Model selection
  modSel(coywtd_fl)
  summary(coywtd_hunt1)
  
  
  ####  Bobcat-Mule Deer Hunting Season  ####
  (bobmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_hunting_UMF, silent = TRUE))
  (bobmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_md_hunting_UMF, silent = TRUE))
  (bobmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  (bobmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  bobmd_fld <- fitList(bobmd_trail, bobmd_hunt, bobmd_pub, bobmd_pubhunt)
  #' Model selection
  modSel(bobmd_fld)
  
  (bobmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  (bobmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_md_hunting_UMF, silent = TRUE))
  (bobmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bob_md_hunting_UMF, silent = TRUE)) 
  (bobmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bob_md_hunting_UMF, silent = TRUE))
  (bobmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_md_hunting_UMF, silent = TRUE)) 
  (bobmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_md_hunting_UMF, silent = TRUE)) #FAIL?
  (bobmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_md_hunting_UMF, silent = TRUE)) 
  (bobmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_md_hunting_UMF, silent = TRUE)) 
  (bobmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_md_hunting_UMF, silent = TRUE))
  (bobmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_md_hunting_UMF, silent = TRUE))
  (bobmd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, bob_md_hunting_UMF, silent = TRUE)) 
  (bobmd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, bob_md_hunting_UMF, silent = TRUE))
  (bobmd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, bob_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  bobmd_fl <- fitList(bobmd_null0, bobmd_null1, bobmd_hab0, bobmd_hab1, bobmd_hunt0, bobmd_hunt1, bobmd_hunt2, bobmd_pub0, bobmd_pub1, bobmd_pub2, bobmd_pubhunt0, bobmd_pubhunt1, bobmd_pubhunt2)
  #' Model selection
  modSel(bobmd_fl)
  summary(bobmd_hunt0) # ok but the hunting activity effect on bobcat sub-model has unusually large effect size...
  
  
  ####  Bobcat-White-tailed Deer Hunting Season  ####
  (bobwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE))
  (bobwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE))
  (bobwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  (bobwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  bobwtd_fld <- fitList(bobwtd_trail, bobwtd_hunt, bobwtd_pub, bobwtd_pubhunt)
  #' Model selection
  modSel(bobwtd_fld)
  
  (bobwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  (bobwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_wtd_hunting_UMF, silent = TRUE))
  (bobwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (bobwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bob_wtd_hunting_UMF, silent = TRUE))
  (bobwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (bobwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_wtd_hunting_UMF, silent = TRUE))
  (bobwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_wtd_hunting_UMF, silent = TRUE)) #FAIL?
  (bobwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (bobwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_wtd_hunting_UMF, silent = TRUE))
  (bobwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_wtd_hunting_UMF, silent = TRUE))
  (bobwtd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (bobwtd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, bob_wtd_hunting_UMF, silent = TRUE)) #FAIL?
  (bobwtd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, bob_wtd_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  bobwtd_fl <- fitList(bobwtd_null0, bobwtd_null1, bobwtd_hab0, bobwtd_hab1, bobwtd_hunt0, bobwtd_hunt1, bobwtd_pub0, bobwtd_pub1, bobwtd_pub2, bobwtd_pubhunt0)
  #' Model selection
  modSel(bobwtd_fl)
  summary(bobwtd_pub1) #even though interaction appears independent? 
  summary(bobwtd_pub0)
  
  
  
  
  
  
  
  
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info

  #'  Function to save occupancy results
  occ_out <- function(mod, spp, season, model) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
        # Model = rep(model, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
      # relocate(Model, .before = Species)
    return(out)
  }
  
  #'  Run each model through function
  #'  Full models
  bob_s1819_occ <- occ_out(bob_s1819_global, "Bobcat", "Summer") #, "Global"
  bob_w1820_occ <- occ_out(bob_w1820_global, "Bobcat", "Winter")
  coug_s1819_occ <- occ_out(coug_s1819_global, "Cougar", "Summer")
  coug_w1820_occ <- occ_out(coug_w1820_global, "Cougar", "Winter")
  coy_s1819_occ <- occ_out(coy_s1819_global, "Coyote", "Summer")
  coy_w1820_occ <- occ_out(coy_w1820_global, "Coyote", "Winter")
  wolf_s1819_occ <- occ_out(wolf_s1819_global2, "Wolf", "Summer")
  wolf_w1820_occ <- occ_out(wolf_w1820_global2, "Wolf", "Winter")
  elk_s1819_occ <- occ_out(elk_s1819_global2, "Elk", "Summer")
  elk_w1820_occ <- occ_out(elk_w1820_global2, "Elk", "Winter")
  md_s1819_occ <- occ_out(md_s1819_global, "Mule Deer", "Summer")
  md_w1820_occ <- occ_out(md_w1820_global, "Mule Deer", "Winter")
  wtd_s1819_occ <- occ_out(wtd_s1819_global2, "White-tailed Deer", "Summer")
  wtd_w1820_occ <- occ_out(wtd_w1820_global2, "White-tailed Deer", "Winter")
  #'  Top models identified by dredging
  # bob_s1819_occt <- occ_out(bob_s1819_top, "Bobcat", "Summer", "Top")
  # bob_w1820_occt <- occ_out(bob_w1820_top, "Bobcat", "Winter", "Top")
  # coug_s1819_occt <- occ_out(coug_s1819_top, "Cougar", "Summer", "Top")
  # coug_w1820_occt <- occ_out(coug_w1820_top, "Cougar", "Winter", "Top")
  # coy_s1819_occt <- occ_out(coy_s1819_top, "Coyote", "Summer", "Top")
  # coy_w1820_occt <- occ_out(coy_w1820_top, "Coyote", "Winter", "Top")
  # wolf_s1819_occt <- occ_out(wolf_s1819_top, "Wolf", "Summer", "Top")
  # wolf_w1820_occt <- occ_out(wolf_w1820_top, "Wolf", "Winter", "Top")
  # elk_s1819_occt <- occ_out(elk_s1819_top, "Elk", "Summer", "Top")
  # elk_w1820_occt <- occ_out(elk_w1820_top, "Elk", "Winter", "Top")
  # md_s1819_occt <- occ_out(md_s1819_top, "Mule Deer", "Summer", "Top")
  # md_w1820_occt <- occ_out(md_w1820_top, "Mule Deer", "Winter", "Top")
  # wtd_s1819_occt <- occ_out(wtd_s1819_top, "White-tailed Deer", "Summer", "Top")
  # wtd_w1820_occt <- occ_out(wtd_w1820_top, "White-tailed Deer", "Winter", "Top")
  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  summer_occ <- rbind(bob_s1819_occ, coug_s1819_occ, coy_s1819_occ, wolf_s1819_occ,
                      elk_s1819_occ, md_s1819_occ, wtd_s1819_occ)
  winter_occ <- rbind(bob_w1820_occ, coug_w1820_occ, coy_w1820_occ, wolf_w1820_occ,
                      elk_w1820_occ, md_w1820_occ, wtd_w1820_occ)
  occ_results <- rbind(summer_occ, winter_occ) %>%
    arrange(Species)
  colnames(occ_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval") #"Model", 
  #' #'  Top models identified by derdging
  #' summer_occ_top <- rbind(bob_s1819_occt, coug_s1819_occt, coy_s1819_occt, wolf_s1819_occt,
  #'                     elk_s1819_occt, md_s1819_occt, wtd_s1819_occt)
  #' winter_occ_top <- rbind(bob_w1820_occt, coug_w1820_occt, coy_w1820_occt, wolf_w1820_occt,
  #'                     elk_w1820_occt, md_w1820_occt, wtd_w1820_occt)
  #' occ_results_top <- rbind(summer_occ_top, winter_occ_top) %>%
  #'   arrange(Species) 
  #' colnames(occ_results_top) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  #' #'  Combine full model and top model results
  #' summer_occ_combo <- rbind(bob_s1819_occ, bob_s1819_occt, coug_s1819_occ, coug_s1819_occt, 
  #'                         coy_s1819_occ, coy_s1819_occt, wolf_s1819_occ, wolf_s1819_occt,
  #'                         elk_s1819_occ, elk_s1819_occt, md_s1819_occ, md_s1819_occt, 
  #'                         wtd_s1819_occ, wtd_s1819_occt)
  #' winter_occ_combo <- rbind(bob_w1820_occ, bob_w1820_occt, coug_w1820_occ, coug_w1820_occt, 
  #'                         coy_w1820_occ, coy_w1820_occt, wolf_w1820_occ, wolf_w1820_occt,
  #'                         elk_w1820_occ, elk_w1820_occt, md_w1820_occ, md_w1820_occt, 
  #'                         wtd_w1820_occ, wtd_w1820_occt)
  #' occ_results_combo <- rbind(summer_occ_combo, winter_occ_combo) %>%
  #'   arrange(Species) 
  #' colnames(occ_results_combo) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are easier to look at
  rounddig <- 2
  results_psi <- occ_results %>%
    mutate(
      Estimate = round(Estimate, rounddig),
      SE = round(SE, rounddig),
      z = round(z, rounddig),
      Pval = round(Pval, rounddig)
    )
  # results_psi_top <- occ_results_top %>%
  #   mutate(
  #     Estimate = round(Estimate, rounddig),
  #     SE = round(SE, rounddig),
  #     z = round(z, rounddig),
  #     Pval = round(Pval, rounddig)
  #   )
  # results_psi_combo <- occ_results_combo %>%
  #   mutate(
  #     Estimate = round(Estimate, rounddig),
  #     SE = round(SE, rounddig),
  #     z = round(z, rounddig),
  #     Pval = round(Pval, rounddig)
  #   )
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_psi_wide <- results_psi %>%  #results_psi_combo
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept Pval"), sep = "_") %>%
    separate("AreaOK", c("AreaOK (SE)", "AreaOK Pval"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev Pval"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope Pval"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix Pval"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass Pval"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub Pval"), sep = "_") %>%
    separate("RoadDensity", c("RoadDensity (SE)", "RoadDensity Pval"), sep = "_") %>%
    # separate("HumanMod", c("HumanMod (SE)", "HumanMod Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
    arrange(match(Season, c("Summer", "Winter")))
  
  #'  Save!
  write.csv(results_psi, paste0("./Outputs/Tables/OccMod_OccProb_Results_NoHM_", Sys.Date(), ".csv"))  #'  KEEP TRACK of whether HumanMod was excluded from models
  write.csv(results_psi_wide, paste0("./Outputs/Tables/OccMod_OccProb_Results_wide_NoHM_", Sys.Date(), ".csv"))
  
 
  #'  Function to save detection results
  det_out <- function(mod, spp, season, model) {
    out <- summary(mod@estimates)$det %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$det),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
        # Model = rep(model, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
      # relocate(Model, .before = Species)
    return(out)
  }
  
  #'  Run each model through detection function
  bob_s1819_det <- det_out(bob_s1819_global, "Bobcat", "Summer") #, "Global"
  bob_w1820_det <- det_out(bob_w1820_global, "Bobcat", "Winter")
  coug_s1819_det <- det_out(coug_s1819_global, "Cougar", "Summer")
  coug_w1820_det <- det_out(coug_w1820_global, "Cougar", "Winter")
  coy_s1819_det <- det_out(coy_s1819_global, "Coyote", "Summer")
  coy_w1820_det <- det_out(coy_w1820_global, "Coyote", "Winter")
  wolf_s1819_det <- det_out(wolf_s1819_global2, "Wolf", "Summer")
  wolf_w1820_det <- det_out(wolf_w1820_global2, "Wolf", "Winter")
  elk_s1819_det <- det_out(elk_s1819_global2, "Elk", "Summer")
  elk_w1820_det <- det_out(elk_w1820_global2, "Elk", "Winter")
  md_s1819_det <- det_out(md_s1819_global, "Mule Deer", "Summer")
  md_w1820_det <- det_out(md_w1820_global, "Mule Deer", "Winter")
  wtd_s1819_det <- det_out(wtd_s1819_global2, "White-tailed Deer", "Summer")
  wtd_w1820_det <- det_out(wtd_w1820_global2, "White-tailed Deer", "Winter")

  # bob_s1819_dett <- det_out(bob_s1819_top, "Bobcat", "Summer", "Top")
  # bob_w1820_dett <- det_out(bob_w1820_top, "Bobcat", "Winter", "Top")
  # coug_s1819_dett <- det_out(coug_s1819_top, "Cougar", "Summer", "Top")
  # coug_w1820_dett <- det_out(coug_w1820_top, "Cougar", "Winter", "Top")
  # coy_s1819_dett <- det_out(coy_s1819_top, "Coyote", "Summer", "Top")
  # coy_w1820_dett <- det_out(coy_w1820_top, "Coyote", "Winter", "Top")
  # wolf_s1819_dett <- det_out(wolf_s1819_top, "Wolf", "Summer", "Top")
  # wolf_w1820_dett <- det_out(wolf_w1820_top, "Wolf", "Winter", "Top")
  # elk_s1819_dett <- det_out(elk_s1819_top, "Elk", "Summer", "Top")
  # elk_w1820_dett <- det_out(elk_w1820_top, "Elk", "Winter", "Top")
  # md_s1819_dett <- det_out(md_s1819_top, "Mule Deer", "Summer", "Top")
  # md_w1820_dett <- det_out(md_w1820_top, "Mule Deer", "Winter", "Top")
  # wtd_s1819_dett <- det_out(wtd_s1819_top, "White-tailed Deer", "Summer", "Top")
  # wtd_w1820_dett <- det_out(wtd_w1820_top, "White-tailed Deer", "Winter", "Top")
  
  #'  Merge into larger data frames for easy comparison
  summer_det <- rbind(bob_s1819_det, coug_s1819_det, coy_s1819_det, wolf_s1819_det,
                      elk_s1819_det, md_s1819_det, wtd_s1819_det)
  winter_det <- rbind(bob_w1820_det, coug_w1820_det, coy_w1820_det, wolf_w1820_det,
                      elk_w1820_det, md_w1820_det, wtd_w1820_det)
  det_results <- rbind(summer_det, winter_det) %>%
    arrange(Species)
  colnames(det_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval") #"Model", 
  
  # summer_det_top <- rbind(bob_s1819_dett, coug_s1819_dett, coy_s1819_dett, wolf_s1819_dett,
  #                     elk_s1819_dett, md_s1819_dett, wtd_s1819_dett)
  # winter_det_top <- rbind(bob_w1820_dett, coug_w1820_dett, coy_w1820_dett, wolf_w1820_dett,
  #                     elk_w1820_dett, md_w1820_dett, wtd_w1820_dett)
  # det_results_top <- rbind(summer_det_top, winter_det_top) %>%
  #   arrange(Species)
  # colnames(det_results_top) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  # 
  # summer_det_combo <- rbind(bob_s1819_det, bob_s1819_dett, coug_s1819_det, coug_s1819_dett, 
  #                           coy_s1819_det, coy_s1819_dett, wolf_s1819_det, wolf_s1819_dett,
  #                           elk_s1819_det, elk_s1819_dett, md_s1819_det, md_s1819_dett, 
  #                           wtd_s1819_det, wtd_s1819_dett)
  # winter_det_combo <- rbind(bob_w1820_det, bob_w1820_dett, coug_w1820_det, coug_w1820_dett, 
  #                           coy_w1820_det, coy_w1820_dett, wolf_w1820_det, wolf_w1820_dett,
  #                           elk_w1820_det, elk_w1820_dett, md_w1820_det, md_w1820_dett, 
  #                           wtd_w1820_det, wtd_w1820_dett)
  # det_results_combo <- rbind(summer_det_combo, winter_det_combo) %>%
  #   arrange(Species)
  # colnames(det_results_combo) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are a little easier to interpret
  results_det <- det_results %>%
    mutate(
      Estimate = round(Estimate, 2),
      SE = round(SE, 2),
      z = round(z, 2),
      Pval = round(Pval, 2)
    )
  # results_det_top <- det_results_top %>%
  #   mutate(
  #     Estimate = round(Estimate, 3),
  #     SE = round(SE, 3),
  #     z = round(z, 3),
  #     Pval = round(Pval, 3),
  #     Parameter = ifelse(Parameter == "Distance:Height", "Height:Distance", Parameter)
  #   )
  # results_det_combo <- det_results_combo %>%
  #   mutate(
  #     Estimate = round(Estimate, 3),
  #     SE = round(SE, 3),
  #     z = round(z, 3),
  #     Pval = round(Pval, 3),
  #     Parameter = ifelse(Parameter == "Distance:Height", "Height:Distance", Parameter)
  #   )
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_det_wide <- results_det %>% #results_det_combo
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) %>%
    #' #'  Bold significant variables- doesn't work if continue manipulating data frame
    #' condformat(.) %>%
    #' rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept Pval"), sep = "_") %>%
    separate("TrailDirt road", c("Road (SE)", "Road Pval"), sep = "_") %>%
    separate("TrailDecommissioned road", c("Decom Road (SE)", "Decom Road Pval"), sep = "_") %>%
    separate("Height", c("Height (SE)", "Height Pval"), sep = "_") %>%
    separate("Temp_smr", c("Summer Temp (SE)", "Smr Temp Pval"), sep = "_") %>%
    separate("Temp_wtr", c("Winter Temp (SE)", "Wtr Temp Pval"), sep = "_") %>%
    separate("Distance", c("Distance (SE)", "Distance Pval"), sep = "_") %>%
    separate("Height:Distance", c("Height*Distance (SE)", "Height*Distance Pval"), sep = "_") %>%
    separate("YearYear2", c("Year2 (SE)", "Year2 Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Wolf", "Coyote", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
    arrange(match(Season, c("Summer", "Winter")))

  #'  Save!
  write.csv(results_det, paste0("./Outputs/Tables/OccMod_DetProb_Results_NoHM_", Sys.Date(), ".csv"))  #'  KEEP TRACK of whether human mod was excluded from models
  write.csv(results_det_wide, paste0("./Outputs/Tables/OccMod_DetProb_Results_NoHM_wide", Sys.Date(), ".csv"))


  #'  Predict probability of occupancy across sites
  mu_occ <- function(mod, species, season) {
    #'  Predict occupancy probability for all camera sties
    occu_mean <- predict(object = mod, type = "state") %>%  # do I provide newdata = sitecovariates if it's the same values as what went into the model?
      #'  Average occupancy probabilities across sites for mean psi
      summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE)
    #'  Predict occupancy probability for all camera sties
    det_mean <- predict(object = mod, type = "det") %>%
      #'  Average occupancy probabilities across sites for mean psi
      summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE) 
    predicted <- as.data.frame(rbind(occu_mean, det_mean))
    colnames(predicted) <- c("Mean", "SE")
    Parameter <- c("Occupancy", "Detection")
    Species <- species
    Season <- season
    predicted <- cbind(predicted, Parameter)
    predicted <- cbind(predicted, Species)
    predicted <- cbind(predicted, Season)
    return(predicted)
  }
  #'  Estimate mean probability of occupancy and detection per species and season
  md_predict_smr <- mu_occ(md_s1819_global, "Mule Deer", "Summer")
  md_predict_wtr <- mu_occ(md_w1820_global, "Mule Deer", "Winter")
  elk_predict_smr <- mu_occ(elk_s1819_global2, "Elk", "Summer")
  elk_predict_wtr <- mu_occ(elk_w1820_global2, "Elk", "Winter")
  wtd_predict_smr <- mu_occ(wtd_s1819_global2, "White-tailed Deer", "Summer")
  wtd_predict_wtr <- mu_occ(wtd_w1820_global2, "White-tailed Deer", "Winter")
  coug_predict_smr <- mu_occ(coug_s1819_global, "Cougar", "Summer")
  coug_predict_wtr <- mu_occ(coug_w1820_global, "Cougar", "Winter")
  wolf_predict_smr <- mu_occ(wolf_s1819_global2, "Wolf", "Summer")
  wolf_predict_wtr <- mu_occ(wolf_w1820_global2, "Wolf", "Winter")
  bob_predict_smr <- mu_occ(bob_s1819_global, "Bobcat", "Summer")
  bob_predict_wtr <- mu_occ(bob_w1820_global, "Bobcat", "Winter")
  coy_predict_smr <- mu_occ(coy_s1819_global, "Coyote", "Summer")
  coy_predict_wtr <- mu_occ(coy_w1820_global, "Coyote", "Winter")

  #'  Make a pretty table
  Mean_tbl <- bind_rows(md_predict_smr, md_predict_wtr, elk_predict_smr, elk_predict_wtr, wtd_predict_smr, 
              wtd_predict_wtr, coug_predict_smr, coug_predict_wtr, wolf_predict_smr, 
              wolf_predict_wtr, bob_predict_smr, bob_predict_wtr, coy_predict_smr, 
              coy_predict_wtr) %>%
    relocate(Species, .before = Mean) %>%
    relocate(Season, .after = Species) %>%
    relocate(Parameter, .after = Season) %>%
    arrange(Parameter, Mean, Species)
  
  #'  Save
  write.csv(Mean_tbl, paste0("./Outputs/Tables/OccMod_Mean_Estimates_NoHM_", Sys.Date(), ".csv"))  #"  Keep track of whether human mod was excluded from analyses

 
  #'  Save workspace
  save.image(file = "./Outputs/OccMod_script_results.RData")
  