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
  (gs_cougmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_graze <- occuMulti(detFormulas_graze, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougmd_fld <- fitList(gs_cougmd_trail, gs_cougmd_graze)
  #' Model selection
  modSel(gs_cougmd_fld)
  
  (gs_cougmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, coug_md_grazing_UMF, silent = TRUE))
  (gs_cougmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, coug_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougmd_fl <- fitList(gs_cougmd_null0, gs_cougmd_null1, gs_cougmd_hab0, gs_cougmd_hab1, gs_cougmd_graze0, gs_cougmd_graze1, gs_cougmd_graze2)
  #' Model selection
  modSel(gs_cougmd_fl)
  summary(gs_cougmd_hab0)
  summary(gs_cougmd_hab1) # f12 not significant
  
  ####  Cougar-ELK Grazing Season  ####
  (gs_cougelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougelk_fld <- fitList(gs_cougelk_trail, gs_cougelk_dgraze)
  #' Model selection
  modSel(gs_cougelk_fld)
  
  (gs_cougelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, coug_elk_grazing_UMF, silent = TRUE))
  (gs_cougelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, coug_elk_grazing_UMF, silent = TRUE)) # fails
  (gs_cougelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, coug_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougelk_fl <- fitList(gs_cougelk_null0, gs_cougelk_null1, gs_cougelk_hab0, gs_cougelk_hab1, gs_cougelk_graze0, gs_cougelk_graze2)
  #' Model selection
  modSel(gs_cougelk_fl)
  summary(gs_cougelk_hab0)
  summary(gs_cougelk_hab1) # f12 not significant
  
  ####  Cougar-White-tailed Deer Grazing Season  ####
  (gs_cougwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougwtd_fld <- fitList(gs_cougwtd_trail, gs_cougwtd_dgraze)
  #' Model selection
  modSel(gs_cougwtd_fld)
  
  (gs_cougwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_wtd_grazing_UMF, silent = TRUE))
  (gs_cougwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougwtd_fl <- fitList(gs_cougwtd_null0, gs_cougwtd_null1, gs_cougwtd_hab0, gs_cougwtd_hab1, gs_cougwtd_graze0, gs_cougwtd_graze1, gs_cougwtd_graze2)
  #' Model selection
  modSel(gs_cougwtd_fl)
  #summary(gs_cougwtd_graze2) # f12 not significant
  summary(gs_cougwtd_graze0)
  summary(gs_cougwtd_graze1) # f12 not significant
  
  ####  Cougar-Moose Grazing Season  ####
  (gs_cougmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougmoose_fld <- fitList(gs_cougmoose_trail, gs_cougmoose_dgraze)
  #' Model selection
  modSel(gs_cougmoose_fld)
  
  (gs_cougmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coug_moose_grazing_UMF, silent = TRUE))
  (gs_cougmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coug_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_cougmoose_fl <- fitList(gs_cougmoose_null0, gs_cougmoose_null1, gs_cougmoose_hab0, gs_cougmoose_hab1, gs_cougmoose_graze0, gs_cougmoose_graze1, gs_cougmoose_graze2)
  #' Model selection
  modSel(gs_cougmoose_fl)
  summary(gs_cougmoose_hab0)
  summary(gs_cougmoose_hab1) # f12 not significant
  summary(gs_cougmoose_graze0)
  
  
  ####  Wolf-Mule Deer Grazing Season  ####
  (gs_wolfmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfmd_fld <- fitList(gs_wolfmd_trail, gs_wolfmd_dgraze)
  #' Model selection
  modSel(gs_wolfmd_fld)
  
  (gs_wolfmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, wolf_md_grazing_UMF, silent = TRUE))
  (gs_wolfmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, wolf_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfmd_fl <- fitList(gs_wolfmd_null0, gs_wolfmd_null1, gs_wolfmd_hab0, gs_wolfmd_hab1, gs_wolfmd_graze0, gs_wolfmd_graze1, gs_wolfmd_graze2)
  #' Model selection
  modSel(gs_wolfmd_fl)
  summary(gs_wolfmd_hab0)
  summary(gs_wolfmd_hab1) # f12 not significant
  
  
  ####  Wolf-Elk Grazing Season  ####
  (gs_wolfelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfelk_fld <- fitList(gs_wolfelk_trail, gs_wolfelk_dgraze)
  #' Model selection
  modSel(gs_wolfelk_fld)
  
  (gs_wolfelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, wolf_elk_grazing_UMF, silent = TRUE))
  (gs_wolfelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, wolf_elk_grazing_UMF, silent = TRUE)) 
  (gs_wolfelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, wolf_elk_grazing_UMF, silent = TRUE)) # FAILS
  #' List of fitted models
  gs_wolfelk_fl <- fitList(gs_wolfelk_null0, gs_wolfelk_null1, gs_wolfelk_hab0, gs_wolfelk_hab1, gs_wolfelk_graze0, gs_wolfelk_graze1)
  #' Model selection
  modSel(gs_wolfelk_fl)
  summary(gs_wolfelk_hab0)
  summary(gs_wolfelk_graze0) # not adding any significant info
  summary(gs_wolfelk_hab1) # f1 not significant
  
  ####  Wolf-White-tailed Deer Grazing Season  ####
  (gs_wolfwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfwtd_fld <- fitList(gs_wolfwtd_trail, gs_wolfwtd_dgraze)
  #' Model selection
  modSel(gs_wolfwtd_fld)
  
  (gs_wolfwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_wtd_grazing_UMF, silent = TRUE))
  (gs_wolfwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfwtd_fl <- fitList(gs_wolfwtd_null0, gs_wolfwtd_null1, gs_wolfwtd_hab0, gs_wolfwtd_hab1, gs_wolfwtd_graze0, gs_wolfwtd_graze1, gs_wolfwtd_graze2)
  #' Model selection
  modSel(gs_wolfwtd_fl)
  summary(gs_wolfwtd_graze0)
  summary(gs_wolfwtd_graze1) # f12 not significant
  
  ####  Wolf-Moose Grazing Season  ####
  (gs_wolfmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfmoose_fld <- fitList(gs_wolfmoose_trail, gs_wolfmoose_dgraze)
  #' Model selection
  modSel(gs_wolfmoose_fld)
  
  (gs_wolfmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, wolf_moose_grazing_UMF, silent = TRUE))
  (gs_wolfmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, wolf_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_wolfmoose_fl <- fitList(gs_wolfmoose_null0, gs_wolfmoose_null1, gs_wolfmoose_hab0, gs_wolfmoose_hab1, gs_wolfmoose_graze0, gs_wolfmoose_graze1, gs_wolfmoose_graze2)
  #' Model selection
  modSel(gs_wolfmoose_fl)
  summary(gs_wolfmoose_graze0)
  summary(gs_wolfmoose_hab0)
  summary(gs_wolfmoose_graze1) # f12 not significant
  
  
  
  ####  Bear-Mule Deer Grazing Season  ####
  (gs_bearmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearmd_fld <- fitList(gs_bearmd_trail, gs_bearmd_dgraze)
  #' Model selection
  modSel(gs_bearmd_fld)
  
  (gs_bearmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bear_md_grazing_UMF, silent = TRUE))
  (gs_bearmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bear_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearmd_fl <- fitList(gs_bearmd_null0, gs_bearmd_null1, gs_bearmd_hab0, gs_bearmd_hab1, gs_bearmd_graze0, gs_bearmd_graze1, gs_bearmd_graze2)
  #' Model selection
  modSel(gs_bearmd_fl)
  summary(gs_bearmd_graze1) # f12 not significant
  summary(gs_bearmd_graze0) # grazing not actually significant
  summary(gs_bearmd_graze2) # f12 and grazing not significant
  summary(gs_bearmd_hab1) # f12 not significant
  summary(gs_bearmd_hab0)
  
  ####  Bear-Elk Grazing Season  ####
  (gs_bearelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearelk_fld <- fitList(gs_bearelk_trail, gs_bearelk_dgraze)
  #' Model selection
  modSel(gs_bearelk_fld)
  
  (gs_bearelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bear_elk_grazing_UMF, silent = TRUE))
  (gs_bearelk_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bear_elk_grazing_UMF, silent = TRUE)) 
  (gs_bearelk_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bear_elk_grazing_UMF, silent = TRUE)) # FAILS
  #' List of fitted models
  gs_bearelk_fl <- fitList(gs_bearelk_null0, gs_bearelk_null1, gs_bearelk_hab0, gs_bearelk_hab1, gs_bearelk_graze0, gs_bearelk_graze1)
  #' Model selection
  modSel(gs_bearelk_fl)
  summary(gs_bearelk_graze1) # f12 not significant
  summary(gs_bearelk_graze0) # grazing not significant
  summary(gs_bearelk_hab1) # f12 not signigicant
  summary(gs_bearelk_hab0)
  
  ####  Bear-White-tailed Deer Grazing Season  ####
  (gs_bearwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearwtd_fld <- fitList(gs_bearwtd_trail, gs_bearwtd_dgraze)
  #' Model selection
  modSel(gs_bearwtd_fld)
  
  (gs_bearwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_wtd_grazing_UMF, silent = TRUE))
  (gs_bearwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearwtd_fl <- fitList(gs_bearwtd_null0, gs_bearwtd_null1, gs_bearwtd_hab0, gs_bearwtd_hab1, gs_bearwtd_graze0, gs_bearwtd_graze1, gs_bearwtd_graze2)
  #' Model selection
  modSel(gs_bearwtd_fl)
  summary(gs_bearwtd_graze0)
  
  
  ####  Bear-Moose Grazing Season  ####
  (gs_bearmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearmoose_fld <- fitList(gs_bearmoose_trail, gs_bearmoose_dgraze)
  #' Model selection
  modSel(gs_bearmoose_fld)
  
  (gs_bearmoose_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bear_moose_grazing_UMF, silent = TRUE))
  (gs_bearmoose_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bear_moose_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bearmoose_fl <- fitList(gs_bearmoose_null0, gs_bearmoose_null1, gs_bearmoose_hab0, gs_bearmoose_hab1, gs_bearmoose_graze0, gs_bearmoose_graze1, gs_bearmoose_graze2)
  #' Model selection
  modSel(gs_bearmoose_fl)
  summary(gs_bearmoose_graze0)
  summary(gs_bearmoose_graze2) # f12 not significant
  summary(gs_bearmoose_graze1) # f12 not significant
  
  
  
  ####  Coyote-Mule Deer Grazing Season  ####
  (gs_coymd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_coymd_fld <- fitList(gs_coymd_trail, gs_coymd_dgraze)
  #' Model selection
  modSel(gs_coymd_fld)
  
  (gs_coymd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_md_grazing_UMF, silent = TRUE))
  (gs_coymd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_coymd_fl <- fitList(gs_coymd_null0, gs_coymd_null1, gs_coymd_hab0, gs_coymd_hab1, gs_coymd_graze0, gs_coymd_graze1, gs_coymd_graze2)
  #' Model selection
  modSel(gs_coymd_fl)
  summary(gs_coymd_graze2)
  
  ####  Coyote-White-tailed Deer Grazing Season  ####
  (gs_coywtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_coywtd_fld <- fitList(gs_coywtd_trail, gs_coywtd_dgraze)
  #' Model selection
  modSel(gs_coywtd_fld)
  
  (gs_coywtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, coy_wtd_grazing_UMF, silent = TRUE))
  (gs_coywtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, coy_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_coywtd_fl <- fitList(gs_coywtd_null0, gs_coywtd_null1, gs_coywtd_hab0, gs_coywtd_hab1, gs_coywtd_graze0, gs_coywtd_graze1, gs_coywtd_graze2)
  #' Model selection
  modSel(gs_coywtd_fl)
  summary(gs_coywtd_graze2) 
  summary(gs_coywtd_graze1)
  
  
  
  ####  Bobcat-Mule Deer Grazing Season  ####
  (gs_bobmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bobmd_fld <- fitList(gs_bobmd_trail, gs_bobmd_dgraze)
  #' Model selection
  modSel(gs_bobmd_fld)
  
  (gs_bobmd_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_graze0 <- occuMulti(detFormulas_trail, occFormulas_graze0, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_graze1 <- occuMulti(detFormulas_trail, occFormulas_graze1, bob_md_grazing_UMF, silent = TRUE))
  (gs_bobmd_graze2 <- occuMulti(detFormulas_trail, occFormulas_graze2, bob_md_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bobmd_fl <- fitList(gs_bobmd_null0, gs_bobmd_null1, gs_bobmd_hab0, gs_bobmd_hab1, gs_bobmd_graze0, gs_bobmd_graze1, gs_bobmd_graze2)
  #' Model selection
  modSel(gs_bobmd_fl)
  summary(gs_bobmd_graze0)
  summary(gs_bobmd_graze2) # f12 not significant
  summary(gs_bobmd_graze1) # f12 not significant
  
  ####  Bobcat-White-tailed Deer Grazing Season  ####
  (gs_bobwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_dgraze <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bobwtd_fld <- fitList(gs_bobwtd_trail, gs_bobwtd_dgraze)
  #' Model selection
  modSel(gs_bobwtd_fld)
  
  (gs_bobwtd_null0 <- occuMulti(detFormulas_graze, occFormulas_null, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_null1 <- occuMulti(detFormulas_graze, occFormulas_null1, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_hab0 <- occuMulti(detFormulas_graze, occFormulas_hab0, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_hab1 <- occuMulti(detFormulas_graze, occFormulas_hab1, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_graze0 <- occuMulti(detFormulas_graze, occFormulas_graze0, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_graze1 <- occuMulti(detFormulas_graze, occFormulas_graze1, bob_wtd_grazing_UMF, silent = TRUE))
  (gs_bobwtd_graze2 <- occuMulti(detFormulas_graze, occFormulas_graze2, bob_wtd_grazing_UMF, silent = TRUE))
  #' List of fitted models
  gs_bobwtd_fl <- fitList(gs_bobwtd_null0, gs_bobwtd_null1, gs_bobwtd_hab0, gs_bobwtd_hab1, gs_bobwtd_graze0, gs_bobwtd_graze1, gs_bobwtd_graze2)
  #' Model selection
  modSel(gs_bobwtd_fl)
  summary(gs_bobwtd_graze2) # f12 and grazing not actually significant
  summary(gs_bobwtd_graze0)
  summary(gs_bobwtd_graze1) # f12 not significant

  
  
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
  # occFormulas_huntpubish0 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
  #                          "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                          "~0")
  # occFormulas_huntpubish1 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
  #                          "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                          "~1")
  # occFormulas_huntpubish2 <- c("~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity", 
  #                          "~Study_Area + Elev + I(Elev^2) + PercForest + HuntingActivity + Public", 
  #                          "~HuntingActivity + Public")
  
  
  ####  Cougar-Mule Deer Hunting Season  ####
  #'  Note: using occupancy model formulas that drop Public covariate from predator sub-model
  (hs_cougmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougmd_fld <- fitList(hs_cougmd_trail, hs_cougmd_hunt, hs_cougmd_pub, hs_cougmd_pubhunt)
  #' Model selection
  modSel(hs_cougmd_fld)
  
  (hs_cougmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_md_hunting_UMF, silent = TRUE)) 
  (hs_cougmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_md_hunting_UMF, silent = TRUE)) 
  (hs_cougmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_md_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_md_hunting_UMF, silent = TRUE)) 
  (hs_cougmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_md_hunting_UMF, silent = TRUE))
  (hs_cougmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougmd_fl <- fitList(hs_cougmd_null0, hs_cougmd_null1, hs_cougmd_hab0, hs_cougmd_hab1, hs_cougmd_hunt0, hs_cougmd_hunt1, hs_cougmd_pub0, hs_cougmd_pub1, hs_cougmd_pub2) 
  #' Model selection
  modSel(hs_cougmd_fl)
  summary(hs_cougmd_pub2) 
  summary(hs_cougmd_hunt0)
  summary(hs_cougmd_hunt1)

  ####  Cougar-Elk Hunting Season  ####
  (hs_cougelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougelk_fld <- fitList(hs_cougelk_trail, hs_cougelk_hunt, hs_cougelk_pub, hs_cougelk_pubhunt)
  #' Model selection
  modSel(hs_cougelk_fld)
  
  (hs_cougelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, coug_elk_hunting_UMF, silent = TRUE)) 
  (hs_cougelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pubish0, coug_elk_hunting_UMF, silent = TRUE))
  (hs_cougelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pubish1, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pubish2, coug_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_cougelk_fl <- fitList(hs_cougelk_null0, hs_cougelk_hab0, hs_cougelk_hunt0, hs_cougelk_pub0) 
  #' Model selection
  modSel(hs_cougelk_fl)
  summary(hs_cougelk_pub0)
  
  ####  Cougar-White-tailed Deer Hunting Season  ####
  (hs_cougwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougwtd_fld <- fitList(hs_cougwtd_trail, hs_cougwtd_hunt, hs_cougwtd_pub, hs_cougwtd_pubhunt)
  #' Model selection
  modSel(hs_cougwtd_fld)
  
  (hs_cougwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_wtd_hunting_UMF, silent = TRUE)) 
  (hs_cougwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_wtd_hunting_UMF, silent = TRUE)) 
  (hs_cougwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_wtd_hunting_UMF, silent = TRUE))
  (hs_cougwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_wtd_hunting_UMF, silent = TRUE)) 
  (hs_cougwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougwtd_fl <- fitList(hs_cougwtd_null0, hs_cougwtd_null1, hs_cougwtd_hab0, hs_cougwtd_hab1, hs_cougwtd_hunt0, hs_cougwtd_hunt1, hs_cougwtd_hunt2, hs_cougwtd_pub0, hs_cougwtd_pub1, hs_cougwtd_pub2) 
  #' Model selection
  modSel(hs_cougwtd_fl)
  summary(hs_cougwtd_pub2)
  summary(hs_cougwtd_pub0)
  
  ####  Cougar-Moose Hunting Season  ####
  (hs_cougmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_pub <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougmoose_fld <- fitList(hs_cougmoose_trail, hs_cougmoose_hunt, hs_cougmoose_pub, hs_cougmoose_pubhunt)
  #' Model selection
  modSel(hs_cougmoose_fld)
  
  (hs_cougmoose_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coug_moose_hunting_UMF, silent = TRUE)) 
  (hs_cougmoose_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coug_moose_hunting_UMF, silent = TRUE)) # FAIL
  (hs_cougmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pubish0, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pubish1, coug_moose_hunting_UMF, silent = TRUE))
  (hs_cougmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pubish2, coug_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_cougmoose_fl <- fitList(hs_cougmoose_null0, hs_cougmoose_null1, hs_cougmoose_hab0, hs_cougmoose_hunt0, hs_cougmoose_pub0, hs_cougmoose_pub1, hs_cougmoose_pub2) 
  #' Model selection
  modSel(hs_cougmoose_fl)
  summary(hs_cougmoose_hunt0)
  
  
  ####  Wolf-Mule Deer Hunting Season  ####
  #'  Note: using detection & occupancy model formulas that drop Public covariate from predator sub-models
  (hs_wolfmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_wolfmd_fld <- fitList(hs_wolfmd_trail, hs_wolfmd_hunt, hs_wolfmd_pub, hs_wolfmd_pubhunt) 
  #' Model selection
  modSel(hs_wolfmd_fld)
  
  (hs_wolfmd_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_md_hunting_UMF, silent = TRUE)) 
  (hs_wolfmd_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_md_hunting_UMF, silent = TRUE)) 
  (hs_wolfmd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_md_hunting_UMF, silent = TRUE))
  (hs_wolfmd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_md_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_wolfmd_fl <- fitList(hs_wolfmd_null0, hs_wolfmd_null1, hs_wolfmd_hab0, hs_wolfmd_hab1, hs_wolfmd_hunt0, hs_wolfmd_hunt1, hs_wolfmd_hunt2, hs_wolfmd_pub0, hs_wolfmd_pub1) 
  #' Model selection
  modSel(hs_wolfmd_fl)
  summary(hs_wolfmd_hunt1)
  summary(hs_wolfmd_hunt2)
  
  
  ####  Wolf-Elk Hunting Season  ####
  (hs_wolfelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_wolfelk_fld <- fitList(hs_wolfelk_trail, hs_wolfelk_hunt, hs_wolfelk_pub, hs_wolfelk_pubhunt) 
  #' Model selection
  modSel(hs_wolfelk_fld)
  
  (hs_wolfelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, wolf_elk_hunting_UMF, silent = TRUE)) 
  (hs_wolfelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, wolf_elk_hunting_UMF, silent = TRUE)) 
  (hs_wolfelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pubish0, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pubish1, wolf_elk_hunting_UMF, silent = TRUE))
  (hs_wolfelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pubish2, wolf_elk_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_wolfelk_fl <- fitList(hs_wolfelk_null0, hs_wolfelk_null1, hs_wolfelk_hab0, hs_wolfelk_hab1, hs_wolfelk_hunt0, hs_wolfelk_hunt1, hs_wolfelk_hunt2, hs_wolfelk_pub0, hs_wolfelk_pub1) 
  #' Model selection
  modSel(hs_wolfelk_fl)
  summary(hs_wolfelk_hunt0)
  
  ####  Wolf-White-tailed Deer Hunting Season  ####
  (hs_wolfwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_wolfwtd_fld <- fitList(hs_wolfwtd_trail, hs_wolfwtd_hunt, hs_wolfwtd_pub, hs_wolfwtd_pubhunt)
  #' Model selection
  modSel(hs_wolfwtd_fld)
  
  (hs_wolfwtd_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_wtd_hunting_UMF, silent = TRUE)) 
  (hs_wolfwtd_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_wtd_hunting_UMF, silent = TRUE)) 
  (hs_wolfwtd_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  (hs_wolfwtd_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_wtd_hunting_UMF, silent = TRUE))  
  (hs_wolfwtd_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_wtd_hunting_UMF, silent = TRUE))
  (hs_wolfwtd_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_wtd_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_wolfwtd_fl <- fitList(hs_wolfwtd_null0, hs_wolfwtd_null1, hs_wolfwtd_hab0, hs_wolfwtd_hab1, hs_wolfwtd_hunt0, hs_wolfwtd_hunt2, hs_wolfwtd_pub0, hs_wolfwtd_pub1) 
  #' Model selection
  modSel(hs_wolfwtd_fl)
  summary(hs_wolfwtd_hunt0)
  
  ####  Wolf-Moose Hunting Season  ####
  (hs_wolfmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_pub <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_pubhunt <- occuMulti(detFormulas_pubhuntish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_wolfmoose_fld <- fitList(hs_wolfmoose_trail, hs_wolfmoose_hunt, hs_wolfmoose_pub, hs_wolfmoose_pubhunt)
  #' Model selection
  modSel(hs_wolfmoose_fld)
  
  (hs_wolfmoose_null0 <- occuMulti(detFormulas_pubish, occFormulas_null, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_null1 <- occuMulti(detFormulas_pubish, occFormulas_null1, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_hab0 <- occuMulti(detFormulas_pubish, occFormulas_hab0, wolf_moose_hunting_UMF, silent = TRUE)) 
  (hs_wolfmoose_hab1 <- occuMulti(detFormulas_pubish, occFormulas_hab1, wolf_moose_hunting_UMF, silent = TRUE)) 
  (hs_wolfmoose_hunt0 <- occuMulti(detFormulas_pubish, occFormulas_hunt0, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_hunt1 <- occuMulti(detFormulas_pubish, occFormulas_hunt1, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_hunt2 <- occuMulti(detFormulas_pubish, occFormulas_hunt2, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_pub0 <- occuMulti(detFormulas_pubish, occFormulas_pubish0, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_pub1 <- occuMulti(detFormulas_pubish, occFormulas_pubish1, wolf_moose_hunting_UMF, silent = TRUE))
  (hs_wolfmoose_pub2 <- occuMulti(detFormulas_pubish, occFormulas_pubish2, wolf_moose_hunting_UMF, silent = TRUE)) # FAIL
  #' List of fitted models
  hs_wolfmoose_fl <- fitList(hs_wolfmoose_null0, hs_wolfmoose_null1, hs_wolfmoose_hab0, hs_wolfmoose_hab1, hs_wolfmoose_hunt0, hs_wolfmoose_hunt1, hs_wolfmoose_hunt2, hs_wolfmoose_pub0, hs_wolfmoose_pub1)
  #' Model selection
  modSel(hs_wolfmoose_fl)
  summary(hs_wolfmoose_hunt2)
  
  
  ####  Bear-Mule Deer Hunting Season  ####
  (hs_bearmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearmd_fld <- fitList(hs_bearmd_trail, hs_bearmd_hunt, hs_bearmd_pub, hs_bearmd_pubhunt)
  #' Model selection
  modSel(hs_bearmd_fld)
  
  (hs_bearmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_md_hunting_UMF, silent = TRUE)) 
  (hs_bearmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_md_hunting_UMF, silent = TRUE)) 
  (hs_bearmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_md_hunting_UMF, silent = TRUE)) 
  (hs_bearmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_md_hunting_UMF, silent = TRUE))
  (hs_bearmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearmd_fl <- fitList(hs_bearmd_null0, hs_bearmd_null1, hs_bearmd_hab0, hs_bearmd_hab1, hs_bearmd_hunt0, hs_bearmd_hunt1, hs_bearmd_hunt2, hs_bearmd_pub0, hs_bearmd_pub1, hs_bearmd_pub2) 
  #' Model selection
  modSel(hs_bearmd_fl)
  summary(hs_bearmd_hab0)
  summary(hs_bearmd_hab1)
  
  ####  Bear-Elk Hunting Season  ####
  (hs_bearelk_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearelk_fld <- fitList(hs_bearelk_trail, hs_bearelk_hunt, hs_bearelk_pub, hs_bearelk_pubhunt)
  #' Model selection
  modSel(hs_bearelk_fld)
  
  (hs_bearelk_null0 <- occuMulti(detFormulas_trail, occFormulas_null, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_null1 <- occuMulti(detFormulas_trail, occFormulas_null1, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_hab0 <- occuMulti(detFormulas_trail, occFormulas_hab0, bear_elk_hunting_UMF, silent = TRUE)) 
  (hs_bearelk_hab1 <- occuMulti(detFormulas_trail, occFormulas_hab1, bear_elk_hunting_UMF, silent = TRUE)) 
  (hs_bearelk_hunt0 <- occuMulti(detFormulas_trail, occFormulas_hunt0, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_hunt1 <- occuMulti(detFormulas_trail, occFormulas_hunt1, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_hunt2 <- occuMulti(detFormulas_trail, occFormulas_hunt2, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_pub0 <- occuMulti(detFormulas_trail, occFormulas_pub0, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_pub1 <- occuMulti(detFormulas_trail, occFormulas_pub1, bear_elk_hunting_UMF, silent = TRUE))
  (hs_bearelk_pub2 <- occuMulti(detFormulas_trail, occFormulas_pub2, bear_elk_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  hs_bearelk_fl <- fitList(hs_bearelk_null0, hs_bearelk_null1, hs_bearelk_hab0, hs_bearelk_hab1, hs_bearelk_hunt0, hs_bearelk_hunt1, hs_bearelk_hunt2, hs_bearelk_pub0, hs_bearelk_pub1)
  #' Model selection
  modSel(hs_bearelk_fl)
  summary(hs_bearelk_pub1) # f12 not significant
  summary(hs_bearelk_hab1) # f12 not significant
  
  
  ####  Bear-White-tailed Deer Hunting Season  ####
  (hs_bearwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearwtd_fld <- fitList(hs_bearwtd_trail, hs_bearwtd_hunt, hs_bearwtd_pub, hs_bearwtd_pubhunt)
  #' Model selection
  modSel(hs_bearwtd_fld)
  
  (hs_bearwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bearwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_wtd_hunting_UMF, silent = TRUE))
  (hs_bearwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_wtd_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearwtd_fl <- fitList(hs_bearwtd_null0, hs_bearwtd_null1, hs_bearwtd_hab0, hs_bearwtd_hab1, hs_bearwtd_hunt0, hs_bearwtd_hunt1, hs_bearwtd_hunt2, hs_bearwtd_pub0, hs_bearwtd_pub1, hs_bearwtd_pub2) 
  #' Model selection
  modSel(hs_bearwtd_fl)
  summary(hs_bearwtd_pub0)
  summary(hs_bearwtd_hab0)
  
  ####  Bear-Moose Hunting Season  ####
  (hs_bearmoose_trail <- occuMulti(detFormulas_trail, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_pub <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bearmoose_fld <- fitList(hs_bearmoose_trail, hs_bearmoose_hunt, hs_bearmoose_pub, hs_bearmoose_pubhunt)
  #' Model selection
  modSel(hs_bearmoose_fld)

  (hs_bearmoose_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bear_moose_hunting_UMF, silent = TRUE)) 
  (hs_bearmoose_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bear_moose_hunting_UMF, silent = TRUE)) 
  (hs_bearmoose_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  (hs_bearmoose_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bear_moose_hunting_UMF, silent = TRUE)) 
  (hs_bearmoose_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bear_moose_hunting_UMF, silent = TRUE))
  (hs_bearmoose_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bear_moose_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  hs_bearmoose_fl <- fitList(hs_bearmoose_null0, hs_bearmoose_null1, hs_bearmoose_hab0, hs_bearmoose_hab1, hs_bearmoose_hunt0, hs_bearmoose_hunt2, hs_bearmoose_pub0, hs_bearmoose_pub1) 
  #' Model selection
  modSel(hs_bearmoose_fl)
  summary(hs_bearmoose_hunt0)
  summary(hs_bearmoose_hab0)
  
  
  ####  Coyote-Mule Deer Hunting Season  ####
  (hs_coymd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  hs_coymd_fld <- fitList(hs_coymd_trail, hs_coymd_hunt, hs_coymd_pub, hs_coymd_pubhunt)
  #' Model selection
  modSel(hs_coymd_fld)
  
  (hs_coymd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_md_hunting_UMF, silent = TRUE)) # FAIL
  (hs_coymd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_md_hunting_UMF, silent = TRUE)) 
  (hs_coymd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_md_hunting_UMF, silent = TRUE))
  (hs_coymd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_md_hunting_UMF, silent = TRUE))
  # (coymd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, coy_md_hunting_UMF, silent = TRUE)) 
  # (coymd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coy_md_hunting_UMF, silent = TRUE))
  # (coymd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coy_md_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  hs_coymd_fl <- fitList(hs_coymd_null0, hs_coymd_null1, hs_coymd_hab0, hs_coymd_hab1, hs_coymd_hunt0, hs_coymd_hunt1, hs_coymd_pub0, hs_coymd_pub1, hs_coymd_pub2)
  #' Model selection
  modSel(hs_coymd_fl)
  summary(hs_coymd_hunt0)
  summary(hs_coymd_hunt1)
  
  ####  Coyote-White-tailed Deer Hunting Season  ####
  (hs_coywtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  hs_coywtd_fld <- fitList(hs_coywtd_trail, hs_coywtd_hunt, hs_coywtd_pub, hs_coywtd_pubhunt)
  #' Model selection
  modSel(hs_coywtd_fld)
  
  (hs_coywtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, coy_wtd_hunting_UMF, silent = TRUE)) # FAIL
  (hs_coywtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, coy_wtd_hunting_UMF, silent = TRUE)) 
  (hs_coywtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, coy_wtd_hunting_UMF, silent = TRUE))
  (hs_coywtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, coy_wtd_hunting_UMF, silent = TRUE))
  # (hs_coywtd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, coy_wtd_hunting_UMF, silent = TRUE)) 
  # (hs_coywtd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, coy_wtd_hunting_UMF, silent = TRUE))
  # (hs_coywtd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, coy_wtd_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  hs_coywtd_fl <- fitList(hs_coywtd_null0, hs_coywtd_null1, hs_coywtd_hab0, hs_coywtd_hab1, hs_coywtd_hunt0, hs_coywtd_hunt1, hs_coywtd_pub0, hs_coywtd_pub1, hs_coywtd_pub2)
  #' Model selection
  modSel(hs_coywtd_fl)
  summary(hs_coywtd_hunt1)
  
  
  ####  Bobcat-Mule Deer Hunting Season  ####
  (hs_bobmd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  hs_bobmd_fld <- fitList(hs_bobmd_trail, hs_bobmd_hunt, hs_bobmd_pub, hs_bobmd_pubhunt)
  #' Model selection
  modSel(hs_bobmd_fld)
  
  (hs_bobmd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_md_hunting_UMF, silent = TRUE)) 
  (hs_bobmd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_md_hunting_UMF, silent = TRUE))
  (hs_bobmd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_md_hunting_UMF, silent = TRUE))
  # (hs_bobmd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, bob_md_hunting_UMF, silent = TRUE)) 
  # (hs_bobmd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, bob_md_hunting_UMF, silent = TRUE))
  # (hs_bobmd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, bob_md_hunting_UMF, silent = TRUE))
  #' List of fitted models
  hs_bobmd_fl <- fitList(hs_bobmd_null0, hs_bobmd_null1, hs_bobmd_hab0, hs_bobmd_hab1, hs_bobmd_hunt0, hs_bobmd_hunt1, hs_bobmd_hunt2, hs_bobmd_pub0, hs_bobmd_pub1, hs_bobmd_pub2)
  #' Model selection
  modSel(hs_bobmd_fl)
  summary(hs_bobmd_hunt0) # ok but the hunting activity effect on bobcat sub-model has unusually large effect size...
  summary(hs_bobmd_hunt1)
  summary(hs_bobmd_pub0)
  
  
  ####  Bobcat-White-tailed Deer Hunting Season  ####
  (hs_bobwtd_trail <- occuMulti(detFormulas_trail, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_hunt <- occuMulti(detFormulas_hunt, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_pub <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_pubhunt <- occuMulti(detFormulas_pubhunt, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  #' List of fitted models
  hs_bobwtd_fld <- fitList(hs_bobwtd_trail, hs_bobwtd_hunt, hs_bobwtd_pub, hs_bobwtd_pubhunt)
  #' Model selection
  modSel(hs_bobwtd_fld)
  
  (hs_bobwtd_null0 <- occuMulti(detFormulas_pub, occFormulas_null, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_null1 <- occuMulti(detFormulas_pub, occFormulas_null1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_hab0 <- occuMulti(detFormulas_pub, occFormulas_hab0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_hab1 <- occuMulti(detFormulas_pub, occFormulas_hab1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_hunt0 <- occuMulti(detFormulas_pub, occFormulas_hunt0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_hunt1 <- occuMulti(detFormulas_pub, occFormulas_hunt1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_hunt2 <- occuMulti(detFormulas_pub, occFormulas_hunt2, bob_wtd_hunting_UMF, silent = TRUE)) # FAIL massive effect sizes and sign switches on intercept
  (hs_bobwtd_pub0 <- occuMulti(detFormulas_pub, occFormulas_pub0, bob_wtd_hunting_UMF, silent = TRUE)) 
  (hs_bobwtd_pub1 <- occuMulti(detFormulas_pub, occFormulas_pub1, bob_wtd_hunting_UMF, silent = TRUE))
  (hs_bobwtd_pub2 <- occuMulti(detFormulas_pub, occFormulas_pub2, bob_wtd_hunting_UMF, silent = TRUE))
  # (hs_bobwtd_pubhunt0 <- occuMulti(detFormulas_pub, occFormulas_huntpubish0, bob_wtd_hunting_UMF, silent = TRUE)) 
  # (hs_bobwtd_pubhunt1 <- occuMulti(detFormulas_pub, occFormulas_huntpubish1, bob_wtd_hunting_UMF, silent = TRUE)) #FAIL?
  # (hs_bobwtd_pubhunt2 <- occuMulti(detFormulas_pub, occFormulas_huntpubish2, bob_wtd_hunting_UMF, silent = TRUE)) #FAIL
  #' List of fitted models
  hs_bobwtd_fl <- fitList(hs_bobwtd_null0, hs_bobwtd_null1, hs_bobwtd_hab0, hs_bobwtd_hab1, hs_bobwtd_hunt0, hs_bobwtd_hunt1, hs_bobwtd_pub0, hs_bobwtd_pub1, hs_bobwtd_pub2)
  #' Model selection
  modSel(hs_bobwtd_fl)
  summary(hs_bobwtd_pub1) # f12 not significant
  summary(hs_bobwtd_pub0)
  
  
  #' Save model outputs in one giant R image
  save.image(file = "./Outputs/MultiSpp_CoOcc_Models.RData")
  
  
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
  occ_cougwtd_grazing <- occ_out(gs_cougwtd_graze0, "Cougar", "White-tailed Deer", "Grazing")
  occ_cougmoose_grazing <- occ_out(gs_cougmoose_hab0, "Cougar", "Moose", "Grazing")
  occ_wolfmd_grazing <- occ_out(gs_wolfmd_hab0, "Wolf", "Mule Deer", "Grazing")
  occ_wolfelk_grazing <- occ_out(gs_wolfelk_hab0, "Wolf", "Elk", "Grazing")
  occ_wolfwtd_grazing <- occ_out(gs_wolfwtd_graze0, "Wolf", "White-tailed Deer", "Grazing")
  occ_wolfmoose_grazing <- occ_out(gs_wolfmoose_graze0, "Wolf", "Moose", "Grazing")
  occ_bearmd_grazing <- occ_out(gs_bearmd_hab0, "Black Bear", "Mule Deer", "Grazing")
  occ_bearelk_grazing <- occ_out(gs_bearelk_hab0, "Black Bear", "Elk", "Grazing")
  occ_bearwtd_grazing <- occ_out(gs_bearwtd_graze0, "Black Bear", "White-tailed Deer", "Grazing")
  occ_bearmoose_grazing <- occ_out(gs_bearmoose_graze0, "Black Bear", "Moose", "Grazing")
  occ_coymd_grazing <- occ_out(gs_coymd_graze2, "Coyote", "Mule Deer", "Grazing")
  occ_coyelk_grazing <- occ_out(gs_coywtd_graze2, "Coyote", "White-tailed Deer", "Grazing")
  occ_bobwtd_grazing <- occ_out(gs_bobmd_graze0, "Bobcat", "Mule Deer Deer", "Grazing")
  occ_bobmoose_grazing <- occ_out(gs_bobwtd_graze0, "Bobcat", "White-tailed Deer", "Grazing")
  #'  Hunting season models
  occ_cougmd_hunting <- occ_out(hs_cougmd_pub2, "Cougar", "Mule Deer", "Hunting")
  occ_cougelk_hunting <- occ_out(hs_cougelk_pub0, "Cougar", "Elk", "Hunting")
  occ_cougwtd_hunting <- occ_out(hs_cougwtd_pub2, "Cougar", "White-tailed Deer", "Hunting")
  occ_cougmoose_hunting <- occ_out(hs_cougmoose_hunt0, "Cougar", "Moose", "Hunting")
  occ_wolfmd_hunting <- occ_out(hs_wolfmd_hunt1, "Wolf", "Mule Deer", "Hunting")
  occ_wolfelk_hunting <- occ_out(hs_wolfelk_hunt0, "Wolf", "Elk", "Hunting")
  occ_wolfwtd_hunting <- occ_out(hs_wolfwtd_hunt0, "Wolf", "White-tailed Deer", "Hunting")
  occ_wolfmoose_hunting <- occ_out(hs_wolfmoose_hunt2, "Wolf", "Moose", "Hunting")
  occ_bearmd_hunting <- occ_out(hs_bearmd_hab0, "Black Bear", "Mule Deer", "Hunting")
  occ_bearelk_hunting <- occ_out(hs_bearelk_pub1, "Black Bear", "Elk", "Hunting")
  occ_bearwtd_hunting <- occ_out(hs_bearwtd_pub0, "Black Bear", "White-tailed Deer", "Hunting")
  occ_bearmoose_hunting <- occ_out(hs_bearmoose_hunt0, "Black Bear", "Moose", "Hunting")
  occ_coymd_hunting <- occ_out(hs_coymd_hunt0, "Coyote", "Mule Deer", "Hunting")
  occ_coyelk_hunting <- occ_out(hs_coywtd_hunt1, "Coyote", "White-tailed Deer", "Hunting")
  occ_bobwtd_hunting <- occ_out(hs_bobmd_hunt0, "Bobcat", "Mule Deer Deer", "Hunting")
  occ_bobmoose_hunting <- occ_out(hs_bobwtd_pub0, "Bobcat", "White-tailed Deer", "Hunting")
  
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
  det_cougwtd_grazing <- det_out(gs_cougwtd_graze0, "Cougar", "White-tailed Deer", "Grazing")
  det_cougmoose_grazing <- det_out(gs_cougmoose_hab0, "Cougar", "Moose", "Grazing")
  det_wolfmd_grazing <- det_out(gs_wolfmd_hab0, "Wolf", "Mule Deer", "Grazing")
  det_wolfelk_grazing <- det_out(gs_wolfelk_hab0, "Wolf", "Elk", "Grazing")
  det_wolfwtd_grazing <- det_out(gs_wolfwtd_graze0, "Wolf", "White-tailed Deer", "Grazing")
  det_wolfmoose_grazing <- det_out(gs_wolfmoose_graze0, "Wolf", "Moose", "Grazing")
  det_bearmd_grazing <- det_out(gs_bearmd_hab0, "Black Bear", "Mule Deer", "Grazing")
  det_bearelk_grazing <- det_out(gs_bearelk_hab0, "Black Bear", "Elk", "Grazing")
  det_bearwtd_grazing <- det_out(gs_bearwtd_graze0, "Black Bear", "White-tailed Deer", "Grazing")
  det_bearmoose_grazing <- det_out(gs_bearmoose_graze0, "Black Bear", "Moose", "Grazing")
  det_coymd_grazing <- det_out(gs_coymd_graze2, "Coyote", "Mule Deer", "Grazing")
  det_coyelk_grazing <- det_out(gs_coywtd_graze2, "Coyote", "White-tailed Deer", "Grazing")
  det_bobwtd_grazing <- det_out(gs_bobmd_graze0, "Bobcat", "Mule Deer Deer", "Grazing")
  det_bobmoose_grazing <- det_out(gs_bobwtd_graze0, "Bobcat", "White-tailed Deer", "Grazing")
  #'  Hunting season models
  det_cougmd_hunting <- det_out(hs_cougmd_pub2, "Cougar", "Mule Deer", "Hunting")
  det_cougelk_hunting <- det_out(hs_cougelk_pub0, "Cougar", "Elk", "Hunting")
  det_cougwtd_hunting <- det_out(hs_cougwtd_pub2, "Cougar", "White-tailed Deer", "Hunting")
  det_cougmoose_hunting <- det_out(hs_cougmoose_hunt0, "Cougar", "Moose", "Hunting")
  det_wolfmd_hunting <- det_out(hs_wolfmd_hunt1, "Wolf", "Mule Deer", "Hunting")
  det_wolfelk_hunting <- det_out(hs_wolfelk_hunt0, "Wolf", "Elk", "Hunting")
  det_wolfwtd_hunting <- det_out(hs_wolfwtd_hunt0, "Wolf", "White-tailed Deer", "Hunting")
  det_wolfmoose_hunting <- det_out(hs_wolfmoose_hunt2, "Wolf", "Moose", "Hunting")
  det_bearmd_hunting <- det_out(hs_bearmd_hab0, "Black Bear", "Mule Deer", "Hunting")
  det_bearelk_hunting <- det_out(hs_bearelk_pub1, "Black Bear", "Elk", "Hunting")
  det_bearwtd_hunting <- det_out(hs_bearwtd_pub0, "Black Bear", "White-tailed Deer", "Hunting")
  det_bearmoose_hunting <- det_out(hs_bearmoose_hunt0, "Black Bear", "Moose", "Hunting")
  det_coymd_hunting <- det_out(hs_coymd_hunt0, "Coyote", "Mule Deer", "Hunting")
  det_coyelk_hunting <- det_out(hs_coywtd_hunt1, "Coyote", "White-tailed Deer", "Hunting")
  det_bobwtd_hunting <- det_out(hs_bobmd_hunt0, "Bobcat", "Mule Deer Deer", "Hunting")
  det_bobmoose_hunting <- det_out(hs_bobwtd_pub0, "Bobcat", "White-tailed Deer", "Hunting")

  
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

  