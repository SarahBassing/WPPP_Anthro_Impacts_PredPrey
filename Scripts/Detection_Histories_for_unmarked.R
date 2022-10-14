  #'  ============================================
  #'  Camera Trap Detection Histories
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2022
  #'  ============================================
  #'  Script to combine species detection data with camera station data (mainly,
  #'  camera operation dates and problems time periods). This uses the camtrapR 
  #'  package to generate species-specific encounter histories that can be used
  #'  for occupancy models.
  #'  
  #'  Combines: 
  #'  "full_camdata_DATE.csv" from Detections_by_Camera_Station.R in 
  #'      WPPP_CameraTrapping.Rproj
  #'     -Contains ALL detections of animals, humans, & vehicles (no empties),
  #'      camera coordinates and deployment covariates (cam height, etc.)
  #'  "All_Camera_Stations_18-21_updated_DATE.csv"
  #'     -Contains camera locations (including updated names/locations when a
  #'     camera was moved), deployment & pull dates, and problem dates
  #'      
  #'  
  #'  1. Drop unnecessary fields and format dates and times correctly
  #'  2. Truncate detection data to desired date range
  #'  3. Create camera operation table (deals with cameras that were temporarily
  #'     inoperable during desired date range)
  #'  4. Create species-specific detection probabilities
  #'  
  #'  Acknowledgments: Script is based on code originally provided by Mitch 
  #'  Parsons & Michael Havrda, UW SEFS.
  #'  ============================================
  
  #' #'  Clean workspace & load libraries
  #' rm(list = ls())
  
  library(camtrapR)
  library(chron)
  library(lubridate)
  library(tidyverse)
  
  #'  Read in data, format, and filter
  #'  Heads up: DON'T format the dates in the cam_stations file yet!
  cam_stationsYr1 <- read.csv("./Data/All_Camera_Stations_18-19_updated_1.21.21.csv")
  cam_stationsYr2 <- read.csv("./Data/All_Camera_Stations_19-20.csv")
  cam_stationsYr3 <- read.csv("./Data/All_Camera_Stations_20-21.csv")
  cam_stations <- rbind(cam_stationsYr1, cam_stationsYr2, cam_stationsYr3)

  #'  Count number of camera stations per year and study area
  ncam_yr1 <- nrow(cam_stationsYr1)
  nNE_yr1 <- nrow(filter(cam_stationsYr1, grepl("NE", CameraLocation)))
  nOK_yr1 <- nrow(filter(cam_stationsYr1, grepl("OK", CameraLocation)))
  ncam_yr2 <- nrow(cam_stationsYr2)
  nNE_yr2 <- nrow(filter(cam_stationsYr2, grepl("NE", CameraLocation)))
  nOK_yr2 <- nrow(filter(cam_stationsYr2, grepl("OK", CameraLocation)))
  ncam_yr3 <- nrow(cam_stationsYr3)
  nNE_yr3 <- nrow(filter(cam_stationsYr3, grepl("NE", CameraLocation)))
  nOK_yr3 <- nrow(filter(cam_stationsYr3, grepl("OK", CameraLocation)))
   
  megadata <- read.csv("./Data/full_camdata18-21_2022-04-27.csv") %>%  #2022-04-14
    dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation", 
                  "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle", 
                  "Species", "HumanActivity", "Count") %>%
    filter(!grepl("Moultrie", CameraLocation)) %>%
    #  Need to have something in the Species column for each detection
    mutate(
      Species = ifelse(Human == "TRUE" | Human == "true", "Human", Species),
      Species = ifelse(Vehicle == "TRUE" | Vehicle == "true", "Vehicle", Species),
      Species = ifelse(Species == "", "NA", Species),
      HumanActivity = ifelse(HumanActivity == "", "NA", HumanActivity)
    ) %>%
    #  Remove rows where no detection occurred but snuck into this data set somehow
    filter(!(Animal == "FALSE" & Human == "FALSE" & Vehicle == "FALSE") | (Animal == "false" & Human == "false" & Vehicle == "false")) %>%
    #'  Remove observations of humans
    filter(!is.na(Species)) %>%
    # filter(Species != "Cattle") %>%
    mutate(
      DateTime = as.POSIXct(DateTime,
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles"),
      Date = as.Date(Date, format = "%Y-%m-%d"),
      Time = chron(times = Time)
    )
  
  #'  Total number of images
  npix <- nrow(megadata)
  

  ####  Extract independent detections for wildlife species  ####
  #'  Create a column identifying whether each image is an "independent" event
  #'  If camera site is diff from previous row then give unique value. If not then...
  #'  If species detected is diff from previous row at same site then give unique value. If not then...
  #'  If DateTime is >30 min from previous DateTime at same site for same species then give unique value. If not then...
  #'  Capture value is the same as that in the previous row.
  dat <- arrange(megadata, CameraLocation, DateTime)
  caps <- c()
  caps[1] <- 1
  for (i in 2:nrow(dat)){
    if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps[i] = i
    else (if (dat$Species[i-1] != dat$Species[i]) caps[i] = i
          else (if (difftime(dat$DateTime[i], dat$DateTime[i-1], units = c("mins")) > 30) caps[i] = i
                else caps[i] = caps[i-1]))
  }

  caps <- as.factor(caps)

  #'  Add new column to larger data set
  capdata <- cbind(as.data.frame(dat), caps)

  #'  Retain only the first image from each unique detection event
  detections <- capdata %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  #'  Keep in mind the 30 minute gap too large for cattle & human detections
  #'  Use Cattle_Hunter_Detections.R to generate DH for these
  
  
  #'  Pick date ranges for livestock-focused analyses & hunter-focused analyses
  #'  Livestock activity
  cows <- filter(detections, Species == "Cattle") %>%
    dplyr::select(Date, Time, CameraLocation) %>%
    mutate(Month = month(Date)) 
  hist(cows$Month); summary(cows)
  #'  Detections range March through Dec, peaking between July and Sept. 
  #'  July through Sept core activity on cameras based on 1st - 3rd quartiles
  
  #'  Human/Hunter activity
  rifle <- filter(detections, HumanActivity == "Hunter Rifle") %>%
    dplyr::select(Date, Time, CameraLocation, Count) %>%
    mutate(Month = month(Date))
  hist(rifle$Month); summary(rifle)
  #'  Oct through Nov core activity on cameras; Massive peak in October
  
  bow <- filter(detections, HumanActivity == "Hunter Bow") %>%
    dplyr::select(Date, Time, CameraLocation, Count) %>%
    mutate(Month = month(Date))
  hist(bow$Month); summary(bow)
  #'  Peak in Sept, very few Oct through Dec
  
  hiker <- filter(detections, HumanActivity == "Hiker") %>%
    dplyr::select(Date, Time, CameraLocation, Count) %>%
    mutate(Month = month(Date))
  hist(hiker$Month); summary(hiker)
  #'  June through Sept core activity on camera based on quartiles
  
  truck <- filter(detections, HumanActivity == "Vehicle Truck Car") %>%
    dplyr::select(Date, Time, CameraLocation, Count) %>%
    mutate(Month = month(Date))
  hist(truck$Month); summary(truck)
  #'  June through Oct core activity based on quartiles, peak in Oct
  
  atv <- filter(detections, HumanActivity == "Vehicle ATV") %>%
    dplyr::select(Date, Time, CameraLocation, Count) %>%
    mutate(Month = month(Date))
  hist(atv$Month); summary(atv)
  #'  June through Oct core activity based on quartiles, peak in Oct
  
  bike <- filter(detections, HumanActivity == "Dirt Bike") %>%
    dplyr::select(Date, Time, CameraLocation, Count) %>%
    mutate(Month = month(Date))
  hist(bike$Month); summary(bike)
  #'  May through Aug core activity based on quartiles, peak in May
  
  #'  Livestock analyses: July through September
  #'  Hunter analyses: October through November
  #'  NE general rifle seasons: 
  #'       2018: 10/13/18 - 10/26/18 and 11/10/18 - 11/19/18
  #'       2019: 10/12/19 - 10/25/19 and 11/9/19 - 11/19/19
  #'       2020: 10/17/20 - 10/30/20 and 11/7/20 - 11/19/20
  #'  OK general rifle seasons: 
  #'       2018: 10/13/18 - 10/23/18 
  #'       2019: 10/12/19 - 10/22/19
  #'       2020: 10/17/20 - 10/27/20
  
  
  #'  Filter dates to specific range 
  #'  Grazing Season 2018: 07/01/2018 - 09/29/2018 (thirteen 7-day sampling periods)
  images_graze2018 <- detections %>%
    filter(Date > "2018-06-30") %>%
    filter(Date < "2018-09-30") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species")  
  #'  Subset by study area
  NE_graze18 <- filter(images_graze2018, grepl("NE", CameraLocation))
  OK_graze18 <- filter(images_graze2018, grepl("OK", CameraLocation))
  #'  Hunting Season 2018: 10/1/2018 - 11/25/2018 (eight 7-day sampling periods)
  images_hunt2018 <- detections %>%
    filter(Date > "2018-09-30") %>%
    filter(Date < "2018-11-26") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species") 
  #'  Subset by study area
  NE_hunt18 <- filter(images_hunt2018, grepl("NE", CameraLocation))
  OK_hunt18 <- filter(images_hunt2018, grepl("OK", CameraLocation))
  #'  Grazing Season 2019: 07/01/2019 - 09/29/2019 (thirteen 7-day sampling periods)
  images_graze2019 <- detections %>%
    filter(Date > "2019-06-30") %>%
    filter(Date < "2019-09-30") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species")  
  #'  Subset by study area
  NE_graze19 <- filter(images_graze2019, grepl("NE", CameraLocation))
  OK_graze19 <- filter(images_graze2019, grepl("OK", CameraLocation))
  #'  Hunting Season 2019: 10/1/2019 - 11/25/2019 (eight 7-day sampling periods)
  images_hunt2019 <- detections %>%
    filter(Date > "2019-09-30") %>%
    filter(Date < "2019-11-26") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species")  
  #'  Subset by study area
  NE_hunt19 <- filter(images_hunt2019, grepl("NE", CameraLocation))
  OK_hunt19 <- filter(images_hunt2019, grepl("OK", CameraLocation))
  #'  Grazing Season 2020: 07/01/2020 - 09/29/2020 (thirteen 7-day sampling periods)
  images_graze2020 <- detections %>%
    filter(Date > "2020-06-30") %>%
    filter(Date < "2020-09-30") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species")  
  #'  Subset by study area
  NE_graze20 <- filter(images_graze2020, grepl("NE", CameraLocation))
  OK_graze20 <- filter(images_graze2020, grepl("OK", CameraLocation))
  #'  Hunting Season 2020: 10/1/2020 - 11/25/2020 (eight 7-day sampling periods)
  images_hunt2020 <- detections %>%
    filter(Date > "2020-09-30") %>%
    filter(Date < "2020-11-26") %>%
    dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "Species")  
  #'  Subset by study area
  NE_hunt20 <- filter(images_hunt2020, grepl("NE", CameraLocation))
  OK_hunt20 <- filter(images_hunt2020, grepl("OK", CameraLocation))
  
  
  #'  Double check these were truncated correctly
  min(images_graze2018$Date); max(images_graze2018$Date)
  min(images_hunt2018$Date); max(images_hunt2018$Date)
  
  min(images_graze2019$Date); max(images_graze2019$Date)
  min(images_hunt2019$Date); max(images_hunt2019$Date)
  
  min(images_graze2020$Date); max(images_graze2020$Date)
  min(images_hunt2020$Date); max(images_hunt2020$Date)
  
  
  #'  Save for data summary and making data available for publication
  CamTrap_Detections <- as.data.frame(rbind(images_graze2018, images_hunt2018, 
                                            images_graze2019, images_hunt2019,
                                            images_graze2020, images_hunt2020)) %>%
    mutate(
      Season = ifelse(Date < "2018-09-30", "Grazing18", "Grazing19"),
      Season = ifelse(Date > "2020-01-01" & Date < "2020-09-30", "Grazing20", Season),
      Season = ifelse(Date > "2018-09-30" & Date < "2018-11-30", "Hunt18", Season),
      Season = ifelse(Date > "2019-09-30" & Date < "2019-11-30", "Hunt19", Season),
      Season = ifelse(Date > "2020-09-30", "Hunt20", Season)
    ) %>%
    filter(Species == "Mule Deer" | Species == "Elk" | Species == "White-tailed Deer" | 
             Species == "Cougar" | Species == "Wolf" | Species == "Black Bear" | 
             Species == "Bobcat" | Species == "Coyote")
  
  # write.csv(CamTrap_Detections, paste0("./Outputs/CamTrap_Detections_", Sys.Date(), ".csv"))
  
  
  #'  Calculate number of trap nights per camera
  #'  Length of deployment
  trapnights_full <- as.numeric(as.Date(cam_stations$Retrieval_date, format = "%m/%d/%Y") - as.Date(cam_stations$Setup_date, format = "%m/%d/%Y"))
  hist(trapnights_full)
  #'  Number of camera trap nights (across all cameras)
  sum(trapnights_full)
  #'  Average number trap nights/camera
  sum(trapnights_full)/length(trapnights_full)

  
  ####  Camera Operation Table  ####
  #'  ------------------------------ 
  #'  Creates a matrix with each camera & dates it deployed
  #'  Define date format at this step, not before!
  #'  1 = operating; 0 = not operating but deployed; NA = not deployed
  camop_problem <- cameraOperation(CTtable = cam_stations,
                                   stationCol = "CameraLocation",
                                   setupCol = "Setup_date",
                                   retrievalCol = "Retrieval_date",
                                   hasProblems = TRUE,
                                   dateFormat = "%m/%d/%Y", # Define date format here!
                                   writecsv = FALSE) 

  probs <- as.data.frame(camop_problem)
  
  #'  Filter data by study area
  NEcams <- filter(cam_stations, grepl("NE", CameraLocation))
  camop_problem_NE <- cameraOperation(CTtable = NEcams,
                                      stationCol = "CameraLocation",
                                      setupCol = "Setup_date",
                                      retrievalCol = "Retrieval_date",
                                      hasProblems = TRUE,
                                      dateFormat = "%m/%d/%Y", 
                                      writecsv = FALSE) 

  OKcams <- filter(cam_stations, grepl("OK", CameraLocation))
  camop_problem_OK <- cameraOperation(CTtable = OKcams,
                                      stationCol = "CameraLocation",
                                      setupCol = "Setup_date",
                                      retrievalCol = "Retrieval_date",
                                      hasProblems = TRUE,
                                      dateFormat = "%m/%d/%Y", 
                                      writecsv = FALSE) 
  
  
  ####  Detection Histories  ####
  #'  ---------------------------
  #'  Function to create season-specific detection histories for each species
  #'  for each season of interest and create a sampling effort matrix/season.
  #'  
  #'  Key arguments:
  #'  -occasionLength: currently using 7 day sampling occasions
  #'  -day1: sampling occasion begins upon station setup date ("station"), 
  #'   first day of survey ("survey"), or a specific date ("2018-06-15")
  #'   FYI this defines start date but NOT end date so DH goes until camera pulls
  #'  -includeEffort: compute # active trap days/station/occasion- effects DH
  #'   if FALSE then when camera is not set up or malfunctioning (NA or 0 in
  #'   camOp) during part of sampling occasion DH will be 0 for that occasion
  #'  -scaleEffort: center & scale effort (I plan to do this later in model)
  #'  -dateAsOccasionNames: if day1 = "survey" then uses 1st and last day of 
  #'   occasion as name for each sampling occasion
  #'  -output: return binary detections or counts of detections; don't want to
  #'   use "count" right now b/c would count each image not independent events
  #'  -occasionStartTime: time of day (full hour) at which to begin occasions
  #'   default is midnight (0)
  #'  -unmarkedMultFrameInput: create input for multi-season occmod in unmarked
  #'   ONLY use if running multi-season models & need to add more info to camop
  #'   
  #'  FYI, cannot have any NAs in the Species column or this doesn't work
  #'  Need to remove columns that extend beyond date range of interest!
  #'  Grazing: July 1 - Sept 29 = 13 weeks
  #'  Hunting: Oct 1 - Nov 26 = 8 weeks 
  
  DH <- function(images, spp, start_date) {
    det_hist <- detectionHistory(recordTable = images,
                                 camOp = camop_problem,
                                 stationCol = "CameraLocation",
                                 speciesCol = "Species",
                                 recordDateTimeCol = "DateTime",
                                 recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                 species = spp,
                                 occasionLength = 7,
                                 day1 = start_date, 
                                 # datesAsOccasionNames = TRUE,
                                 # occasionStartTime = 12, # starts at noon
                                 timeZone = "America/Los_Angeles",
                                 output = "binary",
                                 includeEffort = TRUE,
                                 scaleEffort = FALSE,
                                 # writecsv = TRUE,
                                 outDir = "./Data/Detection_Histories")
    
    return(det_hist)
  }
  #'  Create species and season-specific detection histories
  #'  Year1 cams 1:125 (n=125), Year2 cams 126:242 (n=117), Year3 cams 243:361 (n=119)
  #'  Grazing season has 13 occasions, hunting season has 8 occasions (1 week/occasion)
  
  ####  BOBCATS  ####
  #'  Grazing season
  bob_graze18 <- DH(images_graze2018, "Bobcat", "2018-07-01")
  DH_bob_graze18 <- bob_graze18[[1]][1:125,1:13]   
  bob_graze19 <- DH(images_graze2019, "Bobcat", "2019-07-01")
  DH_bob_graze19 <- bob_graze19[[1]][126:242,1:13] 
  bob_graze20 <- DH(images_graze2020, "Bobcat", "2020-07-01")
  DH_bob_graze20 <- bob_graze20[[1]][243:361,1:13] 
  
  DH_bob_graze1820 <- rbind(DH_bob_graze18, DH_bob_graze19, DH_bob_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  #'  All detections are missing at sites: 10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336
  DH_bob_graze1820 <- DH_bob_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_bob_graze1820_NE <- DH_bob_graze1820[grepl("NE", row.names(DH_bob_graze1820)),]
  DH_bob_graze1820_OK <- DH_bob_graze1820[grepl("OK", row.names(DH_bob_graze1820)),]
  
  #'  Hunting season
  bob_hunt18 <- DH(images_hunt2018, "Bobcat", "2018-10-01")
  DH_bob_hunt18 <- bob_hunt18[[1]][1:125,1:8]  
  bob_hunt19 <- DH(images_hunt2019, "Bobcat", "2019-10-01")
  DH_bob_hunt19 <- bob_hunt19[[1]][126:242,1:8] 
  bob_hunt20 <- DH(images_hunt2020, "Bobcat", "2020-10-01")
  DH_bob_hunt20 <- bob_hunt20[[1]][243:361,1:8]
  
  DH_bob_hunt1820 <- rbind(DH_bob_hunt18, DH_bob_hunt19, DH_bob_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  #'  All detections are missing at sites: 16, 23, 25, 27, 29, 38, 85, 111, 119, 
  #'  128, 129, 144, 146, 156, 162, 216, 274, 282, 283
  DH_bob_hunt1820 <- DH_bob_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_bob_hunt1820_NE <- DH_bob_hunt1820[grepl("NE", row.names(DH_bob_hunt1820)),]
  DH_bob_hunt1820_OK <- DH_bob_hunt1820[grepl("OK", row.names(DH_bob_hunt1820)),]
  
  ####  BLACK BEARS  ####
  #'  Grazing season
  bear_graze18 <- DH(images_graze2018, "Black Bear", "2018-07-01")
  DH_bear_graze18 <- bear_graze18[[1]][1:125,1:13]   
  bear_graze19 <- DH(images_graze2019, "Black Bear", "2019-07-01")
  DH_bear_graze19 <- bear_graze19[[1]][126:242,1:13] 
  bear_graze20 <- DH(images_graze2020, "Black Bear", "2020-07-01")
  DH_bear_graze20 <- bear_graze20[[1]][243:361,1:13] 
  
  DH_bear_graze1820 <- rbind(DH_bear_graze18, DH_bear_graze19, DH_bear_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_bear_graze1820 <- DH_bear_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_bear_graze1820_NE <- DH_bear_graze1820[grepl("NE", row.names(DH_bear_graze1820)),]
  DH_bear_graze1820_OK <- DH_bear_graze1820[grepl("OK", row.names(DH_bear_graze1820)),]
  
  #'  Hunting season
  bear_hunt18 <- DH(images_hunt2018, "Black Bear", "2018-10-01")
  DH_bear_hunt18 <- bear_hunt18[[1]][1:125,1:8]  
  bear_hunt19 <- DH(images_hunt2019, "Black Bear", "2019-10-01")
  DH_bear_hunt19 <- bear_hunt19[[1]][126:242,1:8] 
  bear_hunt20 <- DH(images_hunt2020, "Black Bear", "2020-10-01")
  DH_bear_hunt20 <- bear_hunt20[[1]][243:361,1:8]
  
  DH_bear_hunt1820 <- rbind(DH_bear_hunt18, DH_bear_hunt19, DH_bear_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_bear_hunt1820 <- DH_bear_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_bear_hunt1820_NE <- DH_bear_hunt1820[grepl("NE", row.names(DH_bear_hunt1820)),]
  DH_bear_hunt1820_OK <- DH_bear_hunt1820[grepl("OK", row.names(DH_bear_hunt1820)),]
  
  ####  COUGARS  ####
  #'  Grazing season
  coug_graze18 <- DH(images_graze2018, "Cougar", "2018-07-01")
  DH_coug_graze18 <- coug_graze18[[1]][1:125,1:13]   
  coug_graze19 <- DH(images_graze2019, "Cougar", "2019-07-01")
  DH_coug_graze19 <- coug_graze19[[1]][126:242,1:13] 
  coug_graze20 <- DH(images_graze2020, "Cougar", "2020-07-01")
  DH_coug_graze20 <- coug_graze20[[1]][243:361,1:13] 
  
  DH_coug_graze1820 <- rbind(DH_coug_graze18, DH_coug_graze19, DH_coug_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_coug_graze1820 <- DH_coug_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_coug_graze1820_NE <- DH_coug_graze1820[grepl("NE", row.names(DH_coug_graze1820)),]
  DH_coug_graze1820_OK <- DH_coug_graze1820[grepl("OK", row.names(DH_coug_graze1820)),]
  
  #'  Hunting season
  coug_hunt18 <- DH(images_hunt2018, "Cougar", "2018-10-01")
  DH_coug_hunt18 <- coug_hunt18[[1]][1:125,1:8]  
  coug_hunt19 <- DH(images_hunt2019, "Cougar", "2019-10-01")
  DH_coug_hunt19 <- coug_hunt19[[1]][126:242,1:8] 
  coug_hunt20 <- DH(images_hunt2020, "Cougar", "2020-10-01")
  DH_coug_hunt20 <- coug_hunt20[[1]][243:361,1:8]
  
  DH_coug_hunt1820 <- rbind(DH_coug_hunt18, DH_coug_hunt19, DH_coug_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_coug_hunt1820 <- DH_coug_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_coug_hunt1820_NE <- DH_coug_hunt1820[grepl("NE", row.names(DH_coug_hunt1820)),]
  DH_coug_hunt1820_OK <- DH_coug_hunt1820[grepl("OK", row.names(DH_coug_hunt1820)),]
  
  ####  COYOTES  ####
  #'  Grazing season
  coy_graze18 <- DH(images_graze2018, "Coyote", "2018-07-01")
  DH_coy_graze18 <- coy_graze18[[1]][1:125,1:13]   
  coy_graze19 <- DH(images_graze2019, "Coyote", "2019-07-01")
  DH_coy_graze19 <- coy_graze19[[1]][126:242,1:13] 
  coy_graze20 <- DH(images_graze2020, "Coyote", "2020-07-01")
  DH_coy_graze20 <- coy_graze20[[1]][243:361,1:13] 
  
  DH_coy_graze1820 <- rbind(DH_coy_graze18, DH_coy_graze19, DH_coy_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_coy_graze1820 <- DH_coy_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_coy_graze1820_NE <- DH_coy_graze1820[grepl("NE", row.names(DH_coy_graze1820)),]
  DH_coy_graze1820_OK <- DH_coy_graze1820[grepl("OK", row.names(DH_coy_graze1820)),]
  
  #'  Hunting season
  coy_hunt18 <- DH(images_hunt2018, "Coyote", "2018-10-01")
  DH_coy_hunt18 <- coy_hunt18[[1]][1:125,1:8]  
  coy_hunt19 <- DH(images_hunt2019, "Coyote", "2019-10-01")
  DH_coy_hunt19 <- coy_hunt19[[1]][126:242,1:8] 
  coy_hunt20 <- DH(images_hunt2020, "Coyote", "2020-10-01")
  DH_coy_hunt20 <- coy_hunt20[[1]][243:361,1:8]
  
  DH_coy_hunt1820 <- rbind(DH_coy_hunt18, DH_coy_hunt19, DH_coy_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_coy_hunt1820 <- DH_coy_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_coy_hunt1820_NE <- DH_coy_hunt1820[grepl("NE", row.names(DH_coy_hunt1820)),]
  DH_coy_hunt1820_OK <- DH_coy_hunt1820[grepl("OK", row.names(DH_coy_hunt1820)),]
  
  ####  ELK  ####
  #'  Grazing season
  elk_graze18 <- DH(images_graze2018, "Elk", "2018-07-01")
  DH_elk_graze18 <- elk_graze18[[1]][1:125,1:13]   
  elk_graze19 <- DH(images_graze2019, "Elk", "2019-07-01")
  DH_elk_graze19 <- elk_graze19[[1]][126:242,1:13] 
  elk_graze20 <- DH(images_graze2020, "Elk", "2020-07-01")
  DH_elk_graze20 <- elk_graze20[[1]][243:361,1:13] 
  
  DH_elk_graze1820 <- rbind(DH_elk_graze18, DH_elk_graze19, DH_elk_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_elk_graze1820 <- DH_elk_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_elk_graze1820_NE <- DH_elk_graze1820[grepl("NE", row.names(DH_elk_graze1820)),]
  DH_elk_graze1820_OK <- DH_elk_graze1820[grepl("OK", row.names(DH_elk_graze1820)),]
  
  #'  Hunting season
  elk_hunt18 <- DH(images_hunt2018, "Elk", "2018-10-01")
  DH_elk_hunt18 <- elk_hunt18[[1]][1:125,1:8]  
  elk_hunt19 <- DH(images_hunt2019, "Elk", "2019-10-01")
  DH_elk_hunt19 <- elk_hunt19[[1]][126:242,1:8] 
  elk_hunt20 <- DH(images_hunt2020, "Elk", "2020-10-01")
  DH_elk_hunt20 <- elk_hunt20[[1]][243:361,1:8]
  
  DH_elk_hunt1820 <- rbind(DH_elk_hunt18, DH_elk_hunt19, DH_elk_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_elk_hunt1820 <- DH_elk_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_elk_hunt1820_NE <- DH_elk_hunt1820[grepl("NE", row.names(DH_elk_hunt1820)),]
  DH_elk_hunt1820_OK <- DH_elk_hunt1820[grepl("OK", row.names(DH_elk_hunt1820)),]
  
  ####  MOOSE  ####
  #'  Grazing season
  moose_graze18 <- DH(images_graze2018, "Moose", "2018-07-01")
  DH_moose_graze18 <- moose_graze18[[1]][1:125,1:13]   
  moose_graze19 <- DH(images_graze2019, "Moose", "2019-07-01")
  DH_moose_graze19 <- moose_graze19[[1]][126:242,1:13] 
  moose_graze20 <- DH(images_graze2020, "Moose", "2020-07-01")
  DH_moose_graze20 <- moose_graze20[[1]][243:361,1:13] 
  
  DH_moose_graze1820 <- rbind(DH_moose_graze18, DH_moose_graze19, DH_moose_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_moose_graze1820 <- DH_moose_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_moose_graze1820_NE <- DH_moose_graze1820[grepl("NE", row.names(DH_moose_graze1820)),]
  DH_moose_graze1820_OK <- DH_moose_graze1820[grepl("OK", row.names(DH_moose_graze1820)),]
  
  #'  Hunting season
  moose_hunt18 <- DH(images_hunt2018, "Moose", "2018-10-01")
  DH_moose_hunt18 <- moose_hunt18[[1]][1:125,1:8]  
  moose_hunt19 <- DH(images_hunt2019, "Moose", "2019-10-01")
  DH_moose_hunt19 <- moose_hunt19[[1]][126:242,1:8] 
  moose_hunt20 <- DH(images_hunt2020, "Moose", "2020-10-01")
  DH_moose_hunt20 <- moose_hunt20[[1]][243:361,1:8]
  
  DH_moose_hunt1820 <- rbind(DH_moose_hunt18, DH_moose_hunt19, DH_moose_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_moose_hunt1820 <- DH_moose_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_moose_hunt1820_NE <- DH_moose_hunt1820[grepl("NE", row.names(DH_moose_hunt1820)),]
  DH_moose_hunt1820_OK <- DH_moose_hunt1820[grepl("OK", row.names(DH_moose_hunt1820)),]
  
  ####  MULE DEER  ####
  #'  Grazing season
  md_graze18 <- DH(images_graze2018, "Mule Deer", "2018-07-01")
  DH_md_graze18 <- md_graze18[[1]][1:125,1:13]   
  md_graze19 <- DH(images_graze2019, "Mule Deer", "2019-07-01")
  DH_md_graze19 <- md_graze19[[1]][126:242,1:13] 
  md_graze20 <- DH(images_graze2020, "Mule Deer", "2020-07-01")
  DH_md_graze20 <- md_graze20[[1]][243:361,1:13] 
  
  DH_md_graze1820 <- rbind(DH_md_graze18, DH_md_graze19, DH_md_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_md_graze1820 <- DH_md_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_md_graze1820_NE <- DH_md_graze1820[grepl("NE", row.names(DH_md_graze1820)),]
  DH_md_graze1820_OK <- DH_md_graze1820[grepl("OK", row.names(DH_md_graze1820)),]
  
  #'  Hunting season
  md_hunt18 <- DH(images_hunt2018, "Mule Deer", "2018-10-01")
  DH_md_hunt18 <- md_hunt18[[1]][1:125,1:8]  
  md_hunt19 <- DH(images_hunt2019, "Mule Deer", "2019-10-01")
  DH_md_hunt19 <- md_hunt19[[1]][126:242,1:8] 
  md_hunt20 <- DH(images_hunt2020, "Mule Deer", "2020-10-01")
  DH_md_hunt20 <- md_hunt20[[1]][243:361,1:8]
  
  DH_md_hunt1820 <- rbind(DH_md_hunt18, DH_md_hunt19, DH_md_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_md_hunt1820 <- DH_md_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_md_hunt1820_NE <- DH_md_hunt1820[grepl("NE", row.names(DH_md_hunt1820)),]
  DH_md_hunt1820_OK <- DH_md_hunt1820[grepl("OK", row.names(DH_md_hunt1820)),]
  
  ####  WHITE-TAILED DEER  ####
  #'  Grazing season
  wtd_graze18 <- DH(images_graze2018, "White-tailed Deer", "2018-07-01")
  DH_wtd_graze18 <- wtd_graze18[[1]][1:125,1:13]   
  wtd_graze19 <- DH(images_graze2019, "White-tailed Deer", "2019-07-01")
  DH_wtd_graze19 <- wtd_graze19[[1]][126:242,1:13] 
  wtd_graze20 <- DH(images_graze2020, "White-tailed Deer", "2020-07-01")
  DH_wtd_graze20 <- wtd_graze20[[1]][243:361,1:13] 
  
  DH_wtd_graze1820 <- rbind(DH_wtd_graze18, DH_wtd_graze19, DH_wtd_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_wtd_graze1820 <- DH_wtd_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_wtd_graze1820_NE <- DH_wtd_graze1820[grepl("NE", row.names(DH_wtd_graze1820)),]
  DH_wtd_graze1820_OK <- DH_wtd_graze1820[grepl("OK", row.names(DH_wtd_graze1820)),]
  
  #'  Hunting season
  wtd_hunt18 <- DH(images_hunt2018, "White-tailed Deer", "2018-10-01")
  DH_wtd_hunt18 <- wtd_hunt18[[1]][1:125,1:8]  
  wtd_hunt19 <- DH(images_hunt2019, "White-tailed Deer", "2019-10-01")
  DH_wtd_hunt19 <- wtd_hunt19[[1]][126:242,1:8] 
  wtd_hunt20 <- DH(images_hunt2020, "White-tailed Deer", "2020-10-01")
  DH_wtd_hunt20 <- wtd_hunt20[[1]][243:361,1:8]
  
  DH_wtd_hunt1820 <- rbind(DH_wtd_hunt18, DH_wtd_hunt19, DH_wtd_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_wtd_hunt1820 <- DH_wtd_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_wtd_hunt1820_NE <- DH_wtd_hunt1820[grepl("NE", row.names(DH_wtd_hunt1820)),]
  DH_wtd_hunt1820_OK <- DH_wtd_hunt1820[grepl("OK", row.names(DH_wtd_hunt1820)),]
  
  ####  WOLVES  ####
  #'  Grazing season
  wolf_graze18 <- DH(images_graze2018, "Wolf", "2018-07-01")
  DH_wolf_graze18 <- wolf_graze18[[1]][1:125,1:13]   
  wolf_graze19 <- DH(images_graze2019, "Wolf", "2019-07-01")
  DH_wolf_graze19 <- wolf_graze19[[1]][126:242,1:13] 
  wolf_graze20 <- DH(images_graze2020, "Wolf", "2020-07-01")
  DH_wolf_graze20 <- wolf_graze20[[1]][243:361,1:13]
  
  DH_wolf_graze1820 <- rbind(DH_wolf_graze18, DH_wolf_graze19, DH_wolf_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_wolf_graze1820 <- DH_wolf_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  DH_wolf_graze1820_NE <- DH_wolf_graze1820[grepl("NE", row.names(DH_wolf_graze1820)),]
  DH_wolf_graze1820_OK <- DH_wolf_graze1820[grepl("OK", row.names(DH_wolf_graze1820)),]
  
  #'  Hunting season
  wolf_hunt18 <- DH(images_hunt2018, "Wolf", "2018-10-01")
  DH_wolf_hunt18 <- wolf_hunt18[[1]][1:125,1:8]  
  wolf_hunt19 <- DH(images_hunt2019, "Wolf", "2019-10-01")
  DH_wolf_hunt19 <- wolf_hunt19[[1]][126:242,1:8] 
  wolf_hunt20 <- DH(images_hunt2020, "Wolf", "2020-10-01")
  DH_wolf_hunt20 <- wolf_hunt20[[1]][243:361,1:8]
  
  DH_wolf_hunt1820 <- rbind(DH_wolf_hunt18, DH_wolf_hunt19, DH_wolf_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  DH_wolf_hunt1820 <- DH_wolf_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  DH_wolf_hunt1820_NE <- DH_wolf_hunt1820[grepl("NE", row.names(DH_wolf_hunt1820)),]
  DH_wolf_hunt1820_OK <- DH_wolf_hunt1820[grepl("OK", row.names(DH_wolf_hunt1820)),]
  
  ####  SAMPLING EFFORT  ####
  #'  -------------------
  #'  Save sampling effort for each camera and season
  #'  This is the same for all species so only need to do this once
  Effort_graze18 <- bob_graze18[[2]][1:125,1:13]
  Effort_graze19 <- bob_graze19[[2]][126:242,1:13]
  Effort_graze20 <- bob_graze20[[2]][243:361,1:13]
  Effort_hunt18 <- bob_hunt18[[2]][1:125,1:8]
  Effort_hunt19 <- bob_hunt19[[2]][126:242,1:8]
  Effort_hunt20 <- bob_hunt20[[2]][243:361,1:8]
  
  Effort_graze1820 <- rbind(Effort_graze18, Effort_graze19, Effort_graze20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  Effort_graze1820 <- Effort_graze1820[-c(10, 13, 17, 27, 30, 61, 62, 110, 216, 274, 283, 336),]
  Effort_graze1820_NE <- Effort_graze1820[grepl("NE", row.names(Effort_graze1820)),]
  Effort_graze1820_OK <- Effort_graze1820[grepl("OK", row.names(Effort_graze1820)),]
  
  Effort_hunt1820 <- rbind(Effort_hunt18, Effort_hunt19, Effort_hunt20)
  #'  Remove rows missing detection data for ALL occasions (camera completely inoperable)
  Effort_hunt1820 <- Effort_hunt1820[-c(16, 23, 25, 27, 29, 38, 85, 111, 119, 128, 129, 144, 146, 156, 162, 216, 274, 282, 283),]
  Effort_hunt1820_NE <- Effort_hunt1820[grepl("NE", row.names(Effort_hunt1820)),]
  Effort_hunt1820_OK <- Effort_hunt1820[grepl("OK", row.names(Effort_hunt1820)),]
  
  #'  Summary stats on trap nights (active sites only)
  #'  Total number of trap nights (across camera sites)
  trapnights_graze <- sum(Effort_graze1820, na.rm = T) #27359
  trapnights_hunt <- sum(Effort_hunt1820, na.rm = T)   #18432
  #'  Remove rows of all NAs (inactive camera stations)
  eff_graze1820 <- Effort_graze1820[rowSums(is.na(Effort_graze1820)) != ncol(Effort_graze1820), ]
  eff_hunt1820 <- Effort_hunt1820[rowSums(is.na(Effort_hunt1820)) != ncol(Effort_hunt1820), ]
  nrow(eff_graze1820) #349
  nrow(eff_hunt1820)  #342
  #'  Average number of trap nights per sampling occasion & SE
  mean(eff_graze1820, na.rm = TRUE); sd(eff_graze1820, na.rm = TRUE)/length(eff_graze1820)
  mean(eff_hunt1820, na.rm = TRUE); sd(eff_hunt1820, na.rm = TRUE)/length(eff_hunt1820)
  #'  Average number, range, & SE of trap nights per season
  summary(rowSums(eff_graze1820, na.rm = TRUE)); sd(rowSums(eff_graze1820, na.rm = TRUE))/nrow(eff_graze1820) #length(eff_smr1819)
  summary(rowSums(eff_hunt1820, na.rm = TRUE)); sd(rowSums(eff_hunt1820, na.rm = TRUE))/nrow(eff_hunt1820) #length(eff_wtr1820)
  #'  Plot and visualize distribution vs mean (red, dashed) and median (blue, solid)
  hist(rowSums(eff_graze1820, na.rm = TRUE), breaks = 50, xlim = c(0, 100), xlab = "Number of Trap Nights", main = "Number of trap nights across camera sites, \nSummer 2018 & 2019")
  abline(v = mean(rowSums(eff_graze1820, na.rm = TRUE)), col = "red", lty = 2)
  abline(v = median(rowSums(eff_graze1820, na.rm = TRUE)), col = "blue", lty = 1, lwd = 2)
  hist(rowSums(eff_hunt1820, na.rm = TRUE), breaks = 50, xlim = c(0, 100), xlab = "Number of Trap Nights", main = "Number of trap nights across camera sites, \nWinter 2018/2019 & 2019/2020")
  abline(v = mean(rowSums(eff_hunt1820, na.rm = TRUE)), col = "red", lty = 2)
  abline(v = median(rowSums(eff_hunt1820, na.rm = TRUE)), col = "blue", lty = 1, lwd = 2)
  #'  Total number of trap nights 
  sum(rowSums(eff_graze1820, na.rm = TRUE)) 
  sum(rowSums(eff_hunt1820, na.rm = TRUE))
  #'  Total number of active cameras per season and by study area
  CameraLocation <- rownames(eff_graze1820)
  eff_graze1820 <- as.data.frame(cbind(CameraLocation, eff_graze1820))
  nrow(eff_graze1820)
  nrow(eff_graze1820[grepl("NE", CameraLocation),])
  nrow(eff_graze1820[grepl("OK", CameraLocation),])
  CameraLocation <- rownames(eff_hunt1820)
  eff_hunt1820 <- as.data.frame(cbind(CameraLocation, eff_hunt1820))
  nrow(eff_hunt1820)
  nrow(eff_hunt1820[grepl("NE", CameraLocation),])
  nrow(eff_hunt1820[grepl("OK", CameraLocation),])
  
  
  #'  How many sampling occasions had incomplete sampling effort
  #'  Ignoring sampling occasions when camera was inactive
  loweff_graze <- sum(Effort_graze1820 < 7, na.rm = TRUE)
  loweff_hunt <- sum(Effort_hunt1820 < 7, na.rm = TRUE)
  #'  Number of sampling occasions total
  nocc_graze <- length(Effort_graze1820[!is.na(Effort_graze1820)])
  nocc_hunt <- length(Effort_hunt1820[!is.na(Effort_hunt1820)])
  #'  What percent of sampling occasions had incomplete sampling effort?
  #'  If less than 5%, no going to worry about accounting for sampling effort
  #'  in detection process
  loweff_graze/nocc_graze
  loweff_hunt/nocc_hunt
  
  
  ####  Summary Stats  ####
  #'  -----------------
  #'  Number of independent detections per season from camera traps
  ndet <- CamTrap_Detections %>%
    group_by(Season, Species) %>%
    summarise(n = n()) %>%
    ungroup()
  summary_dets <- group_by(ndet, Season) %>% 
    summarize(mu_dets = mean(n), sd = sd(n), se_dets = sd(n)/sqrt(n())) %>% 
    ungroup() 
  mean_dets <- ndet %>%
    mutate(Anthro = ifelse(Season == "Grazing18" | Season == "Grazing19" | Season == "Grazing20", "Grazing", "Hunting")) %>%
    group_by(Anthro) %>%
    summarize(mu_des = mean(n), sd = sd(n), se_dets = sd(n)/sqrt(n())) %>%
    ungroup()
  
  #'  Number of detections per species per season
  ndet_by_spp <- CamTrap_Detections %>%
    mutate(
      Season2 = ifelse(grepl("Grazing", Season), "Grazing", "Hunting")
    ) %>%
    dplyr::select(CameraLocation, Species, Season, Season2) %>%
    group_by(Species, Season2) %>%
    summarise(ndet = n()) %>%
    ungroup()
  colnames(ndet_by_spp) <- c("Species", "Season", "Independent detections")
  
  #'  Percent of cameras where a species was detected
  #'  Based on number of active cameras in each study area and season
  ncams_graze <- nrow(eff_graze1820)
  NEcams_graze <- nrow(eff_graze1820[grepl("NE", eff_graze1820$CameraLocation),])
  OKcams_graze <- nrow(eff_graze1820[grepl("OK", eff_graze1820$CameraLocation),])
  ncams_hunt <- nrow(eff_hunt1820)
  NEcams_hunt <- nrow(eff_hunt1820[grepl("NE", eff_hunt1820$CameraLocation),])
  OKcams_hunt <- nrow(eff_hunt1820[grepl("OK", eff_hunt1820$CameraLocation),])
  
  perc_cams <- CamTrap_Detections %>%
    mutate(
      Season2 = ifelse(grepl("Grazing", Season), "Grazing", "Hunting")
    ) %>%
    dplyr::select(CameraLocation, Species, Season, Season2) %>%
    group_by(Species, Season, CameraLocation) %>%
    filter(row_number(CameraLocation) == 1) %>%
    ungroup() %>%
    group_by(Species, Season2) %>%
    summarise(ncams = n()) %>%
    ungroup() %>%
    mutate(
      propcams = ifelse(Season2 == "Grazing", ncams/ncams_graze, ncams/ncams_hunt),
      propcams = round(propcams, 2)
      #' #'  If I want to only focus on NE cameras for elk
      #' propcams = ifelse(Species == "Elk" & Season2 == "Grazing", ncams/NEcams_graze, propcams),
      #' propcams = ifelse(Species == "Elk" & Season2 == "Hunting", ncams/NEcams_hunt, propcams)
    ) %>%
    dplyr::select(-ncams)
  colnames(perc_cams) <- c("Species", "Season", "Percent of cameras")
  
  detection_summary_data <- full_join(ndet_by_spp, perc_cams, by = c("Species", "Season"))

  #'  Save
  # write.csv(detection_summary_data, file = "./Outputs/Detection_Summary_Table.csv")
  
  
  
  