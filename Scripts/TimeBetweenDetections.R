  #'  =======================================
  #'  Time-between-detections analysis
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing, Cameron Ho, Hunter Hicks
  #'  July 2022
  #'  =======================================
  #'  Calculate times-between-detections of predators and prey at sites where 
  #'  cattle/hunters are and are not detected. Also calculate times-between-
  #'  detections of wildlife and cattle/hunters. Analyzed whether times-between-
  #'  detections differ as a result of recent cattle/hunter present.
  #'  ================================
  
  #'  Load packages
  library(data.table)
  library(lubridate)
  library(chron)
  library(tidyverse)
  
  #'  Read in and format data
  megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdata18-21_2022-04-27.csv") %>%
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
    #'  Remove observations that are still NA
    filter(!is.na(Species)) %>%
    mutate(
      DateTime = as.POSIXct(DateTime,
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles"),
      Date = as.Date(Date, format = "%Y-%m-%d"),
      Time = chron(times = Time)
    )
  
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
  
  #'  Add column identifying predators, prey, cattle, hunters, and other
  capdata <- capdata %>%
    mutate(Category = ifelse(Species == "Bobcat" | Species == "Black Bear" | 
                                     Species == "Cougar" | Species == "Coyote" | 
                                     Species == "Lynx" | Species == "Wolf", "Predator", "Other"),
           Category = ifelse(Species == "Elk" | Species == "Mule Deer" | 
                                     Species == "Moose" | Species == "White-tailed Deer", "Prey", Category),
           Category = ifelse(Species == "Cattle", "Cattle", Category),
           Category = ifelse(HumanActivity == "Hunter Bow" & !is.na(HumanActivity), "Hunter", Category),
           Category = ifelse(HumanActivity == "Hunter Rifle" & !is.na(HumanActivity), "Hunter", Category)) 
  
 
  #'  Filter predator, cattle, and hunter data to the last image of each unique detection event
  lastpredator <- capdata[capdata$Category == "Predator",] %>% 
    group_by(caps) %>% 
    slice_tail() %>%
    ungroup()
  lastcow <- capdata[capdata$Category == "Cattle",] %>% 
    group_by(caps) %>% 
    slice_tail() %>%
    ungroup()
  lasthunter <- capdata[capdata$Category == "Hunter",] %>% 
    group_by(caps) %>% 
    slice_tail() %>%
    ungroup()
  
  #'  Filter predator, prey, and other data to the first image from each unique detection event
  firstprey <- capdata[capdata$Category == "Prey",] %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  firstpredator <- capdata[capdata$Category == "Predator",] %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  firstother <- capdata[capdata$Category == "Other",] %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  
  #'  Merge and filter to specific date ranges of interest
  #'  Want to include detections of ALL specie and humans to account for non-focal
  #'  species that are detected between species of interest, but doesn't matter if
  #'  its the first or last image of their unique detection events here
  grazing_season_dets <- rbind(lastpredator, firstprey, lastcow, lasthunter, firstother) %>% 
    arrange(CameraLocation, DateTime) %>%
    filter(Date > "2018-06-30" & Date < "2018-09-30" | Date > "2019-06-30" & Date < "2019-09-30" | 
             Date > "2020-06-30" & Date < "2020-09-30")
  hunting_season_dets <- rbind(lastpredator, firstprey, lastcow, lasthunter, firstother) %>% 
    arrange(CameraLocation, DateTime) %>%
    filter(Date > "2018-09-30" & Date < "2018-11-26" | Date > "2019-09-30" & Date < "2019-11-26" | 
             Date > "2020-09-30" & Date < "2020-11-26")
  #'  Use first image of predator detection in relation to cattle and hunters
  grazing_predator_dets <- rbind(firstpredator, lastcow, firstprey, lasthunter, firstother) %>%
    arrange(CameraLocation, DateTime) %>%
    filter(Date > "2018-06-30" & Date < "2018-09-30" | Date > "2019-06-30" & Date < "2019-09-30" | 
             Date > "2020-06-30" & Date < "2020-09-30") 
  hunting_predator_dets <- rbind(firstpredator, lasthunter, lastcow, firstprey, firstother) %>%
    arrange(CameraLocation, DateTime) %>%
    filter(Date > "2018-09-30" & Date < "2018-11-26" | Date > "2019-09-30" & Date < "2019-11-26" | 
             Date > "2020-09-30" & Date < "2020-11-26")
  
  #'  Group multiple detection events of same category (but of different species) 
  #'  when they occur sequentially, then reduce to a single observation (e.g., 
  #'  we only care about the LAST of the last predator detections in a series of 
  #'  predator detections).
  thin_dat <- function(dets) {
    dat <- arrange(dets, CameraLocation, DateTime)
    caps_new <- c()
    caps_new[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps_new[i] = i
      else(if (dat$Category[i-1] != dat$Category[i]) caps_new[i] = i
           else caps_new[i] = caps_new[i-1])
    }
    
    caps_new <- as.factor(caps_new)
    
    #'  Add new column to larger data set
    capdata <- cbind(as.data.frame(dat), caps_new)
    
    #'  Remove all extra detections when multiple detections of same category occur in a row
    firstpreyspp <- capdata[capdata$Category == "Prey",] %>%
      group_by(caps_new) %>% 
      slice(1L) %>%
      ungroup()
    lasteverythingelse <- capdata[capdata$Category != "Prey",] %>%
      group_by(caps_new) %>% 
      slice_tail() %>%
      ungroup()
    #'  Combin into final dataset
    dets <- rbind(firstpreyspp, lasteverythingelse) %>%
      arrange(CameraLocation, DateTime)
    return(dets)
  }
  grazing_season_dets_thin <- thin_dat(grazing_season_dets)
  hunting_season_dets_thin <- thin_dat(hunting_season_dets)
  grazing_predator_dets_thin <- thin_dat(grazing_predator_dets)
  hunting_predator_dets_thin <- thin_dat(hunting_predator_dets)
  
  #'  Function to reduce detections to just series of a spp1 detection followed
  #'  by a spp2 detection (e.g., predator then prey, but no other, cattle, or hunter 
  #'  detections in between)
  spppair_dat <- function(dets, spp1, spp2) {
    #'  Assign same ID to all detection events from the same camera
    dat <- arrange(dets, CameraLocation, DateTime)
    cam <- c()
    cam[1] <- 1
    for (i in 2:nrow(dat)){
      if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) cam[i] = i
      else cam[i] = cam[i-1]
    }
    
    #'  Identify images where a spp1 is detected right before a spp2
    spp_caps1 <- c()
    spp_caps1[1] <- "N"
    for (i in 2:nrow(dat)){
      if (dat$Category[i-1] == spp1 & dat$Category[i] == spp2) spp_caps1[i-1] = "Y"
      else spp_caps1[i-1] = "N"
    }
    #'  Add "N" to very end of vector so the length matches number of rows in dets
    spp_caps1[nrow(dat)] <- "N"
    
    #'  Identify images where a spp2 is detected right after a spp1
    spp_caps2 <- c()
    spp_caps2[1] <- "N"
    for (i in 2:nrow(dat)){
      if (dat$Category[i-1] == spp1 & dat$Category[i] == spp2) spp_caps2[i] = "Y"
      else spp_caps2[i] = "N"
    }

    #'  Identify the species pair for easier filtering later on
    spp1spp2 <- c()
    spp1spp2[1] <- NA
    for(i in 2:nrow(dat)){
      if (dat$Category[i-1] == spp1 & dat$Category[i] == spp2) spp1spp2[i] = paste0(dat$Species[i-1], "_", dat$Species[i])
      else spp1spp2[i] =  NA
    }
    
    #'  Add new columns to larger data set and filter
    spp_new1 <- as.factor(spp_caps1)
    spp_new2 <- as.factor(spp_caps2)
    spp_pair <- as.factor(spp1spp2)
    capdata <- cbind(as.data.frame(dat), cam, spp_new1, spp_new2, spp_pair) %>%
      #'  Filter detection events to just situations where spp2 follows spp1
      filter(spp_new1 == "Y" | spp_new2 == "Y")
    
    return(capdata)
  }
  resp2pred_graze <- spppair_dat(grazing_season_dets_thin, spp1 = "Predator", spp2 = "Prey")
  resp2pred_hunt <- spppair_dat(hunting_season_dets_thin, spp1 = "Predator", spp2 = "Prey")
  predresp2cattle <- spppair_dat(grazing_predator_dets_thin, spp1 = "Cattle", spp2 = "Predator")
  preyresp2cattle <- spppair_dat(grazing_predator_dets_thin, spp1 = "Cattle", spp2 = "Prey")
  predresp2hunt <- spppair_dat(hunting_predator_dets_thin, spp1 = "Hunter", spp2 = "Predator")
  preyresp2hunt <- spppair_dat(hunting_predator_dets_thin, spp1 = "Hunter", spp2 = "Prey")
  
  
  ####  Calculate times between detection events  ####
  #'  --------------------------------------------
  #'  Function to calculate time between detection events of two focal species
  #'  Data structured so only last image of spp1 and first image of spp2 per
  #'  detection event are included in data frame.
  tbd <- function(detection_data, spp1, unittime) {
    #'  Create empty vector to be filled
    detection_data$TimeSinceLastDet <- c()
    #'  Fill first element of the vector to get it started 
    detection_data$TimeSinceLastDet[1] <- 0
    #'  Loop throug each row to calculate elapsed time since previous detection
    for (i in 2:nrow(detection_data)){
      #'  If previous detection was spp1, set time to 0
      if (detection_data$Category[i-1] == spp1) detection_data$TimeSinceLastDet[i] = 0
      #'  If current detection is spp2 and follows detection of spp1, calcualte
      #'  the difference in time from previous detection to current detection
      if (detection_data$Category[i] != spp1) detection_data$TimeSinceLastDet[i] = difftime(detection_data$DateTime[i], detection_data$DateTime[i-1], units = unittime)
    }
    return(detection_data)
  }
  #'  Calculate time between detections for different pairs of species of interest
  #'  spp1 should be the species detected first, unittime is the unit of time 
  #'  to make calculations in (options are: "sec", "min", "hour", "day")
  tbd_pred.prey_graze <- tbd(resp2pred_graze, spp1 = "Predator", unittime = "min")
  tbd_pred.prey_hunt <- tbd(resp2pred_hunt, spp1 = "Predator", unittime = "min")
  tbd_cattle.pred <- tbd(predresp2cattle, spp1 = "Cattle", unittime = "min")
  tbd_cattle.prey <- tbd(preyresp2cattle, spp1 = "Cattle", unittime = "min")
  tbd_hunter.pred <- tbd(predresp2hunt, spp1 = "Hunter", unittime = "min")
  tbd_hunter.prey <- tbd(preyresp2hunt, spp1 = "Hunter", unittime = "min")

  #'  Split data by species pairs of interest
  unique(tbd_pred.prey_graze$spp_pair)
  tbd_coug.md_graze <- filter(tbd_pred.prey_graze, spp_pair == "Cougar_Mule Deer")
  tbd_coug.elk_graze <- filter(tbd_pred.prey_graze, spp_pair == "Cougar_Elk")
  tbd_coug.wtd_graze <- filter(tbd_pred.prey_graze, spp_pair == "Cougar_White-tailed Deer")
  tbd_coug.moose_graze <- filter(tbd_pred.prey_graze, spp_pair == "Cougar_Moose")
  tbd_wolf.md_graze <- filter(tbd_pred.prey_graze, spp_pair == "Wolf_Mule Deer")
  tbd_wolf.elk_graze <- filter(tbd_pred.prey_graze, spp_pair == "Wolf_Elk")
  tbd_wolf.wtd_graze <- filter(tbd_pred.prey_graze, spp_pair == "Wolf_White-tailed Deer")
  tbd_wolf.moose_graze <- filter(tbd_pred.prey_graze, spp_pair == "Wolf_Moose")
  tbd_bear.md_graze <- filter(tbd_pred.prey_graze, spp_pair == "Black Bear_Mule Deer")
  tbd_bear.elk_graze <- filter(tbd_pred.prey_graze, spp_pair == "Black Bear_Elk")
  tbd_bear.wtd_graze <- filter(tbd_pred.prey_graze, spp_pair == "Black Bear_White-tailed Deer")
  tbd_bear.moose_graze <- filter(tbd_pred.prey_graze, spp_pair == "Black Bear_Moose")
  tbd_bob.md_graze <- filter(tbd_pred.prey_graze, spp_pair == "Bobcat_Mule Deer")
  tbd_bob.wtd_graze <- filter(tbd_pred.prey_graze, spp_pair == "Bobcat_White-tailed Deer")
  tbd_coy.md_graze <- filter(tbd_pred.prey_graze, spp_pair == "Coyote_Mule Deer")
  tbd_coy.wtd_graze <- filter(tbd_pred.prey_graze, spp_pair == "Coyote_White-tailed Deer")
  
  unique(tbd_pred.prey_hunt$spp_pair)
  tbd_coug.md_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Cougar_Mule Deer")
  tbd_coug.elk_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Cougar_Elk")
  tbd_coug.wtd_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Cougar_White-tailed Deer")
  tbd_coug.moose_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Cougar_Moose")
  tbd_wolf.md_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Wolf_Mule Deer")
  tbd_wolf.elk_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Wolf_Elk")
  tbd_wolf.wtd_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Wolf_White-tailed Deer")
  tbd_wolf.moose_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Wolf_Moose")
  tbd_bear.md_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Black Bear_Mule Deer")
  tbd_bear.elk_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Black Bear_Elk")
  tbd_bear.wtd_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Black Bear_White-tailed Deer")
  tbd_bear.moose_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Black Bear_Moose")
  tbd_bob.md_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Bobcat_Mule Deer")
  tbd_bob.wtd_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Bobcat_White-tailed Deer")
  tbd_coy.md_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Coyote_Mule Deer")
  tbd_coy.wtd_hunt <- filter(tbd_pred.prey_hunt, spp_pair == "Coyote_White-tailed Deer")
  
  unique(tbd_cattle.pred$spp_pair)
  tbd_coug.cattle <- filter(tbd_cattle.pred, spp_pair == "Cattle_Cougar")
  tbd_wolf.cattle <- filter(tbd_cattle.pred, spp_pair == "Cattle_Wolf")
  tbd_bear.cattle <- filter(tbd_cattle.pred, spp_pair == "Cattle_Black Bear")
  tbd_bob.cattle <- filter(tbd_cattle.pred, spp_pair == "Cattle_Bobcat")
  tbd_coy.cattle <- filter(tbd_cattle.pred, spp_pair == "Cattle_Coyote")
  tbd_lynx.cattle <- filter(tbd_cattle.pred, spp_pair == "Cattle_Lynx")
  unique(tbd_cattle.prey$spp_pair)
  tbd_md.cattle <- filter(tbd_cattle.prey, spp_pair == "Cattle_Mule Deer")
  tbd_elk.cattle <- filter(tbd_cattle.prey, spp_pair == "Cattle_Elk")
  tbd_wtd.cattle <- filter(tbd_cattle.prey, spp_pair == "Cattle_White-tailed Deer")
  tbd_moose.cattle <- filter(tbd_cattle.prey, spp_pair == "Cattle_Moose")
  
  unique(tbd_hunter.pred$spp_pair)
  tbd_coug.hunter <- filter(tbd_hunter.pred, spp_pair == "Human_Cougar")
  tbd_wolf.hunter <- filter(tbd_hunter.pred, spp_pair == "Human_Wolf")
  tbd_bear.hunter <- filter(tbd_hunter.pred, spp_pair == "Human_Black Bear")
  tbd_bob.hunter <- filter(tbd_hunter.pred, spp_pair == "Human_Bobcat")
  tbd_coy.hunter <- filter(tbd_hunter.pred, spp_pair == "Human_Coyote")
  unique(tbd_hunter.prey$spp_pair)
  tbd_md.hunter <- filter(tbd_hunter.prey, spp_pair == "Human_Mule Deer")
  tbd_elk.hunter <- filter(tbd_hunter.prey, spp_pair == "Human_Elk")
  tbd_wtd.hunter <- filter(tbd_hunter.prey, spp_pair == "Human_White-tailed Deer")
  tbd_moose.hunter <- filter(tbd_hunter.prey, spp_pair == "Human_Moose")
  
  
  ####  NEXT- run some models - get binary covariates pulled together (Y/N for grazing, hunting, public land)
  ####  Permutation test approach or same analysis but prey-prey and predator-predator?
  ####  Think about other ways to compare relationships, eg. maybe only care about 
  ####  tbd for pred and prey, even if something else (other category) is detected
  ####  btwn them
  ####  Do I split out by species combos or just lump all species together (e.g.,
  ####  time between detection of any predator and any prey)?
  
  
  
    
  
  
  