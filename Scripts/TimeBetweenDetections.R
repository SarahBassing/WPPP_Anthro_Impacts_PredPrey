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
  
  #' #'  Filter data to grazing and hunting seasons
  #' grazing_dat <- filter(capdata, Date > "2018-06-30" & Date < "2018-09-30" | 
  #'                         Date > "2019-06-30" & Date < "2019-09-30" | 
  #'                         Date > "2020-06-30" & Date < "2020-09-30")
  #' hunting_dat <- filter(capdata, Date > "2018-09-30" & Date < "2018-11-26" | 
  #'                         Date > "2019-09-30" & Date < "2019-11-26" | 
  #'                         Date > "2020-09-30" & Date < "2020-11-26")
  
 
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
  cat_dat <- function(dets) {
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
  grazing_season_dets_thin <- cat_dat(grazing_season_dets)
  hunting_season_dets_thin <- cat_dat(hunting_season_dets)
  grazing_predator_dets_thin <- cat_dat(grazing_predator_dets)
  hunting_predator_dets_thin <- cat_dat(hunting_predator_dets)
  
  detect_data <- grazing_season_dets_thin
  
  #'  Function to calculate time between detection events of focal species
  #'  CANNOT have multiple observations of same category in a row so MUST reduce
  #'  detection data down to the very last or very first detection of a category,
  #'  regardless if different species within a category are detected in a row
  #'  (e.g., predator [coyote], predator [bear], predator [bear], then prey [deer] 
  #'  reduces down to --> predator [bear] then prey [deer]). We care about the
  #'  time between detections of the LAST predator and FIRST prey species only.
  tbd <- function(detect_data, spp1) {
    detect_data$TimeSinceSpp1 <- NA
    #'  Calculate time between most recent spp1 detection & responding spp's detection
    camera <- unique(detect_data$CameraLocation)
    for(i in 1:length(camera)){
      #'  Select which camera to process
      tmp_camera <- detect_data[detect_data$CameraLocation == camera[i],]
      #'  Create empty object to hold detection event
      MostRecentSpp1Det <- NA

      #'  Loop through each image for a camera
      for(j in 1:nrow(tmp_camera)){
        #'  Select which image to compare to most recent detection of spp1
        tmp_image <- detect_data[j,]

        #'  If the image is of spp1, mark it as the most recent spp1 observation
        if(tmp_image$Category == spp1){
          MostRecentSpp1Det <- tmp_image
        }

        #'  If the image is of a responding species (e.g., prey)
        #'  Suppress a warning created by is.na() when MostRecentSpp1Det is a row and not NA
        suppressWarnings(if(tmp_image$Category != spp1){
          #'  Mark MostRecentSpp1Det as NA if there has not been a spp1 detected yet
          if(is.na(MostRecentSpp1Det) == TRUE){detect_data$TimeSinceSpp1[j] <- NA}
          #'  Otherwise, calculate difference between latest spp1 detection event
          #'  and responding species detection
          else(detect_data$TimeSinceSpp1[j] <- difftime(tmp_image$DateTime, MostRecentSpp1Det$DateTime, units="mins"))
        })
      }
    }
    return(detect_data)
  }
  #'  Calculate time between detection of a predator and every other category 
  #'  until the next predator detection in the grazing and hunting seasons
  response2predator_graze <- tbd(grazing_season_dets_thin, spp1 = "Predator")
  response2predator_hunt <- tbd(hunting_season_dets_thin, spp1 = "Predator")
  #'  Calculate time between detection of cattle/hunters and every other category 
  #'  until the next cattle/hunter detection in the grazing and hunting seasons
  response2cattle_graze <- tbd(grazing_predator_dets_thin, spp1 = "Cattle")
  response2hunter_hunt <- tbd(hunting_predator_dets_thin, spp1 = "Hunter")
  
  ####  NEXT- reduce these to just the first detection after the species of interest
  ####  Then run some models?
  ####  Think about other ways to compare relationships, eg. maybe only care about 
  ####  tbd for pred and prey, even if something else (other category) is detected
  ####  btwn them, or only care about tbd of apex pred and prey, even if other pred 
  ####  detected btwn them???
  
  
  
    
  
  
  