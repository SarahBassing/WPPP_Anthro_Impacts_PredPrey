  #'  =======================================================
  #'  Temporal overlap patterns in response to cattle/hunters
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  June 2022
  #'  =======================================================
  #'  Script to estimate temporal overlap between predators and prey during the
  #'  livestock grazing and hunting seasons. 
  #'  =======================================================
  
  #'  Load libraries
  library(lubridate)
  library(chron)
  library(overlap)
  library(sp)
  library(tidyverse)
  
  #'  Load and format detection data
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
      HumanActivity = ifelse(HumanActivity == "", "NA", HumanActivity),
      HumanActivity = ifelse(HumanActivity == "Hunter Rifle", "Hunter", HumanActivity),
      HumanActivity = ifelse(HumanActivity == "Hunter Bow", "Hunter", HumanActivity),
      HumanActivity = ifelse(HumanActivity == "Vehicle Truck Car", "Vehicle", HumanActivity),
      HumanActivity = ifelse(HumanActivity == "Vehicle ATV", "Vehicle", HumanActivity)
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
      Time = chron(times = Time),
      radTime = ((hour(DateTime)*60 + minute(DateTime) + second(DateTime)/(60))/1440)*2*pi,
    )
  
  xy <- SpatialPoints(cbind(megadata$Camera_Long, megadata$Camera_Lat), 
                      proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  #'  Convert radians to sun times to account for seasonal changes in sunlight 
  sunTime <- sunTime(megadata$radTime, Dates = megadata$DateTime, Coords = xy)
  head(sunTime)
  
  #'  Add to larger data set
  megadata$sunTime <- sunTime
  
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
  
  #'  Retain only the last image from each unique detection event
  last_det <- capdata %>%
    group_by(caps) %>%
    slice_tail(n = 1) %>%
    ungroup()
  #'  Keep in mind the 30 minute gap too large for cattle & human detections
  #'  Use Cattle_Hunter_Detections.R to generate DH for these
  
  #'  Filter data to desired date ranges
  grazing_filter <- function(dets) {
    cattle2018 <- dets %>%
      filter(Date > "2018-06-30") %>%
      filter(Date < "2018-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity")  
    cattle2019 <- dets %>%
      filter(Date > "2019-06-30") %>%
      filter(Date < "2019-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity")  
    cattle2020 <- dets %>%
      filter(Date > "2020-06-30") %>%
      filter(Date < "2020-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity")  
    detections <- rbind(cattle2018, cattle2019, cattle2020)
    return(detections)
  }
  grazing_first <- grazing_filter(detections)
  grazing_last <- grazing_filter(last_det)
  
  #'  Hunting season
  hunting_filter <- function(dets) {
    hunters2018 <- dets %>%
      filter(Date > "2018-09-30") %>%
      filter(Date < "2018-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity") 
    hunters2019 <- dets %>%
      filter(Date > "2019-09-30") %>%
      filter(Date < "2019-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity")  
    hunters2020 <- dets %>%
      filter(Date > "2020-09-30") %>%
      filter(Date < "2020-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity")  
    last_dets <- rbind(hunters2018, hunters2019, hunters2020)
    return(last_dets)
  }
  hunting_first <- hunting_filter(detections)
  hunting_last <- hunting_filter(last_det)
  
  est_overlap <- function(spp1, spp2, spp3, name1, name2, name3, nboot, i) {
    #'  Create logical vectors (T/F) indicating whether spp1 & spp2 were detected 
    #'  at the same site and reduce detection events to just those cameras --> These 
    #'  species need to spatially overlap for any temporal overlap to be meaningful
    both.present <- spp1$CameraLocation %in% spp2$CameraLocation
    spp1_dat <- cbind(spp1, both.present) 
    spp1_dat <- spp1_dat[spp1_dat$both.present == T,]
    both.present <- spp2$CameraLocation %in% spp1$CameraLocation
    spp2_dat <- cbind(spp2, both.present)
    spp2_dat <- spp2_dat[spp2_dat$both.present == T,]
    #'  Double check the same number of cameras are included in each
    length(unique(spp1_dat$CameraLocation)); length(unique(spp2_dat$CameraLocation))
    
    #'  Create a logical vector (T/F) indicating whether spp3 was detected at each
    #'  site where spp1 was detected (based on cameras where spp1 & spp2 were detected)
    spp3.present <- spp1_dat$CameraLocation %in% spp3$CameraLocation
    spp1_dat <- cbind(spp1_dat, spp3.present)
    #'  Split out data into camera locations where both are present vs spp3 absent
    spp1_spp3.present <- spp1_dat[spp1_dat$spp3.present == T,]
    spp1_spp3.absent <- spp1_dat[spp1_dat$spp3.present == F,]
    
    #'  Create a logical vector (T/F) indicating whether spp3 was detected at each
    #'  site where spp2 was detected (based on cameras where spp1 & spp2 were detected)
    spp3.present <- spp2_dat$CameraLocation %in% spp3$CameraLocation
    spp2_dat <- cbind(spp2_dat, spp3.present)
    #'  Split out data into camera locations where both are present vs spp3 absent
    spp2_spp3.present <- spp2_dat[spp2_dat$spp3.present == T,]
    spp2_spp3.absent <- spp2_dat[spp2_dat$spp3.present == F,]
    
    #'  Review sample size per species- smaller sample will determine which coefficient
    #'  of overlap estimator to use (delta1 for small samples [<50 detection events], 
    #'  delta4 for larger samples [>50 detection events])
    print(nrow(spp1_spp3.present)); print(nrow(spp2_spp3.present))
    print(nrow(spp1_spp3.absent)); print(nrow(spp2_spp3.absent))
    
    #'  Visualize general temporal activity with density plots
    densityPlot(spp1$sunTime,rug = F, col = "blue", main = paste0("Density Plot of ", name1, " Daily Activity"))
    densityPlot(spp2$sunTime,rug = F, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    densityPlot(spp3$sunTime,rug = F, col = "blue", main = paste0("Density Plot of ", name3, " Daily Activity"))
  
    #'  Visualize temporal overlap
    overlapPlot(spp1$sunTime, spp2$sunTime, rug = F, linet = c(1, 1), 
                linec = c("red", "blue"), linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) and ", name2, " (blue)"))
    
    overlapPlot(spp1$sunTime, spp3$sunTime, rug = F, linet = c(1, 1), 
                linec = c("red", "yellow"), linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) and ", name3, " (green)"))
    
    overlapPlot(spp2$sunTime, spp3$sunTime, rug = F, linet = c(1, 1), 
                linec = c("blue", "yellow"), linew = c(2, 2), main = paste0("Overlap Plots of ", name2, " (blue) and ", name3, " (green)"))
    
    
    overlapPlot(spp1_spp3.present$sunTime, spp2_spp3.present$sunTime, rug = F, linet = c(1, 1), 
                linec = c("red", "blue"), linew = c(2, 2), #ylim = c(0, 0.08), 
                main = paste0("Overlap Plots of ", name1, " (red) and ", name2, " (blue) \nwhen ", name3, " are Present"))
    
    overlapPlot(spp1_spp3.absent$sunTime, spp2_spp3.absent$sunTime, rug = F, linet = c(1, 1), 
                linec = c("red", "blue"), linew = c(2, 2),  #ylim = c(0, 0.08),
                main = paste0("Overlap Plots of ", name1, " (red) and ", name2, " (blue) \nwhen ", name3, " are Absent"))
    
    #'  Calculate coefficient of overlap
    dhats_spp1.spp2.spp3 <- overlapEst(A = spp1_spp3.present$sunTime, B = spp2_spp3.present$sunTime)
    print(dhats_spp1.spp2.spp3)
    dhats_spp1.spp2.NOspp3 <- overlapEst(A = spp1_spp3.absent$sunTime, B = spp2_spp3.absent$sunTime)
    print(dhats_spp1.spp2.NOspp3)
    
    #'  Bootstrap to estimate standard errors for when spp3 is present
    spp1.spp3.boot <- resample(spp1_spp3.present$sunTime, nboot)
    spp2.spp3.boot <- resample(spp2_spp3.present$sunTime, nboot)
    dhat.spp3.present.boot <- bootEst(spp1.spp3.boot, spp2.spp3.boot, nboot)
    print(head(dhat.spp3.present.boot))
    
    #'  Bootstrap to estimate standard errors for when spp3 is absent
    spp1.NOspp3.boot <- resample(spp1_spp3.absent$sunTime, nboot)
    spp2.NOspp3.boot <- resample(spp2_spp3.absent$sunTime, nboot)
    dhat.spp3.absent.boot <- bootEst(spp1.NOspp3.boot, spp2.NOspp3.boot, nboot)
    print(head(dhat.spp3.absent.boot))
    
    #'  Pull out vector of bootstrapped estimates for one of the coefficient of
    #'  overlap metrics (1 = delta1, 2 = delta4)
    dh4_spp3.present_CI <- as.vector(dhat.spp3.present.boot[,i])
    dh4_spp3.absent_CI <- as.vector(dhat.spp3.absent.boot[,i])
    
    #'  Bootstrap 95% confidence intervals
    spp3.present_ci <- bootCI(dhats_spp1.spp2.spp3[i], dh4_spp3.present_CI, conf = 0.95)
    print(spp3.present_ci)
    spp3.absent_ci <- bootCI(dhats_spp1.spp2.NOspp3[i], dh4_spp3.absent_CI, conf = 0.95)
    print(spp3.absent_ci)
    
    #'  Save as a giant list
    overlap_list <- list(dhats_spp1.spp2.spp3, dhats_spp1.spp2.NOspp3, 
                         dhat.spp3.present.boot, dhat.spp3.absent.boot, 
                         spp3.present_ci, spp3.absent_ci)
    return(overlap_list)
    }
  coug_md_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                   spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Cougar", name2 = "Mule Deer", 
                                   name3 = "Hunters", nboot = 1000, i = 1)
  coug_elk_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                    spp2 = filter(hunting_first, Species == "Elk"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Cougar", name2 = "Elk", 
                                    name3 = "Hunters", nboot = 1000, i = 1)
  coug_wtd_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                    spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Cougar", name2 = "White-tailed Deer", 
                                    name3 = "Hunters", nboot = 1000, i = 1)
  coug_moose_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                      spp2 = filter(hunting_first, Species == "Moose"), 
                                      spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                      name1 = "Cougar", name2 = "Moose", 
                                      name3 = "Hunters", nboot = 1000, i = 1)
  wolf_md_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                   spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Wolf", name2 = "Mule Deer", 
                                   name3 = "Hunters", nboot = 1000, i = 1)
  wolf_elk_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                    spp2 = filter(hunting_first, Species == "Elk"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Wolf", name2 = "Elk", 
                                    name3 = "Hunters", nboot = 1000, i = 1)
  wolf_wtd_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                    spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Wolf", name2 = "White-tailed Deer", 
                                    name3 = "Hunters", nboot = 1000, i = 1)
  wolf_moose_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                      spp2 = filter(hunting_first, Species == "Moose"), 
                                      spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                      name1 = "Wolf", name2 = "Moose", 
                                      name3 = "Hunters", nboot = 1000, i = 2)
  bear_md_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                  spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                  spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                  name1 = "Black Bear", name2 = "Mule Deer", 
                                  name3 = "Hunters", nboot = 1000, i = 2)
  bear_elk_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                   spp2 = filter(hunting_first, Species == "Elk"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Black Bear", name2 = "Elk", 
                                   name3 = "Hunters", nboot = 1000, i = 2)
  bear_wtd_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                    spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Black Bear", name2 = "White-tailed Deer", 
                                    name3 = "Hunters", nboot = 1000, i = 2)
  bear_moose_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                   spp2 = filter(hunting_first, Species == "Moose"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Black Bear", name2 = "Moose", 
                                   name3 = "Hunters", nboot = 1000, i = 2)
  bob_md_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                  spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                  spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                  name1 = "Bobcat", name2 = "Mule Deer", 
                                  name3 = "Hunters", nboot = 1000, i = 1)
  bob_wtd_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                   spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Bobcat", name2 = "White-tailed Deer", 
                                   name3 = "Hunters", nboot = 1000, i = 1)
  coy_md_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                  spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                  spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                  name1 = "Coyote", name2 = "Mule Deer", 
                                  name3 = "Hunters", nboot = 1000, i = 2)
  coy_wtd_hunt_over <- est_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                   spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Coyote", name2 = "White-tailed Deer", 
                                   name3 = "Hunters", nboot = 1000, i = 2)
  
  pred_prey_hunt_overlap <- list(coug_md_hunt_over, coug_elk_hunt_over, coug_wtd_hunt_over, coug_moose_hunt_over,
                                 wolf_md_hunt_over, wolf_elk_hunt_over, wolf_wtd_hunt_over, wolf_moose_hunt_over,
                                 bear_md_hunt_over, bear_elk_hunt_over, bear_wtd_hunt_over, bear_moose_hunt_over,
                                 bob_md_hunt_over, bob_elk_hunt_over, coy_wtd_hunt_over, coy_moose_hunt_over)
  
  save(pred_prey_hunt_overlap, file = "./Temporal Overlap/pred-prey_hunt_overlap.RData")
  save.image(file = "./Temporal Overlap/temporal_overlap_analyses_workspace.RData")
  
  # figure out how to change Dhat based on sample size... e.g. bobcat definitely need Dhat1
  
  
    ####  THINGS TO THINK ABOUT: Do I include vehicles in "hunter presence"? otherwise not that many cameras will have hunting activity
  
  #'  Unique, independent detections
  coyote <- filter(hunting_first, Species == "Coyote")
  wtd <- filter(hunting_first, Species == "White-tailed Deer")
  hunt <- filter(hunting_first, HumanActivity == "Hunter")
  
  #'  Create a logical vector (T/F) for spp1 locations indicating human detections
  hunt.present <- coyote$CameraLocation %in% hunt$CameraLocation
  coyote <- cbind(coyote, hunt.present)
  coy.present <- coyote[coyote$hunt.present == T,]
  coy.absent <- coyote[coyote$hunt.present == F,]
  
  hunt.present <- wtd$CameraLocation %in% hunt$CameraLocation
  wtd <- cbind(wtd, hunt.present)
  wtd.present <- wtd[wtd$hunt.present == T,]
  wtd.absent <- wtd[wtd$hunt.present == F,]
  
  #'  Review sample size per species- smaller sample will determine which coefficient
  #'  of overlap estimator to use (delta1 for small samples, delta4 for larger samples)
  nrow(coy)
  nrow(coy.present)
  nrow(coy.absent)
  nrow(wtd)
  nrow(wtd.present)
  nrow(wtd.absent)
  
  #'  Density plots of temporal activity patterns
  densityPlot(coyote$sunTime,rug = F, col = "blue")
  densityPlot(wtd$sunTime,rug = F,col = "blue")
  densityPlot(hunt$sunTime,rug = F,col = "blue")
  
  
  #'  Overlap plots
  overlapPlot(coyote$sunTime, wtd$sunTime, rug = F, linet = c(1, 1), 
              linec = c("red", "blue"), linew = c(2, 2))
  
  overlapPlot(coy.present$sunTime, wtd.present$sunTime, rug = F, linet = c(1, 1), 
              linec = c("red", "blue"), linew = c(2, 2))
  
  overlapPlot(coy.absent$sunTime, wtd.absent$sunTime, rug = F, linet = c(1, 1), 
              linec = c("red", "blue"), linew = c(2, 2))
  
  
  #'  Coefficient of overlap
  dhats.coy.wtd <- overlapEst(A = coyote$sunTime, B = wtd$sunTime)
  dhats.coy.wtd
  
  #'  Bootstrap to estimate standard errors
  coy.boot <- resample(coyote$sunTime, 1000)
  wtd.boot <- resample(wtd$sunTime, 1000)
  dhat.boot <- bootEst(coy.boot, wtd.boot, 1000)
  
  dh4_CI <- as.vector(dhat.boot[,2])
  ci <- bootCI(dhats.coy.wtd[2], dh4_CI, conf = 0.95)
  
  