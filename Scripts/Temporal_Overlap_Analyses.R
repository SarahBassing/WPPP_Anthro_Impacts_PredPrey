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
  library(circular)
  library(ggplot2)
  library(sp)
  library(tidyverse)
  
  #'  Load and format detection data
  megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdata18-21_2022-04-27.csv") %>%  
    dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation", 
                  "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle", 
                  "Species", "HumanActivity", "Count", "Land_Mgnt", "Land_Owner") %>%
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
      HumanActivity = ifelse(HumanActivity == "Vehicle ATV", "Vehicle", HumanActivity),
      #'  Identify if cameras were on public (1) or private (0) land
      Public1 = ifelse(Land_Mgnt != "Private", 1, 0),
      #'  Even though timberland is private, it's generally open to public hunting
      #'  so considering it public for these purposes
      Public1 = ifelse(Land_Mgnt == "Private" & Land_Owner == "Private timber", 1, Public1)
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
  stations_data <- read.csv("./Data/stations_data.csv") %>%
    dplyr::select(c("CameraLocation", "Year", "Study_Area", "PublicGrazing"))
  megadata <- right_join(megadata, stations_data, by = "CameraLocation")
  
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
  detections_OK <- detections %>%
    filter(grepl("OK", CameraLocation))
  
  #'  Retain only the last image from each unique detection event
  last_det <- capdata %>%
    group_by(caps) %>%
    slice_tail(n = 1) %>%
    ungroup()
  last_det_OK <- last_det %>%
    filter(grepl("OK", CameraLocation))
  
  #'  Filter data to desired date ranges
  grazing_filter <- function(dets) {
    cattle2018 <- dets %>%
      filter(Date > "2018-06-30") %>%
      filter(Date < "2018-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity", "PublicGrazing")  
    cattle2019 <- dets %>%
      filter(Date > "2019-06-30") %>%
      filter(Date < "2019-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity", "PublicGrazing")  
    cattle2020 <- dets %>%
      filter(Date > "2020-06-30") %>%
      filter(Date < "2020-09-30") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity", "PublicGrazing")  
    detections <- rbind(cattle2018, cattle2019, cattle2020)
    return(detections)
  }
  grazing_first <- grazing_filter(detections)
  grazing_last <- grazing_filter(last_det)
  grazing_first_OK <- grazing_filter(detections_OK)
  grazing_last_OK <- grazing_filter(last_det_OK)
  
  #'  Hunting season
  hunting_filter <- function(dets) {
    hunters2018 <- dets %>%
      filter(Date > "2018-09-30") %>%
      filter(Date < "2018-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity", "Public1") 
    hunters2019 <- dets %>%
      filter(Date > "2019-09-30") %>%
      filter(Date < "2019-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity", "Public1")  
    hunters2020 <- dets %>%
      filter(Date > "2020-09-30") %>%
      filter(Date < "2020-11-26") %>%
      dplyr::select("File", "CameraLocation", "DateTime", "Date", "Time", "radTime", "sunTime", "Species", "HumanActivity", "Public1")  
    last_dets <- rbind(hunters2018, hunters2019, hunters2020)
    return(last_dets)
  }
  hunting_first <- hunting_filter(detections)
  hunting_last <- hunting_filter(last_det)
  
  ####  Overlap analyses for predator-prey activity ####
  #'  ----------------------------------------------
  #'  Function to estimate temporal overlap between predators (spp1) and prey (spp2)
  #'  at camera sites where cattle/humans (spp3) are present (detected) vs absent (not detected).
  #'  Using spp3 to also represent cameras on public vs private land in this function.
  pred_prey_overlap <- function(spp1, spp2, spp3, name1, name2, name3, nboot, dhat) { #i
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
    #'  of overlap estimator to use (Meredith & Ridout 2017 suggest using delta1 
    #'  for small samples [<50 detection events] and delta4 for larger samples 
    #'  [>50 detection events]).
    ndet_spp1_spp3.present <- nrow(spp1_spp3.present)
    ndet_spp2_spp3.present <- nrow(spp2_spp3.present)
    ndet_spp1_spp3.absent <- nrow(spp1_spp3.absent)
    ndet_spp2_spp3.absent <- nrow(spp2_spp3.absent)
    print(ndet_spp1_spp3.present); print(ndet_spp2_spp3.present)
    print(ndet_spp1_spp3.absent); print(ndet_spp2_spp3.absent)
    
    #'  Visualize general temporal activity with density plots
    densityPlot(spp1$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name1, " Daily Activity"))
    densityPlot(spp2$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    densityPlot(spp3$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name3, " Daily Activity"))

    #'  Visualize temporal overlap
    overlapPlot(spp1$sunTime, spp2$sunTime, rug = F, linet = c(1, 1),
                linec = c("red", "blue"), linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) and ", name2, " (blue)"))

    overlapPlot(spp1$sunTime, spp3$sunTime, rug = F, linet = c(1, 1),
                linec = c("red", "yellow"), linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) and ", name3, " (green)"))

    overlapPlot(spp2$sunTime, spp3$sunTime, rug = F, linet = c(1, 1),
                linec = c("blue", "yellow"), linew = c(2, 2), main = paste0("Overlap Plots of ", name2, " (blue) and ", name3, " (green)"))


    overlapPlot(spp1_spp3.present$sunTime, spp2_spp3.present$sunTime, rug = F, linet = c(1, 1),
                linec = c("red", "blue"), linew = c(2, 2),
                main = paste0("Overlap Plots of ", name1, " (red) and ", name2, " (blue) \nwhen ", name3, " are Present"))

    overlapPlot(spp1_spp3.absent$sunTime, spp2_spp3.absent$sunTime, rug = F, linet = c(1, 1),
                linec = c("red", "blue"), linew = c(2, 2),
                main = paste0("Overlap Plots of ", name1, " (red) and ", name2, " (blue) \nwhen ", name3, " are Absent"))

    #'  Calculate coefficient of overlap
    dhats_spp1.spp2.spp3 <- overlapEst(A = spp1_spp3.present$sunTime,
                                       B = spp2_spp3.present$sunTime, type = dhat)
    dhats_spp1.spp2.NOspp3 <- overlapEst(A = spp1_spp3.absent$sunTime,
                                         B = spp2_spp3.absent$sunTime, type = dhat)

    #'  Bootstrap to estimate standard errors
    #'  FYI: smooth = TRUE is default and allows bootstrap to randomly sample from
    #'  a distribution of times that have a wider range than the original sample
    #'  (see pg. 5 in Overlap package vignette for details)
    spp12.spp3.boot <- bootstrap(spp1_spp3.present$sunTime, spp2_spp3.present$sunTime,
                                 nboot, smooth = TRUE, type = dhat)
    spp12.NOspp3.boot <- bootstrap(spp1_spp3.absent$sunTime, spp2_spp3.absent$sunTime,
                                 nboot, smooth = TRUE, type = dhat)
    #'  Bootstrap mean will be a little different then detla coefficient due to
    #'  bootstrap bias (BSmean - delta) that needs to be accounted for in 95% CIs
    BSmean.present <- mean(spp12.spp3.boot)
    BSmean.absent <- mean(spp12.NOspp3.boot)

    #'  Bootstrap 95% Confidence Intervals
    #'  norm0 uses the standard deviation of bootstrap results to calculate CI (delta +/- 1.96*SDboot)
    #'  basic0 takes the 2.5% and 97.5% percentiles and adjusts based on BS bias (percentile - BSbias)
    #'  If sampling distribution is normal, norm0 and basic0 should be similar;
    #'  if sampling distribution is skewed (i.e., if delta is close to 0 or 1) then
    #'  basic0 is the better estimator - using basic0 b/c should be good in either situation
    #'  Using bootCIlogit instead of bootCI so that bias corrections are done on
    #'  the logit scale, then backtransformed. Without this, 95% CIs can fall
    #'  outside (0, 1) interval. See Overlap vignette for more details.
    spp3.present_CI <- bootCIlogit(dhats_spp1.spp2.spp3, spp12.spp3.boot) #[i]
    spp3.absent_CI <- bootCIlogit(dhats_spp1.spp2.NOspp3, spp12.NOspp3.boot) #[i]

    #'  Print results
    #'  Effect of spp3 being present
    print("Overlap coefficients when spp3 is present"); print(dhats_spp1.spp2.spp3)
    print("Bootstrap mean"); print(BSmean.present)
    print("Bootstrap 95% CI"); print(spp3.present_CI)

    #'  Effect of spp3 being absent
    print("Overlap coefficients when spp3 is present"); print(dhats_spp1.spp2.NOspp3)
    print("Bootstrap mean"); print(BSmean.absent)
    print("Bootstrap 95% CI"); print(spp3.absent_CI)

    #'  Save as a giant list
    overlap_list <- list(dhats_spp1.spp2.spp3, dhats_spp1.spp2.NOspp3,
                         spp12.spp3.boot, spp12.NOspp3.boot,
                         spp3.present_CI, spp3.absent_CI, ndet_spp1_spp3.present,
                         ndet_spp2_spp3.present, ndet_spp1_spp3.absent, ndet_spp2_spp3.absent)
    names(overlap_list) <- c("dhat_spp3.present", "dhat_spp3.absent", "dhat_spp3.present_boot",
                             "dhat_spp3.absent_boot", "spp3.present_CI", "spp3.absent_CI",
                             "ndet_spp1_spp3.present", "ndet_spp2_spp3.present",
                             "ndet_spp1_spp3.absent", "ndet_spp2_spp3.absent")

    return(overlap_list)
  }
  ####  Predator-Prey Overlap Grazing Season  ####
  #'  Estimate temporal overlap between predators and prey when cattle are/aren't detected
  #'  Focusing on only OK study area since big difference in number of cameras 
  #'  with cattle in NE vs OK, pooling across study areas can bias results
  coug_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                         name1 = "Cougar", name2 = "Mule Deer", 
                                         name3 = "Cattle", nboot = 10000, dhat = "Dhat1") #i = 1
  # coug_elk_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
  #                                         spp2 = filter(grazing_first_OK, Species == "Elk"), 
  #                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                         name1 = "Cougar", name2 = "Elk", 
  #                                         name3 = "Cattle", nboot = 10000, dhat = "Dhat1")  
  coug_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                          name1 = "Cougar", name2 = "White-tailed Deer", 
                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat1")   # cattle absent sample size for coug SO small not really useful
  coug_moose_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"), 
                                            spp2 = filter(grazing_first_OK, Species == "Moose"), 
                                            spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                            name1 = "Cougar", name2 = "Moose", 
                                            name3 = "Cattle", nboot = 10000, dhat = "Dhat1")  
  # wolf_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
  #                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                         name1 = "Wolf", name2 = "Mule Deer", 
  #                                         name3 = "Cattle", nboot = 10000, dhat = "Dhat1") 
  # wolf_elk_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "Elk"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Wolf", name2 = "Elk", 
  #                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat1") 
  # wolf_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Wolf", name2 = "White-tailed Deer", 
  #                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat1") 
  # wolf_moose_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                            spp2 = filter(grazing_first_OK, Species == "Moose"), 
  #                                            spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                            name1 = "Wolf", name2 = "Moose", 
  #                                            name3 = "Cattle", nboot = 10000, dhat = "Dhat1")  
  bear_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                          spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                          name1 = "Black bear", name2 = "Mule Deer", 
                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat4") 
  # bear_elk_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "Elk"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Black bear", name2 = "Elk", 
  #                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat1")    # so few elk detections it's not really worth it
  bear_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                           spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                           spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                           name1 = "Black bear", name2 = "White-tailed Deer", 
                                           name3 = "Cattle", nboot = 10000, dhat = "Dhat4")   # bear-wtd: cattle present all >> 50, cattle absent bear = 50
  bear_moose_graze_over_dhat1 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                             spp2 = filter(grazing_first_OK, Species == "Moose"), 
                                             spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                             name1 = "Black bear", name2 = "Moose", 
                                             name3 = "Cattle", nboot = 10000, dhat = "Dhat1")   # bear-moose: cattle present all >50, cattle absent moose <50
  bear_moose_graze_over_dhat4 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                             spp2 = filter(grazing_first_OK, Species == "Moose"), 
                                             spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                             name1 = "Black bear", name2 = "Moose", 
                                             name3 = "Cattle", nboot = 10000, dhat = "Dhat4")   # bear-moose: cattle present all >50, cattle absent moose <50
  bob_md_graze_over_dhat1 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
                                          spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                          name1 = "Bobcat", name2 = "Mule Deer", 
                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat1")   # bob-md: cattle present bobcat >> 50, cattle absent bobcat <50
  bob_md_graze_over_dhat4 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                         name1 = "Bobcat", name2 = "Mule Deer", 
                                         name3 = "Cattle", nboot = 10000, dhat = "Dhat4")   # bob-md: cattle present bobcat >> 50, cattle absent bobcat <50
  # bob_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
  #                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
  #                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                          name1 = "Bobcat", name2 = "White-tailed Deer", 
  #                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat1") 
  coy_md_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"), 
                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                         spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                         name1 = "Coyote", name2 = "Mule Deer", 
                                         name3 = "Cattle", nboot = 10000, dhat = "Dhat4") 
  coy_wtd_graze_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"), 
                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                          spp3 = filter(grazing_first_OK, Species == "Cattle"), 
                                          name1 = "Coyote", name2 = "White-tailed Deer", 
                                          name3 = "Cattle", nboot = 10000, dhat = "Dhat4") 
  
  #'  Save all overlap results
  #'  Keep in mind bear-wtd and bobcat-md repeat with different dhat estimates
  #'  due to differences in sample size when cattle are/are not detected
  pred_prey_graze_overlap <- list(coug_md_graze_over, coug_wtd_graze_over, coug_moose_graze_over, #coug_elk_graze_over, 
                                  #wolf_md_graze_over, wolf_wtd_graze_over, wolf_moose_graze_over, #wolf_elk_graze_over, 
                                  bear_md_graze_over, bear_wtd_graze_over, #bear_elk_graze_over,
                                  bear_moose_graze_over_dhat1, bear_moose_graze_over_dhat4,  
                                  bob_md_graze_over_dhat1, bob_md_graze_over_dhat4, #bob_wtd_graze_over, 
                                  coy_md_graze_over, coy_wtd_graze_over) 
  save(pred_prey_graze_overlap, file = paste0("./Outputs/Temporal Overlap/pred_prey_graze_overlap_OK_", Sys.Date(), ".RData"))
  
  ####  Predator-Prey Overlap Grazing Allotments  ####
  #'  Estimate temporal overlap between predators and prey on grazing allotments 
  #'  vs off grazing allotments
  coug_md_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"),
                                         spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                         name1 = "Cougar", name2 = "Mule Deer",
                                         name3 = "Allotments", nboot = 10000, dhat = "Dhat1") #, i = 1
  # coug_elk_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
  #                                         spp2 = filter(grazing_first_OK, Species == "Elk"),
  #                                         spp3 = filter(grazing_first_OK, PublicGrazing == 1),
  #                                         name1 = "Cougar", name2 = "Elk",
  #                                         name3 = "Allotments", nboot = 10000, dhat = "Dhat1")
  coug_wtd_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
                                           spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"),
                                           spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                           name1 = "Cougar", name2 = "White-tailed Deer",
                                           name3 = "Allotments", nboot = 10000, dhat = "Dhat1")  # almost no cougar detections where cattle absent
  coug_moose_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
                                           spp2 = filter(grazing_first_OK, Species == "Moose"),
                                           spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                           name1 = "Cougar", name2 = "Moose",
                                           name3 = "Allotments", nboot = 10000, dhat = "Dhat1")
  wolf_md_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"),
                                          spp2 = filter(grazing_first_OK, Species == "Mule Deer"),
                                          spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                          name1 = "Wolf", name2 = "Mule Deer",
                                          name3 = "Allotments", nboot = 10000, dhat = "Dhat1") #, i = 1
  # wolf_elk_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"),
  #                                         spp2 = filter(grazing_first_OK, Species == "Elk"),
  #                                         spp3 = filter(grazing_first_OK, PublicGrazing == 1),
  #                                         name1 = "Wolf", name2 = "Elk",
  #                                         name3 = "Allotments", nboot = 10000, dhat = "Dhat1")
  # wolf_wtd_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"),
  #                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"),
  #                                          spp3 = filter(grazing_first_OK, PublicGrazing == 1),
  #                                          name1 = "Wolf", name2 = "White-tailed Deer",
  #                                          name3 = "Allotments", nboot = 10000, dhat = "Dhat1")
  wolf_moose_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"),
                                             spp2 = filter(grazing_first_OK, Species == "Moose"),
                                             spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                             name1 = "Wolf", name2 = "Moose",
                                             name3 = "Allotments", nboot = 10000, dhat = "Dhat1")
  bear_md_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"),
                                          spp2 = filter(grazing_first_OK, Species == "Mule Deer"),
                                          spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                          name1 = "Black bear", name2 = "Mule Deer",
                                          name3 = "Allotments", nboot = 10000, dhat = "Dhat4") #, i = 1
  # bear_elk_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"),
  #                                         spp2 = filter(grazing_first_OK, Species == "Elk"),
  #                                         spp3 = filter(grazing_first_OK, PublicGrazing == 1),
  #                                         name1 = "Wolf", name2 = "Elk",
  #                                         name3 = "Allotments", nboot = 10000, dhat = "Dhat1")
  bear_wtd_allot_over_dhat1 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"),
                                           spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"),
                                           spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                           name1 = "Black Bear", name2 = "White-tailed Deer",
                                           name3 = "Allotments", nboot = 100, dhat = "Dhat1")    # bear-wtd: cattle present all >>50, cattle absent bear <50
  bear_wtd_allot_over_dhat4 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"),
                                           spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"),
                                           spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                           name1 = "Black Bear", name2 = "White-tailed Deer",
                                           name3 = "Allotments", nboot = 10000, dhat = "Dhat4")  # bear-wtd: cattle present all >>50, cattle absent bear <50
  bear_moose_allot_over_dhat1 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"),
                                             spp2 = filter(grazing_first_OK, Species == "Moose"),
                                             spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                             name1 = "Black Bear", name2 = "Moose",
                                             name3 = "Allotments", nboot = 10000, dhat = "Dhat1")  # bear-moose: cattle present all >>50, cattle absent all <50
  bear_moose_allot_over_dhat4 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"),
                                             spp2 = filter(grazing_first_OK, Species == "Moose"),
                                             spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                             name1 = "Black Bear", name2 = "Moose",
                                             name3 = "Allotments", nboot = 10000, dhat = "Dhat4")  # bear-moose: cattle present all >>50, cattle absent all <50
  bob_md_allot_over_dhat1 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"),
                                          spp2 = filter(grazing_first_OK, Species == "Mule Deer"),
                                          spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                          name1 = "Bobcat", name2 = "Mule Deer",
                                          name3 = "Allotments", nboot = 10000, dhat = "Dhat1")       # bob-md: cattle present all >>50, cattle absent bobcat <50
  bob_md_allot_over_dhat4 <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"),
                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"),
                                         spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                         name1 = "Bobcat", name2 = "Mule Deer",
                                         name3 = "Allotments", nboot = 10000, dhat = "Dhat4")      # bob-md: cattle present all >>50, cattle absent bobcat <50
  # bob_wtd_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"),
  #                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"),
  #                                          spp3 = filter(grazing_first_OK, PublicGrazing == 1),
  #                                          name1 = "Bobcat", name2 = "White-tailed Deer",
  #                                          name3 = "Allotments", nboot = 10000, dhat = "Dhat1") 
  coy_md_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"),
                                         spp2 = filter(grazing_first_OK, Species == "Mule Deer"),
                                         spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                         name1 = "Coyote", name2 = "Mule Deer",
                                         name3 = "Allotments", nboot = 10000, dhat = "Dhat4") #, i = 1
  coy_wtd_allot_over <- pred_prey_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"),
                                          spp2 = filter(grazing_first_OK, Species == "White-tailed Deer"),
                                          spp3 = filter(grazing_first_OK, PublicGrazing == 1),
                                          name1 = "Coyote", name2 = "White-tailed Deer",
                                          name3 = "Allotments", nboot = 10000, dhat = "Dhat4")
  
  #'  Save all overlap results
  #'  Keep in mind bear-wtd and bobcat-md repeat with different dhat estimates
  #'  due to differences in sample size when cattle are/are not detected
  pred_prey_allot_overlap <- list(coug_md_allot_over, coug_wtd_allot_over, coug_moose_allot_over, #coug_elk_allot_allot, 
                                  #wolf_md_allot_over, wolf_moose_allot_over, #wolf_elk_allot_over, wolf_wtd_allot_over,
                                  bear_md_allot_over, bear_wtd_allot_over_dhat1, #bear_elk_allot_over,
                                  bear_wtd_allot_over_dhat4, bear_moose_allot_over_dhat1,  
                                  bear_moose_allot_over_dhat4, bob_md_allot_over_dhat1, #bob_wtd_allot_over,
                                  bob_md_allot_over_dhat4, coy_md_allot_over, coy_wtd_allot_over)  
  save(pred_prey_allot_overlap, file = paste0("./Outputs/Temporal Overlap/pred_prey_allot_overlap_OK_", Sys.Date(), ".RData"))
  
  
  
  ####  Predator-Prey Overlap Hunting Season ####
  #'  Estimate temporal overlap between predators and prey when hunters are/aren't detected
  coug_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                   spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Cougar", name2 = "Mule Deer", 
                                   name3 = "Hunters", nboot = 10000, dhat = "Dhat1") #, i = 1
  coug_elk_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                    spp2 = filter(hunting_first, Species == "Elk"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Cougar", name2 = "Elk", 
                                    name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  coug_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                    spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Cougar", name2 = "White-tailed Deer", 
                                    name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  coug_moose_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                      spp2 = filter(hunting_first, Species == "Moose"), 
                                      spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                      name1 = "Cougar", name2 = "Moose", 
                                      name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  wolf_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"),
                                   spp2 = filter(hunting_first, Species == "Mule Deer"),
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"),
                                   name1 = "Wolf", name2 = "Mule Deer",
                                   name3 = "Hunters", nboot = 10000, dhat = "Dhat1")        #  So few detections of either spp without hunters, not worth the comparison
  # wolf_elk_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"),
  #                                   spp2 = filter(hunting_first, Species == "Elk"),
  #                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"),
  #                                   name1 = "Wolf", name2 = "Elk",
  #                                   name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  wolf_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                    spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Wolf", name2 = "White-tailed Deer", 
                                    name3 = "Hunters", nboot = 10000, dhat = "Dhat1")       #  Very few wolf detections without hunters, not worth the comparison
  wolf_moose_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                      spp2 = filter(hunting_first, Species == "Moose"), 
                                      spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                      name1 = "Wolf", name2 = "Moose", 
                                      name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  bear_md_hunt_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                  spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                  spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                  name1 = "Black Bear", name2 = "Mule Deer", 
                                  name3 = "Hunters", nboot = 10000, dhat = "Dhat1")        # bear-md: Hunter present bear <50, hunter absent all >50
  bear_md_hunt_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                         spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                         spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                         name1 = "Black Bear", name2 = "Mule Deer", 
                                         name3 = "Hunters", nboot = 10000, dhat = "Dhat4") # bear-md: Hunter present bear <50, hunter absent all >50
  bear_elk_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                   spp2 = filter(hunting_first, Species == "Elk"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Black Bear", name2 = "Elk", 
                                   name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  bear_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                    spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                    spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                    name1 = "Black Bear", name2 = "White-tailed Deer", 
                                    name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  bear_moose_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                   spp2 = filter(hunting_first, Species == "Moose"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Black Bear", name2 = "Moose", 
                                   name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  bob_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                  spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                  spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                  name1 = "Bobcat", name2 = "Mule Deer", 
                                  name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  bob_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                   spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Bobcat", name2 = "White-tailed Deer", 
                                   name3 = "Hunters", nboot = 10000, dhat = "Dhat1")
  coy_md_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                  spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                  spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                  name1 = "Coyote", name2 = "Mule Deer", 
                                  name3 = "Hunters", nboot = 10000, dhat = "Dhat4")
  coy_wtd_hunt_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                   spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                   spp3 = filter(hunting_first, HumanActivity == "Hunter"), 
                                   name1 = "Coyote", name2 = "White-tailed Deer", 
                                   name3 = "Hunters", nboot = 10000, dhat = "Dhat4")
  
  #'  Bundle and save!
  pred_prey_hunt_overlap <- list(coug_md_hunt_over, coug_elk_hunt_over, coug_wtd_hunt_over, coug_moose_hunt_over,
                                 wolf_md_hunt_over, wolf_wtd_hunt_over, wolf_moose_hunt_over, #wolf_elk_hunt_over, 
                                 bear_md_hunt_over_dhat1, bear_md_hunt_over_dhat4, 
                                 bear_elk_hunt_over, bear_wtd_hunt_over, bear_moose_hunt_over,
                                 bob_md_hunt_over, bob_wtd_hunt_over, coy_md_hunt_over, coy_wtd_hunt_over)
  save(pred_prey_hunt_overlap, file = paste0("./Outputs/Temporal Overlap/pred_prey_hunt_overlap_", Sys.Date(), ".RData"))
  
  
  ####  Predator-Prey Overlap Public vs Private Land  ####
  #'  Estimate temporal overlap between predators and prey on public vs private land
  # coug_md_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"),
  #                                        spp2 = filter(hunting_first, Species == "Mule Deer"),
  #                                        spp3 = filter(hunting_first, Public1 == 1),
  #                                        name1 = "Cougar", name2 = "Mule Deer",
  #                                        name3 = "Public lands", nboot = 10000, dhat = "Dhat1") #, i = 1
  coug_elk_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                          spp2 = filter(hunting_first, Species == "Elk"), 
                                          spp3 = filter(hunting_first, Public1 == 1), 
                                          name1 = "Cougar", name2 = "Elk", 
                                          name3 = "Public lands", nboot = 10000, dhat = "Dhat1")
  coug_wtd_public_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                          spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                          spp3 = filter(hunting_first, Public1 == 1), 
                                          name1 = "Cougar", name2 = "White-tailed Deer", 
                                          name3 = "Public lands", nboot = 10000, dhat = "Dhat1")    # coug-wtd: public all >50, private cougar <50
  coug_wtd_public_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                            spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                            spp3 = filter(hunting_first, Public1 == 1), 
                                            name1 = "Cougar", name2 = "White-tailed Deer", 
                                            name3 = "Public lands", nboot = 10000, dhat = "Dhat4")  # coug-wtd: public all >50, private cougar <50
  coug_moose_public_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                            spp2 = filter(hunting_first, Species == "Moose"), 
                                            spp3 = filter(hunting_first, Public1 == 1), 
                                            name1 = "Cougar", name2 = "Moose", 
                                            name3 = "Public lands", nboot = 10000, dhat = "Dhat1")    # coug-moose: public all >50, private all <50
  coug_moose_public_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Cougar"), 
                                              spp2 = filter(hunting_first, Species == "Moose"), 
                                              spp3 = filter(hunting_first, Public1 == 1), 
                                              name1 = "Cougar", name2 = "Moose", 
                                              name3 = "Public lands", nboot = 10000, dhat = "Dhat4")  # coug-moose: public all >50, private all <50
  # wolf_md_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
  #                                        spp2 = filter(hunting_first, Species == "Mule Deer"), 
  #                                        spp3 = filter(hunting_first, Public1 == 1), 
  #                                        name1 = "Wolf", name2 = "Mule Deer", 
  #                                        name3 = "Public lands", nboot = 10000, dhat = "Dhat1")
  # wolf_elk_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
  #                                         spp2 = filter(hunting_first, Species == "Elk"), 
  #                                         spp3 = filter(hunting_first, Public1 == 1), 
  #                                         name1 = "Wolf", name2 = "Elk", 
  #                                         name3 = "Public lands", nboot = 10000, dhat = "Dhat1")
  # wolf_wtd_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
  #                                         spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
  #                                         spp3 = filter(hunting_first, Public1 == 1), 
  #                                         name1 = "Wolf", name2 = "White-tailed Deer", 
  #                                         name3 = "Public lands", nboot = 10000, dhat = "Dhat1")
  # wolf_moose_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
  #                                           spp2 = filter(hunting_first, Species == "Moose"), 
  #                                           spp3 = filter(hunting_first, Public1 == 1), 
  #                                           name1 = "Wolf", name2 = "Moose", 
  #                                           name3 = "Public lands", nboot = 10000, dhat = "Dhat1")
  bear_md_public_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                         spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                         spp3 = filter(hunting_first, Public1 == 1), 
                                         name1 = "Black Bear", name2 = "Mule Deer", 
                                         name3 = "Public lands", nboot = 10000, dhat = "Dhat1")    # bear-md: public all >50, private all <50
  bear_md_public_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                                 spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                                 spp3 = filter(hunting_first, Public1 == 1), 
                                                 name1 = "Black Bear", name2 = "Mule Deer", 
                                                 name3 = "Public lands", nboot = 10000, dhat = "Dhat4")    # bear-md: public all >50, private all <50
  bear_elk_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                          spp2 = filter(hunting_first, Species == "Elk"), 
                                          spp3 = filter(hunting_first, Public1 == 1), 
                                          name1 = "Black Bear", name2 = "Elk", 
                                          name3 = "Public lands", nboot = 10000, dhat = "Dhat1")
  bear_wtd_public_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                          spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                          spp3 = filter(hunting_first, Public1 == 1), 
                                          name1 = "Black Bear", name2 = "White-tailed Deer", 
                                          name3 = "Public lands", nboot = 10000, dhat = "Dhat1")  # bear-wtd: public all >50, private bear <50
  bear_wtd_public_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                                  spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                                  spp3 = filter(hunting_first, Public1 == 1), 
                                                  name1 = "Black Bear", name2 = "White-tailed Deer", 
                                                  name3 = "Public lands", nboot = 10000, dhat = "Dhat4")  # bear-wtd: public all >50, private bear <50
  # bear_moose_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Black Bear"),
  #                                           spp2 = filter(hunting_first, Species == "Moose"),
  #                                           spp3 = filter(hunting_first, Public1 == 1),
  #                                           name1 = "Black Bear", name2 = "Moose",
  #                                           name3 = "Public lands", nboot = 100, dhat = "Dhat1")
  bob_md_public_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                        spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                        spp3 = filter(hunting_first, Public1 == 1), 
                                        name1 = "Bobcat", name2 = "Mule Deer", 
                                        name3 = "Public lands", nboot = 10000, dhat = "Dhat1")  # bob-md: public all >50, private all <50
  bob_md_public_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                          spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                          spp3 = filter(hunting_first, Public1 == 1), 
                                          name1 = "Bobcat", name2 = "Mule Deer", 
                                          name3 = "Public lands", nboot = 10000, dhat = "Dhat4")  # bob-md: public all >50, private all <50
  bob_wtd_public_over_dhat1 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                         spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                         spp3 = filter(hunting_first, Public1 == 1), 
                                         name1 = "Bobcat", name2 = "White-tailed Deer", 
                                         name3 = "Public lands", nboot = 10000, dhat = "Dhat1")   # bob-wtd: public all >50, private bobcat <50
  bob_wtd_public_over_dhat4 <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                           spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                           spp3 = filter(hunting_first, Public1 == 1), 
                                           name1 = "Bobcat", name2 = "White-tailed Deer", 
                                           name3 = "Public lands", nboot = 10000, dhat = "Dhat4") # bob-wtd: public all >50, private bobcat <50
  coy_md_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                        spp2 = filter(hunting_first, Species == "Mule Deer"), 
                                        spp3 = filter(hunting_first, Public1 == 1), 
                                        name1 = "Coyote", name2 = "Mule Deer", 
                                        name3 = "Public lands", nboot = 10000, dhat = "Dhat4")
  coy_wtd_public_over <- pred_prey_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                         spp2 = filter(hunting_first, Species == "White-tailed Deer"), 
                                         spp3 = filter(hunting_first, Public1 == 1), 
                                         name1 = "Coyote", name2 = "White-tailed Deer", 
                                         name3 = "Public lands", nboot = 10000, dhat = "Dhat4")
  
  #'  Save
  pred_prey_public_overlap <- list(coug_elk_public_over, coug_wtd_public_over_dhat1, 
                                   coug_wtd_public_over_dhat4, coug_moose_public_over_dhat1, coug_moose_public_over_dhat4,#coug_md_public_over,
                                   #wolf_md_public_over, wolf_wtd_public_over, wolf_moose_public_over, #wolf_elk_public_over, 
                                   bear_md_public_over_dhat1, bear_md_public_over_dhat4, 
                                   bear_elk_public_over, bear_wtd_public_over_dhat1, bear_wtd_public_over_dhat4, #bear_moose_public_over,
                                   bob_md_public_over_dhat1, bob_md_public_over_dhat4, 
                                   bob_wtd_public_over_dhat1, bob_wtd_public_over_dhat4, 
                                   coy_md_public_over, coy_wtd_public_over)
  save(pred_prey_public_overlap, file = paste0("./Outputs/Temporal Overlap/pred_prey_public_overlap_", Sys.Date(), ".RData"))

  
  ####  Result tables & figures  ####
  #'  ---------------------------
  #load("./Outputs/Temporal Overlap/pred_prey_graze_overlap_2022-07-01.RData")
  load("./Outputs/Temporal Overlap/pred_prey_graze_overlap_OK_2022-08-02.RData")
  load("./Outputs/Temporal Overlap/pred_prey_allot_overlap_OK_2022-08-02.RData")
  load("./Outputs/Temporal Overlap/pred_prey_hunt_overlap_2022-08-02.RData")
  load("./Outputs/Temporal Overlap/pred_prey_public_overlap_2022-08-04.RData")
  
  #'  Create results tables from overlap estimates
  #'  Indexing [2,i] gives norm0 CIs, [4,i] gives basic0 CIs
  results_table <- function(overlap_out, spp1, spp2, spp3) {
    dhat_spp3.present <- round(overlap_out[[1]], 2)
    spp3.present_lci <- round(overlap_out[[5]][4,1], 2)
    spp3.present_uci <- round(overlap_out[[5]][4,2], 2)
    dhat_spp3.absent <- round(overlap_out[[2]], 2)
    spp3.absent_lci <- round(overlap_out[[6]][4,1], 2)
    spp3.absent_uci <- round(overlap_out[[6]][4,2], 2)
    pair <- paste0(spp1, "-", spp2)
    spp <- c(pair, pair)
    predator <- c(spp1, spp1)
    prey <- c(spp2, spp2)
    activity <- c("Detected", "Not detected")
    Dhat <- c(dhat_spp3.present, dhat_spp3.absent)
    l95 <- c(spp3.present_lci, spp3.absent_lci)
    u95 <- c(spp3.present_uci, spp3.absent_uci)
    ndet_predator <- c(overlap_out[[7]], overlap_out[[9]])
    ndet_prey <- c(overlap_out[[8]], overlap_out[[10]])
    df <- as.data.frame(cbind(spp, predator, prey, activity, Dhat, l95, u95, 
                              ndet_predator, ndet_prey))
    rownames(df) <- NULL
    names(df)[names(df) == "spp"] <- "Species.pair"
    names(df)[names(df) == "activity"] <- paste0(spp3, ".activity")
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  Grazing results: Okanogan study area only
  #'  When using two different overlap estimators for same species pairing (e.g.,
  #'  when presence vs absence sample sizes requires different overlap estimators), 
  #'  need to filter to the appropriate estimator given sample size- 
  #'  dhat1 for when cattle are NOT DETECTED, dhat4 for when cattle are DETECTED
  #'  Note: overlap estimate when cattle are detected is always first row,
  #'  overlap estimate when cattle are not detected is always second row 
  coug_md_graze_out <- results_table(pred_prey_graze_overlap[[1]], spp1 = "Cougar", spp2 = "Mule Deer", spp3 = "Grazing")
  coug_wtd_graze_out <- results_table(pred_prey_graze_overlap[[2]], spp1 = "Cougar", spp2 = "White-tailed Deer", spp3 = "Grazing")
  coug_moose_graze_out <- results_table(pred_prey_graze_overlap[[3]], spp1 = "Cougar", spp2 = "Moose", spp3 = "Grazing")
  bear_md_graze_out <- results_table(pred_prey_graze_overlap[[4]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "Grazing")
  bear_wtd_graze_out <- results_table(pred_prey_graze_overlap[[5]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "Grazing")
  #'  bear-moose: cattle present both n>50, cattle absent moose n<50
  bear_moose_graze_out1 <- results_table(pred_prey_graze_overlap[[6]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Grazing")
  bear_moose_graze_out1 <- bear_moose_graze_out1[2,]
  bear_moose_graze_out4 <- results_table(pred_prey_graze_overlap[[7]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Grazing")
  bear_moose_graze_out4 <- bear_moose_graze_out4[1,]
  #'  bob-md: cattle present both n>50, cattle absent bobcat n<50
  bob_md_graze_out1 <- results_table(pred_prey_graze_overlap[[8]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Grazing")
  bob_md_graze_out1 <- bob_md_graze_out1[2,]
  bob_md_graze_out4 <- results_table(pred_prey_graze_overlap[[9]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Grazing")
  bob_md_graze_out4 <- bob_md_graze_out4[1,]
  coy_md_graze_out <- results_table(pred_prey_graze_overlap[[10]], spp1 = "Coyote", spp2 = "Mule Deer", spp3 = "Grazing")
  coy_wtd_graze_out <- results_table(pred_prey_graze_overlap[[11]], spp1 = "Coyote", spp2 = "White-tailed Deer", spp3 = "Grazing")
  
  cattle_overlap_tbl <- rbind(coug_md_graze_out, coug_wtd_graze_out, coug_moose_graze_out,
                              bear_md_graze_out, bear_wtd_graze_out, bear_moose_graze_out4, 
                              bear_moose_graze_out1, bob_md_graze_out4, bob_md_graze_out1, 
                              coy_md_graze_out, coy_wtd_graze_out)
  # write.csv(cattle_overlap_tbl, file = paste0("./Outputs/Temporal Overlap/pred-prey_graze_overlap_", Sys.Date(), ".csv"))
  
  #'  Grazing Allotment results: Okanogan study area only
  coug_md_allot_out <- results_table(pred_prey_allot_overlap[[1]], spp1 = "Cougar", spp2 = "Mule Deer", spp3 = "Allotment")
  coug_wtd_allot_out <- results_table(pred_prey_allot_overlap[[2]], spp1 = "Cougar", spp2 = "White-tailed Deer", spp3 = "Allotment")
  coug_moose_allot_out <- results_table(pred_prey_allot_overlap[[3]], spp1 = "Cougar", spp2 = "Moose", spp3 = "Allotment")
  bear_md_allot_out <- results_table(pred_prey_allot_overlap[[4]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "Allotment")
  #'  bear-wtd: cattle present both n>50, cattle absent bear n<50
  bear_wtd_allot_out1 <- results_table(pred_prey_allot_overlap[[5]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "Allotment")
  bear_wtd_allot_out1 <- bear_wtd_allot_out1[2,]
  bear_wtd_allot_out4 <- results_table(pred_prey_allot_overlap[[6]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "Allotment")
  bear_wtd_allot_out4 <- bear_wtd_allot_out4[1,]
  #'  bear-moose: cattle present both n>50, cattle absent both n<50
  bear_moose_allot_out1 <- results_table(pred_prey_allot_overlap[[7]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Allotment")
  bear_moose_allot_out1 <- bear_moose_allot_out1[2,]
  bear_moose_allot_out4 <- results_table(pred_prey_allot_overlap[[8]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Allotment")
  bear_moose_allot_out4 <- bear_moose_allot_out4[1,]
  #' bob-md: cattle present both n>50, cattle absent bobcat n<50
  bob_md_allot_out1 <- results_table(pred_prey_allot_overlap[[9]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Allotment")
  bob_md_allot_out1 <- bob_md_allot_out1[2,]
  bob_md_allot_out4 <- results_table(pred_prey_allot_overlap[[10]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Allotment")
  bob_md_allot_out4 <- bob_md_allot_out4[1,]
  coy_md_allot_out <- results_table(pred_prey_allot_overlap[[11]], spp1 = "Coyote", spp2 = "Mule Deer", spp3 = "Allotment")
  coy_wtd_allot_out <- results_table(pred_prey_allot_overlap[[12]], spp1 = "Coyote", spp2 = "White-tailed Deer", spp3 = "Allotment")
  
  allotment_overlap_tbl <- rbind(coug_md_allot_out, coug_wtd_allot_out, coug_moose_allot_out,
                                 bear_md_allot_out, bear_wtd_allot_out4, bear_wtd_allot_out1, 
                                 bear_moose_allot_out4, bear_moose_allot_out1, bob_md_allot_out4, 
                                 bob_md_allot_out1, coy_md_allot_out, coy_wtd_allot_out) %>%
    mutate(Public_Grazing = ifelse(Allotment.activity == "Detected", "Permitted", "Not permitted")) %>%
    dplyr::select(-Allotment.activity)
  # write.csv(allotment_overlap_tbl, file = paste0("./Outputs/Temporal Overlap/pred-prey_allot_overlap_", Sys.Date(), ".csv"))
  
  
  #'  Hunting results: NE & OK study areas
  coug_md_hunt_out <- results_table(pred_prey_hunt_overlap[[1]], spp1 = "Cougar", spp2 = "Mule Deer", spp3 = "Hunter")
  coug_elk_hunt_out <- results_table(pred_prey_hunt_overlap[[2]], spp1 = "Cougar", spp2 = "Elk", spp3 = "Hunter")
  coug_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[3]], spp1 = "Cougar", spp2 = "White-tailed Deer", spp3 = "Hunter")
  coug_moose_hunt_out <- results_table(pred_prey_hunt_overlap[[4]], spp1 = "Cougar", spp2 = "Moose", spp3 = "Hunter")
  wolf_md_hunt_out <- results_table(pred_prey_hunt_overlap[[5]], spp1 = "Wolf", spp2 = "Mule Deer", spp3 = "Hunter")
  wolf_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[6]], spp1 = "Wolf", spp2 = "White-tailed Deer", spp3 = "Hunter")
  wolf_moose_hunt_out <- results_table(pred_prey_hunt_overlap[[7]], spp1 = "Wolf", spp2 = "Moose", spp3 = "Hunter")
  # bear-md: hunter present bear n<50, hunter absent both n>50
  bear_md_hunt_out1 <- results_table(pred_prey_hunt_overlap[[8]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "Hunter")
  bear_md_hunt_out1 <- bear_md_hunt_out1[1,]
  bear_md_hunt_out4 <- results_table(pred_prey_hunt_overlap[[9]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "Hunter")
  bear_md_hunt_out4 <- bear_md_hunt_out4[2,]
  bear_elk_hunt_out <- results_table(pred_prey_hunt_overlap[[10]], spp1 = "Black bear", spp2 = "Elk", spp3 = "Hunter")
  bear_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[11]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "Hunter")
  bear_moose_hunt_out <- results_table(pred_prey_hunt_overlap[[12]], spp1 = "Black bear", spp2 = "Moose", spp3 = "Hunter")
  bob_md_hunt_out <- results_table(pred_prey_hunt_overlap[[13]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "Hunter")
  bob_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[14]], spp1 = "Bobcat", spp2 = "White-tailed Deer", spp3 = "Hunter")
  coy_md_hunt_out <- results_table(pred_prey_hunt_overlap[[15]], spp1 = "Coyote", spp2 = "Mule Deer", spp3 = "Hunter")
  coy_wtd_hunt_out <- results_table(pred_prey_hunt_overlap[[16]], spp1 = "Coyote", spp2 = "White-tailed Deer", spp3 = "Hunter")
  
  hunter_overlap_tbl <- rbind(coug_md_hunt_out, coug_elk_hunt_out, coug_wtd_hunt_out, coug_moose_hunt_out,  
                              #wolf_md_hunt_out, wolf_wtd_hunt_out, wolf_moose_hunt_out, 
                              bear_md_hunt_out4, bear_md_hunt_out1, 
                              bear_elk_hunt_out, bear_wtd_hunt_out, bear_moose_hunt_out, 
                              bob_md_hunt_out, bob_wtd_hunt_out, coy_md_hunt_out, coy_wtd_hunt_out)
  # write.csv(hunter_overlap_tbl, file = paste0("./Outputs/Temporal Overlap/pred-prey_hunt_overlap_", Sys.Date(), ".csv"))
  
  #'  Public vs private lands results: NE & OK study areas
  coug_elk_public_out <- results_table(pred_prey_public_overlap[[1]], spp1 = "Cougar", spp2 = "Elk", spp3 = "PublicLand")
  # coug-wtd: public both n>50, private cougar n<50
  coug_wtd_public_out1 <- results_table(pred_prey_public_overlap[[2]], spp1 = "Cougar", spp2 = "White-tailed Deer", spp3 = "PublicLand")
  coug_wtd_public_out1 <- coug_wtd_public_out1[2,]
  coug_wtd_public_out4 <- results_table(pred_prey_public_overlap[[3]], spp1 = "Cougar", spp2 = "White-tailed Deer", spp3 = "PublicLand")
  coug_wtd_public_out4 <- coug_wtd_public_out4[1,]
  # coug-moose: public both n>50, private both n<50
  coug_moose_public_out1 <- results_table(pred_prey_public_overlap[[4]], spp1 = "Cougar", spp2 = "Moose", spp3 = "PublicLand")
  coug_moose_public_out1 <- coug_moose_public_out1[2,]
  coug_moose_public_out4 <- results_table(pred_prey_public_overlap[[5]], spp1 = "Cougar", spp2 = "Moose", spp3 = "PublicLand")
  coug_moose_public_out4 <- coug_moose_public_out4[1,]
  # bear-md: public both n>50, private both n<50
  bear_md_public_out1 <- results_table(pred_prey_public_overlap[[6]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "PublicLand")
  bear_md_public_out1 <- bear_md_public_out1[2,]
  bear_md_public_out4 <- results_table(pred_prey_public_overlap[[7]], spp1 = "Black bear", spp2 = "Mule Deer", spp3 = "PublicLand")
  bear_md_public_out4 <- bear_md_public_out4[1,]
  bear_elk_public_out <- results_table(pred_prey_public_overlap[[8]], spp1 = "Black bear", spp2 = "Elk", spp3 = "PublicLand")
  # bear-wtd: public both n>50, private bear n<50
  bear_wtd_public_out1 <- results_table(pred_prey_public_overlap[[9]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "PublicLand")
  bear_wtd_public_out1 <- bear_wtd_public_out1[2,]
  bear_wtd_public_out4 <- results_table(pred_prey_public_overlap[[10]], spp1 = "Black bear", spp2 = "White-tailed Deer", spp3 = "PublicLand")
  bear_wtd_public_out4 <- bear_wtd_public_out4[1,]
  # bob-md: public both n>50, private both n<50
  bob_md_public_out1 <- results_table(pred_prey_public_overlap[[11]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "PublicLand")
  bob_md_public_out1 <- bob_md_public_out1[2,]
  bob_md_public_out4 <- results_table(pred_prey_public_overlap[[12]], spp1 = "Bobcat", spp2 = "Mule Deer", spp3 = "PublicLand")
  bob_md_public_out4 <- bob_md_public_out4[1,]
  # bob-wtd: public both n>50, private bobcat n<50
  bob_wtd_public_out1 <- results_table(pred_prey_public_overlap[[13]], spp1 = "Bobcat", spp2 = "White-tailed Deer", spp3 = "PublicLand")
  bob_wtd_public_out1 <- bob_wtd_public_out1[2,]
  bob_wtd_public_out4 <- results_table(pred_prey_public_overlap[[14]], spp1 = "Bobcat", spp2 = "White-tailed Deer", spp3 = "PublicLand")
  bob_wtd_public_out4 <- bob_wtd_public_out4[1,]
  coy_md_public_out <- results_table(pred_prey_public_overlap[[15]], spp1 = "Coyote", spp2 = "Mule Deer", spp3 = "PublicLand")
  coy_wtd_public_out <- results_table(pred_prey_public_overlap[[16]], spp1 = "Coyote", spp2 = "White-tailed Deer", spp3 = "PublicLand")
  
  public_overlap_tbl <- rbind(coug_elk_public_out, coug_wtd_public_out4, coug_wtd_public_out1, 
                              coug_moose_public_out4, coug_moose_public_out1, bear_md_public_out4, 
                              bear_md_public_out1, bear_elk_public_out, bear_wtd_public_out4,
                              bear_wtd_public_out1, bob_md_public_out4, bob_md_public_out1, bob_wtd_public_out4, 
                              bob_wtd_public_out1, coy_md_public_out, coy_wtd_public_out) %>%
    mutate(PublicLand.activity = ifelse(PublicLand.activity == "Detected", "Public/Timber", "Private"))
  # write.csv(public_overlap_tbl, file = paste0("./Outputs/Temporal Overlap/pred-prey_public_overlap_", Sys.Date(), ".csv"))
  
  
  ####  Figures for visualization  ####
  #'  -----------------------------
  #' #'  Split up data into species-specific groups 
  #' coug_cattle <- filter(cattle_overlap_tbl, grepl("Cougar", `Species.pair`))
  #' wolf_cattle <- filter(cattle_overlap_tbl, grepl("Wolf", `Species.pair`))
  #' bear_cattle <- filter(cattle_overlap_tbl, grepl("Black bear", `Species.pair`))
  #' bob_cattle <- filter(cattle_overlap_tbl, grepl("Bobcat", `Species.pair`))
  #' coy_cattle <- filter(cattle_overlap_tbl, grepl("Coyote", `Species.pair`))
  #' 
  #' #'  Plot it
  #' ggplot(coug_cattle, aes(x = `Species.pair`, y = Dhat, group = Dhat)) +  
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = `Species.pair`), width = 0.2, position = position_dodge(width=0.3)) +
  #'   geom_point(stat = 'identity', aes(col = `Species.pair`, shape = `Grazing.activity`), size = 3.5,position = position_dodge(width=0.3)) + 
  #'   ylim(0,1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + #theme(legend.position = "none") +
  #'   guides(color = "none") +
  #'   ggtitle("Coefficient of overlap when grazing is and is not detected")
  #' ggplot(wolf_cattle, aes(x = `Species.pair`, y = Dhat, group = Dhat)) +  
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = `Species.pair`), width = 0.2, position = position_dodge(width=0.3)) +
  #'   geom_point(stat = 'identity', aes(col = `Species.pair`, shape = `Grazing.activity`), size = 3.5,position = position_dodge(width=0.3)) + 
  #'   ylim(0,1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + #theme(legend.position = "none") +
  #'   guides(color = "none") +
  #'   ggtitle("Coefficient of overlap when grazing is and is not detected")
  #' ggplot(bear_cattle, aes(x = `Species.pair`, y = Dhat, group = Dhat)) +  
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = `Species.pair`), width = 0.2, position = position_dodge(width=0.3)) +
  #'   geom_point(stat = 'identity', aes(col = `Species.pair`, shape = `Grazing.activity`), size = 3.5,position = position_dodge(width=0.3)) + 
  #'   ylim(0,1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + #theme(legend.position = "none") +
  #'   guides(color = "none") +
  #'   ggtitle("Coefficient of overlap when grazing is and is not detected")
  #' ggplot(bob_cattle, aes(x = `Species.pair`, y = Dhat, group = Dhat)) +  
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = `Species.pair`), width = 0.2, position = position_dodge(width=0.3)) +
  #'   geom_point(stat = 'identity', aes(col = `Species.pair`, shape = `Grazing.activity`), size = 3.5,position = position_dodge(width=0.3)) + 
  #'   ylim(0,1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + #theme(legend.position = "none") +
  #'   guides(color = "none") +
  #'   ggtitle("Coefficient of overlap when grazing is and is not detected")
  #' ggplot(coy_cattle, aes(x = `Species.pair`, y = Dhat, group = l95)) +  #grouping by l95 b/c dhat identical for coy-wtd
  #'   geom_errorbar(aes(ymin = l95, ymax = u95, col = `Species.pair`), width = 0.2, position = position_dodge(width=0.3)) +
  #'   geom_point(stat = 'identity', aes(col = `Species.pair`, shape = `Grazing.activity`), size = 3.5,position = position_dodge(width=0.3)) + 
  #'   ylim(0,1) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + #theme(legend.position = "none") +
  #'   guides(color = "none") +
  #'   ggtitle("Coefficient of overlap when grazing is and is not detected")
  
  #'  Make one single facet_grid plot by grouped by predator species
  overlap_grazing_effect <- ggplot(cattle_overlap_tbl, aes(x = `Species.pair`, y = Dhat, group = Grazing.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = Grazing.activity), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Grazing activity")) + 
    ggtitle("Coefficient of overlap when cattle activity is and is not detected on camera") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  overlap_grazing_effect
  ggsave(overlap_grazing_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_Grazing_Effect_Plot_OKonly.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  #'  Make one single facet_grid plot by grouped by predator species
  allotment_overlap_tbl$Public_Grazing <- factor(allotment_overlap_tbl$Public_Grazing, levels = c("Permitted", "Not permitted"))
  overlap_allot_effect <- ggplot(allotment_overlap_tbl, aes(x = `Species.pair`, y = Dhat, group = Public_Grazing)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = Public_Grazing), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Public grazing")) + 
    ggtitle("Coefficient of overlap on and off public grazing allotments") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  overlap_allot_effect
  ggsave(overlap_allot_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_Allotment_Effect_Plot_OKonly.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  #'  Make one single facet_grid plot by grouped by predator species
  overlap_hunting_effect <- ggplot(hunter_overlap_tbl, aes(x = `Species.pair`, y = Dhat, group = Hunter.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = Hunter.activity), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Hunter activity")) + 
    ggtitle("Coefficient of overlap when hunters are and are not detected") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  overlap_hunting_effect
  ggsave(overlap_hunting_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_Hunter_Effect_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  #'  Make one single facet_grid plot by grouped by predator species
  overlap_public_effect <- ggplot(public_overlap_tbl, aes(x = `Species.pair`, y = Dhat, group = PublicLand.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = PublicLand.activity), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Property Ownership")) + 
    ggtitle("Coefficient of overlap on public vs private land") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  overlap_public_effect
  ggsave(overlap_public_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_PublicVPrivate_Effect_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  
  ####  Single species temporal overlap  ####
  #'  -----------------------------------
  #'  Function to estimate differences in temporal activity for a species at camera  
  #'  sites where cattle/humans are present (detected) vs absent (not detected).
  #'  Using spp3 to also represent cameras on public vs private land in this function.
  spp_overlap <- function(spp1, spp2, name1, name2, nboot, dhat) { #, i
    #'  Create logical vectors (T/F) indicating which cameras spp1 was detected 
    #'  at with and without spp2
    both.present <- spp1$CameraLocation %in% spp2$CameraLocation
    spp1_dat <- cbind(spp1, both.present) 
    #'  Split out data into camera locations where both are present vs spp2 absent
    spp1_spp2.present <- spp1_dat[spp1_dat$both.present == T,]
    spp1_spp2.absent <- spp1_dat[spp1_dat$both.present == F,]
    
    #'  Review sample size per species- smaller sample will determine which coefficient
    #'  of overlap estimator to use (delta1 for small samples [<50 detection events], 
    #'  delta4 for larger samples [>75 detection events])
    ndet_spp2.present <- nrow(spp1_spp2.present)
    ndet_spp2.absent <- nrow(spp1_spp2.absent)
    print(ndet_spp2.present); print(ndet_spp2.absent)
    
    #'  Visualize general temporal activity with density plots
    densityPlot(spp1$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name1, " Daily Activity"))
    densityPlot(spp2$sunTime, rug = T, col = "blue", main = paste0("Density Plot of ", name2, " Daily Activity"))
    
    #'  Visualize temporal overlap
    overlapPlot(spp1$sunTime, spp2$sunTime, rug = F, linet = c(1, 1), 
                linec = c("red", "blue"), linew = c(2, 2), main = paste0("Overlap Plots of ", name1, " (red) and ", name2, " (blue)"))
    
    overlapPlot(spp1_spp2.present$sunTime, spp1_spp2.absent$sunTime, rug = F, linet = c(1, 1), 
                linec = c("red", "blue"), linew = c(2, 2),  
                main = paste0("Overlap Plots of ", name1, " \nwhen ", name2, " are Present (red) and Absente (blue)"))
    
    #'  Calculate coefficient of overlap
    dhats_spp1.spp2 <- overlapEst(A = spp1_spp2.present$sunTime, 
                                  B = spp1_spp2.absent$sunTime, type = dhat) 
    
    #'  Bootstrap to estimate standard errors
    #'  FYI: smooth = TRUE is default and allows bootstrap to randomly sample from 
    #'  a distribution of times that have a wider range than the original sample
    #'  (see pg. 5 in Overlap package vignette for details) 
    spp1.spp2.boot <- bootstrap(spp1_spp2.present$sunTime, spp1_spp2.absent$sunTime, 
                                 nboot, smooth = TRUE, type = dhat)  
    #'  Bootstrap mean will be a little different then detla coefficient due to
    #'  bootstrap bias (BSmean - delta) that needs to be accounted for in 95% CIs
    BSmean <- mean(spp1.spp2.boot)
    
    #'  Bootstrap 95% Confidence Intervals
    #'  norm0 uses the standard deviation of bootstrap results to calculate CI (delta +/- 1.96*SDboot)
    #'  basic0 takes the 2.5% and 97.5% percentiles and adjusts based on BS bias (percentile - BSbias)
    #'  If sampling distribution is normal, norm0 and basic0 should be similar;
    #'  if sampling distribution is skewed (i.e., if delta is close to 0 or 1) then
    #'  basic0 is the better estimator
    #'  Using bootCIlogit instead of bootCI so that bias corrections are done on
    #'  the logit scale, then backtransformed. Without this, 95% CIs can fall
    #'  outside (0, 1) interval. See Overlap vignette for more details.
    CI <- bootCIlogit(dhats_spp1.spp2, spp1.spp2.boot) #dhats_spp1.spp2[i]
    
    watson.two.test
    
    #'  Print results
    #'  Effect of spp2 being present
    print("Overlap coefficients when spp2 is present"); print(dhats_spp1.spp2)
    print("Bootstrap mean"); print(BSmean)
    print("Bootstrap 95% CI"); print(CI)
    
    #'  Save as a giant list
    overlap_list <- list(dhats_spp1.spp2, spp1.spp2.boot, CI, ndet_spp2.present, ndet_spp2.absent)
    names(overlap_list) <- c("dhats_spp1.spp2", "spp1.spp2.boot", "CI", "ndet_spp2.present", "ndet_spp2.absent")
    
    return(overlap_list)
  }
  #'  Estimate temporal overlap for a species when cattle are/aren't detected
  #'  THINK ABOUT focusing on only OK study area since big difference in number of
  #'  cameras with cattle in NE vs OK, pooling across study areas could be confounding
  #'  any apparent temporal patterns
  ####  Single-Species Overlap Grazing Season  ####
  coug_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
                                 spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                 name1 = "Cougar", name2 = "Cattle", nboot = 10000, dhat = "Dhat1") #i = 1
  # wolf_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
  #                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
  #                                name1 = "Wolf", name2 = "Cattle", nboot = 10000, dhat = "Dhat1") # only 2 observations when cattle are present
  bear_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                 spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                 name1 = "Black bear", name2 = "Cattle", nboot = 10000, dhat = "Dhat4")
  bob_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
                                 spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                 name1 = "Bobcat", name2 = "Cattle", nboot = 10000, dhat = "Dhat4")
  coy_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"), 
                                 spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                 name1 = "Coyote", name2 = "Cattle", nboot = 10000, dhat = "Dhat4")
  md_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Mule Deer"), 
                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                name1 = "Mule deer", name2 = "Cattle", nboot = 10000, dhat = "Dhat4")
  # elk_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Elk"), 
  #                               spp2 = filter(grazing_first_OK, Species == "Cattle"), 
  #                               name1 = "Elk", name2 = "Cattle", nboot = 10000, dhat = "Dhat1")  # only 4 observations where cattle absent - OK study area!
  wtd_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                name1 = "White-tailed deer", name2 = "Cattle", nboot = 10000, dhat = "Dhat4")
  moose_graze_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Moose"), 
                                spp2 = filter(grazing_first_OK, Species == "Cattle"), 
                                name1 = "Moose", name2 = "Cattle", nboot = 10000, dhat = "Dhat4")
  ####  Single-Species Overlap Grazing Season on/off Allotments  ####
  coug_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Cougar"),
                                 spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                                 name1 = "Cougar", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat1") #i = 1
  wolf_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Wolf"), 
                                 spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                                 name1 = "Wolf", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat1") # only 9 observations when off allotments
  bear_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Black Bear"), 
                                 spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                                 name1 = "Black bear", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat4")
  bob_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Bobcat"), 
                                spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                                name1 = "Bobcat", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat1") 
  coy_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Coyote"), 
                                spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                                name1 = "Coyote", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat4")
  md_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Mule Deer"), 
                               spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                               name1 = "Mule deer", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat4")
  # elk_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Elk"), 
  #                               spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
  #                               name1 = "Elk", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat1")  # only 1 observation where cattle absent - OK study area!
  wtd_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "White-tailed Deer"), 
                                spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                                name1 = "White-tailed deer", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat4")
  moose_allot_over <- spp_overlap(spp1 = filter(grazing_first_OK, Species == "Moose"), 
                                  spp2 = filter(grazing_first_OK, PublicGrazing == 1), 
                                  name1 = "Moose", name2 = "Grazing Allotments", nboot = 10000, dhat = "Dhat1")
  ####  Single-Species Overlap Hunting Season  ####
  #'  Estimate temporal overlap for a species when hunters are/are not present
  coug_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Cougar"),
                                 spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                 name1 = "Cougar", name2 = "Hunters", nboot = 10000, dhat = "Dhat4") #i = 1
  wolf_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Wolf"), 
                                 spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                 name1 = "Wolf", name2 = "Hunters", nboot = 10000, dhat = "Dhat1")
  bear_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                 spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                 name1 = "Black bear", name2 = "Hunters", nboot = 10000, dhat = "Dhat4")
  bob_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                                spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                name1 = "Bobcat", name2 = "Hunters", nboot = 10000, dhat = "Dhat4")
  coy_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                                spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                name1 = "Coyote", name2 = "Hunters", nboot = 10000, dhat = "Dhat4")
  md_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Mule Deer"), 
                               spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                               name1 = "Mule deer", name2 = "Hunters", nboot = 10000, dhat = "Dhat4")
  elk_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Elk"), 
                                spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                name1 = "Elk", name2 = "Hunters", nboot = 10000, dhat = "Dhat1")
  wtd_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "White-tailed Deer"), 
                                spp2 = filter(hunting_first, HumanActivity == "Hunter"), 
                                name1 = "White-tailed deer", name2 = "Hunters", nboot = 10000, dhat = "Dhat4")
  moose_hunt_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Moose"), 
                                  spp2 = filter(hunting_first,HumanActivity == "Hunter"), 
                                  name1 = "Moose", name2 = "Hunters", nboot = 10000, dhat = "Dhat4")
  ####  Single-Species Overlap Public vs Private Land  ####
  #'  Estimate temporal overlap for a species on private vs public land
  coug_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Cougar"),
                                spp2 = filter(hunting_first, Public1 == 1), 
                                name1 = "Cougar", name2 = "Public lands", nboot = 10000, dhat = "Dhat1") #i = 1
  # wolf_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Wolf"),
  #                               spp2 = filter(hunting_first, Public1 == 1),
  #                               name1 = "Wolf", name2 = "Public lands", nboot = 10000, dhat = "Dhat1")  # 0 detections on private
  bear_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Black Bear"), 
                                spp2 = filter(hunting_first, Public1 == 1), 
                                name1 = "Black bear", name2 = "Public lands", nboot = 10000, dhat = "Dhat1")
  bob_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Bobcat"), 
                               spp2 = filter(hunting_first, Public1 == 1), 
                               name1 = "Bobcat", name2 = "Public lands", nboot = 10000, dhat = "Dhat1")
  coy_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Coyote"), 
                               spp2 = filter(hunting_first, Public1 == 1), 
                               name1 = "Coyote", name2 = "Public lands", nboot = 10000, dhat = "Dhat4")
  md_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Mule Deer"), 
                              spp2 = filter(hunting_first, Public1 == 1), 
                              name1 = "Mule deer", name2 = "Public lands", nboot = 10000, dhat = "Dhat4")
  elk_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Elk"), 
                               spp2 = filter(hunting_first, Public1 == 1), 
                               name1 = "Elk", name2 = "Public lands", nboot = 10000, dhat = "Dhat1")
  wtd_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "White-tailed Deer"), 
                               spp2 = filter(hunting_first, Public1 == 1), 
                               name1 = "White-tailed deer", name2 = "Public lands", nboot = 10000, dhat = "Dhat4")
  moose_public_over <- spp_overlap(spp1 = filter(hunting_first, Species == "Moose"), 
                                 spp2 = filter(hunting_first, Public1 == 1), 
                                 name1 = "Moose", name2 = "Public lands", nboot = 10000, dhat = "Dhat1")
  
  graze_overlap <- list(coug_graze_over, bear_graze_over, bob_graze_over, #wolf_graze_over, elk_graze_over, 
                        coy_graze_over, md_graze_over, wtd_graze_over, moose_graze_over)
  save(graze_overlap, file = paste0("./Outputs/Temporal Overlap/grazing_effect_overlap_OK_", Sys.Date(), ".RData"))
  
  allot_overlap <- list(coug_allot_over, wolf_allot_over, bear_allot_over, bob_allot_over, 
                        coy_allot_over, md_allot_over, wtd_allot_over, moose_allot_over) #elk_allot_over, 
  save(allot_overlap, file = paste0("./Outputs/Temporal Overlap/allotment_effect_overlap_OK_", Sys.Date(), ".RData"))

  hunt_overlap <- list(coug_hunt_over, wolf_hunt_over, bear_hunt_over, bob_hunt_over, 
                       coy_hunt_over, md_hunt_over, elk_hunt_over, wtd_hunt_over, moose_hunt_over)
  save(hunt_overlap, file = paste0("./Outputs/Temporal Overlap/hunter_effect_overlap_", Sys.Date(), ".RData"))
  
  public_overlap <- list(coug_public_over, bear_public_over, bob_public_over, #wolf_public_over, 
                       coy_public_over, md_public_over, elk_public_over, wtd_public_over, moose_public_over)
  save(public_overlap, file = paste0("./Outputs/Temporal Overlap/public_private_effect_overlap_", Sys.Date(), ".RData"))
  
  
  ####  Figures for visualization  ####
  #'  -----------------------------
  load("./Outputs/Temporal Overlap/grazing_effect_overlap_OK_2022-08-04.RData")
  load("./Outputs/Temporal Overlap/allotment_effect_overlap_OK_2022-08-05.RData")
  load("./Outputs/Temporal Overlap/hunter_effect_overlap_2022-08-05.RData")
  load("./Outputs/Temporal Overlap/public_private_effect_overlap_2022-08-05.RData")
  
  
  #'  Create results tables from overlap estimates
  results_table <- function(overlap_out, spp1) {
    Dhat <- round(overlap_out[[1]], 2)
    l95 <- round(overlap_out[[3]][2,1], 2)
    u95 <- round(overlap_out[[3]][2,2], 2)
    Species <- spp1
    ndet_present <- overlap_out[[4]]
    ndet_absent <- overlap_out[[5]]
    df <- as.data.frame(cbind(Species, Dhat, l95, u95, ndet_present, ndet_absent))
    rownames(df) <- NULL
    df <- mutate(df, Dhat = as.numeric(Dhat),
                 l95 = as.numeric(l95),
                 u95 = as.numeric(u95))
    return(df)
  }
  #'  Coefficient of overlap for each species when grazing is and is not detected
  coug_graze_out <- results_table(graze_overlap[[1]], spp1 = "Cougar")
  # wolf_graze_out <- results_table(graze_overlap[[2]], spp1 = "Wolf")
  bear_graze_out <- results_table(graze_overlap[[2]], spp1 = "Black bear")
  bob_graze_out <- results_table(graze_overlap[[3]], spp1 = "Bobcat")
  coy_graze_out <- results_table(graze_overlap[[4]], spp1 = "Coyote")
  md_graze_out <- results_table(graze_overlap[[5]], spp1 = "Mule deer")
  # elk_graze_out <- results_table(graze_overlap[[7]], spp1 = "Elk")
  wtd_graze_out <- results_table(graze_overlap[[6]], spp1 = "White-tailed deer")
  moose_graze_out <- results_table(graze_overlap[[7]], spp1 = "Moose")
  #'  Coefficient of overlap for each species on and off grazing allotments
  coug_allot_out <- results_table(allot_overlap[[1]], spp1 = "Cougar")
  wolf_allot_out <- results_table(allot_overlap[[2]], spp1 = "Wolf")
  bear_allot_out <- results_table(allot_overlap[[3]], spp1 = "Black bear")
  bob_allot_out <- results_table(allot_overlap[[4]], spp1 = "Bobcat")
  coy_allot_out <- results_table(allot_overlap[[5]], spp1 = "Coyote")
  md_allot_out <- results_table(allot_overlap[[6]], spp1 = "Mule deer")
  # elk_allot_out <- results_table(allot_overlap[[7]], spp1 = "Elk")
  wtd_allot_out <- results_table(allot_overlap[[7]], spp1 = "White-tailed deer")
  moose_allot_out <- results_table(allot_overlap[[8]], spp1 = "Moose")
  #'  Coefficient of overlap for each species when hunters are and are not detected
  coug_hunt_out <- results_table(hunt_overlap[[1]], spp1 = "Cougar")
  wolf_hunt_out <- results_table(hunt_overlap[[2]], spp1 = "Wolf")
  bear_hunt_out <- results_table(hunt_overlap[[3]], spp1 = "Black bear")
  bob_hunt_out <- results_table(hunt_overlap[[4]], spp1 = "Bobcat")
  coy_hunt_out <- results_table(hunt_overlap[[5]], spp1 = "Coyote")
  md_hunt_out <- results_table(hunt_overlap[[6]], spp1 = "Mule deer")
  elk_hunt_out <- results_table(hunt_overlap[[7]], spp1 = "Elk")
  wtd_hunt_out <- results_table(hunt_overlap[[8]], spp1 = "White-tailed deer")
  moose_hunt_out <- results_table(hunt_overlap[[9]], spp1 = "Moose")
  #'  Coefficient of overlap for each species on public vs private lands
  coug_public_out <- results_table(public_overlap[[1]], spp1 = "Cougar")
  #wolf_public_out <- results_table(public_overlap[[2]], spp1 = "Wolf")
  bear_public_out <- results_table(public_overlap[[2]], spp1 = "Black bear")
  bob_public_out <- results_table(public_overlap[[3]], spp1 = "Bobcat")
  coy_public_out <- results_table(public_overlap[[4]], spp1 = "Coyote")
  md_public_out <- results_table(public_overlap[[5]], spp1 = "Mule deer")
  elk_public_out <- results_table(public_overlap[[6]], spp1 = "Elk")
  wtd_public_out <- results_table(public_overlap[[7]], spp1 = "White-tailed deer")
  moose_public_out <- results_table(public_overlap[[8]], spp1 = "Moose")
  
  grazing_effects <- rbind(coug_graze_out, bear_graze_out, bob_graze_out, #wolf_graze_out, elk_graze_out, 
                           coy_graze_out, md_graze_out, wtd_graze_out, moose_graze_out)
  write.csv(grazing_effects, file = paste0("./Outputs/Temporal Overlap/graze_effect_overlap_tbl_OK_", Sys.Date(), ".csv"))
  allot_effects <- rbind(coug_allot_out, wolf_allot_out, bear_allot_out, bob_allot_out,
                         coy_allot_out, md_allot_out, elk_allot_out, wtd_allot_out, moose_allot_out)
  write.csv(allot_effects, file = paste0("./Outputs/Temporal Overlap/allotment_effect_overlap_tbl_OK_", Sys.Date(), ".csv"))
  hunter_effects <- rbind(coug_hunt_out, wolf_hunt_out, bear_hunt_out, bob_hunt_out,
                          coy_hunt_out, md_hunt_out, elk_hunt_out, wtd_hunt_out, moose_hunt_out)
  write.csv(hunter_effects, file = paste0("./Outputs/Temporal Overlap/hunter_effect_overlap_tbl_", Sys.Date(), ".csv"))
  public_effects <- rbind(coug_public_out, bear_public_out, bob_public_out, #wolf_public_out, 
                          coy_public_out, md_public_out, elk_public_out, wtd_public_out, moose_public_out) 
  write.csv(public_effects, file = paste0("./Outputs/Temporal Overlap/public_v_private_effect_overlap_tbl_", Sys.Date(), ".csv"))
  
  
  ####  Plot coefficient of overlap estimates for all species  ####
  #'  ---------------------------------------------------------
  #'  Plot coefficient of overlap estimates for each species
  # spp_overlap_grazing_plot <- ggplot(grazing_effects, aes(x = Species, y = Dhat)) +  
  #   geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0.2) +
  #   geom_point(stat = 'identity', aes(col = Species), size = 3.5) + 
  #   ylim(0,1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  #   guides(color = guide_legend(title = "Species")) + 
  #   ggtitle("Coefficient of overlap when grazing is and is not detected")
  # #ggsave(spp_overlap_grazing_plot, filename = "./Outputs/Temporal Overlap/Overlap_Grazing_Effect_Spp_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  # #ggsave(spp_overlap_grazing_plot, filename = "./Outputs/Temporal Overlap/Overlap_Grazing_Effect_Spp_Plot_OK.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  # 
  # spp_overlap_hunter_plot <- ggplot(hunter_effects, aes(x = Species, y = Dhat)) +  
  #   geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0.2) +
  #   geom_point(stat = 'identity', aes(col = Species), size = 3.5) + 
  #   ylim(0,1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  #   guides(color = guide_legend(title = "Species")) + 
  #   ggtitle("Coefficient of overlap when hunters are and are not detected")
  # #ggsave(spp_overlap_hunter_plot, filename = "./Outputs/Temporal Overlap/Overlap_Hunter_Effect_Spp_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  # 
  # spp_overlap_public_plot <- ggplot(public_effects, aes(x = Species, y = Dhat)) +  
  #   geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0.2) +
  #   geom_point(stat = 'identity', aes(col = Species), size = 3.5) + 
  #   ylim(0,1) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  #   guides(color = guide_legend(title = "Species")) + 
  #   ggtitle("Coefficient of overlap on public/timber vs private lands")
  # # ggsave(spp_overlap_public_plot, filename = "./Outputs/Temporal Overlap/Overlap_PublicVPrivate_Effect_Spp_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  
  #'  Make one single facet_grid plot by grouped by predator species
  spp_overlap_grazing_effect <- ggplot(grazing_effects, aes(x = `Species.pair`, y = Dhat, group = Grazing.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = Grazing.activity), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Grazing activity")) + 
    ggtitle("Coefficient of overlap when cattle activity is and is not detected on camera") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  spp_overlap_grazing_effect
  ggsave(spp_overlap_grazing_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_Grazing_Effect_Plot_OKonly_Spp_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  #'  Make one single facet_grid plot by grouped by predator species
  allotment_overlap_tbl$Public_Grazing <- factor(allot_effects$Public_Grazing, levels = c("Permitted", "Not permitted"))
  spp_overlap_allot_effect <- ggplot(allot_effects, aes(x = `Species.pair`, y = Dhat, group = Public_Grazing)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = Public_Grazing), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Public grazing")) + 
    ggtitle("Coefficient of overlap on and off public grazing allotments") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  spp_overlap_allot_effect
  ggsave(spp_overlap_allot_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_Allotment_Effect_Plot_OKonly_Spp_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  #'  Make one single facet_grid plot by grouped by predator species
  spp_overlap_hunting_effect <- ggplot(hunter_effects, aes(x = `Species.pair`, y = Dhat, group = Hunter.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = Hunter.activity), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Hunter activity")) + 
    ggtitle("Coefficient of overlap when hunters are and are not detected") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  spp_overlap_hunting_effect
  ggsave(spp_overlap_hunting_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_Hunter_Effect_Plot_Spp_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  #'  Make one single facet_grid plot by grouped by predator species
  spp_overlap_public_effect <- ggplot(public_effects, aes(x = `Species.pair`, y = Dhat, group = PublicLand.activity)) +   
    geom_errorbar(aes(ymin = l95, ymax = u95, col = predator), width = 0.3, position = position_dodge(width = 0.4)) +
    geom_point(stat = 'identity', aes(col = predator, shape = PublicLand.activity), size = 2.75, position = position_dodge(width = 0.4)) + 
    ylim(0,1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    guides(color = "none", shape = guide_legend(title = "Property Ownership")) + 
    ggtitle("Coefficient of overlap on public vs private land") +
    xlab("Species pairing") + ylab("Coefficient of overlap (Dhat)") +
    facet_grid(~predator, scales = "free", space = "free") 
  spp_overlap_public_effect
  ggsave(spp_overlap_public_effect, filename = "./Outputs/Temporal Overlap/Figures/Overlap_PublicVPrivate_Effect_Plot_Spp_Plot.tiff", width = 9, height = 7, dpi = 600, units = "in", device='tiff')
  
  