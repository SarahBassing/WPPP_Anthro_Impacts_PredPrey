  #'  Cattle grazing activity metrics and correlation analysis
  #'  Cameron Ho
  #'  Sarah Bassing 
  #'  Quantitative Ecology Lab
  #'  March 2022
  #'  ============================================================
  #'  Script examining possible cattle grazing activity metrics from camera trap 
  #'  data on public lands to use as covariates representing impact of grazing 
  #'  activity on wildlife occurrence & activity patterns.
  #'  Notes:
  #'  Requires sourcing functions from the 'Covariate Functions.R' Script to calculate 
  #'  different covariate metrics:
  #'    Covariate 1 - Sum of individual images in each time step of interest (daily/weekly)
  #'    Covariate 2 - Sum of number of unique detention events in each time step
  #'    Covariate 3 - Sum of maximum number of individuals detected in a single 
  #'    image per detection event in each time step
  #'    Covariate 4 - Sum of minutes across unique detection events in each time step
  #'  ============================================================
  
  #'  packages
  library(data.table)
  library(lubridate)
  library(tidyverse)
  
  #'  reading in complete cow dataset
  cow_data <- read.csv("./Data/All_cattle_detections.csv")

  
# Decommissioned code to pull and combine multiple camera csv  -------------
  #'  reading in cow data from multiple csv and combining them
  #cow_data <- list.files(path = "./Data/Cow count csv files", pattern = "*.csv", full.names = T) %>%
    #lapply(read.csv) %>%
    #bind_rows() %>%
    #dplyr::select(-X)
  
  #' fixing date format and making sure DateTime column is filled for all images
  #cow_data$Date <- dmy(cow_data$Date)
  #cow_data$DateTime <- as.POSIXct(paste(cow_data$Date, cow_data$Time), format="%Y-%m-%d %H:%M:%S")
  
  #' pulling only cow data
  #cow_data <- as.data.frame(cow_data[cow_data$Species == "Cattle",])
  
# -------------------------------------------------------------------------

  
  #' read in the functions that calculate the counts/duration for covariates
  source("./Scripts/Covariate Functions.R")

  
  #' finding unique detections
  cow_data <- uniq(cow_data, 5)
  
  
  # CATTLE ----------------------------------------------------------------
  #' vector of which cameras will be included in the subset
  #set.seed(7729) #seed for sampling Camera
  #cow_cam_subset <- sample(unique(cow_data$CameraLocation), round(length(unique(cow_data$CameraLocation))*.2))
  
  #' paring cow_data down to just the subsampled locations
  #cow_data_subset <- subset(cow_data, cow_data$CameraLocation %in% cow_cam_subset)
  
  #' Based on random sample of cameras where individual cows were counted
  #' compiling the data for each covariate
  week_cow_cov1 <-  do.call(rbind, cov1(cow_data, "weeks", "cow"))
  day_cow_cov1 <- do.call(rbind, cov1(cow_data, "days", "cow"))
  
  week_cow_cov2 <- do.call(rbind, cov2(cow_data, "weeks", "cow", m = 5))
  day_cow_cov2 <- do.call(rbind, cov2(cow_data, "days", "cow", m = 5))
  
  week_cow_cov3 <- do.call(rbind, cov3(cow_data, "weeks", "cow", m = 5))
  day_cow_cov3 <- do.call(rbind, cov3(cow_data, "days", "cow", m = 5))
  
  week_cow_cov4 <- do.call(rbind, cov4(cow_data, "weeks", "cow", m = 5))
  day_cow_cov4 <- do.call(rbind, cov4(cow_data, "days", "cow", m = 5))
  
  
  # creating dataframes with the Counts/Durations to be correlated
  week_cow_cor_df <- cbind(week_cow_cov1, week_cow_cov2$n_Detections, week_cow_cov3$Count, 
                           week_cow_cov4$Duration)
  colnames(week_cow_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                 "Cov1_n_Images", "Cov2_n_Detections", 
                                 "Cov3_max_Individuals", "Cov4_Duration") 
  day_cow_cor_df <- cbind(day_cow_cov1, day_cow_cov2$n_Detections, day_cow_cov3$Count, 
                          day_cow_cov4$Duration)
  colnames(day_cow_cor_df) <- c("CameraLocation", "Date", "Cov1_n_Images", 
                                "Cov2_n_Detections", "Cov3_max_Individuals", 
                                "Cov4_Duration") 
  
  # Weeks
  cor(week_cow_cor_df[,4:7]) 
  # m = 5 minutes: r = 0.95 - 0.98 (everything is highly correlated)
  # m = 10 minutes: r = 0.92 - 0.98
  # m = 30 minutes: r = 0.89 - 0.97
  # m = 60 minutes: r = 0.81 - 0.96
  # correlation declines slowly as m increases for weekly cattle grazing activity
  
  
  # Days
  cor(day_cow_cor_df[3:6])
  # m = 5 minutes: r = 0.89 - 0.95 (everything highly correlated)
  # m = 10 minutes: r = 0.84 - 0.95
  # m = 30 minutes: r = 0.74 - 0.92
  # m = 60 minutes: r = 0.68 - 0.91
  # correlation declines noticeably as m increases for daily cattle grazing activity
  
  
  
  
  ####  ALL CATTLE DETECTIONS ####
  #'  ----------------------------
  
  #'  Read in all detection data
  # det <- read.csv("./Data/Bassing_AllDetections18-21_2022-04-03.csv")
  megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdata18-21_2022-04-14.csv") %>% 
    dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation", 
                  "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle", 
                  "Species", "HumanActivity", "Count") %>%
    filter(!grepl("Moultrie", CameraLocation))
  moo <- megadata %>%
    filter(Species == "Cattle") 
  
  # write.csv(moo, "./Data/All_cattle_detections.csv")
  
  #' compiling the data for each covariate
  week_cow_npix <-  do.call(rbind, cov1(moo, "weeks", "cow"))
  day_cow_npix <- do.call(rbind, cov1(moo, "days", "cow"))
  
  week_cow_ndet <- do.call(rbind, cov2(moo, "weeks", "cow", m = 5))
  day_cow_ndet <- do.call(rbind, cov2(moo, "days", "cow", m = 5))
  
  week_cow_nmax <- do.call(rbind, cov3(moo, "weeks", "cow", m = 5))
  day_cow_nmax <- do.call(rbind, cov3(moo, "days", "cow", m = 5))
  
  week_cow_ntime <- do.call(rbind, cov4(moo, "weeks", "cow", m = 5))
  day_cow_ntime <- do.call(rbind, cov4(moo, "days", "cow", m = 5))
  
  
  # Create data frames with the sums/times
  week_cow_df <- week_cow_npix %>%
    full_join(week_cow_ndet, by = c("CameraLocation", "StartDate", "EndDate")) %>%
    full_join(week_cow_nmax, by = c("CameraLocation", "StartDate", "EndDate")) %>%
    full_join(week_cow_ntime, by = c("CameraLocation", "StartDate", "EndDate"))
  day_cow_df <- day_cow_npix %>%
    full_join(day_cow_ndet, by = c("CameraLocation", "Date")) %>%
    full_join(day_cow_nmax, by = c("CameraLocation", "Date")) %>%
    full_join(day_cow_ntime, by = c("CameraLocation", "Date"))
  
  #'  Save
  write.csv(week_cow_df, file = "./Outputs/CamTrap_Activity/weekly_cattle_activity.csv")
  write.csv(day_cow_df, file = "./Outputs/CamTrap_Activity/daily_cattle_activity.csv")
  
  
  #'  Check correlation among variables for all cameras
  #'  Remember max_individuals does not include real counts so this metric is not
  #'  very useful
  #'  Weeks
  cor(week_cow_cor_df[,4:7]) 
  #'  Days
  cor(day_cow_cor_df[3:6])
  
  
  
  
  # week_cow_df <- cbind(week_cow_npix, week_cow_ndet$Count, week_cow_nmax$Count, 
  #                      week_cow_ntime$Duration)
  # colnames(week_cow_df) <- c("CameraLocation", "StartDate", "EndDate", "n_images", 
  #                            "n_detections", "n_cow_max", "n_minutes") 
  # day_cow_df <- cbind(day_cow_npix, day_cow_ndet$Count, day_cow_nmax$Count, 
  #                     day_cow_ntime$Duration)
  # colnames(day_cow_df) <- c("CameraLocation", "Date", "n_images", "n_detections", 
  #                           "n_cow_max", "n_minutes") 
  
