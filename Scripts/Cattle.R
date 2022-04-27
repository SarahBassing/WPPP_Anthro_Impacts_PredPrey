  #'  Cattle grazing activity metrics and correlation analysis
  #'  Cameron Ho
  #'  Sarah Bassing 
  #'  Quantitative Ecology Lab
  #'  March 2022
  #'  ============================================================
  #'  Script examines & generates cattle grazing activity metrics from camera  
  #'  data on public lands to use as covariates representing impact of grazing 
  #'  activity on wildlife occurrence & activity patterns. Goal is to determine 
  #'  which activity metrics are most correlated with the actual count of cattle 
  #'  detected in each image within a given time period. Therefore, comparing 
  #'  the sum of the maximum number of uniquely identifiable individuals per 
  #'  detection event, i.e., covariate #3 - the maximum number of cows in a single 
  #'  image, to other detection metrics. Then generates these metrics for entire
  #'  data set.
  #'  
  #'  Requires sourcing functions from the 'Covariate Functions.R' Script to calculate 
  #'  different covariate metrics:
  #'    Covariate 1 - Sum of individual images in each time step of interest (daily/weekly)
  #'    Covariate 2 - Sum of number of unique detention events in each time step
  #'    Covariate 3 - Sum of maximum number of individuals detected in a single 
  #'    image per detection event in each time step
  #'    Covariate 4 - Sum of minutes across unique detection events in each time step
  #'  ============================================================
  
  #'  Load packages
  library(data.table)
  library(lubridate)
  library(tidyverse)
  
  #'  Source functions that calculate the counts/duration for covariates
  source("./Scripts/Covariate Functions.R")
  
  #'  Read in complete cow dataset
  #' megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdata18-21_2022-04-27.csv") %>% #2022-04-14
  #'   dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation",
  #'                 "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle",
  #'                 "Species", "HumanActivity", "Count") %>%
  #'   filter(!grepl("Moultrie", CameraLocation))
  #' moo <- megadata %>%
  #'   filter(Species == "Cattle")
  #' write.csv(moo, "./Data/All_cattle_detections.csv")
  #' cow_data <- read.csv("./Data/All_cattle_detections.csv")
  #' 
  #' #'  Identify unique detection events
  #' cow_data <- uniq(cow_data, 5)
  #' uniq_df <- write.csv(cow_data, "./Data/uniq_detections.csv")
  cow_data <- read.csv("./Data/uniq_detections.csv")
  moo <- cow_data
  
  
  ####  Correlation test among cow count data & other activity metrics  ####  
  #'  ------------------------------------------------------------------
  
  #'  Read in cow data from multiple csv and combining them
  #'  These data include actual counts of each cow detected in an image
  cow_data <- list.files(path = "./Data/Cow count csv files", pattern = "*.csv", full.names = T) %>%
    lapply(read.csv) %>%
    bind_rows() %>%
    dplyr::select(-X)

  #'  Fix date format and make sure DateTime column is filled for all images
  cow_data$Date <- dmy(cow_data$Date)
  cow_data$DateTime <- as.POSIXct(paste(cow_data$Date, cow_data$Time), format="%Y-%m-%d %H:%M:%S")

  #'  Pull only cow data
  cow_data <- as.data.frame(cow_data[cow_data$Species == "Cattle",])
  
  #'  Generate cattle activity metrics using 4 methods
  #'  1. Number of cattle images per week (or day)
  week_cow_cov1 <-  do.call(rbind, cov1(cow_data, "weeks", "cow"))
  day_cow_cov1 <- do.call(rbind, cov1(cow_data, "days", "cow"))
  
  #'  2. Number of unique detections per week (or day)
  week_cow_cov2 <- do.call(rbind, cov2(cow_data, "weeks", "cow"))
  day_cow_cov2 <- do.call(rbind, cov2(cow_data, "days", "cow"))
  
  #'  3. Maximum number of cows detected in a single image per detection event
  #'  per week (or day)
  week_cow_cov3 <- do.call(rbind, cov3(cow_data, "weeks", "cow"))
  day_cow_cov3 <- do.call(rbind, cov3(cow_data, "days", "cow"))
  
  #'  4. Amount of time cattle spent in front of camera per week (or day)
  week_cow_cov4 <- do.call(rbind, cov4(cow_data, "weeks", "cow"))
  day_cow_cov4 <- do.call(rbind, cov4(cow_data, "days", "cow"))
  
  
  #'  Create dfs with the Counts/Durations
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
  
  #'  Run correlation test among four metrics
  #'  Weekly cattle activity
  cor(week_cow_cor_df[,4:7]) 
  #'  m = 5 minutes: r = 0.95 - 0.98 (everything is highly correlated)
  #'  m = 10 minutes: r = 0.92 - 0.98
  #'  m = 30 minutes: r = 0.89 - 0.97
  #'  m = 60 minutes: r = 0.81 - 0.96
  #'  Correlation declines slowly as m increases for weekly cattle grazing activity
  
  
  #'  Daily cattle activity
  cor(day_cow_cor_df[3:6])
  #'  m = 5 minutes: r = 0.89 - 0.95 (everything highly correlated)
  #'  m = 10 minutes: r = 0.84 - 0.95
  #'  m = 30 minutes: r = 0.74 - 0.92
  #'  m = 60 minutes: r = 0.68 - 0.91
  #'  Correlation declines noticeably as m increases for daily cattle grazing activity
  
  
  
  
  ####  Generate cattle metrics for ALL cattle detections ####
  #'  ----------------------------------------------------
  
  #' #'  Read in all detection data
  #' megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdata18-21_2022-04-14.csv") %>%
  #'   dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation",
  #'                 "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle",
  #'                 "Species", "HumanActivity", "Count") %>%
  #'   filter(!grepl("Moultrie", CameraLocation))
  #' moo <- megadata %>%
  #'   filter(Species == "Cattle")
  #' write.csv(moo, "./Data/All_cattle_detections.csv")
  
  #'  Data already run through first_uniq function to identify detection events
  #'  This speeds up the code below
  cow_data <- read.csv("./Data/uniq_detections.csv")
  moo <- cow_data
  
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
  #'  Remember max_individuals does not include real counts so correlations with
  #'  this variable are meaningless
  #'  Weeks
  cor(week_cow_df[,4:7]) 
  #'  Days
  cor(day_cow_df[3:6])
  
  
