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
  #'    Covariate 1 - Raw Counts of Images in each time step of interest (daily/weekly)
  #'    Covariate 2 - Number of Unique Detentions in each time step of interest
  #'    Covariate 3 - Sum of Max Count from each Unique Detection in each time step 
  #'    Covariate 4 - Number of Minutes of hunter activity in each time step of interest
  #'    Covariate 5 - Time Between Unique Detection of Cattle and Detection of Prey 
  #'    Species (INCOMPLETE)
  #'  ============================================================
  
  #'  packages
  library(data.table)
  library(tidyverse)
  
  #'  reading in cow data from multiple csv and combining them
  cow_data <- list.files(path = "./Data/Cow count csv files", pattern = "*.csv", full.names = T) %>%
    lapply(read.csv) %>%
    bind_rows() %>%
    dplyr::select(-X)
  
  #' fixing date format and making sure DateTime column is filled for all images
  cow_data$Date <- dmy(cow_data$Date)
  cow_data$DateTime <- as.POSIXct(paste(cow_data$Date, cow_data$Time), format="%Y-%m-%d %H:%M:%S")
  
  #' pulling only cow data
  cow_data <- as.data.frame(cow_data[cow_data$Species == "Cattle",])
  
  
  #' read in the functions that calculate the counts/duration for covariates
  source("./Scripts/Covariate Functions.R")
  
  
  # CATTLE ----------------------------------------------------------------
  #' vector of which cameras will be included in the subset
  #set.seed(7729) #seed for sampling Camera
  #cow_cam_subset <- sample(unique(cow_data$CameraLocation), round(length(unique(cow_data$CameraLocation))*.2))
  
  #' paring cow_data down to just the subsampled locations
  #cow_data_subset <- subset(cow_data, cow_data$CameraLocation %in% cow_cam_subset)
  
  
  #' compiling the data for each covariate
  week_cow_cov1 <-  do.call(rbind, cov1(cow_data, "weeks", "cow"))
  day_cow_cov1 <- do.call(rbind, cov1(cow_data, "days", "cow"))
  
  week_cow_cov2 <- do.call(rbind, cov2(cow_data, "weeks", "cow"))
  day_cow_cov2 <- do.call(rbind, cov2(cow_data, "days", "cow"))
  
  week_cow_cov3 <- do.call(rbind, cov3(cow_data, "weeks", "cow"))
  day_cow_cov3 <- do.call(rbind, cov3(cow_data, "days", "cow"))
  
  week_cow_cov4 <- do.call(rbind, cov4(cow_data, "weeks", "cow"))
  day_cow_cov4 <- do.call(rbind, cov4(cow_data, "days", "cow"))
  
  
  # creating dataframes with the Counts/Durations to be correlated
  week_cow_cor_df <- cbind(week_cow_cov1, week_cow_cov2$Count, week_cow_cov3$Count, 
                           week_cow_cov4$Duration)
  colnames(week_cow_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                 "Cov1", "Cov2", "Cov3", "Cov4") 
  day_cow_cor_df <- cbind(day_cow_cov1, day_cow_cov2$Count, day_cow_cov3$Count, 
                          day_cow_cov4$Duration)
  colnames(day_cow_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", 
                                "Cov4") 
  
  # Weeks
  cor(week_cow_cor_df[,4:7])
  # everything is highly correlated (r = 0.95 - 0.98)
  
  
  # Days
  cor(day_cow_cor_df[3:6])
  # everything highly correlated (r = 0.89 - 0.95)
  
  