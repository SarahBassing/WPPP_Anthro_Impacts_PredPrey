  #'  Hunter activity metrics and correlation analysis
  #'  Cameron Ho
  #'  Sarah Bassing 
  #'  Quantitative Ecology Lab
  #'  March 2022
  #'  ============================================================
  #'  Script examining possible hunter activity metrics (and potentially inclusion 
  #'  of other human activities) from camera trap data to use as covariates 
  #'  representing impact of hunter activity on wildlife occurrence & activity patterns.
  #'  Notes:
  #'  Requires sourcing functions from the 'Covariate Functions.R' Script to calculate 
  #'  different covariate metrics:
  #'    Covariate 1 - Raw Counts of Images in each time step of interest (daily/weekly)
  #'    Covariate 2 - Number of Unique Detentions in each time step of interest
  #'    Covariate 3 - Sum of Max Count from each Unique Detection in each time step 
  #'    Covariate 4 - Number of Minutes of hunter activity in each time step of interest
  #'  ============================================================
  
  #'  Packages
  library(data.table)
  
  
  #'  reading in human data
  human_data <- read.csv("./Data/FallHuman_detections_2022-03-15.csv")
  
  
  #'  Splitting into other human activities
  hunter_data <- human_data[human_data$HumanActivity %like% "Hunter",]
  rifle_data <- human_data[human_data$HumanActivity == "Hunter Rifle",]
  bow_data <- human_data[human_data$HumanActivity == "Hunter Bow",]
  vehicle_data <- human_data[human_data$Vehicle == T,]
  car_data <- human_data[human_data$HumanActivity == "Vehicle Truck Car",]
  atv_data <- human_data[human_data$HumanActivity == "Vehicle ATV",]
  hike_data <- human_data[human_data$HumanActivity == "Hiker",]
  
  
  #' read in the functions that calculate the counts/duration for covariates
  source("./Scripts/Covariate Functions.R")
  
  
  
  #' ALL HUNT ----------------------------------------------------------------
  #' Rifle and bow hunters combined
  
  #' compiling the data for each covariate (ALL HUNT)
  week_hunt_cov1 <-  do.call(rbind, cov1(hunter_data, "weeks", "human"))
  day_hunt_cov1 <- do.call(rbind, cov1(hunter_data, "days", "human"))
  
  week_hunt_cov2 <- do.call(rbind, cov2(hunter_data, "weeks", "human"))
  day_hunt_cov2 <- do.call(rbind, cov2(hunter_data, "days", "human"))
  
  week_hunt_cov3 <- do.call(rbind, cov3(hunter_data, "weeks", "human"))
  day_hunt_cov3 <- do.call(rbind, cov3(hunter_data, "days", "human"))
  
  week_hunt_cov4 <- do.call(rbind, cov4(hunter_data, "weeks", "human"))
  day_hunt_cov4 <- do.call(rbind, cov4(hunter_data, "days", "human"))
  
  
  ## Pearson's Correlations --------------------------------------------------
  #' creating dataframes with the Counts/Durations to be correlated
  week_hunt_cor_df <- cbind(week_hunt_cov1, week_hunt_cov2$Count, 
                            week_hunt_cov3$Count, week_hunt_cov4$Duration)
  colnames(week_hunt_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                  "Cov1", "Cov2", "Cov3", "Cov4") 
  day_hunt_cor_df <- cbind(day_hunt_cov1, day_hunt_cov2$Count, day_hunt_cov3$Count, 
                           day_hunt_cov4$Duration)
  colnames(day_hunt_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", 
                                 "Cov4") 
  
  # Weeks
  cor(week_hunt_cor_df[,4:7])
  #cov1, cov2, cov3 are all highly correlated (r = 0.88 - 0.96)
  #cov4 is most correlated with cov1 (0.69), much lower with other covs (~0.52)
  
  # Days
  cor(day_hunt_cor_df[3:6])
  #cov1, cov2, cov3 are all highly correlated (r = 0.86 - 0.94)
  #cov4 is most correlated with cov1 (0.67), much lower with other covs (~0.47)
  
  # RIFLE HUNT --------------------------------------------------------------
  
  #' compiling the data for each covariate (RIFLE HUNT)
  week_rifle_cov1 <-  do.call(rbind, cov1(rifle_data, "weeks", "human"))
  day_rifle_cov1 <- do.call(rbind, cov1(rifle_data, "days", "human"))
  
  week_rifle_cov2 <- do.call(rbind, cov2(rifle_data, "weeks", "human"))
  day_rifle_cov2 <- do.call(rbind, cov2(rifle_data, "days", "human"))
  
  week_rifle_cov3 <- do.call(rbind, cov3(rifle_data, "weeks", "human"))
  day_rifle_cov3 <- do.call(rbind, cov3(rifle_data, "days", "human"))
  
  week_rifle_cov4 <- do.call(rbind, cov4(rifle_data, "weeks", "human"))
  day_rifle_cov4 <- do.call(rbind, cov4(rifle_data, "days", "human"))
  
  
  ## Pearson's Correlations --------------------------------------------------
  # creating dataframes with the Counts/Durations to be correlated
  week_rifle_cor_df <- cbind(week_rifle_cov1, week_rifle_cov2$Count, 
                             week_rifle_cov3$Count, week_rifle_cov4$Duration)
  colnames(week_rifle_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                   "Cov1", "Cov2", "Cov3", "Cov4") 
  day_rifle_cor_df <- cbind(day_rifle_cov1, day_rifle_cov2$Count, day_rifle_cov3$Count, 
                            day_rifle_cov4$Duration)
  colnames(day_rifle_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", 
                                  "Cov4") 
  
  # Weeks
  cor(week_rifle_cor_df[,4:7])
  #cov1, cov2, cov3 are all highly correlated (r = 0.88 - 0.96)
  #cov4 is most correlated with cov1 (0.69), much lower with other covs (~0.53)
  
  # Days
  cor(day_rifle_cor_df[3:6])
  #cov1, cov2, cov3 are all highly correlated (r = 0.85 - 0.94)
  #cov4 is most correlated with cov1 (0.67), much lower with other covs (~0.47)
  
  
  # BOW HUNT ----------------------------------------------------------------
  
  #' compiling the data for each covariate (BOW HUNT)
  week_bow_cov1 <-  do.call(rbind, cov1(bow_data, "weeks", "human"))
  day_bow_cov1 <- do.call(rbind, cov1(bow_data, "days", "human"))
  
  week_bow_cov2 <- do.call(rbind, cov2(bow_data, "weeks", "human"))
  day_bow_cov2 <- do.call(rbind, cov2(bow_data, "days", "human"))
  
  week_bow_cov3 <- do.call(rbind, cov3(bow_data, "weeks", "human"))
  day_bow_cov3 <- do.call(rbind, cov3(bow_data, "days", "human"))
  
  week_bow_cov4 <- do.call(rbind, cov4(bow_data, "weeks", "human"))
  day_bow_cov4 <- do.call(rbind, cov4(bow_data, "days", "human"))
  
  
  ## Pearson's Correlations --------------------------------------------------
  # creating dataframes with the Counts/Durations to be correlated
  week_bow_cor_df <- cbind(week_bow_cov1, week_bow_cov2$Count, week_bow_cov3$Count, 
                           week_bow_cov4$Duration)
  colnames(week_bow_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                 "Cov1", "Cov2", "Cov3", "Cov4") 
  day_bow_cor_df <- cbind(day_bow_cov1, day_bow_cov2$Count, day_bow_cov3$Count, 
                          day_bow_cov4$Duration)
  colnames(day_bow_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", 
                                "Cov4") 
  
  # Weeks
  cor(week_bow_cor_df[,4:7]) 
  # everything is highly correlated (r = 0.84 - 0.97)
  
  # Days
  cor(day_bow_cor_df[3:6])
  # everything is highly correlated (r = 0.80 - 0.95)
  
  
  # VEHICLES ----------------------------------------------------------------
  #' This is done with Vehicle = T, so includes trucks, ATVs, dirt bikes, & snowmobiles
  
  #' compiling the data for each covariate
  week_vehicle_cov1 <-  do.call(rbind, cov1(vehicle_data, "weeks", "human"))
  day_vehicle_cov1 <- do.call(rbind, cov1(vehicle_data, "days", "human"))
  
  week_vehicle_cov2 <- do.call(rbind, cov2(vehicle_data, "weeks", "human"))
  day_vehicle_cov2 <- do.call(rbind, cov2(vehicle_data, "days", "human"))
  
  week_vehicle_cov3 <- do.call(rbind, cov3(vehicle_data, "weeks", "human"))
  day_vehicle_cov3 <- do.call(rbind, cov3(vehicle_data, "days", "human"))
  
  week_vehicle_cov4 <- do.call(rbind, cov4(vehicle_data, "weeks", "human"))
  day_vehicle_cov4 <- do.call(rbind, cov4(vehicle_data, "days", "human"))
  
  
  ## Pearson's Correlations --------------------------------------------------
  # creating dataframes with the Counts/Durations to be correlated
  week_vehicle_cor_df <- cbind(week_vehicle_cov1, week_vehicle_cov2$Count, 
                               week_vehicle_cov3$Count, week_vehicle_cov4$Duration)
  colnames(week_vehicle_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                     "Cov1", "Cov2", "Cov3", "Cov4") 
  day_vehicle_cor_df <- cbind(day_vehicle_cov1, day_vehicle_cov2$Count, 
                              day_vehicle_cov3$Count, day_vehicle_cov4$Duration)
  colnames(day_vehicle_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", 
                                    "Cov3", "Cov4") 
  
  # Weeks
  cor(week_vehicle_cor_df[,4:7])
  # everything is highly correlated (r = 0.81 - 0.99)
  
  # Days
  cor(day_vehicle_cor_df[3:6])
  # everything is highly correlated (r = 0.69 - 0.99)
  
  
  # CARS ----------------------------------------------------------------
  #'  Truck/cars only!
  
  #' compiling the data for each covariate
  week_car_cov1 <-  do.call(rbind, cov1(car_data, "weeks", "human"))
  day_car_cov1 <- do.call(rbind, cov1(car_data, "days", "human"))
  
  week_car_cov2 <- do.call(rbind, cov2(car_data, "weeks", "human"))
  day_car_cov2 <- do.call(rbind, cov2(car_data, "days", "human"))
  
  week_car_cov3 <- do.call(rbind, cov3(car_data, "weeks", "human"))
  day_car_cov3 <- do.call(rbind, cov3(car_data, "days", "human"))
  
  week_car_cov4 <- do.call(rbind, cov4(car_data, "weeks", "human"))
  day_car_cov4 <- do.call(rbind, cov4(car_data, "days", "human"))
  
  
  ## Pearson's Correlations --------------------------------------------------
  # creating dataframes with the Counts/Durations to be correlated
  week_car_cor_df <- cbind(week_car_cov1, week_car_cov2$Count, week_car_cov3$Count, 
                           week_car_cov4$Duration)
  colnames(week_car_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                 "Cov1", "Cov2", "Cov3", "Cov4") 
  day_car_cor_df <- cbind(day_car_cov1, day_car_cov2$Count, day_car_cov3$Count, 
                          day_car_cov4$Duration)
  colnames(day_car_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", 
                                "Cov4") 
  
  # Weeks
  cor(week_car_cor_df[,4:7])
  # everything is highly correlated (r = 0.82 - 0.99)
  
  # Days
  cor(day_car_cor_df[3:6])
  # everything is highly correlated (r = 0.71 - 0.99)
  
  
  # ATV ----------------------------------------------------------------
  
  #' compiling the data for each covariate
  week_atv_cov1 <-  do.call(rbind, cov1(atv_data, "weeks", "human"))
  day_atv_cov1 <- do.call(rbind, cov1(atv_data, "days", "human"))
  
  week_atv_cov2 <- do.call(rbind, cov2(atv_data, "weeks", "human"))
  day_atv_cov2 <- do.call(rbind, cov2(atv_data, "days", "human"))
  
  week_atv_cov3 <- do.call(rbind, cov3(atv_data, "weeks", "human"))
  day_atv_cov3 <- do.call(rbind, cov3(atv_data, "days", "human"))
  
  week_atv_cov4 <- do.call(rbind, cov4(atv_data, "weeks", "human"))
  day_atv_cov4 <- do.call(rbind, cov4(atv_data, "days", "human"))
  
  
  ## Pearson's Correlations --------------------------------------------------
  # creating dataframes with the Counts/Durations to be correlated
  week_atv_cor_df <- cbind(week_atv_cov1, week_atv_cov2$Count, week_atv_cov3$Count, 
                           week_atv_cov4$Duration)
  colnames(week_atv_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                 "Cov1", "Cov2", "Cov3", "Cov4") 
  day_atv_cor_df <- cbind(day_atv_cov1, day_atv_cov2$Count, day_atv_cov3$Count, 
                          day_atv_cov4$Duration)
  colnames(day_atv_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", 
                                "Cov4") 
  
  # Weeks
  cor(week_atv_cor_df[,4:7])
  # cov2 and cov3 highly correlated (.99)
  # cov1 not correlated with cov2 or cov3 (r = 0.46- 0.49)
  # cov1 moderately correlated with cov4 (0.61)
  
  
  # Days
  cor(day_atv_cor_df[3:6])
  # cov2 and cov3 highly correlated (.99)
  # cov1 not correlated with cov2, cov3, or cov4 (r = 0.39 - 0.58)
  
  
  # HIKER ----------------------------------------------------------------
  
  #' compiling the data for each covariate
  week_hike_cov1 <-  do.call(rbind, cov1(hike_data, "weeks", "human"))
  day_hike_cov1 <- do.call(rbind, cov1(hike_data, "days", "human"))
  
  week_hike_cov2 <- do.call(rbind, cov2(hike_data, "weeks", "human"))
  day_hike_cov2 <- do.call(rbind, cov2(hike_data, "days", "human"))
  
  week_hike_cov3 <- do.call(rbind, cov3(hike_data, "weeks", "human"))
  day_hike_cov3 <- do.call(rbind, cov3(hike_data, "days", "human"))
  
  week_hike_cov4 <- do.call(rbind, cov4(hike_data, "weeks", "human"))
  day_hike_cov4 <- do.call(rbind, cov4(hike_data, "days", "human"))
  
  
  ## Pearson's Correlations --------------------------------------------------
  # creating dataframes with the Counts/Durations to be correlated
  week_hike_cor_df <- cbind(week_hike_cov1, week_hike_cov2$Count, week_hike_cov3$Count, 
                            week_hike_cov4$Duration)
  colnames(week_hike_cor_df) <- c("CameraLocation", "StartDate", "EndDate", 
                                  "Cov1", "Cov2", "Cov3", "Cov4") 
  day_hike_cor_df <- cbind(day_hike_cov1, day_hike_cov2$Count, day_hike_cov3$Count, 
                           day_hike_cov4$Duration)
  colnames(day_hike_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", 
                                 "Cov4") 
  
  # Weeks
  cor(week_hike_cor_df[,4:7])
  # everything is highly correlated (r = 0.74 - 0.99)
  
  
  # Days
  cor(day_hike_cor_df[3:6])
  # cov1 and cov3 are correlated (0.65), cov1 and cov2 not correlated (0.52)
  # cov2 and cov3 are highly correlated (0.96)
  # cov4 is correlated with cov1 - cov3 (r = 0.61 - 0.93)
