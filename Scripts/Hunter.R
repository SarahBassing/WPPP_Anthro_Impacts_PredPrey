#'  Cameron Ho
#'  Sarah Bassing Quantitative Ecology Lab
#'  Examining various possible covariate options for hunter activity (and potentially inclusion of other human activities) 
#'     impact on ungulate abundance/behavior
#'  Requires the functions from the 'Covariate Functions' Script

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
atv_list <- human_data[human_data$HumanActivity == "Vehicle ATV",]
hike_list <- human_data[human_data$HumanActivity == "Hiker",]


source("./Scripts/Covariate Functions.R")



# ALL HUNT ----------------------------------------------------------------

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
# creating dataframes with the Counts/Durations to be correlated
week_hunt_cor_df <- cbind(week_hunt_cov1, week_hunt_cov2$Count, week_hunt_cov3$Count, week_hunt_cov4$Duration)
colnames(week_hunt_cor_df) <- c("CameraLocation", "StartDate", "EndDate", "Cov1", "Cov2", "Cov3", "Cov4") 
day_hunt_cor_df <- cbind(day_hunt_cov1, day_hunt_cov2$Count, day_hunt_cov3$Count, day_hunt_cov4$Duration)
colnames(day_hunt_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", "Cov4") 

# Weeks
cor(week_hunt_cor_df[,4:7])
#cov1, cov2, cov3 are all correlated pretty highly (~.9)
#cov4 is most similar to cov1 (0.69)

# Days
cor(day_hunt_cor_df[3:6])
#cov1, cov2, cov3 are all correlated pretty highly (~.9)
#cov4 is most similar to cov1 (0.67)

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
week_rifle_cor_df <- cbind(week_rifle_cov1, week_rifle_cov2$Count, week_rifle_cov3$Count, week_rifle_cov4$Duration)
colnames(week_rifle_cor_df) <- c("CameraLocation", "StartDate", "EndDate", "Cov1", "Cov2", "Cov3", "Cov4") 
day_rifle_cor_df <- cbind(day_rifle_cov1, day_rifle_cov2$Count, day_rifle_cov3$Count, day_rifle_cov4$Duration)
colnames(day_rifle_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", "Cov4") 

# Weeks
cor(week_rifle_cor_df[,4:7])
#cov1, cov2, cov3 are all correlated pretty highly (~.9)
#cov4 is most similar to cov1 (0.69)

# Days
cor(day_rifle_cor_df[3:6])
#cov1, cov2, cov3 are all correlated pretty highly (~.9)
#cov4 is most similar to cov1 (0.67)


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
week_bow_cor_df <- cbind(week_bow_cov1, week_bow_cov2$Count, week_bow_cov3$Count, week_bow_cov4$Duration)
colnames(week_bow_cor_df) <- c("CameraLocation", "StartDate", "EndDate", "Cov1", "Cov2", "Cov3", "Cov4") 
day_bow_cor_df <- cbind(day_bow_cov1, day_bow_cov2$Count, day_bow_cov3$Count, day_bow_cov4$Duration)
colnames(day_bow_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", "Cov4") 

# Weeks
cor(week_bow_cor_df[,4:7]) 
#everything is super correlated (~.9)

# Days
cor(day_bow_cor_df[3:6])
#everything is super correlated (~.9)


# VEHICLES ----------------------------------------------------------------
#' this is done with Vehicle=T, so would also include dirt bikes and snowmobile

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
week_vehicle_cor_df <- cbind(week_vehicle_cov1, week_vehicle_cov2$Count, week_vehicle_cov3$Count, week_vehicle_cov4$Duration)
colnames(week_vehicle_cor_df) <- c("CameraLocation", "StartDate", "EndDate", "Cov1", "Cov2", "Cov3", "Cov4") 
day_vehicle_cor_df <- cbind(day_vehicle_cov1, day_vehicle_cov2$Count, day_vehicle_cov3$Count, day_vehicle_cov4$Duration)
colnames(day_vehicle_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", "Cov4") 

# Weeks
cor(week_vehicle_cor_df[,4:7])
#everything is super correlated (~.9)

# Days
cor(day_vehicle_cor_df[3:6])
#everything is pretty correlated (~.7)


# CARS ----------------------------------------------------------------

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
week_car_cor_df <- cbind(week_car_cov1, week_car_cov2$Count, week_car_cov3$Count, week_car_cov4$Duration)
colnames(week_car_cor_df) <- c("CameraLocation", "StartDate", "EndDate", "Cov1", "Cov2", "Cov3", "Cov4") 
day_car_cor_df <- cbind(day_car_cov1, day_car_cov2$Count, day_car_cov3$Count, day_car_cov4$Duration)
colnames(day_car_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", "Cov4") 

# Weeks
cor(week_car_cor_df[,4:7])
#everything is super correlated (~.85)

# Days
cor(day_car_cor_df[3:6])
#everything is pretty correlated (~.75)


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
week_atv_cor_df <- cbind(week_atv_cov1, week_atv_cov2$Count, week_atv_cov3$Count, week_atv_cov4$Duration)
colnames(week_atv_cor_df) <- c("CameraLocation", "StartDate", "EndDate", "Cov1", "Cov2", "Cov3", "Cov4") 
day_atv_cor_df <- cbind(day_atv_cov1, day_atv_cov2$Count, day_atv_cov3$Count, day_atv_cov4$Duration)
colnames(day_atv_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", "Cov4") 

# Weeks
cor(week_atv_cor_df[,4:7])
#cov2 and cov3 highly correlated (.99)


# Days
cor(day_atv_cor_df[3:6])
#cov2 and cov3 highly correlated (.99)


# HIKER ----------------------------------------------------------------

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
week_atv_cor_df <- cbind(week_atv_cov1, week_atv_cov2$Count, week_atv_cov3$Count, week_atv_cov4$Duration)
colnames(week_atv_cor_df) <- c("CameraLocation", "StartDate", "EndDate", "Cov1", "Cov2", "Cov3", "Cov4") 
day_atv_cor_df <- cbind(day_atv_cov1, day_atv_cov2$Count, day_atv_cov3$Count, day_atv_cov4$Duration)
colnames(day_atv_cor_df) <- c("CameraLocation", "Date", "Cov1", "Cov2", "Cov3", "Cov4") 

# Weeks
cor(week_atv_cor_df[,4:7])
#cov2 and cov3 highly correlated (.99)


# Days
cor(day_atv_cor_df[3:6])
#cov2 and cov3 highly correlated (.99)
