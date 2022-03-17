#'  Cameron Ho
#'  Sarah Bassing Quantitative Ecology Lab
#'  Examining various possible covariate options for cattle grazing impact on ungulate abundance/behavior

setwd("G:/My Drive/1 Volunteers/Side projects for interns/Hunter-Cattle-Activity") #Sarah
setwd("C:/Users/billy/Desktop/Unprocessed Images/Hunter-Cattle Covariate Examination (R)") #Cameron


#'  packages
library(data.table) 


#'  reading in human data
cow_data <- read.csv("")


# ALL HUNT ----------------------------------------------------------------
#' vector of which cameras will be included in the subset
set.seed(7729) #seed for sampling Camera
hunt_cam_subset <- sample(unique(hunter_data$CameraLocation), round(length(unique(hunter_data$CameraLocation))*.2))

#' compiling the data for each covariate (ALL HUNT)
week_hunt_cov1 <-  cov1(hunter_data, "weeks", "human")
day_hunt_cov1 <- cov1(hunter_data, "days", "human")

week_hunt_cov2 <- cov2(hunter_data, "weeks", "human")
day_hunt_cov2 <- cov2(hunter_data, "days", "human")

week_hunt_cov3 <- cov32(hunter_data, "weeks", "human")
day_hunt_cov3 <- cov3(hunter_data, "days", "human")

week_hunt_cov4 <- cov4(hunter_data, "weeks", "human")
day_hunt_cov4 <- cov4(hunter_data, "days", "human")


#' subsetting
# Camera
week_hunt_cov1_sub <- week_hunt_cov1[names(week_hunt_cov1) %in% hunt_cam_subset]
week_hunt_cov2_sub <- week_hunt_cov2[names(week_hunt_cov2) %in% hunt_cam_subset]
week_hunt_cov3_sub <- week_hunt_cov3[names(week_hunt_cov3) %in% hunt_cam_subset]
week_hunt_cov4_sub <- week_hunt_cov4[names(week_hunt_cov4) %in% hunt_cam_subset]
day_hunt_cov1_sub <- day_hunt_cov1[names(day_hunt_cov1) %in% hunt_cam_subset]
day_hunt_cov2_sub <- day_hunt_cov2[names(day_hunt_cov2) %in% hunt_cam_subset]
day_hunt_cov3_sub <- day_hunt_cov3[names(day_hunt_cov3) %in% hunt_cam_subset]
day_hunt_cov4_sub <- day_hunt_cov4[names(day_hunt_cov4) %in% hunt_cam_subset]
# Steps
week_hunt_cov1_sub <- subsample_steps(week_hunt_cov1_sub, seed = 4971)
week_hunt_cov2_sub <- subsample_steps(week_hunt_cov2_sub, seed = 4971)
week_hunt_cov3_sub <- subsample_steps(week_hunt_cov3_sub, seed = 4971)
week_hunt_cov4_sub <- subsample_steps(week_hunt_cov4_sub, seed = 4971)
day_hunt_cov1_sub <- subsample_steps(day_hunt_cov1_sub, seed = 4971)
day_hunt_cov2_sub <- subsample_steps(day_hunt_cov2_sub, seed = 4971)
day_hunt_cov3_sub <- subsample_steps(day_hunt_cov3_sub, seed = 4971)
day_hunt_cov4_sub <- subsample_steps(day_hunt_cov4_sub, seed = 4971)