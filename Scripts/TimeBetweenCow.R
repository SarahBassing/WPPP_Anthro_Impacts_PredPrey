#'  Cameron Ho
#'  Sarah Bassing Quantitative Ecology Lab
#'  Creates data table that shows the time between a cow detection and a prey species detection
#'  Requires the functions from the 'Covariate Functions' Script

  
#'  packages
library(data.table)
library(tidyverse)


#'  reading in cow data from multiple csv and combining them
data <- list.files(path = "./Data/Cow count csv files", pattern = "*.csv", full.names = T) %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  dplyr::select(-X)


#'  fixing date format and making sure DateTime column is filled for all images
data$Date <- dmy(data$Date)
data$DateTime <- as.POSIXct(paste(data$Date, data$Time), format="%Y-%m-%d %H:%M:%S")


#'  read in the functions that calculate the counts/duration for covariates
source("./Scripts/Covariate Functions.R")


#'  pulling first image in unique detections from prey species
first_wtd <- first_uniq(data[data$Species == "White-tailed Deer",])
first_mule <- first_uniq(data[data$Species == "Mule Deer",])
first_elk <- first_uniq(data[data$Species == "Elk",])
first_moose <- first_uniq(data[data$Species == "Moose",])
first_hare <- first_uniq(data[data$Species == "Snowshoe Hare",])
first_rabbitspp <- first_uniq(data[data$Species == "Rabbit Spp",])
first_unkdeer <- first_uniq(data[data$Species == "Unknown Deer",])

#'  pulling last image in unique detections from cows
last_cow <- last_uniq(data[data$Species == "Cattle",])


#'  combine dataframes from cow and prey
detect_data <- as.data.frame(do.call(rbind, list(first_wtd, first_mule, first_elk, first_moose, 
                                   first_hare, first_rabbitspp, first_unkdeer, last_cow)))


#'  sorting by camera and date/time
detect_data <- arrange(detect_data, CameraLocation, DateTime)


#'  loop calculating the time between the most recent cow detection and prey detection
camera <- unique(detect_data$CameraLocation)
for(i in 1:length(camera)){
  # selecting which camera to process
  tmp_camera <- detect_data[detect_data$CameraLocation == camera[i],]
  
  MostRecentCow <- NA
  
  # loop that goes through each image for of a camera
  for(j in 1:nrow(tmp_camera)){
    #selecting which image to compare
    tmp_image <- detect_data[j,]
    
    # if the image is a cow, mark it as the most recent cow sighting
    if(tmp_image$Species == "Cattle"){
      MostRecentCow <- tmp_image
    }
    
    # if the image is a prey species
    # suppressing a warning created by is.na() when MostRecentCow is a row and not NA
    suppressWarnings(if(tmp_image$Species != "Cattle"){
      # marking TimeSinceCow as NA if there has not been a cow detected yet
      # **!!! change this if you want another indicator that a cow has not yet been detected!!!**
      if(is.na(MostRecentCow) == TRUE){detect_data$TimeSinceCow[j] <- NA}
      # calculating differences between prey detection and latest cow
      else(detect_data$TimeSinceCow[j] <- difftime(tmp_image$DateTime, MostRecentCow$DateTime, units="mins"))
    })
  }
}
