#'  Cameron Ho
#'  Sarah Bassing Quantitative Ecology Lab
#'  Functions to manipulate data sets used to assess correlation of covariates

#'  ==========================
#'  NOTES 
#'  Currently using 5 minute interval to identify unique detection events
#'  ==========================

#' Required Packages -------------------------------------------------------
library(lubridate)
library(data.table)
library(tidyverse)


#' First Image from Unique Detections
first_uniq <- function(data){  
  
  #'  Create a column identifying whether each image is an "independent" event
  use_data <- data %>%
    arrange(CameraLocation, DateTime) 

  #'  Empty vector to be filled
  caps <- c()

  #'  Fill first element of the vector to get it started (1st detection event)
  caps[1] <- 1

  #'  Giant for loop to run through the rest of the dataset
  #'  1. If camera site is diff from previous row then give unique value. If not then...
  #'  2. If species detected is diff from previous row at same site then give unique value. If not then...
  #'  3. If DateTime is >30 min from previous DateTime at same site for same species then give unique value. If not then...
  #'  4. Capture value is the same as that in the previous row
  for (i in 2:nrow(data)){
    if (use_data$CameraLocation[i-1] != use_data$CameraLocation[i]) caps[i] = i
    else (if (use_data$Species[i-1] != use_data$Species[i]) caps[i] = i
          else (if (difftime(use_data$DateTime[i], use_data$DateTime[i-1], units =
                             c("mins")) > 5) caps[i] = i        #change the number based on the time between detections that is unique
                else caps[i] = caps[i-1]))
  } # close loop
  #'  Format the caps vector as a factor
  caps <- as.factor(caps)

  #'  Add new column to larger data set- this can then be filtered different ways to pull out specific observations for unique detection events
  capdata <- cbind(as.data.frame(use_data), caps)
  
  #'  Filter unique detections to pull out the first image of every detection event
  firstdetection <- capdata %>% 
    group_by(caps) %>% 
    slice(1L) %>%
    ungroup()
  
  return(firstdetection)
}




#' Max Count from Unique Detections
max_uniq <- function(data){  
  
  #'  Create a column identifying whether each image is an "independent" event
  use_data <- data %>%
    arrange(CameraLocation, DateTime) 

  #'  Empty vector to be filled
  caps <- c()

  #'  Fill first element of the vector to get it started (1st detection event)
  caps[1] <- 1

  #'  Giant for loop to run through the rest of the dataset
  #'  1. If camera site is diff from previous row then give unique value. If not then...
  #'  2. If species detected is diff from previous row at same site then give unique value. If not then...
  #'  3. If DateTime is >30 min from previous DateTime at same site for same species then give unique value. If not then...
  #'  4. Capture value is the same as that in the previous row
  for (i in 2:nrow(use_data)){
    if (use_data$CameraLocation[i-1] != use_data$CameraLocation[i]) caps[i] = i
    else (if (use_data$Species[i-1] != use_data$Species[i]) caps[i] = i
          else (if (difftime(use_data$DateTime[i], use_data$DateTime[i-1], units =
                             c("mins")) > 5) caps[i] = i        #change the number based on the time between detections that is unique
                else caps[i] = caps[i-1]))
  } # close loop
  #'  Format the caps vector as a factor
  caps <- as.factor(caps)

  #'  Add new column to larger data set- this can then be filtered different ways to pull out specific observations for unique detection events
  capdata <- cbind(as.data.frame(use_data), caps)
  
  # Filter unique detections to pull out the max count from every detection event
  maxdetection <- capdata %>% 
    group_by(caps) %>% 
    slice_max(Count) %>%
    slice_head(n=1) %>%
    ungroup()
  
  return(maxdetection)
}





#' Unique Detections
uniq <- function(data){

  #'  Create a column identifying whether each image is an "independent" event
  use_data <- data %>%
    arrange(CameraLocation, DateTime) 

  #'  Empty vector to be filled
  caps <- c()

  #'  Fill first element of the vector to get it started (1st detection event)
  caps[1] <- 1

  #'  Giant for loop to run through the rest of the dataset
  #'  1. If camera site is diff from previous row then give unique value. If not   #’  then...
  #'  2. If species detected is diff from previous row at same site then give unique #’  value. If not then...
  #'  3. If DateTime is >30 min from previous DateTime at same site for same species #’  then give unique value. If not then...
  #'  4. Capture value is the same as that in the previous row
  for (i in 2:nrow(use_data)){
    if (use_data$CameraLocation[i-1] != use_data$CameraLocation[i]) caps[i] = i
    else (if (use_data$Species[i-1] != use_data$Species[i]) caps[i] = i
          else (if (difftime(use_data$DateTime[i], use_data$DateTime[i-1], units =
                             c("mins")) > 5) caps[i] = i        #change the number based on the time between detections that is unique
                else caps[i] = caps[i-1]))
  } # close loop
  #'  Format the caps vector as a factor
  caps <- as.factor(caps)

  #'  Add new column to larger data set- this can then be filtered different ways #’  to pull out specific observations for unique detection events
  capdata <- cbind(as.data.frame(use_data), caps)

  return(capdata)
}





#' Covariate 1 - Raw Counts of Images in each Timestep
# Careful with the input data set since the function will not discriminate by Species/HumanAct
# time = "weeks", "days"
# interest = "human", "cow"
cov1 <- function(data, time, interest){
  
  locations <- unique(data$CameraLocation)
  
  #' creating list that will be filled for each camera
  cov1_list <- vector("list", length = length(locations))
  names(cov1_list) <- locations
  
  
  #'  loop that outputs a list of data frames with all camera locations
  #'  each data frames has the raw count in each timestep
  #'  1(i). goes through each camera location 
  #'  2(j). goes trough each timestep
  for(i in 1:length(locations)){
    
    #temporary data set for camera location i
    tmp_data <- data[data$CameraLocation == locations[i],]
    
    #creating a fixed time sequence depending on 'interest'
    year <- year(tmp_data[1,]$Date)
    if(interest == "human"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
    
    if(interest == "cow"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
    
    
    #temporary data frame to put
    #if statements for difference in data frame setup between days and weeks
    if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                             StartDate = time_seq[1:(length(time_seq)-1)], 
                                             EndDate = time_seq[2:length(time_seq)]-1, 
                                             Count = 0)}
    if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i], Date = time_seq , Count = 0)}
    
    #loop to calculate raw count per timestep
    for(j in 1:(length(time_seq)-1)){
      
      tmp_df[j, "Count"] <- nrow(tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                                          & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),])
    }
    cov1_list[[i]] <- tmp_df
  }
  return(cov1_list)
}





#' Covariate 2 -  Number of Unique Detentions in each Timestep
# Careful with the input data set since the function will not discriminate by Species/HumanAct
# time = "weeks", "days"
# interest = "human", "cow"
cov2 <- function(data, time, interest){
  
  first_uniq_data <- first_uniq(data)
  
  locations <- unique(first_uniq_data$CameraLocation)
  
  #' creating list that will be filled for each camera
  cov2_list <- vector("list", length = length(locations))
  names(cov2_list) <- locations
  
  
  #'  loops that outputs a list of data frames with all camera locations
  #'  1(i). goes through each camera location 
  #'  2(j). goes trough each timestep
  #'  each data frames has the raw count in each timestep
  for(i in 1:length(locations)){
    
    #temporary data set for camera location i
    tmp_data <- first_uniq_data[first_uniq_data$CameraLocation == locations[i],]
    
    #creating a fixed time sequence depending on 'interest'
    year <- year(tmp_data[1,]$Date)
    if(interest == "human"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
    
    if(interest == "cow"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}

    
    
    #temporary data frame to put
    #if statements for difference in data frame setup between days and weeks
    if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                             StartDate = time_seq[1:(length(time_seq)-1)], 
                                             EndDate = time_seq[2:length(time_seq)]-1, 
                                             Count = 0)}
    if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i],Date = time_seq , Count = 0)}
    
    #loop to calculate raw count per timestep
    for(j in 1:(length(time_seq)-1)){
      
      tmp_df[j, "Count"] <- nrow(tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                                          & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),])
    }
    cov2_list[[i]] <- tmp_df
  }
  return(cov2_list)
}





#' Covariate 3 -  Sum of Max Count from each Unique Detection in each Timestep
# Careful with the input data set since the function will not discriminate by Species/HumanAct
# time = "weeks", "days"
# interest = "human", "cow"
cov3 <- function(data, time, interest){
  
  max_uniq_data <- max_uniq(data)
  
  locations <- unique(max_uniq_data$CameraLocation)
  
  #' creating list that will be filled for each camera
  cov3_list <- vector("list", length = length(locations))
  names(cov3_list) <- locations
  
  
  #'  loops that outputs a list of data frames with all camera locations
  #'  1(i). goes through each camera location 
  #'  2(j). goes trough each timestep
  #'  each data frames has the raw count in each timestep
  for(i in 1:length(locations)){
    
    #temporary data set for camera location i
    tmp_data <- max_uniq_data[max_uniq_data$CameraLocation == locations[i],]
    
    #creating a fixed time sequence depending on 'interest'
    year <- year(tmp_data[1,]$Date)
    if(interest == "human"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
    
    if(interest == "cow"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
    
    
    
    #temporary data frame to put
    #if statements for difference in data frame setup between days and weeks
    if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                             StartDate = time_seq[1:(length(time_seq)-1)], 
                                             EndDate = time_seq[2:length(time_seq)]-1, 
                                             Count = 0)}
    if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i],Date = time_seq , Count = 0)}
    
    #loop to calculate raw count per timestep
    for(j in 1:(length(time_seq)-1)){

      tmp_df[j, "Count"]  <- sum(tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                                          & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),]$Count)
    }
    cov3_list[[i]] <- tmp_df
  }
  return(cov3_list)
}





#' Covariate 4 -  Number of Minutes in each Timestep
# Careful with the input data set since the function will not discriminate by Species/HumanAct
# time = "weeks", "days"
# interest = "human", "cow"
cov4 <- function(data, time, interest){
  
  uniq_data <- uniq(data)
  
  locations <- unique(uniq_data$CameraLocation)
  
  #' creating list that will be filled for each camera
  cov4_list <- vector("list", length = length(locations))
  names(cov4_list) <- locations
  
  
  #'  loops that outputs a list of data frames with all camera locations
  #'  1(i). goes through each camera location 
  #'  2(j). goes trough each timestep
  #'  3(k). goes through each unique detection calculating the durations
  #'  each data frames has the raw count in each timestep
  for(i in 1:length(locations)){
    
    #temporary data set for camera location i
    tmp_data <- uniq_data[uniq_data$CameraLocation == locations[i],]
    
    #creating a fixed time sequence depending on 'interest'
    year <- year(tmp_data[1,]$Date)
    if(interest == "human"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
    
    if(interest == "cow"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
    
    
    
    #temporary data frame to put
    #if statements for difference in data frame setup between days and weeks
    if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                             StartDate = time_seq[1:(length(time_seq)-1)], 
                                             EndDate = time_seq[2:length(time_seq)]-1, 
                                             Duration = 0)}
    if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i],Date = time_seq , Duration = 0)}
    
    for(j in 1:(length(time_seq)-1)){
      
      tmp_timestep <- tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                                      & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),]
      
      detections <- as.vector(unique(tmp_timestep$caps))
      tmp_duration <- vector()
      
      #if there are no detection in the timestep set duration to 0
      if(nrow(tmp_timestep) != 0){
        for(k in 1:length(detections)){
          
          tmp_detection <- tmp_timestep[tmp_timestep$caps == detections[k],]
          
          #can change time metric
          tmp_duration[k] <- as.numeric(difftime(tmp_detection[nrow(tmp_detection),]$DateTime, tmp_detection[1,]$DateTime, units="secs"))
        }
        tmp_df[j, "Duration"] <- sum(tmp_duration)
      } else(tmp_df[j, "Duration"] <- 0)

    }
    cov4_list[[i]] <- tmp_df
  }
  return(cov4_list)
}





#' Covariate 5 -  Average Duration between Unique Detections in each Timestep (INCOMPLETE)
# Careful with the input data set since the function will not discriminate by Species/HumanAct
# time = "weeks", "days"
# interest = "human", "cow"
cov5 <- function(data, time, interest){
  
  uniq_data <- uniq(data)
  
  locations <- unique(uniq_data$CameraLocation)
  
  #' creating list that will be filled for each camera
  cov5_list <- vector("list", length = length(locations))
  names(cov5_list) <- locations
  
  
  #'  loops that outputs a list of data frames with all camera locations
  #'  1(i). goes through each camera location 
  #'  2(j). goes trough each timestep
  #'  3(k). goes through each unique detection calculating the durations
  #'  each data frames has the raw count in each timestep
  for(i in 1:length(locations)){
    
    #temporary data set for camera location i
    tmp_data <- uniq_data[uniq_data$CameraLocation == locations[i],]
    
    #creating a fixed time sequence depending on 'interest'
    year <- year(tmp_data[1,]$Date)
    if(interest == "human"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
    
    if(interest == "cow"){
      if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
      else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
    
    
    
    #temporary data frame to put
    #if statements for difference in data frame setup between days and weeks
    if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                             StartDate = time_seq[1:(length(time_seq)-1)], 
                                             EndDate = time_seq[2:length(time_seq)]-1, 
                                             Duration = 0)}
    if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i],Date = time_seq , Duration = 0)}
    
    for(j in 1:(length(time_seq)-1)){
      
      tmp_timestep <- tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                               & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),]
      
      detections <- as.vector(unique(tmp_timestep$caps))
      tmp_duration <- vector()
      
      
      for(k in 1:(length(detections)-1)){
        
        unique1 <- tmp_timestep[tmp_timestep$caps == detections[k],][sum(tmp_timestep$caps == detections[k]),]
        unique2 <- tmp_timestep[tmp_timestep$caps == detections[k+1],][1,]
        
        #can change time metric
        tmp_duration[k] <- as.numeric(difftime(tmp_detection[nrow(tmp_detection),]$DateTime, tmp_detection[1,]$DateTime, units="secs"))
      }
      
      tmp_df[j, "Duration"] <- mean(tmp_duration)
      
    }
    cov5_list[[i]] <- tmp_df
  }
  return(cov5_list)
}





#' Subsample Days/Weeks
subsample_steps <- function(data, percent = .2, seed){
  if(is.na(seed) == T){set.seed(seed)}

  step_subset <- lapply(data, function(x){sample(x$StartDate, size = round(nrow(x)*percent))})
  
  #applying the subset steps to the camera data frames
  subset_list <- list()
  for(i in 1:length(data)){
    tmp_df <- data[[i]]
    subset_list[[i]] <- tmp_df[tmp_df$StartDate %in% step_subset[[i]],]
  }
  
  names(subset_list) <- names(data)
  
  return(do.call(rbind, subset_list))
}

