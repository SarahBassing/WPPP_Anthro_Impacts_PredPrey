  #'  Functions to generate different metrics of activity at camera sites
  #'  Cameron Ho
  #'  Sarah Bassing
  #'  Quantitative Ecology Lab
  #'  March 2022
  #'  ==========================
  #'  NOTES:
  #'  Adjust temporal interval to identify unique detection events with m (in 
  #'  minutes) in if-else statements in the first_uniq, last_uniq, max_uniq, 
  #'  and uniq functions.
  #'  ==========================
  
  #' Required Packages 
  library(lubridate)
  library(data.table)
  library(tidyverse)
  library(progress)
  
  
  ####  Functions to identify different types of detections  ####
  #'  -------------------------------------------------------
  
  #' First Image from each Unique Detections
  first_uniq <- function(data, m){  
    
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
    #'  3. If DateTime is >m min from previous DateTime at same site for same species then give unique value. If not then...
    #'  4. Capture value is the same as that in the previous row
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Setting unique detections" , total = nrow(data)-1)
    for (i in 2:nrow(data)){
      if (use_data$CameraLocation[i-1] != use_data$CameraLocation[i]) caps[i] = i
      else (if (use_data$Species[i-1] != use_data$Species[i]) caps[i] = i
            else (if (difftime(use_data$DateTime[i], use_data$DateTime[i-1], units =
                               c("mins")) > m) caps[i] = i        
                  else caps[i] = caps[i-1]))
      
      #progress bar tick
      progbar$tick()
    } # close loop
    progbar$close
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
  

  #' Last Image from Unique Detections
  last_uniq <- function(data, m){  
    
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
    #'  3. If DateTime is >m min from previous DateTime at same site for same species then give unique value. If not then...
    #'  4. Capture value is the same as that in the previous row
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Setting unique detections" , total = nrow(data)-1)
    for (i in 2:nrow(data)){
      if (use_data$CameraLocation[i-1] != use_data$CameraLocation[i]) caps[i] = i
      else (if (use_data$Species[i-1] != use_data$Species[i]) caps[i] = i
            else (if (difftime(use_data$DateTime[i], use_data$DateTime[i-1], units =
                               c("mins")) > m) caps[i] = i        
                  else caps[i] = caps[i-1]))
      
      #progress bar tick
      progbar$tick()
    } # close loop
    progbar$close
    #'  Format the caps vector as a factor
    caps <- as.factor(caps)
    
    #'  Add new column to larger data set- this can then be filtered different ways to pull out specific observations for unique detection events
    capdata <- cbind(as.data.frame(use_data), caps)
    
    #'  Filter unique detections to pull out the last image of every detection event
    lastdetection <- capdata %>% 
      group_by(caps) %>% 
      slice_tail(n = 1) %>%
      ungroup()
    
    return(lastdetection)
  }
  

  #' Maximum Count of individuals in an image from Unique Detections
  max_uniq <- function(data, m){  
    
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
    #'  3. If DateTime is >m min from previous DateTime at same site for same species then give unique value. If not then...
    #'  4. Capture value is the same as that in the previous row
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Setting unique detections" , total = nrow(data)-1)
    for (i in 2:nrow(use_data)){
      if (use_data$CameraLocation[i-1] != use_data$CameraLocation[i]) caps[i] = i
      else (if (use_data$Species[i-1] != use_data$Species[i]) caps[i] = i
            else (if (difftime(use_data$DateTime[i], use_data$DateTime[i-1], units =
                               c("mins")) > m) caps[i] = i        
                  else caps[i] = caps[i-1]))
      
      #progress bar tick
      progbar$tick()
    } # close loop
    progbar$close
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
  uniq <- function(data, m){
  
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
    #'  3. If DateTime is >m min from previous DateTime at same site for same species then give unique value. If not then...
    #'  4. Capture value is the same as that in the previous row
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Setting unique detections" , total = nrow(data)-1)
    for (i in 2:nrow(use_data)){
      if (use_data$CameraLocation[i-1] != use_data$CameraLocation[i]) caps[i] = i
      else (if (use_data$Species[i-1] != use_data$Species[i]) caps[i] = i
            else (if (difftime(use_data$DateTime[i], use_data$DateTime[i-1], units =
                               c("mins")) > m) caps[i] = i        
                  else caps[i] = caps[i-1]))
      
      #progress bar tick
      progbar$tick()
    } # close loop
    progbar$close
    #'  Format the caps vector as a factor
    caps <- as.factor(caps)
  
    #'  Add new column to larger data set- this can then be filtered different ways 
    #'  to pull out specific observations for unique detection events
    capdata <- cbind(as.data.frame(use_data), caps)
  
    return(capdata)
  }
  
  
  
  
  ####  Functions to generate activity covariates  ####
  #'  ---------------------------------------------
  
  #'  Covariate 1 - Sum of individual images in each time step of interest (daily/weekly)
  #'  Careful with the input data set since the function will not discriminate by Species/HumanAct
  #'  time = "weeks", "days"
  #'  interest = "human", "cow"
  cov1 <- function(data, time, interest){
    
    #'  Identify unique camera location names
    locations <- unique(data$CameraLocation)
    
    #' creating list that will be filled for each camera
    cov1_list <- vector("list", length = length(locations))
    names(cov1_list) <- locations
    
    
    #'  loop that outputs a list of data frames with all camera locations
    #'  each data frames has the raw count in each timestep
    #'  1(i). goes through each camera location 
    #'  2(j). goes trough each timestep
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Calculating covariate 1 values" , total = length(locations))
    for(i in 1:length(locations)){
      
      #temporary data set for camera location i
      tmp_data <- data[data$CameraLocation == locations[i],] %>%
        arrange(CameraLocation, DateTime)
        
      
      #creating a fixed time sequence depending on 'interest'
      #requiring a check to see if last row in tmp_data has the same year as first row
      #if the data runs past a year, will add extra rows to account for another year
      year <- year(tmp_data[1,]$Date)
      year_check <- tmp_data[nrow(tmp_data),]$Date
      if(interest == "human"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
      if(year_check > paste0(year+1,"-7-31") & interest == "human"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-30")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-31")), by=time)))}
      
      
      if(interest == "cow"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
      if(year != year(year_check) & interest == "cow"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-11-04")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-10-31")), by=time)))}
      
      #temporary data frame to put
      #if statements for difference in data frame setup between days and weeks
      if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                               StartDate = time_seq[1:(length(time_seq)-1)], 
                                               EndDate = time_seq[2:length(time_seq)]-1, 
                                               n_Images = 0)}
      if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i], Date = time_seq , n_Images = 0)}
      
      # loop to sum total number of images per time step
      for(j in 1:(length(time_seq)-1)){
        
        tmp_df[j, "n_Images"] <- nrow(tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                                            & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),])
      }
      cov1_list[[i]] <- tmp_df
      
      #progress bar tick
      progbar$tick()
    }
    progbar$close
    
    return(cov1_list)
  }
  

  #'  Covariate 2 -  Sum number of unique detention events in each time step
  #'  Careful with the input data set since the function will not discriminate by Species/HumanAct
  #'  time = "weeks", "days"
  #'  interest = "human", "cow"
  #'  m = interval used to define unique detection event
  cov2 <- function(data, time, interest, m){
    
    #'  Filter unique detections to pull out the first image of every detection event
    #'  FOR COW: Using output of uniq function as input 
    #'  FOR HUNT: Using first_uniq function
    if(interest == "cow"){
    first_uniq_data <- data %>% 
      group_by(caps) %>% 
      slice(1L) %>%
      ungroup()} else(first_uniq_data <- uniq_first(data, m))
    
    #'  Identify unique camera location names
    locations <- unique(first_uniq_data$CameraLocation)
    
    #' Create list that will be filled for each camera
    cov2_list <- vector("list", length = length(locations))
    names(cov2_list) <- locations
    
    
    #'  loops that outputs a list of data frames with all camera locations
    #'  1(i). goes through each camera location 
    #'  2(j). goes trough each timestep
    #'  each data frames has the raw count in each timestep
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Calculating covariate 2 values" , total = length(locations))
    for(i in 1:length(locations)){
      
      #temporary data set for camera location i
      tmp_data <- first_uniq_data[first_uniq_data$CameraLocation == locations[i],]
      
      #creating a fixed time sequence depending on 'interest'
      #requiring a check to see if last row in tmp_data has the same year as first row
      #if the data runs past a year, will add extra rows to account for another year
      year <- year(tmp_data[1,]$Date)
      year_check <- tmp_data[nrow(tmp_data),]$Date
      if(interest == "human"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
      if(year_check > paste0(year+1,"-7-31") & interest == "human"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-30")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-31")), by=time)))}
      
      
      if(interest == "cow"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
      if(year != year(year_check) & interest == "cow"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-11-04")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-10-31")), by=time)))}
      
      
      #temporary data frame to put
      #if statements for difference in data frame setup between days and weeks
      if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                               StartDate = time_seq[1:(length(time_seq)-1)], 
                                               EndDate = time_seq[2:length(time_seq)]-1, 
                                               n_Detections = 0)}
      if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i],Date = time_seq , n_Detections = 0)}
      
      # loop to sum total number of unique detection events per time step
      for(j in 1:(length(time_seq)-1)){
        
        tmp_df[j, "n_Detections"] <- nrow(tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                                            & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),])
      }
      cov2_list[[i]] <- tmp_df
      
      #progress bar tick
      progbar$tick()
    }
    progbar$close
    
    return(cov2_list)
  }
  

  
  #'  Covariate 3 -  Sum of maximum number of individuals detected in single image 
  #'                 from each unique detection event in each time step
  #'  Careful with the input data set since the function will not discriminate by Species/HumanAct
  #'  time = "weeks", "days"
  #'  interest = "human", "cow"
  #'  m = interval used to define unique detection event
  cov3 <- function(data, time, interest, m){
    
    #'  Filter unique detections to pull out the max count from every detection event
    #'  FOR COW: Using output of uniq function as input 
    #'  FOR HUNT: Using first_uniq function
    if(interest == "cow"){
      max_uniq_data <- data %>% 
        group_by(caps) %>% 
        slice_max(Count) %>%
        slice_head(n=1) %>%
        ungroup()} else(max_uniq_data <- max_uniq(data, m))

    
    #'  Identify unique camera location names
    locations <- unique(max_uniq_data$CameraLocation)
    
    #' creating list that will be filled for each camera
    cov3_list <- vector("list", length = length(locations))
    names(cov3_list) <- locations
    
    
    #'  loops that outputs a list of data frames with all camera locations
    #'  1(i). goes through each camera location 
    #'  2(j). goes trough each timestep
    #'  each data frames has the raw count in each timestep
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Calculating covariate 3 values" , total = length(locations))
    for(i in 1:length(locations)){
      
      #temporary data set for camera location i
      tmp_data <- max_uniq_data[max_uniq_data$CameraLocation == locations[i],]
      
      #creating a fixed time sequence depending on 'interest'
      #requiring a check to see if last row in tmp_data has the same year as first row
      #if the data runs past a year, will add extra rows to account for another year
      year <- year(tmp_data[1,]$Date)
      year_check <- tmp_data[nrow(tmp_data),]$Date
      if(interest == "human"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
      if(year_check > paste0(year+1,"-7-31") & interest == "human"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-30")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-31")), by=time)))}
      
      
      if(interest == "cow"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
      if(year != year(year_check) & interest == "cow"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-11-04")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-10-31")), by=time)))}
      
      
      
      #temporary data frame to put
      #if statements for difference in data frame setup between days and weeks
      if(time == "weeks"){tmp_df <- data.frame(CameraLocation = locations[i],
                                               StartDate = time_seq[1:(length(time_seq)-1)], 
                                               EndDate = time_seq[2:length(time_seq)]-1, 
                                               Count = 0)}
      if(time == "days"){tmp_df <- data.frame(CameraLocation = locations[i],Date = time_seq , Count = 0)}
      
      # loop to sum maximum number of individuals detected in a single image during 
      # each unique detection event per time step
      for(j in 1:(length(time_seq)-1)){
  
        tmp_df[j, "Count"]  <- sum(tmp_data[as.Date(tmp_data$DateTime) >= as.Date(time_seq[j]) 
                                            & as.Date(tmp_data$DateTime) < as.Date(time_seq[j+1]),]$Count)
      }
      cov3_list[[i]] <- tmp_df
      
      #progress bar tick
      progbar$tick()
    }
    progbar$close
    
    return(cov3_list)
  }
  
  
  
  #'  Covariate 4 -  Sum of minutes across unique detection events in each time step
  #'  Careful with the input data set since the function will not discriminate by Species/HumanAct
  #'  time = "weeks", "days"
  #'  interest = "human", "cow"
  #'  m = interval used to define unique detection event
  cov4 <- function(data, time, interest, m){
    
    #'  FOR COW: Using output of uniq function as input 
    #'  FOR HUNT: Using first_uniq function
    if(interest == "cow"){
      uniq_data <- data} else(uniq_data <- uniq(data, m))

    
    #'  Identify unique camera location names
    locations <- unique(uniq_data$CameraLocation)
    
    #' creating list that will be filled for each camera
    cov4_list <- vector("list", length = length(locations))
    names(cov4_list) <- locations
    
    
    #'  loops that outputs a list of data frames with all camera locations
    #'  1(i). goes through each camera location 
    #'  2(j). goes trough each timestep
    #'  3(k). goes through each unique detection calculating the durations
    #'  each data frames has the raw count in each timestep
    progbar <- progress_bar$new(format = "[:bar] :percent in :elapsedfull. Calculating covariate 4 values" , total = length(locations))
    for(i in 1:length(locations)){
      
      #temporary data set for camera location i
      tmp_data <- uniq_data[uniq_data$CameraLocation == locations[i],]
      
      #creating a fixed time sequence depending on 'interest'
      #requiring a check to see if last row in tmp_data has the same year as first row
      #if the data runs past a year, will add extra rows to account for another year
      year <- year(tmp_data[1,]$Date)
      year_check <- tmp_data[nrow(tmp_data),]$Date
      if(interest == "human"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-30")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-8-01")), to=as.Date(paste0(year+1,"-1-31")), by=time))}
      if(year_check > paste0(year+1,"-7-31") & interest == "human"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-30")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-8-01")), to=as.Date(paste0(year(year_check)+1,"-1-31")), by=time)))}
      
      
      if(interest == "cow"){
        if(time == "weeks"){time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-11-04")), by=time)}
        else(time_seq <- seq(from=as.Date(paste0(year,"-4-01")), to=as.Date(paste0(year,"-10-31")), by=time))}
      if(year != year(year_check) & interest == "cow"){
        if(time == "weeks"){time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-11-04")), by=time))}
        else(time_seq <- c(time_seq, seq(from=as.Date(paste0(year(year_check),"-4-01")), to=as.Date(paste0(year(year_check),"-10-31")), by=time)))}
      
      
      
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
        
        #if there are no detection in the time step set duration to 0
        if(nrow(tmp_timestep) != 0){
          for(k in 1:length(detections)){
            
            tmp_detection <- tmp_timestep[tmp_timestep$caps == detections[k],]
            
            #can change time metric
            tmp_duration[k] <- as.numeric(difftime(tmp_detection[nrow(tmp_detection),]$DateTime, tmp_detection[1,]$DateTime, units="secs"))
          }
          #'  Sum amount of time passed within each detection event
          tmp_df[j, "Duration"] <- sum(tmp_duration)
        } else(tmp_df[j, "Duration"] <- 0)
  
      }
      cov4_list[[i]] <- tmp_df
      
      #progress bar tick
      progbar$tick()
    }
    progbar$close
    
    return(cov4_list)
  }

  