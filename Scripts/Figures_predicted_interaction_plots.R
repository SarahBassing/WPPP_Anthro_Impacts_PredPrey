  #'  ================================
  #'  Co-Occurrence plots 
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  June 2022
  #'  ================================
  #'  Plot species interactions using multi-species occupancy model results.
  #'  ================================
    
  #'  Load libraries
  library(unmarked)
  library(ggplot2)
  library(tidyverse)
  
  #'  Read in models
  load(file = "./Outputs/MultiSpp_CoOcc_Models.RData")
  
  #'  Predict effect of elevation on species occupancy
  r <- range(stations$Elev)
  x <- seq(r[1], r[2], length.out = 100)
  x_scaled <- (x - mean(stations$Elev)) / (sd(stations$Elev))
  nd <- data.frame(Elev = x_scaled, GrazingActivity = 0, PercForest = 0, Public = 1, Study_Area = 1)
  #'  Marginal occupancy
  pred.coy <- predict(coy.md.cow_grz, type = "state", species = "coyote", 
                      newdata = nd, se.fit = TRUE, nsims = 10^5) %>%
    mutate(Species = "Coyote")
  pred.md <- predict(coy.md.cow_grz, type = "state", species = "mule_deer", 
                      newdata = nd, se.fit = TRUE, nsims = 10^5) %>%
    mutate(Species = "Mule Deer")
  pred.cattle <- predict(coy.md.cow_grz, type = "state", species = "cattle", 
                     newdata = nd, se.fit = TRUE, nsims = 10^5) %>%
    mutate(Species = "Cattle")
  
  #'  Conditional occupancy when all covariates held at their mean (logit scale)
  nd_mu <- data.frame(Elev = 0, GrazingActivity = 0, PercForest = 0, Public = 1, Study_Area = 1)
  coy.none <- predict(coy.md.cow_grz, type = "state", species = "coyote", 
                      cond = c("-mule_deer", "-cattle"), newdata = nd_mu, se.fit = TRUE, nsims = 10^5)
  coy.md <- predict(coy.md.cow_grz, type = "state", species = "coyote", 
                      cond = c("mule_deer", "-cattle"), newdata = nd_mu, se.fit = TRUE, nsims = 10^5)
  coy.cattle <- predict(coy.md.cow_grz, type = "state", species = "coyote", 
                    cond = c("-mule_deer", "cattle"), newdata = nd_mu, se.fit = TRUE, nsims = 10^5)
  coy.both <- predict(coy.md.cow_grz, type = "state", species = "coyote", 
                        cond = c("mule_deer", "cattle"), newdata = nd_mu, se.fit = TRUE, nsims = 10^5)
  round(cond.occ <- rbind("Neither" = coy.none, "Mule Deer" = coy.md, "Cattle" = coy.cattle,
                          "Both" = coy.both), 2)
  
  #'  Interaction parameter in response to grazing activity
  r <- range(stations$GrazingActivity)
  x <- seq(r[1], r[2], length.out = 100)
  x_scaled <- (x - mean(stations$GrazingActivity)) / (sd(stations$GrazingActivity))
  nd <- data.frame(Elev = 0, GrazingActivity = x_scaled, PercForest = 0, Public = 1, Study_Area = 1)
  pr_nomd <- predict(coy.md.cow_grz, type = "state", species = "coyote", 
                     cond = "-mule_deer", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "Mule deer absent") 
  pr_md <- predict(coy.md.cow_grz, type = "state", species = "coyote", 
                   cond = "mule_deer", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "Mule deer present")
  pr_nocoy <- predict(coy.md.cow_grz, type = "state", species = "mule_deer", 
                     cond = "-coyote", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "Coyote absent")
  pr_coy <- predict(coy.md.cow_grz, type = "state", species = "mule_deer", 
                   cond = "coyote", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "Coyote present")
  
  pr_nowolf <- predict(wolf.wtd.cow_grz, type = "state", species = "wtd", 
                     cond = "-wolf", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "Wolf absent") 
  pr_wolf <- predict(wolf.wtd.cow_grz, type = "state", species = "wtd", 
                   cond = "wolf", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "Wolf present")
  pr_nowtd <- predict(wolf.wtd.cow_grz, type = "state", species = "wolf", 
                      cond = "-wtd", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "White-tailed deer absent")
  pr_wtd <- predict(wolf.wtd.cow_grz, type = "state", species = "wolf", 
                    cond = "wtd", newdata = nd, se.fit = TRUE, nsims = 10^6) %>%
    mutate(Interacting_Species = "White-tailed deer present")
  
  
  grazing <- as.data.frame(x)
  grazing <- rbind(grazing, grazing)
  nd_rep <- rbind(nd, nd)
  nd_rep <- cbind(nd_rep, grazing)
  
  coy_coocc <- rbind(pr_nomd, pr_md)
  coy_coocc <- cbind(nd_rep, coy_coocc)
  md_coocc <- rbind(pr_nocoy, pr_coy)
  md_coocc <- cbind(nd_rep, md_coocc)
  
  wtd_coocc <- rbind(pr_nowolf, pr_wolf)
  wtd_coocc <- cbind(nd_rep, wtd_coocc)
  wolf_coocc <- rbind(pr_nowtd, pr_wtd)
  wolf_coocc <- cbind(nd_rep, wolf_coocc)
  
  
  ggplot(coy_coocc, aes(x = x, y = Predicted, group = Interacting_Species, colour = Interacting_Species)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("#CC6666", "#9999CC")) +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted coyote occupancy") +
    ggtitle("Effect of cattle grazing on coyote occupancy with and without \nmule deer presence")
    
  ggplot(md_coocc, aes(x = x, y = Predicted, group = Interacting_Species, colour = Interacting_Species)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("#CC6666", "#9999CC")) +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted mule deer occupancy") +
    ggtitle("Effect of cattle grazing on mule deer occupancy with and without \ncoyotes present")
  
  ggplot(wolf_coocc, aes(x = x, y = Predicted, group = Interacting_Species, colour = Interacting_Species)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("#CC6666", "#9999CC")) +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted wolf occupancy") +
    ggtitle("Effect of cattle grazing on wolf occupancy with and without \nwhite-tailed deer presence")
  
  ggplot(wtd_coocc, aes(x = x, y = Predicted, group = Interacting_Species, colour = Interacting_Species)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("#CC6666", "#9999CC")) +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted white-tailed deer occupancy") +
    ggtitle("Effect of cattle grazing on white-tailed deer occupancy with \nand without wolves present")
  
 
  
    