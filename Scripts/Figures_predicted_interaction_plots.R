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
  
  scale_cov <- function(cov) {
    #'  Identify range of covariate of interest
    r <- range(cov)
    #'  Create sequence of values starting and ending with range values
    x <- seq(r[1], r[2], length.out = 100)
    #'  Scale values based on covariate mean and standard deviation
    x_scaled <- (x - mean(cov)) / (sd(cov))
    #'  Make single data frame of scaled and unscaled covariate
    x_cov <- as.data.frame(cbind(x, x_scaled))
    
    return(x_cov)
  }
  scaled_elev <- scale_cov(stations$Elev)
  scaled_forest <- scale_cov(stations$PercForest)
  scaled_graze <- scale_cov(stations$GrazingActivity)
  scaled_hunt <- scale_cov(stations$HuntingActivity)
  
  #'  Function to predict marginal probability from each occupancy model
  marginal_occ_g <- function(mod, elev, act, forest, pub, area, spp1, spp2, cov) {
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, GrazingActivity = act, PercForest = forest,
                         Public = pub, Study_Area = area)
    #'  Predict marginal occupancy for each species of interest using top model
    marg_predator <- predict(mod, type = "state", species = spp1, newdata = cov_df,
                             se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Cov = cov) %>% bind_cols(cov_df)
    marg_prey <- predict(mod, type = "state", species = spp2, newdata = cov_df,
                             se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Cov = cov) %>% bind_cols(cov_df)
    marg_cattle <- predict(mod, type = "state", species = "cattle", newdata = cov_df,
                             se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = "Cattle",
             Cov = cov) %>% bind_cols(cov_df)
    
    #'  Create one large data frame with all marginal probabilities
    marg_df <- rbind(marg_predator, marg_prey, marg_cattle)
    
    return(marg_df)

  }
  #'  Predict marginal probabilities for different combinations of species
  #'  NOTE: BE SURE TO USE CORRECT TOP MODEL! 
  #'  Think through which covariates you predict across (scaled_cov, public, study area)
  marg_coug_md_cattle <- marginal_occ_g(coug.md.cow_hab0, elev = scaled_elev[,2], act = 0, 
                                      forest = 0, pub = 1, area = 1, spp1 = "cougar", 
                                      spp2 = "muledeer", cov = scaled_elev[,1])
  marg_coug_elk_cattle <- marginal_occ_g(coug.elk.cow_hab0, elev = scaled_elev[,2], act = 0, 
                                       forest = 0, pub = 1, area = 0, spp1 = "cougar", 
                                       spp2 = "elk", cov = scaled_elev[,1])
  marg_coug_wtd_cattle <- marginal_occ_g(coug.wtd.cow_grz, elev = scaled_elev[,2], act = 0, 
                                       forest = 0, pub = 1, area = 0, spp1 = "cougar", 
                                       spp2 = "wtd", cov = scaled_elev[,1])
  marg_coug_moose_cattle <- marginal_occ_g(coug.moose.cow_grz, elev = scaled_elev[,2], act = 0, 
                                         forest = 0, pub = 1, area = 0, spp1 = "cougar", 
                                         spp2 = "moose", cov = scaled_elev[,1])
  marg_wolf_md_cattle <- marginal_occ_g(wolf.md.cow_grz, elev = scaled_elev[,2], act = 0, 
                                      forest = 0, pub = 1, area = 1, spp1 = "wolf", 
                                      spp2 = "muledeer", cov = scaled_elev[,1])
  marg_wolf_elk_cattle <- marginal_occ_g(wolf.elk.cow_hab0, elev = scaled_elev[,2], act = 0, 
                                       forest = 0, pub = 1, area = 0, spp1 = "wolf", 
                                       spp2 = "elk", cov = scaled_elev[,1])
  marg_wolf_wtd_cattle <- marginal_occ_g(wolf.wtd.cow_grz, elev = scaled_elev[,2], act = 0, 
                                       forest = 0, pub = 1, area = 0, spp1 = "wolf", 
                                       spp2 = "wtd", cov = scaled_elev[,1])
  marg_wolf_moose_cattle <- marginal_occ_g(wolf.moose.cow_hab0, elev = scaled_elev[,2], act = 0, 
                                         forest = 0, pub = 1, area = 0, spp1 = "wolf", 
                                         spp2 = "moose", cov = scaled_elev[,1])
  marg_bear_md_cattle <- marginal_occ_g(bear.md.cow_hab0, elev = scaled_elev[,2], act = 0, 
                                     forest = 0, pub = 1, area = 1, spp1 = "blackbear", 
                                     spp2 = "muledeer", cov = scaled_elev[,1])
  marg_bear_elk_cattle <- marginal_occ_g(bear.elk.cow_hab0, elev = scaled_elev[,2], act = 0, 
                                       forest = 0, pub = 1, area = 0, spp1 = "blackbear", 
                                       spp2 = "elk", cov = scaled_elev[,1])
  marg_bear_wtd_cattle <- marginal_occ_g(bear.wtd.cow_grz, elev = scaled_elev[,2], act = 0, 
                                      forest = 0, pub = 1, area = 0, spp1 = "blackbear", 
                                      spp2 = "wtd", cov = scaled_elev[,1])
  marg_bear_moose_cattle <- marginal_occ_g(bear.moose.cow_grzpub, elev = scaled_elev[,2], act = 0, 
                                         forest = 0, pub = 1, area = 0, spp1 = "blackbear", 
                                         spp2 = "moose", cov = scaled_elev[,1])
  marg_bob_md_cattle <- marginal_occ_g(bob.md.cow_grz, elev = scaled_elev[,2], act = 0, 
                                      forest = 0, pub = 1, area = 1, spp1 = "bobcat", 
                                      spp2 = "mule_deer", cov = scaled_elev[,1])
  marg_bob_wtd_cattle <- marginal_occ_g(bob.wtd.cow_grz, elev = scaled_elev[,2], act = 0, 
                                      forest = 0, pub = 1, area = 0, spp1 = "bobcat", 
                                      spp2 = "wtd", cov = scaled_elev[,1])
  marg_coy_md_cattle <- marginal_occ_g(coy.md.cow_grz, elev = scaled_elev[,2], act = 0, 
                                     forest = 0, pub = 1, area = 1, spp1 = "coyote", 
                                     spp2 = "mule_deer", cov = scaled_elev[,1])
  marg_coy_wtd_cattle <- marginal_occ_g(coy.wtd.cow_grzpub, elev = scaled_elev[,2], act = 0, 
                                      forest = 0, pub = 1, area = 0, spp1 = "coyote", 
                                      spp2 = "wtd", cov = scaled_elev[,1])
  
  
  
  
  ggplot(coy.md.cattle.predictions, aes(x = x, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    ylim(0, 1) +
    xlab("Elevation (m)") + 
    ylab("Marginal occupancy") +
    ggtitle("Marginal occupancy for coyotes, mule deer, and cattle \nas elevation changes")
  
  ggplot(wolf.wtd.cattle.predictions, aes(x = x, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    ylim(0, 1) +
    xlab("Elevation (m)") + 
    ylab("Marginal occupancy") +
    ggtitle("Marginal occupancy for wolves, white-tailed deer, and cattle \nas elevation changes")
  
  
  conditional_occu <- function(mod, pub, area, spp1, spp2, spp3, name1, name2) {
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = 0, GrazingActivity = 0, PercForest = 0,
                         Public = pub, Study_Area = area)
    #'  Create characters for each species that include a "-", necessary for cond argument
    no_spp2 <- paste0("-",spp2); no_spp3 <- paste0("-",spp3)
    spp1_none <- predict(mod, type = "state", species = spp1, cond = c(no_spp2, no_spp3), 
                         newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    spp1_spp2 <- predict(mod, type = "state", species = spp1, cond = c(spp2, no_spp3), 
                         newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    spp1_spp3 <- predict(mod, type = "state", species = spp1, cond = c(no_spp2, spp3), 
                         newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    spp1_both <- predict(mod, type = "state", species = spp1, cond = c(spp2, spp3), 
                         newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    cond.occ <- round(rbind("Neither" = spp1_none, name1 = spp1_spp2, name2 = spp1_spp3,
                            "Both" = spp1_both))
    
    return(cond.occ)
  }
  cond_coy_md_cattle <- conditional_occu(coy.md.cow_grz, pub = 1, area = 1, spp1 = "coyote",
                                         spp2 = "mule_deer", spp3 = "cattle", name1 = "Mule Deer", name2 = "Cattle")
  
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
  
  
  
  
  #'  Function to predict species interactions in response to grazing activity
  spp_interactions_g <- function(mod, elev, act, forest, pub, area, spp1, spp2, cov) {
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, GrazingActivity = act, PercForest = forest,
                         Public = pub, Study_Area = area)
    #'  Create characters for each species that include a "-", necessary for cond argument
    no_spp2 <- paste0("-",spp2); no_spp1 <- paste0("-",spp1)
    
    #'  Predict conditional occupancy when spp2 is absent
    spp2_absent <- predict(mod, type = "state", species = spp1, cond = no_spp2,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp2 is present
    spp2_present <- predict(mod, type = "state", species = spp1, cond = spp2,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Predict conditional occupancy when spp1 is absent
    spp1_absent <- predict(mod, type = "state", species = spp2, cond = no_spp1,
                           newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " absent"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    #'  Predict conditional occupancy when spp1 is present
    spp1_present <- predict(mod, type = "state", species = spp2, cond = spp1,
                            newdata = cov_df, se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " present"),
             Cov = cov) %>% 
      bind_cols(cov_df)
    
    #'  Create one large data frame with all marginal probabilities
    sppX_df <- rbind(spp2_absent, spp2_present, spp1_absent, spp1_present)
    
    return(sppX_df)
  }
  sppX_coy_md_cattle <- spp_interactions_g(coy.md.cow_grz, elev = 0, act = scaled_graze[,2], 
                                     forest = 0, pub = 1, area = 1, spp1 = "coyote", 
                                     spp2 = "mule_deer", cov = scaled_graze[,1])
  sppX_coy_wtd_cattle <- spp_interactions_g(coy.wtd.cow_grzpub, elev = 0, act = scaled_graze[,2], 
                                      forest = 0, pub = 1, area = 0, spp1 = "coyote", 
                                      spp2 = "wtd", cov = scaled_graze[,1])
  
  
  
  
  sppX_coy_wtd <- sppX_coy_wtd_cattle[sppX_coy_wtd_cattle$Species == "coyote",]
  ggplot(sppX_coy_wtd, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("#CC6666", "#9999CC")) +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted coyote occupancy") +
    ggtitle("Effect of cattle grazing on coyote occupancy with and without \nmule deer presence")
  sppX_wtd_coy <- sppX_coy_wtd_cattle[sppX_coy_wtd_cattle$Species == "wtd",]
  ggplot(sppX_wtd_coy, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
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
  
 
  
    