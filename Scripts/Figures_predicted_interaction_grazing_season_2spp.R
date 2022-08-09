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
  library(khroma)
  
  #'  Read in models
  load(file = "./Outputs/MultiSpp_CoOcc_Models_2spp.RData")
  
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
  
  ####  Marginal occupancy probability  ####
  #'  Function to predict marginal probability from each occupancy model
  marginal_occ_g <- function(mod, elev, act, forest, pub, area, spp1, spp2, cov) {
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled or categorical) 
    cov_df <- data.frame(Elev = elev, GrazingActivity = act, PercForest = forest,
                         PublicGrazing = pub, Study_Area = area)
    # cov_df$PublicGrazing <- factor(cov_df$PublicGrazing, levels = c(0, 1))
    # cov_df$Study_Area <- factor(cov_df$Study_Area, levels = c(0, 1))
    #'  Predict marginal occupancy for each species of interest using top model
    marg_predator <- predict(mod, type = "state", species = spp1, newdata = cov_df,
                             se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp1,
             Cov = cov) %>% bind_cols(cov_df)
    marg_prey <- predict(mod, type = "state", species = spp2, newdata = cov_df,
                             se.fit = TRUE, nsims = 10^5) %>%
      mutate(Species = spp2,
             Cov = cov) %>% bind_cols(cov_df)
    
    #'  Create one large data frame with all marginal probabilities
    marg_df <- rbind(marg_predator, marg_prey)
    
    return(marg_df)

  }
  #'  Predict marginal probabilities for different combinations of species
  #'  Think through which covariates you predict across (scaled_cov, public, study area)
  #'  Currently setting public = 1 => on grazing allotments b/c this is generally
  #'  where livestock activity is occurring;
  #'  Setting area = 1 => Okanogan b/c most cattle activity in OK study area  
  #'  except for elk where almost all elk data are from NE (set to 0)
  # marg_coug_md_gs <- marginal_occ_g(gs_cougmd_global, elev = 0, act = scaled_graze[,2],
  #                                     forest = 0, pub = 1, area = 1, spp1 = "cougar",
  #                                     spp2 = "muledeer", cov = scaled_graze[,1])
  # marg_coug_elk_gs <- marginal_occ_g(gs_cougelk_global, elev = 0, act = scaled_graze[,2],
  #                                      forest = 0, pub = 1, area = 0, spp1 = "cougar",
  #                                      spp2 = "elk", cov = scaled_graze[,1])
  marg_coug_wtd_gs <- marginal_occ_g(gs_cougwtd_global, elev = 0, act = scaled_graze[,2],
                                       forest = 0, pub = 1, area = 1, spp1 = "cougar", 
                                       spp2 = "wtd", cov = scaled_graze[,1])
  # marg_coug_moose_gs <- marginal_occ_g(gs_cougmoose_global, elev = 0, act = scaled_graze[,2],
  #                                        forest = 0, pub = 1, area = 1, spp1 = "cougar",
  #                                        spp2 = "moose", cov = scaled_graze[,1])
  # marg_wolf_md_gs <- marginal_occ_g(gs_wolfmd_global, elev = 0, act = scaled_graze[,2],
  #                                     forest = 0, pub = 1, area = 1, spp1 = "wolf",
  #                                     spp2 = "muledeer", cov = scaled_graze[,1])
  # marg_wolf_elk_gs <- marginal_occ_g(gs_wolfelk_global, elev = 0, act = scaled_graze[,2],
  #                                      forest = 0, pub = 0, area = 0, spp1 = "wolf",
  #                                      spp2 = "elk", cov = scaled_graze[,1])
  marg_wolf_wtd_gs <- marginal_occ_g(gs_wolfwtd_global, elev = 0, act = scaled_graze[,2], 
                                       forest = 0, pub = 1, area = 1, spp1 = "wolf", 
                                       spp2 = "wtd", cov = scaled_graze[,1])
  # marg_wolf_moose_gs <- marginal_occ_g(gs_wolfmoose_global, elev = 0, act = scaled_graze[,2],
  #                                        forest = 0, pub = 1, area = 1, spp1 = "wolf",
  #                                        spp2 = "moose", cov = scaled_graze[,1])
  # marg_bear_md_gs <- marginal_occ_g(gs_bearmd_global, elev = 0, act = scaled_graze[,2],
  #                                    forest = 0, pub = 1, area = 1, spp1 = "blackbear",
  #                                    spp2 = "muledeer", cov = scaled_graze[,1])
  # marg_bear_elk_gs <- marginal_occ_g(gs_bearelk_global, elev = 0, act = scaled_graze[,2],
  #                                      forest = 0, pub = 0, area = 0, spp1 = "blackbear",
  #                                      spp2 = "elk", cov = scaled_graze[,1])
  # marg_bear_wtd_gs <- marginal_occ_g(gs_bearwtd_global, elev = 0, act = scaled_graze[,2],
  #                                     forest = 0, pub = 1, area = 1, spp1 = "blackbear",
  #                                     spp2 = "wtd", cov = scaled_graze[,1])
  # marg_bear_moose_gs <- marginal_occ_g(gs_bearmoose_global, elev = 0, act = scaled_graze[,2],
  #                                        forest = 0, pub = 1, area = 1, spp1 = "blackbear",
  #                                        spp2 = "moose", cov = scaled_graze[,1])
  marg_bob_md_gs <- marginal_occ_g(gs_bobmd_global, elev = 0, act = scaled_graze[,2], 
                                      forest = 0, pub = 1, area = 1, spp1 = "bobcat", 
                                      spp2 = "mule_deer", cov = scaled_graze[,1])
  marg_bob_wtd_gs <- marginal_occ_g(gs_bobwtd_global, elev = 0, act = scaled_graze[,2], 
                                      forest = 0, pub = 1, area = 1, spp1 = "bobcat", 
                                      spp2 = "wtd", cov = scaled_graze[,1])
  marg_coy_md_gs <- marginal_occ_g(gs_coymd_global, elev = 0, act = scaled_graze[,2], 
                                     forest = 0, pub = 1, area = 1, spp1 = "coyote", 
                                     spp2 = "mule_deer", cov = scaled_graze[,1])
  marg_coy_wtd_gs <- marginal_occ_g(gs_coywtd_global, elev = 0, act = scaled_graze[,2], 
                                      forest = 0, pub = 1, area = 1, spp1 = "coyote", 
                                      spp2 = "wtd", cov = scaled_graze[,1])
  
  #'  Plot marginal occupancy for each species against changes in cattle activity
  #'  [only for pairings where one or both species were affected by grazing]
  MargOcc_coug_wtd_graze <- ggplot(marg_coug_wtd_gs, aes(x = Cov, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Species), linetype = 0, alpha =0.5) +
    ylim(0, 1) +
    xlab("Cattle activity (mean daily detections)") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Marginal occupancy for cougar and white-tailed deer \nas cattle activity changes")
  MargOcc_coug_wtd_graze
  
  MargOcc_wolf_wtd_graze <- ggplot(marg_wolf_wtd_gs, aes(x = Cov, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Species), linetype = 0, alpha =0.5) +
    ylim(0, 1) +
    xlab("Cattle activity (mean daily detections)") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Marginal occupancy for wolves and white-tailed deer \nas cattle activity changes")
  MargOcc_wolf_wtd_graze
  
  MargOcc_bob_md_graze <- ggplot(marg_bob_md_gs, aes(x = Cov, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Species), linetype = 0, alpha =0.5) +
    ylim(0, 1) +
    xlab("Cattle activity (mean daily detections)") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Marginal occupancy for bobcats and mule deer \nas cattle activity changes")
  MargOcc_bob_md_graze
  
  MargOcc_bob_wtd_graze <- ggplot(marg_bob_wtd_gs, aes(x = Cov, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Species), linetype = 0, alpha =0.5) +
    ylim(0, 1) +
    xlab("Cattle activity (mean daily detections)") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Marginal occupancy for bobcats and white-tailed deer \nas cattle activity changes")
  MargOcc_bob_wtd_graze
  
  MargOcc_coy_md_graze <- ggplot(marg_coy_md_gs, aes(x = Cov, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Species), linetype = 0, alpha =0.5) +
    ylim(0, 1) +
    xlab("Cattle activity (mean daily detections)") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Marginal occupancy for coyotes and mule deer \nas cattle activity changes")
  MargOcc_coy_md_graze
  
  MargOcc_coy_wtd_graze <- ggplot(marg_coy_wtd_gs, aes(x = Cov, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Species), linetype = 0, alpha =0.5) +
    ylim(0, 1) +
    xlab("Cattle activity (mean daily detections)") + 
    ylab("Marginal occupancy probability") +
    scale_fill_discrete(labels = c("Coyote", "White-tailed \ndeer")) +
    scale_color_discrete(labels = c("Coyote", "White-tailed \ndeer")) +
    ggtitle("Marginal occupancy probabilities in response to cattle activity")
  MargOcc_coy_wtd_graze
  
  ggsave("./Outputs/Figures/MargOcc_coug_wtd_graze.tiff", MargOcc_coug_wtd_graze, 
         units = "in", width = 5, height = 4, dpi = 600, device = 'tiff') #, compression = 'lzw'
  ggsave("./Outputs/Figures/MargOcc_wolf_wtd_graze.tiff", MargOcc_wolf_wtd_graze, 
         units = "in", width = 5, height = 4, dpi = 600, device = 'tiff') #, compression = 'lzw'
  ggsave("./Outputs/Figures/MargOcc_bob_md_graze.tiff", MargOcc_bob_md_graze, 
         units = "in", width = 5, height = 4, dpi = 600, device = 'tiff') #, compression = 'lzw'
  ggsave("./Outputs/Figures/MargOcc_bob_wtd_graze.tiff", MargOcc_bob_wtd_graze, 
         units = "in", width = 5, height = 4, dpi = 600, device = 'tiff') #, compression = 'lzw'
  ggsave("./Outputs/Figures/MargOcc_coy_md_graze.tiff", MargOcc_coy_wtd_graze, 
         units = "in", width = 5, height = 4, dpi = 600, device = 'tiff') #, compression = 'lzw'
  ggsave("./Outputs/Figures/MargOcc_coy_wtd_graze.tiff", MargOcc_coy_wtd_graze, 
         units = "in", width = 5, height = 4, dpi = 600, device = 'tiff') #, compression = 'lzw'
  
  #'  Choose colorblind-friendly scheme
  bright <- colour("bright")
  bright(5)
  
  #'  Merge predator and prey responses together for facet wrap plots
  pred_marg_data <- rbind(marg_bob_wtd_gs, marg_coug_wtd_gs) %>% #, marg_wolf_wtd_gs
    #'  Remove wtd data from 2 predator data sets since marginal wtd occ is basically
    #'  the same across all pred-prey models
    filter(!Species == "wtd")
  marg_data <- rbind(pred_marg_data, marg_coy_wtd_gs)
  newlabs <- c("bobcat" = "Bobcat", "cougar" = "Cougar", "coyote" = "Coyote", "wtd" = "White-tailed deer") #"wolf" = "Wolf", 
  #newlabs <- c("coyote" = "Coyote", "wtd" = "White-tailed deer")
  marginal_graze_facet <- ggplot(marg_data, aes(x = Cov, y = Predicted, group = Species, colour = Species)) +
    geom_line(size = 1) +
      scale_colour_bright() +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Species), linetype = 0, alpha =0.5) +
      scale_fill_bright() +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs), nrow = 2) + 
    theme(legend.position = "none") +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Effect of cattle activity on species-specific marginal occupancy")
  marginal_graze_facet
  
  ggsave("./Outputs/Figures/MargOcc_grazing_effect.tiff", marginal_graze_facet, 
         units = "in", width = 6, height = 4, dpi = 300, device = 'tiff') #, compression = 'lzw'
  
  
  ####  Allotment effect on marginal occupancy probabilities  ####
  #'  Function to calculate mean marginal occupancy probability for each species
  #'  when they are on public grazing allotments vs not
  mean_marginal_occu_g <- function(mod, area, spp1, spp2) {
    #'  Create data frame while holding all covs at their mean (0 when scaled) 
    #'  and desired category (0 or 1)
    cov_df0 <- data.frame(Elev = 0, GrazingActivity = 0, PercForest = 0,
                         PublicGrazing = 0, Study_Area = area)
    cov_df1 <- data.frame(Elev = 0, GrazingActivity = 0, PercForest = 0,
                          PublicGrazing = 1, Study_Area = area)
    #'  Predict marginal occupancy based on data set
    spp1_noallot <- predict(mod, type = "state", species = spp1, newdata = cov_df0,  
                            se.fit = TRUE, nsims = 10^5)
    spp1_allot <- predict(mod, type = "state", species = spp1, newdata = cov_df1, 
                          se.fit = TRUE, nsims = 10^5)
    spp2_noallot <- predict(mod, type = "state", species = spp2, newdata = cov_df0,  
                            se.fit = TRUE, nsims = 10^5)
    spp2_allot <- predict(mod, type = "state", species = spp2, newdata = cov_df1, 
                          se.fit = TRUE, nsims = 10^5)
    marg.occ <- round(rbind("Spp1_noallotment" = spp1_noallot, "Spp2_noallotment" = spp2_noallot, 
                            "Spp1_allotment" = spp1_allot,"Spp2_allotment" = spp2_allot), 2)
    Allotment_use <- c("No", "No", "Yes", "Yes")
    Species <- c(spp1, spp2, spp1, spp2)
    marg.occ <- cbind(Allotment_use, Species, marg.occ)
    return(marg.occ)
  }
  (marg_coug_elk_allot <- mean_marginal_occu_g(gs_cougelk_global, area = 0,
                                               spp1 = "cougar", spp2 = "elk"))
  (marg_coug_moose_allot <- mean_marginal_occu_g(gs_cougmoose_global, area = 1,
                                               spp1 = "cougar", spp2 = "moose"))
  (marg_wolf_moose_allot <- mean_marginal_occu_g(gs_wolfmoose_global, area = 1,
                                                 spp1 = "wolf", spp2 = "moose"))
  
  (marg_bear_elk_allot <- mean_marginal_occu_g(gs_bearelk_global, area = 0,
                                               spp1 = "blackbear", spp2 = "elk"))
  (marg_bear_moose_allot <- mean_marginal_occu_g(gs_bearmoose_global, area = 1,
                                                 spp1 = "blackbear", spp2 = "moose"))
 
  marg_coug_elk_allot$Species <- factor(marg_coug_elk_allot$Species, levels = c("cougar","elk"))
  ggplot(marg_coug_elk_allot, aes(x = Allotment_use, y = Predicted, group = Species, colour = Species)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_color_manual(values = c("cougar" = "#D55E00", "elk" = "#0072B2"), labels = c("Cougar", "Elk")) +
    ylim(0, 1) +
    xlab("Public grazing allotment use") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Effect public grazing allotments on \ncougar and elk occurrence")
  marg_coug_moose_allot$Species <- factor(marg_coug_moose_allot$Species, levels = c("cougar","moose"))
  MargOcc_coug_moose_allot <- ggplot(marg_coug_moose_allot, aes(x = Allotment_use, y = Predicted, group = Species, colour = Species)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_color_manual(values = c("cougar" = "#D55E00", "moose" = "#0072B2"), labels = c("Coguar", "Moose")) +
    ylim(0, 1) +
    xlab("Public grazing allotment use") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Effect public grazing allotments on \ncougar and moose occurrence")
  marg_wolf_moose_allot$Species <- factor(marg_wolf_moose_allot$Species, levels = c("wolf","moose"))
  MargOcc_wolf_moose_allot <- ggplot(marg_wolf_moose_allot, aes(x = Allotment_use, y = Predicted, group = Species, colour = Species)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_color_manual(values = c("wolf" = "#D55E00", "moose" = "#0072B2"), labels = c("Wolf", "Moose")) +
    ylim(0, 1) +
    xlab("Public grazing allotment use") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Effect public grazing allotments on \nwolf and moose occurrence")
  marg_bear_elk_allot$Species <- factor(marg_bear_elk_allot$Species, levels = c("blackbear","elk"))
  ggplot(marg_bear_elk_allot, aes(x = Allotment_use, y = Predicted, group = Species, colour = Species)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_color_manual(values = c("blackbear" = "#D55E00", "elk" = "#0072B2"), labels = c("Black bear", "Elk")) +
    ylim(0, 1) +
    xlab("Public grazing allotment use") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Effect public grazing allotments on \nblack bear and elk occurrence")
  marg_bear_moose_allot$Species <- factor(marg_bear_moose_allot$Species, levels = c("blackbear","moose"))
  MargOcc_bear_moose_allot <- ggplot(marg_bear_moose_allot, aes(x = Allotment_use, y = Predicted, group = Species, colour = Species)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_color_manual(values = c("blackbear" = "#D55E00", "moose" = "#0072B2"), labels = c("Black bear", "Moose")) +
    ylim(0, 1) +
    xlab("Public grazing allotment use") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Effect public grazing allotments on \nblack bear and moose occurrence")
  
  ggsave("./Outputs/Figures/MargOcc_coug_moose_allot.tiff", MargOcc_coug_moose_allot, 
         units = "in", width = 5, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/MargOcc_wolf_moose_allot.tiff", MargOcc_wolf_moose_allot, 
         units = "in", width = 5, height = 4, dpi = 300, device = 'tiff') #, compression = 'lzw'
  ggsave("./Outputs/Figures/MargOcc_bear_moose_allot.tiff", MargOcc_bear_moose_allot, 
         units = "in", width = 5, height = 4, dpi = 300, device = 'tiff') #, compression = 'lzw'
  
  #'  Merge predator and prey responses together for facet wrap plots
  marg_coug_elk_allot$Species <- factor(marg_coug_elk_allot$Species, levels = c("cougar","elk"))
  marg_coug_moose_allot$Species <- factor(marg_coug_moose_allot$Species, levels = c("cougar","moose"))
  marg_wolf_moose_allot$Species <- factor(marg_wolf_moose_allot$Species, levels = c("wolf","moose"))
  marg_bear_elk_allot$Species <- factor(marg_bear_elk_allot$Species, levels = c("blackbear","elk"))
  marg_bear_moose_allot$Species <- factor(marg_bear_moose_allot$Species, levels = c("blackbear","moose"))
  marg_coug_elk_allot$Interaction <- "cougar-elk"
  marg_coug_moose_allot$Interaction <- "cougar-moose"
  marg_wolf_moose_allot$Interaction <- "wolf-moose"
  marg_bear_elk_allot$Interaction <- "bear-elk"
  marg_bear_moose_allot$Interaction <- "bear-moose"
  
  # marg_allot_data <- rbind(marg_coug_elk_allot, marg_coug_moose_allot, marg_wolf_moose_allot, marg_bear_elk_allot, marg_bear_moose_allot)
  # marg_allot_data$Species <- factor(marg_allot_data$Species, levels = c("blackbear", "cougar", "wolf", "elk", "moose"))
  # newlabs <- c("bear-elk" = "Black bear-Elk", "bear-moose" = "Black bear-Moose", "cougar-elk" = "Cougar-Elk", "cougar-moose" = "Cougar-Moose", "wolf-moose" = "Wolf-Moose")
  marg_allot_data <- rbind(marg_coug_moose_allot, marg_wolf_moose_allot, marg_bear_moose_allot)
  marg_allot_data$Species <- factor(marg_allot_data$Species, levels = c("blackbear", "cougar", "wolf", "moose"))
  newlabs <- c("bear-moose" = "Black bear-Moose", "cougar-moose" = "Cougar-Moose", "wolf-moose" = "Wolf-Moose")
  MargOcc_allot_facet <- ggplot(marg_allot_data, aes(x = Allotment_use, y = Predicted, group = Species, colour = Species)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_colour_bright(labels = c("Black bear", "Cougar", "Wolf", "Moose")) + #"Elk", "Moose"
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Interaction, scales = "free_y", labeller = as_labeller(newlabs), nrow = 1) + #nrow = 3
    #'  Move legend to empty position in facet, only if nrow = 3 above
    #theme(legend.position = c(1, 0), legend.justification = c(1.5, 0)) +
    #theme(legend.position = "bottom") +
    xlab("On public grazing allotments") + 
    ylab("Marginal occupancy probability") +
    ggtitle("Effect public grazing allotments on marginal occupancy")
  MargOcc_allot_facet
  
  ggsave("./Outputs/Figures/MargOcc_allot_effect.tiff", MargOcc_allot_facet, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  # width = 5.5, height = 6, if nrow = 3
  
  ####  Mean conditional occupancy probabilities  ####
  #'  Function to calculate mean conditional occupancy probability for each species
  #'  with and without the presence of the other
  conditional_occu_g <- function(mod, pub, area, spp1, spp2) { #, spp3
    #'  Create data frame while holding all covs at their mean (0 when scaled) 
    #'  and desired category (0 or 1)
    cov_df <- data.frame(Elev = 0, GrazingActivity = 0, PercForest = 0,
                         PublicGrazing = pub, Study_Area = area)
    #Public = factor(pub, levels = c(0, 1)), Study_Area = factor(area, levels = c(0, 1)))
    #'  Create characters for each species that include a "-", necessary for cond 
    #'  argument in predict when species is not present
    no_spp2 <- paste0("-",spp2)
    spp1_only <- predict(mod, type = "state", species = spp1, cond = no_spp2, 
                         newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    spp1_spp2 <- predict(mod, type = "state", species = spp1, cond = spp2, 
                         newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    cond.occ <- round(rbind("Spp2 absent" = spp1_only, "Spp2 present" = spp1_spp2, 2))
    # no_spp3 <- paste0("-",spp3)
    # spp1_none <- predict(mod, type = "state", species = spp1, cond = c(no_spp2, no_spp3), 
    #                      newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    # spp1_spp2 <- predict(mod, type = "state", species = spp1, cond = c(spp2, no_spp3), 
    #                      newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    # spp1_spp3 <- predict(mod, type = "state", species = spp1, cond = c(no_spp2, spp3), 
    #                      newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    # spp1_both <- predict(mod, type = "state", species = spp1, cond = c(spp2, spp3), 
    #                      newdata = cov_df, se.fit = TRUE, nsims = 10^5)
    # cond.occ <- round(rbind("Neither" = spp1_none, "Spp2 present" = spp1_spp2, 
    #                         "Spp3 present" = spp1_spp3,"Both" = spp1_both), 2)
    
    return(cond.occ)
  }
  #'  Estimate conditional occupancy for pairings of species
  (cond_coug_wtd_gs <- conditional_occu_g(gs_cougwtd_global, pub = 1, area = 1, 
                                          spp1 = "cougar", spp2 = "wtd"))
  (cond_wolf_wtd_gs <- conditional_occu_g(gs_wolfwtd_global, pub = 1, area = 1, 
                                          spp1 = "cougar", spp2 = "wtd"))
  (cond_bob_md_gs <- conditional_occu_g(gs_bobmd_global, pub = 1, area = 1, 
                                        spp1 = "coyote", spp2 = "mule_deer")) 
  (cond_bob_wtd_gs <- conditional_occu_g(gs_bobwtd_global, pub = 1, area = 1, 
                                         spp1 = "coyote", spp2 = "wtd"))
  (cond_coy_md_gs <- conditional_occu_g(gs_coymd_global, pub = 1, area = 1, 
                                        spp1 = "coyote", spp2 = "mule_deer")) 
  (cond_coy_wtd_gs <- conditional_occu_g(gs_coywtd_global, pub = 1, area = 1, 
                                         spp1 = "coyote", spp2 = "wtd")) 
  
  
  ####  Covariate effects on species interaction term  ####
  #'  Function to predict species interactions in response to covariate of interest
  spp_interactions_g <- function(mod, elev, act, forest, pub, area, spp1, spp2, cov) {
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, GrazingActivity = act, PercForest = forest,
                         PublicGrazing = pub, Study_Area = area)
    #Public = factor(pub, levels = c(0, 1)), Study_Area = factor(area, levels = c(0, 1)))
    #'  Create characters for each species that include a "-", necessary for cond 
    #'  argument in predict when species is not present
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
  #'  Predict cond. occupancy for each pairwise interaction over range of cattle activity
  #'  NOTE: setting all other continuous variables to 0 (their mean value) and 
  #'  setting categorical variables as either 0 or 1
  #'  For elk, study area = 0 (NE) because there were so few detections in OK
  #'  For md, wtd, and moose, study area = 1 (OK) b/c that is where the bulk of
  #'  cattle activity occurred and there were sufficient detections of these species
  #'  "pub" variable set to active grazing allotments (1)
  # sppX_coug_md_cattle <- spp_interactions_g(gs_cougmd_global, elev = 0, act = scaled_graze[,2], 
  #                                          forest = 0, pub = 1, area = 1, spp1 = "cougar", 
  #                                          spp2 = "muledeer", cov = scaled_graze[,1])
  # sppX_coug_elk_cattle <- spp_interactions_g(gs_cougelk_global, elev = 0, act = scaled_graze[,2], 
  #                                            forest = 0, pub = 1, area = 0, spp1 = "cougar", 
  #                                            spp2 = "elk", cov = scaled_graze[,1])
  sppX_coug_wtd_cattle <- spp_interactions_g(gs_cougwtd_global, elev = 0, act = scaled_graze[,2], 
                                            forest = 0, pub = 1, area = 1, spp1 = "cougar", 
                                            spp2 = "wtd", cov = scaled_graze[,1])
  # sppX_coug_moose_cattle <- spp_interactions_g(gs_cougmoose_global, elev = 0, act = scaled_graze[,2], 
  #                                            forest = 0, pub = 1, area = 1, spp1 = "cougar", 
  #                                            spp2 = "moose", cov = scaled_graze[,1])
  # sppX_wolf_md_cattle <- spp_interactions_g(gs_wolfmd_global, elev = 0, act = scaled_graze[,2], 
  #                                           forest = 0, pub = 1, area = 1, spp1 = "wolf", 
  #                                           spp2 = "muledeer", cov = scaled_graze[,1])
  # sppX_wolf_elk_cattle <- spp_interactions_g(gs_wolfelk_global, elev = 0, act = scaled_graze[,2], 
  #                                            forest = 0, pub = 1, area = 0, spp1 = "wolf", 
  #                                            spp2 = "elk", cov = scaled_graze[,1])
  sppX_wolf_wtd_cattle <- spp_interactions_g(gs_wolfwtd_global, elev = 0, act = scaled_graze[,2], 
                                             forest = 0, pub = 1, area = 1, spp1 = "wolf", 
                                             spp2 = "wtd", cov = scaled_graze[,1])
  # sppX_wolf_moose_cattle <- spp_interactions_g(gs_wolfmoose_global, elev = 0, act = scaled_graze[,2], 
  #                                              forest = 0, pub = 1, area = 1, spp1 = "wolf", 
  #                                              spp2 = "moose", cov = scaled_graze[,1])
  # sppX_bear_md_cattle <- spp_interactions_g(gs_bearmd_global, elev = 0, act = scaled_graze[,2], 
  #                                           forest = 0, pub = 1, area = 1, spp1 = "blackbear", 
  #                                           spp2 = "muledeer", cov = scaled_graze[,1])
  # sppX_bear_elk_cattle <- spp_interactions_g(gs_bearelk_global, elev = 0, act = scaled_graze[,2], 
  #                                            forest = 0, pub = 1, area = 0, spp1 = "blackbear", 
  #                                            spp2 = "elk", cov = scaled_graze[,1])
  # sppX_bear_wtd_cattle <- spp_interactions_g(gs_bearwtd_global, elev = 0, act = scaled_graze[,2], 
  #                                            forest = 0, pub = 1, area = 1, spp1 = "blackbear", 
  #                                            spp2 = "wtd", cov = scaled_graze[,1])
  # sppX_bear_moose_cattle <- spp_interactions_g(gs_bearmoose_global, elev = 0, act = scaled_graze[,2], 
  #                                              forest = 0, pub = 1, area = 1, spp1 = "blackbear", 
  #                                              spp2 = "moose", cov = scaled_graze[,1])
  sppX_bob_md_cattle <- spp_interactions_g(gs_bobmd_global, elev = 0, act = scaled_graze[,2], 
                                           forest = 0, pub = 1, area = 1, spp1 = "bobcat", 
                                           spp2 = "mule_deer", cov = scaled_graze[,1])
  sppX_bob_wtd_cattle <- spp_interactions_g(gs_bobwtd_global, elev = 0, act = scaled_graze[,2], 
                                            forest = 0, pub = 1, area = 1, spp1 = "bobcat", 
                                            spp2 = "wtd", cov = scaled_graze[,1])
  sppX_coy_md_cattle <- spp_interactions_g(gs_coymd_global, elev = 0, act = scaled_graze[,2], 
                                     forest = 0, pub = 1, area = 1, spp1 = "coyote", 
                                     spp2 = "mule_deer", cov = scaled_graze[,1])
  sppX_coy_wtd_cattle <- spp_interactions_g(gs_coywtd_global, elev = 0, act = scaled_graze[,2], 
                                      forest = 0, pub = 1, area = 1, spp1 = "coyote", 
                                      spp2 = "wtd", cov = scaled_graze[,1])

  ####  Visualize conditional occupancy  ####
  #'  Plot each pairing from the perspective of the predator, then the prey
  #'  Only focus on pairings where the top model contained an interaction term
  sppX_coug_wtd <- sppX_coug_wtd_cattle[sppX_coug_wtd_cattle$Species == "cougar",]
  OccX_coug_wtd_graze <- ggplot(sppX_coug_wtd, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted cougar occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_coug_wtd_graze
  sppX_wtd_coug <- sppX_coug_wtd_cattle[sppX_coug_wtd_cattle$Species == "wtd",]
  OccX_wtd_coug_graze <- ggplot(sppX_wtd_coug, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("cougar absent" = "#CC6666", "cougar present" = "#9999CC"), labels = c("Absent", "Present"), name = "Cougar") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("cougar absent" = "#CC6666", "cougar present" = "#9999CC"), labels = c("Absent", "Present"), name = "Cougar") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted white-tailed deer occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_wtd_coug_graze
  
  sppX_wolf_wtd <- sppX_wolf_wtd_cattle[sppX_wolf_wtd_cattle$Species == "wolf",]
  OccX_wolf_wtd_graze <- ggplot(sppX_wolf_wtd, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") +
    ylab("Predicted wolf occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_wolf_wtd_graze
  sppX_wtd_wolf <- sppX_wolf_wtd_cattle[sppX_wolf_wtd_cattle$Species == "wtd",]
  OccX_wtd_wolf_graze <- ggplot(sppX_wtd_wolf, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("wolf absent" = "#CC6666", "wolf present" = "#9999CC"), labels = c("Absent", "Present"), name = "Wolves") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("wolf absent" = "#CC6666", "wolf present" = "#9999CC"), labels = c("Absent", "Present"), name = "Wolves") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") +
    ylab("Predicted white-tailed deer occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_wtd_wolf_graze
  
  sppX_bob_wtd <- sppX_bob_wtd_cattle[sppX_bob_wtd_cattle$Species == "bobcat",]
  OccX_bob_wtd_graze <- ggplot(sppX_bob_wtd, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted bobcat occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_bob_wtd_graze
  sppX_wtd_bob <- sppX_bob_wtd_cattle[sppX_bob_wtd_cattle$Species == "wtd",]
  OccX_wtd_bob_graze <- ggplot(sppX_wtd_bob, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("bobcat absent" = "#CC6666", "bobcat present" = "#9999CC"), labels = c("Absent", "Present"), name = "Bobcat") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("bobcat absent" = "#CC6666", "bobcat present" = "#9999CC"), labels = c("Absent", "Present"), name = "Bobcat") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted white-tailed deer occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_wtd_bob_graze
  
  sppX_coy_md <- sppX_coy_md_cattle[sppX_coy_md_cattle$Species == "coyote",]
  OccX_coy_md_graze <- ggplot(sppX_coy_md, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("mule_deer absent" = "#CC6666", "mule_deer present" = "#9999CC"), labels = c("Absent", "Present"), name = "Mule deer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("mule_deer absent" = "#CC6666", "mule_deer present" = "#9999CC"), labels = c("Absent", "Present"), name = "Mule deer") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted coyote occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_coy_md_graze
  sppX_md_coy <- sppX_coy_md_cattle[sppX_coy_md_cattle$Species == "mule_deer",]
  OccX_md_coy_graze <- ggplot(sppX_md_coy, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("coyote absent" = "#CC6666", "coyote present" = "#9999CC"), labels = c("Absent", "Present"), name = "Coyote") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("coyote absent" = "#CC6666", "coyote present" = "#9999CC"), labels = c("Absent", "Present"), name = "Coyote") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted mule deer occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_md_coy_graze
  
  sppX_coy_wtd <- sppX_coy_wtd_cattle[sppX_coy_wtd_cattle$Species == "coyote",]
  OccX_coy_wtd_graze <- ggplot(sppX_coy_wtd, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed \ndeer") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted coyote occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_coy_wtd_graze
  sppX_wtd_coy <- sppX_coy_wtd_cattle[sppX_coy_wtd_cattle$Species == "wtd",]
  OccX_wtd_coy_graze <- ggplot(sppX_wtd_coy, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("coyote absent" = "#CC6666", "coyote present" = "#9999CC"), labels = c("Absent", "Present"), name = "Coyote") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("coyote absent" = "#CC6666", "coyote present" = "#9999CC"), labels = c("Absent", "Present"), name = "Coyote") +
    ylim(0, 1) +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Predicted white-tailed deer occupancy") +
    ggtitle("Effect of cattle activity on species interactions")
  OccX_wtd_coy_graze
  
  ggsave("./Outputs/Figures/OccX_coug_wtd_graze.tiff", OccX_coug_wtd_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_wtd_coug_graze.tiff", OccX_wtd_coug_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_wolf_wtd_graze.tiff", OccX_wolf_wtd_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_wtd_wolf_graze.tiff", OccX_wtd_wolf_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_bob_wtd_graze.tiff", OccX_bob_wtd_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_wtd_bob_graze.tiff", OccX_wtd_bob_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_coy_md_graze.tiff", OccX_coy_md_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_md_coy_graze.tiff", OccX_md_coy_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_coy_wtd_graze.tiff", OccX_coy_wtd_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_wtd_coy_graze.tiff", OccX_wtd_coy_graze, 
         units = "in", width = 6, height = 4, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  
  #'  Merge predator and prey responses together for facet wrap plots
  pred_data <- rbind(sppX_coug_wtd, sppX_bob_wtd, sppX_coy_wtd, sppX_wolf_wtd)
  newlabs <- c("bobcat" = "Predator: Bobcat", "cougar" = "Predator: Cougar", "coyote" = "Predator: Coyote", "wolf" = "Predator: Wolf")
  pred_wtd_graze_facet <- ggplot(pred_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed deer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed deer") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of cattle activity on predator co-occurrence with white-tailed deer")
  
  prey_data <- rbind(sppX_wtd_bob, sppX_wtd_coug, sppX_wtd_coy, sppX_wtd_wolf) %>%
    mutate(Predator = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction))
  newlabs <- c("bobcat" = "Predator: Bobcat", "cougar" = "Predator: Cougar", "coyote" = "Predator: Coyote", "wolf" = "Predator: Wolf")
  wtd_pred_graze_facet <- ggplot(prey_data, aes(x = Cov, y = Predicted, group = InterX, colour = InterX)) +
    geom_line(size = 1) +
    scale_color_manual(values=c("absent" = "#CC6666", "present" = "#9999CC"), labels = c("Absent", "Present", ""), name = "Predator") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = InterX), linetype = 0, alpha =0.5) +
    scale_fill_manual(values=c("absent" = "#CC6666", "present" = "#9999CC"), labels = c("Absent", "Present", ""), name = "Predator") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Predator, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional white-tailed deer occupancy") +
    ggtitle("Effect of cattle activity on white-tailed deer co-occurrence with predators")
  wtd_pred_graze_facet
  
  ggsave("./Outputs/Figures/OccX_pred_wtd_graze.tiff", pred_wtd_graze_facet, 
         units = "in", width = 6.5, height = 5, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_wtd_pred_graze.tiff", wtd_pred_graze_facet, 
         units = "in", width = 6.5, height = 5, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  
  
  ####  Allotment effect on conditional occupancy probabilities  ####
  #'  Function to calculate conditional occupancy probabilities for each species
  #'  when they are on public grazing allotments vs not
  allot_conditional_occu_g <- function(mod, area, spp1, spp2, n) {
    #'  Create data frame while holding all covs at their mean (0 when scaled) 
    #'  and desired category (0 or 1)
    cov_df0 <- data.frame(Elev = 0, GrazingActivity = 0, PercForest = 0,
                          PublicGrazing = 0, Study_Area = area)
    cov_df1 <- data.frame(Elev = 0, GrazingActivity = 0, PercForest = 0,
                          PublicGrazing = 1, Study_Area = area)
    #'  Create characters for each species that include a "-", necessary for cond 
    #'  argument in predict when species is not present
    no_spp2 <- paste0("-",spp2); no_spp1 <- paste0("-",spp1)
    
    #'  Predict conditional occupancy when spp2 is absent
    spp2_absent0 <- predict(mod, type = "state", species = spp1, cond = no_spp2,
                           newdata = cov_df0, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " absent"),
             `Interacting species` = "Absent",
             Public_Grazing = "Not permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df0)
    spp2_absent1 <- predict(mod, type = "state", species = spp1, cond = no_spp2,
                            newdata = cov_df1, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " absent"),
             `Interacting species` = "Absent",
             Public_Grazing = "Permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df1)
    #'  Predict conditional occupancy when spp2 is present
    spp2_present0 <- predict(mod, type = "state", species = spp1, cond = spp2,
                            newdata = cov_df0, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " present"),
             `Interacting species` = "Present",
             Public_Grazing = "Not permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df0)
    spp2_present1 <- predict(mod, type = "state", species = spp1, cond = spp2,
                            newdata = cov_df1, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp1,
             Interaction = paste0(spp2, " present"),
             `Interacting species` = "Present",
             Public_Grazing = "Permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df1)
    
    #'  Predict conditional occupancy when spp1 is absent
    spp1_absent0 <- predict(mod, type = "state", species = spp2, cond = no_spp1,
                           newdata = cov_df0, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " absent"),
             `Interacting species` = "Absent",
             Public_Grazing = "Not permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df0)
    spp1_absent1 <- predict(mod, type = "state", species = spp2, cond = no_spp1,
                            newdata = cov_df1, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " absent"),
             `Interacting species` = "Absent",
             Public_Grazing = "Permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df1)
    #'  Predict conditional occupancy when spp1 is present
    spp1_present0 <- predict(mod, type = "state", species = spp2, cond = spp1,
                            newdata = cov_df0, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " present"),
             `Interacting species` = "Present",
             Public_Grazing = "Not permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df0)
    spp1_present1 <- predict(mod, type = "state", species = spp2, cond = spp1,
                             newdata = cov_df1, se.fit = TRUE, nsims = n) %>%
      mutate(Species = spp2,
             Interaction = paste0(spp1, " present"),
             `Interacting species` = "Present",
             Public_Grazing = "Permitted",
             SpeciesPairing = paste0(spp1, "-", spp2)) %>% 
      bind_cols(cov_df1)
    
    #'  Create one large data frame with all marginal probabilities
    sppX_allot_df <- rbind(spp2_absent0, spp2_absent1, spp2_present0, spp2_present1, spp1_absent0, spp1_absent1, spp1_present0, spp1_present1) %>%
      mutate(Predicted = round(Predicted, 2),
             SE = round(SE, 2),
             lower = round(lower, 2),
             upper = round(upper, 2))
    
    return(sppX_allot_df)
  }
  # sppX_coug_md_allot <- allot_conditional_occu_g(gs_cougmd_global, area = 1, spp1 = "cougar", spp2 = "muledeer", n = 10^5)
  sppX_coug_elk_allot <- allot_conditional_occu_g(gs_cougelk_global, area = 0, spp1 = "cougar", spp2 = "elk", n = 10^5)
  sppX_coug_elk_allot$Species <- factor(sppX_coug_elk_allot$Species, levels = c("cougar", "elk"))
  # sppX_coug_wtd_allot <- allot_conditional_occu_g(gs_cougwtd_global, area = 1, spp1 = "cougar", spp2 = "wtd", n = 10^5)
  sppX_coug_moose_allot <- allot_conditional_occu_g(gs_cougmoose_global, area = 1, spp1 = "cougar", spp2 = "moose", n = 10^5)
  # sppX_wolf_md_allot <- allot_conditional_occu_g(gs_wolfmd_global, area = 1, spp1 = "wolf", spp2 = "md", n = 10^5)
  # sppX_wolf_elk_allot <- allot_conditional_occu_g(gs_wolfelk_global, area = 0, spp1 = "wolf", spp2 = "elk", n = 10^5)
  # sppX_wolf_wtd_allot <- allot_conditional_occu_g(gs_wolfwtd_global, area = 1, spp1 = "wolf", spp2 = "wtd", n = 10^5)
  sppX_wolf_moose_allot <- allot_conditional_occu_g(gs_wolfmoose_global, area = 1, spp1 = "wolf", spp2 = "moose", n = 10^5)
  sppX_wolf_moose_allot$Species <- factor(sppX_wolf_moose_allot$Species, levels = c("wolf", "moose"))
  # sppX_bear_md_allot <- allot_conditional_occu_g(gs_bearmd_global, area = 1, spp1 = "blackbear", spp2 = "md", n = 10^5)
  sppX_bear_elk_allot <- allot_conditional_occu_g(gs_bearelk_global, area = 0, spp1 = "blackbear", spp2 = "elk", n = 10^5)
  # sppX_bear_wtd_allot <- allot_conditional_occu_g(gs_bearwtd_global, area = 1, spp1 = "blackbear", spp2 = "wtd", n = 10^5)
  sppX_bear_moose_allot <- allot_conditional_occu_g(gs_bearmoose_global, area = 1, spp1 = "blackbear", spp2 = "moose", n = 10^5)
  # sppX_bob_md_allot <- allot_conditional_occu_g(gs_bobmd_global, area = 1, spp1 = "bobcat", spp2 = "md", n = 10^5)
  # sppX_bob_wtd_allot <- allot_conditional_occu_g(gs_bobwtd_global, area = 1, spp1 = "bobcat", spp2 = "wtd", n = 10^5)
  # sppX_coy_md_allot <- allot_conditional_occu_g(gs_coymd_global, area = 1, spp1 = "coyote", spp2 = "md", n = 10^5)
  # sppX_coy_wtd_allot <- allot_conditional_occu_g(gs_coywtd_global, area = 1, spp1 = "coyote", spp2 = "wtd", n = 10^5)
  
    
  ####  Visualize Allotment Effect on Conditional Occupancy  ####
  sppX_wolf_moose_allot$Species <- factor(sppX_wolf_moose_allot$Species, levels = c("wolf","moose"))
  ggplot(sppX_wolf_moose_allot, aes(x = Public_Grazing, y = Predicted, group = Species, colour = Species)) +
    geom_point(size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_color_manual(values = c("wolf" = "#D55E00", "moose" = "#0072B2"), labels = c("Wolf", "Moose")) +
    ylim(0, 1) +
    xlab("Public grazing") + 
    ylab("Conditional occupancy probability") +
    ggtitle("Effect public grazing allotments on \nwolf and moose occurrence")
  sppX_coug_elk_allot$Species <- factor(sppX_coug_elk_allot$Species, levels = c("cougar","elk"))
  ggplot(sppX_coug_elk_allot, aes(x = Public_Grazing, y = Predicted, group = Interaction, colour = Species)) + 
    geom_point(aes(shape = `Interacting species`), size = 2, position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = "dodge") +
    scale_colour_bright(labels = c("cougar" = "Cougar", "elk" = "Elk"), name = "Species") +
    ylim(0, 1) +
    theme(legend.position="bottom") +
    xlab("Public grazing") + 
    ylab("Conditional occupancy probability") +
    ggtitle("Effect public grazing allotments on cougar and elk occurrence")
  
  ####  ONLY SIGNIFICANT RELATIONSHIPS FOR PAPER  ####
  #'  Conditional occupancy with the effect of cattle activity
  sppX_coug_wtd <- sppX_coug_wtd_cattle[sppX_coug_wtd_cattle$Species == "cougar",]
  sppX_wtd_coug <- sppX_coug_wtd_cattle[sppX_coug_wtd_cattle$Species == "wtd",]
  sppX_coy_wtd <- sppX_coy_wtd_cattle[sppX_coy_wtd_cattle$Species == "coyote",]
  sppX_wtd_coy <- sppX_coy_wtd_cattle[sppX_coy_wtd_cattle$Species == "wtd",]
  
  #'  Merge predator and prey responses together for facet wrap plots
  pred_data <- rbind(sppX_coug_wtd, sppX_coy_wtd) %>%
    mutate(Predator = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction))
  newlabs <- c("cougar" = "Cougar", "coyote" = "Coyote") 
  coug_coy_wtd_graze_facet <- ggplot(pred_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    #scale_color_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed deer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    #scale_fill_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed deer") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of cattle activity on predator co-occurrence with white-tailed deer")
  coug_coy_wtd_graze_facet
  
  prey_data <- rbind(sppX_wtd_coug, sppX_wtd_coy) %>%
    mutate(Predator = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           spp1spp2 = paste(Species, Predator))
  prey_data$Interaction <- factor(prey_data$Interaction, levels = c("cougar absent", "cougar present", "coyote absent", "coyote present"))
  newlabs <- c("wtd cougar" = "White-tailed Deer", "wtd coyote" = "White-tailed Deer")
  wtd_coug_coy_graze_facet <- ggplot(prey_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) + 
    geom_line(size = 1) +
    scale_colour_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present",
                                   "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    # scale_color_manual(values=c("absent" = "#CC6666", "present" = "#9999CC"), labels = c("Absent", "Present", ""), name = "Predator") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present",
                                 "coyote absent" = "Coyote absent", "coyote present" = "Coyote present"), name = "Species Interaction") +
    # scale_fill_manual(values=c("absent" = "#CC6666", "present" = "#9999CC"), labels = c("Absent", "Present", ""), name = "Predator") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~spp1spp2, scales = "free_y", labeller = as_labeller(newlabs)) + 
    theme(legend.position="bottom") +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of cattle activity on white-tailed deer co-occurrence with predators")
  wtd_coug_coy_graze_facet
  
  ggsave("./Outputs/Figures/OccX_coug_coy_wtd_graze.tiff", coug_coy_wtd_graze_facet, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  ggsave("./Outputs/Figures/OccX_wtd_coug_coy_graze.tiff", wtd_coug_coy_graze_facet, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
    
  #'  Conditional occupancy in response to grazing allotments
  sppX_allot <- rbind(sppX_coug_elk_allot, sppX_coug_moose_allot, sppX_wolf_moose_allot, sppX_bear_elk_allot, sppX_bear_moose_allot)
  sppX_allot$Species <- factor(sppX_allot$Species, levels = c("blackbear", "cougar", "wolf", "elk", "moose"))
  newlabs <- c("cougar-elk" = "Cougar-Elk", "cougar-moose" = "Cougar-Moose", "wolf-moose" = "Wolf-Moose",
               "blackbear-elk" = "Black bear-Elk", "blackbear-moose" = "Black bear-Moose")
  
  pred_prey_allot_facet <- ggplot(sppX_allot, aes(x = Public_Grazing, y = Predicted, group = Interaction, colour = Species)) + 
    geom_point(aes(shape = `Interacting species`), size = 2.5, position = position_dodge(width = 0.4)) +
    geom_errorbar(width = 0.2, aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
    scale_colour_bright() +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~SpeciesPairing, scales = "free_y", labeller = as_labeller(newlabs)) + 
    theme(legend.position = c(1, -.05),
          legend.justification = c(1.2, -0.05)) + 
    xlab("Public grazing") + 
    ylab("Conditional occupancy probability") +
    ggtitle("Effect public grazing allotments on predator-prey co-occurrence")
  pred_prey_allot_facet
  ggsave("./Outputs/Figures/OccX_bear_coug_wolf_elk_moose_allot.tiff", pred_prey_allot_facet, 
         units = "in", width = 7, height = 6, dpi = 600, device = 'tiff', ) #, compression = 'lzw'
  
  