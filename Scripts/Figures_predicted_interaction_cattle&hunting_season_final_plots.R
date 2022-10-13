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
  library(patchwork)
  
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
  
  
  #'  ----------------------------------------------------------
  ####  Covariate effects on cattle - species interaction term  ####
  #'  ----------------------------------------------------------
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
  
  #'  Save these since they take so long to generate!
  sppX_cattle_list <- list(sppX_coug_wtd_cattle, sppX_wolf_wtd_cattle, sppX_bob_md_cattle, 
                           sppX_bob_wtd_cattle, sppX_coy_md_cattle, sppX_coy_wtd_cattle)
  names(sppX_cattle_list) <- c("sppX_coug_wtd_cattle", "sppX_wolf_wtd_cattle", "sppX_bob_md_cattle", 
                               "sppX_bob_wtd_cattle", "sppX_coy_md_cattle", "sppX_coy_wtd_cattle")
  save(sppX_cattle_list, file = "./Outputs/sppX_cattle_for_visualizing.RData")
  
  
  
  #'  -----------------------------------------------------------
  ####  Covariate effects on hunter - species interaction term  ####
  #'  -----------------------------------------------------------
  #'  Function to predict species interactions in response to covariate of interest
  spp_interactions_h <- function(mod, elev, act, forest, pub, area, spp1, spp2, cov) {
    #'  Create data frame using the scaled covariate of interest while holding
    #'  all others at their mean (0 when scaled) or desired category (0 or 1)
    cov_df <- data.frame(Elev = elev, HuntingActivity = act, PercForest = forest,
                         Public = pub, Study_Area = area)
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
  #'  For elk, moose, & wtd study area = 0 (NE) vs md study area = 1 (OK) because 
  #'  higher species-specific densities in the respective study areas
  #'  For all species, pub = 1
  sppX_coug_md_hunter <- spp_interactions_h(hs_cougmd_global, elev = 0, act = scaled_hunt[,2],
                                            forest = 0, pub = 1, area = 1, spp1 = "cougar",
                                            spp2 = "muledeer", cov = scaled_hunt[,1])
  # sppX_coug_elk_hunter <- spp_interactions_h(hs_cougelk_global, elev = 0, act = scaled_hunt[,2], 
  #                                            forest = 0, pub = 1, area = 0, spp1 = "cougar", 
  #                                            spp2 = "elk", cov = scaled_hunt[,1])
  sppX_coug_wtd_hunter <- spp_interactions_h(hs_cougwtd_global, elev = 0, act = scaled_hunt[,2],
                                             forest = 0, pub = 1, area = 0, spp1 = "cougar",
                                             spp2 = "wtd", cov = scaled_hunt[,1])
  # sppX_coug_moose_hunter <- spp_interactions_h(hs_cougmoose_global, elev = 0, act = scaled_hunt[,2], 
  #                                            forest = 0, pub = 1, area = 0, spp1 = "cougar", 
  #                                            spp2 = "moose", cov = scaled_hunt[,1])
  # sppX_wolf_md_hunter <- spp_interactions_h(hs_wolfmd_global, elev = 0, act = scaled_hunt[,2],
  #                                          forest = 0, pub = 1, area = 1, spp1 = "wolf",
  #                                          spp2 = "muledeer", cov = scaled_hunt[,1])
  # sppX_wolf_elk_hunter <- spp_interactions_h(hs_wolfelk_global, elev = 0, act = scaled_hunt[,2],
  #                                            forest = 0, pub = 1, area = 0, spp1 = "wolf",
  #                                            spp2 = "elk", cov = scaled_hunt[,1])
  # sppX_wolf_wtd_hunter <- spp_interactions_h(hs_wolfwtd_global, elev = 0, act = scaled_hunt[,2],
  #                                           forest = 0, pub = 1, area = 0, spp1 = "wolf",
  #                                           spp2 = "wtd", cov = scaled_hunt[,1])
  sppX_wolf_moose_hunter <- spp_interactions_h(hs_wolfmoose_global, elev = 0, act = scaled_hunt[,2], 
                                               forest = 0, pub = 1, area = 0, spp1 = "wolf", 
                                               spp2 = "moose", cov = scaled_hunt[,1])
  # sppX_bear_md_hunter <- spp_interactions_h(hs_bearmd_global, elev = 0, act = scaled_hunt[,2], 
  #                                           forest = 0, pub = 1, area = 1, spp1 = "blackbear", 
  #                                           spp2 = "muledeer", cov = scaled_hunt[,1])
  # sppX_bear_elk_hunter <- spp_interactions_h(hs_bearelk_global, elev = 0, act = scaled_hunt[,2], 
  #                                            forest = 0, pub = 1, area = 0, spp1 = "blackbear", 
  #                                            spp2 = "elk", cov = scaled_hunt[,1])
  # sppX_bear_wtd_hunter <- spp_interactions_h(hs_bearwtd_global, elev = 0, act = scaled_hunt[,2], 
  #                                            forest = 0, pub = 1, area = 0, spp1 = "blackbear", 
  #                                            spp2 = "wtd", cov = scaled_hunt[,1])
  # sppX_bear_moose_hunter <- spp_interactions_h(hs_bearmoose_global, elev = 0, act = scaled_hunt[,2], 
  #                                              forest = 0, pub = 1, area = 0, spp1 = "blackbear", 
  #                                              spp2 = "moose", cov = scaled_hunt[,1])
  # sppX_bob_md_hunter <- spp_interactions_h(hs_bobmd_global, elev = 0, act = scaled_hunt[,2], 
  #                                          forest = 0, pub = 1, area = 1, spp1 = "bobcat", 
  #                                          spp2 = "mule_deer", cov = scaled_hunt[,1])
  sppX_bob_wtd_hunter <- spp_interactions_h(hs_bobwtd_global, elev = 0, act = scaled_hunt[,2],
                                            forest = 0, pub = 1, area = 0, spp1 = "bobcat",
                                            spp2 = "wtd", cov = scaled_hunt[,1])
  # sppX_coy_md_hunter <- spp_interactions_h(hs_coymd_global, elev = 0, act = scaled_hunt[,2], 
  #                                          forest = 0, pub = 1, area = 1, spp1 = "coyote", 
  #                                          spp2 = "mule_deer", cov = scaled_hunt[,1])
  sppX_coy_wtd_hunter <- spp_interactions_h(hs_coywtd_global, elev = 0, act = scaled_hunt[,2],
                                            forest = 0, pub = 0, area = 0, spp1 = "coyote",
                                            spp2 = "wtd", cov = scaled_hunt[,1]) #note this is estimated for on private land!
  
  
  #'  Save these since they take so long to generate!
  sppX_hunt_list <- list(sppX_coug_md_hunter, sppX_coug_wtd_hunter, sppX_wolf_moose_hunter, 
                         sppX_bob_wtd_hunter, sppX_coy_wtd_hunter)
  names(sppX_hunt_list) <- c("sppX_coug_md_hunter", "sppX_coug_wtd_hunter", "sppX_wolf_moose_hunter", 
                             "sppX_bob_wtd_hunter", "sppX_coy_wtd_hunter")
  save(sppX_hunt_list, file = "./Outputs/sppX_hunter_for_visualizing.RData")
  
  
  #'  ------------------------------------
  ####  Visualize conditional occupancy  ####
  #'  ------------------------------------
  #'  Load predicted conditional occupancy
  load("./Outputs/sppX_cattle_for_visualizing.RData")
  load("./Outputs/sppX_hunter_for_visualizing.RData")
  
  #'  Cougar-wtd co-occurrence with cattle activity
  coug_wtd_data <- sppX_cattle_list$sppX_coug_wtd_cattle %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "cougar absent", "cougar present")))
  newlabs <- c("cougar" = "Cougar", "wtd" = "White-tailed deer") 
  coug_wtd_graze_facet <- ggplot(coug_wtd_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    #scale_color_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed deer") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("cougar absent" = "Cougar absent", "cougar present" = "Cougar present",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    #scale_fill_manual(values=c("wtd absent" = "#CC6666", "wtd present" = "#9999CC"), labels = c("Absent", "Present"), name = "White-tailed deer") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of cattle activity on cougar - white-tailed deer co-occurrence")
  coug_wtd_graze_facet
  
  ggsave("./Outputs/Figures/OccX_coug_wtd_graze.tiff", coug_wtd_graze_facet, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw') 
  
  #'  Coyote-wtd co-occurrence with cattle activity
  coy_wtd_data <- sppX_cattle_list$sppX_coy_wtd_cattle %>% 
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "coyote absent", "coyote present")))
  newlabs <- c("coyote" = "Coyote", "wtd" = "White-tailed deer") 
  coy_wtd_graze_facet <- ggplot(coy_wtd_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("coyote absent" = "Coyote absent", "coyote present" = "Coyote present", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("coyote absent" = "Coyote absent", "coyote present" = "Coyote present",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Cattle grazing activity (cattle detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of cattle activity on coyote - white-tailed deer co-occurrence")
  coy_wtd_graze_facet
  
  ggsave("./Outputs/Figures/OccX_coy_wtd_graze.tiff", coy_wtd_graze_facet, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw') 
  
  #'  Coyote-wtd co-occurrence with hunter activity
  coy_wtd_data <- sppX_hunt_list$sppX_coy_wtd_hunter %>% #rbind(sppX_coy_wtd, sppX_wtd_coy) %>%
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("wtd absent", "wtd present", "coyote absent", "coyote present")))
  newlabs <- c("coyote" = "Coyote", "wtd" = "White-tailed deer") 
  coy_wtd_hunt_facet <- ggplot(coy_wtd_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("coyote absent" = "Coyote absent", "coyote present" = "Coyote present", 
                                   "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("coyote absent" = "Coyote absent", "coyote present" = "Coyote present",
                                 "wtd absent" = "White-tailed \ndeer absent", "wtd present" = "White-tailed \ndeer present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Hunter activity (hunter detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of hunter activity on coyote - white-tailed deer co-occurrence")
  coy_wtd_hunt_facet
  
  ggsave("./Outputs/Figures/OccX_coy_wtd_hunt.tiff", coy_wtd_hunt_facet, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  #'  Wolf-moose co-occurrence with hunter activity
  wolf_moose_data <- sppX_hunt_list$sppX_wolf_moose_hunter %>% #rbind(sppX_wolf_moose, sppX_moose_wolf) %>%
    mutate(InterXSpp = gsub( " .*$", "", Interaction ),
           InterX = gsub(".* ", "", Interaction),
           Interaction = factor(Interaction, levels = c("moose absent", "moose present", "wolf absent", "wolf present")),
           Species = factor(Species, levels = c("wolf", "moose")))
  newlabs <- c("wolf" = "Wolf", "moose" = "Moose") 
  wolf_moose_hunt_facet <- ggplot(wolf_moose_data, aes(x = Cov, y = Predicted, group = Interaction, colour = Interaction)) +
    geom_line(size = 1) +
    scale_colour_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                   "moose absent" = "Moose absent", "moose present" = "Moose present"), name = "Species Interaction") +
    geom_ribbon(aes(ymin=lower, ymax = upper, fill = Interaction), linetype = 0, alpha =0.5) +
    scale_fill_bright(labels = c("wolf absent" = "Wolf absent", "wolf present" = "Wolf present", 
                                 "moose absent" = "Moose absent", "moose present" = "Moose present"), name = "Species Interaction") +
    ylim(0, 1) +
    theme_bw() +
    facet_wrap(~Species, scales = "free_y", labeller = as_labeller(newlabs)) +
    theme(legend.position="bottom") +
    xlab("Hunter activity (hunter detections/day)") + 
    ylab("Conditional occupancy") +
    ggtitle("Effect of hunter activity on wolf - moose co-occurrence")
  wolf_moose_hunt_facet
  
  ggsave("./Outputs/Figures/OccX_wolf_moose_hunt.tiff", wolf_moose_hunt_facet, 
         units = "in", width = 7, height = 5, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  
  coocc_patwork <- coug_wtd_graze_facet + coy_wtd_graze_facet + wolf_moose_hunt_facet + 
    coy_wtd_hunt_facet + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'a')
  coocc_patwork
  
  ggsave("./Outputs/Figures/OccX_pred-prey_cattl&hunter.tiff", coocc_patwork, 
         units = "in", width = 15, height = 10, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  