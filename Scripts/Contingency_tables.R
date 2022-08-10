  #'  ==================================
  #'  Contingency tables
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  August 2022
  #'  ==================================
  #'  Generate contingency tables to review sample size of detections where 
  #'  species are detected together 
  #'  ==================================
      
  #'  Load libraries
  library(tidyverse)
  
  #'  Source data
  source("./Scripts/Detection_Histories_for_unmarked.R")
  source("./Scripts/Cattle_Hunter_Activity.R")
  source("./Scripts/Data_Formatting_2SppX_OccMods.R")
  
  #'  row sum across each detection history, pair sums for each site, build table
  contingency_table <- function(DH_pred, DH_prey, pred, prey, season) {
    #'  Reduce detection histories to single 1 or 0 at each site (naive occupancy)
    pred_det <- apply(DH_pred, 1, max, na.rm = T)
    prey_det <- apply(DH_prey, 1, max, na.rm = T)
    
    #'  Join in a single data frame and reformat so 0/1 represent categories of
    #'  predator/prey presence and absence
    df <- as.data.frame(cbind(pred_det, prey_det)) %>%
      mutate(pred_det = ifelse(pred_det == 1, "predator", "no predator"),
             prey_det = ifelse(prey_det == 1, "prey", "no prey"))
    
    #'  Create 2-way contingency table
    tbl <- table(df)
    tbl <- as.data.frame(tbl)
    
    #'  Add species specific data
    tbl$Predator <- pred
    tbl$Prey <- prey
    tbl$Season <- season
    
    #'  Reorganize
    tbl <- tbl %>%
      relocate(Predator, .before = "pred_det") %>%
      relocate(Prey, .after = "Predator") %>%
      relocate(Season, .after = "Prey")
    colnames(tbl) <- c("Predator", "Prey", "Season", "Predator detection", "Prey detection", "Frequency")
    
    return(tbl)
  }
  ####  Grazing Season Tables  ####
  coug_md_graze_tbl <- contingency_table(DH_coug_graze1820, DH_md_graze1820, pred = "Cougar", prey = "Mule deer", season = "Grazing")
  coug_elk_graze_tbl <- contingency_table(DH_coug_graze1820, DH_elk_graze1820, pred = "Cougar", prey = "Elk", season = "Grazing")
  coug_wtd_graze_tbl <- contingency_table(DH_coug_graze1820, DH_wtd_graze1820, pred = "Cougar", prey = "White-tailed deer", season = "Grazing")
  coug_moose_graze_tbl <- contingency_table(DH_coug_graze1820, DH_moose_graze1820, pred = "Cougar", prey = "Moose", season = "Grazing")
  coug_ungulate_graze_tbl <- rbind(coug_md_graze_tbl, coug_elk_graze_tbl, coug_wtd_graze_tbl, coug_moose_graze_tbl)
  
  wolf_md_graze_tbl <- contingency_table(DH_wolf_graze1820, DH_md_graze1820, pred = "Wolf", prey = "Mule deer", season = "Grazing")
  wolf_elk_graze_tbl <- contingency_table(DH_wolf_graze1820, DH_elk_graze1820, pred = "Wolf", prey = "Elk", season = "Grazing")
  wolf_wtd_graze_tbl <- contingency_table(DH_wolf_graze1820, DH_wtd_graze1820, pred = "Wolf", prey = "White-tailed deer", season = "Grazing")
  wolf_moose_graze_tbl <- contingency_table(DH_wolf_graze1820, DH_moose_graze1820, pred = "Wolf", prey = "Moose", season = "Grazing")
  wolf_ungulate_graze_tbl <- rbind(wolf_md_graze_tbl, wolf_elk_graze_tbl, wolf_wtd_graze_tbl, wolf_moose_graze_tbl)
  
  bear_md_graze_tbl <- contingency_table(DH_bear_graze1820, DH_md_graze1820, pred = "Black bear", prey = "Mule deer", season = "Grazing")
  bear_elk_graze_tbl <- contingency_table(DH_bear_graze1820, DH_elk_graze1820, pred = "Black bear", prey = "Elk", season = "Grazing")
  bear_wtd_graze_tbl <- contingency_table(DH_bear_graze1820, DH_wtd_graze1820, pred = "Black bear", prey = "White-tailed deer", season = "Grazing")
  bear_moose_graze_tbl <- contingency_table(DH_bear_graze1820, DH_moose_graze1820, pred = "Black bear", prey = "Moose", season = "Grazing")
  bear_ungulate_graze_tbl <- rbind(bear_md_graze_tbl, bear_elk_graze_tbl, bear_wtd_graze_tbl, bear_moose_graze_tbl)
  
  bob_md_graze_tbl <- contingency_table(DH_bob_graze1820, DH_md_graze1820, pred = "Bobcat", prey = "Mule deer", season = "Grazing")
  bob_wtd_graze_tbl <- contingency_table(DH_bob_graze1820, DH_wtd_graze1820, pred = "Bobcat", prey = "White-tailed deer", season = "Grazing")
  bob_ungulate_graze_tbl <- rbind(bob_md_graze_tbl, bob_wtd_graze_tbl)
  
  coy_md_graze_tbl <- contingency_table(DH_coy_graze1820, DH_md_graze1820, pred = "Coyote", prey = "Mule deer", season = "Grazing")
  coy_wtd_graze_tbl <- contingency_table(DH_coy_graze1820, DH_wtd_graze1820, pred = "Coyote", prey = "White-tailed deer", season = "Grazing")
  coy_ungulate_graze_tbl <- rbind(coy_md_graze_tbl, coy_wtd_graze_tbl)
  
  grazing_contingency_tbl <- rbind(coug_ungulate_graze_tbl, wolf_ungulate_graze_tbl, bear_ungulate_graze_tbl,
                                   bob_ungulate_graze_tbl, coy_ungulate_graze_tbl)
  
  ####  Hunting Season Tables  ####
  coug_md_hunt_tbl <- contingency_table(DH_coug_hunt1820, DH_md_hunt1820, pred = "Cougar", prey = "Mule deer", season = "Hunting")
  coug_elk_hunt_tbl <- contingency_table(DH_coug_hunt1820, DH_elk_hunt1820, pred = "Cougar", prey = "Elk", season = "Hunting")
  coug_wtd_hunt_tbl <- contingency_table(DH_coug_hunt1820, DH_wtd_hunt1820, pred = "Cougar", prey = "White-tailed deer", season = "Hunting")
  coug_moose_hunt_tbl <- contingency_table(DH_coug_hunt1820, DH_moose_hunt1820, pred = "Cougar", prey = "Moose", season = "Hunting")
  coug_ungulate_hunt_tbl <- rbind(coug_md_hunt_tbl, coug_elk_hunt_tbl, coug_wtd_hunt_tbl, coug_moose_hunt_tbl)
  
  wolf_md_hunt_tbl <- contingency_table(DH_wolf_hunt1820, DH_md_hunt1820, pred = "Wolf", prey = "Mule deer", season = "Hunting")
  wolf_elk_hunt_tbl <- contingency_table(DH_wolf_hunt1820, DH_elk_hunt1820, pred = "Wolf", prey = "Elk", season = "Hunting")
  wolf_wtd_hunt_tbl <- contingency_table(DH_wolf_hunt1820, DH_wtd_hunt1820, pred = "Wolf", prey = "White-tailed deer", season = "Hunting")
  wolf_moose_hunt_tbl <- contingency_table(DH_wolf_hunt1820, DH_moose_hunt1820, pred = "Wolf", prey = "Moose", season = "Hunting")
  wolf_ungulate_hunt_tbl <- rbind(wolf_md_hunt_tbl, wolf_elk_hunt_tbl, wolf_wtd_hunt_tbl, wolf_moose_hunt_tbl)
  
  bear_md_hunt_tbl <- contingency_table(DH_bear_hunt1820, DH_md_hunt1820, pred = "Black bear", prey = "Mule deer", season = "Hunting")
  bear_elk_hunt_tbl <- contingency_table(DH_bear_hunt1820, DH_elk_hunt1820, pred = "Black bear", prey = "Elk", season = "Hunting")
  bear_wtd_hunt_tbl <- contingency_table(DH_bear_hunt1820, DH_wtd_hunt1820, pred = "Black bear", prey = "White-tailed deer", season = "Hunting")
  bear_moose_hunt_tbl <- contingency_table(DH_bear_hunt1820, DH_moose_hunt1820, pred = "Black bear", prey = "Moose", season = "Hunting")
  bear_ungulate_hunt_tbl <- rbind(bear_md_hunt_tbl, bear_elk_hunt_tbl, bear_wtd_hunt_tbl, bear_moose_hunt_tbl)
  
  bob_md_hunt_tbl <- contingency_table(DH_bob_hunt1820, DH_md_hunt1820, pred = "Bobcat", prey = "Mule deer", season = "Hunting")
  bob_wtd_hunt_tbl <- contingency_table(DH_bob_hunt1820, DH_wtd_hunt1820, pred = "Bobcat", prey = "White-tailed deer", season = "Hunting")
  bob_ungulate_hunt_tbl <- rbind(bob_md_hunt_tbl, bob_wtd_hunt_tbl)
  
  coy_md_hunt_tbl <- contingency_table(DH_coy_hunt1820, DH_md_hunt1820, pred = "Coyote", prey = "Mule deer", season = "Hunting")
  coy_wtd_hunt_tbl <- contingency_table(DH_coy_hunt1820, DH_wtd_hunt1820, pred = "Coyote", prey = "White-tailed deer", season = "Hunting")
  coy_ungulate_hunt_tbl <- rbind(coy_md_hunt_tbl, coy_wtd_hunt_tbl)
  
  hunting_contingency_tbl <- rbind(coug_ungulate_hunt_tbl, wolf_ungulate_hunt_tbl, bear_ungulate_hunt_tbl,
                                   bob_ungulate_hunt_tbl, coy_ungulate_hunt_tbl)
  
  detection_contingency_tbl <- rbind(grazing_contingency_tbl, hunting_contingency_tbl) %>%
    arrange(Season, Predator, Prey)
  write.csv(detection_contingency_tbl, file = paste0("./Outputs/Tables/Contingency_tables_", Sys.Date(), ".csv"))  
  
  