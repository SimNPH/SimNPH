library(SimNPH)

save_folder <- paste0("data/results/merged_", format(Sys.Date(), "%Y-%m-%d"), "/")
dir.create(save_folder)

source("scripts/present_data_common.R")

parameter_cols <- c("recruitment", "n_pat", "hazard_ctrl", "censoring_prop",
                    "delay", "effect_size_ph", "crossing", "hr_before",
                    "prog_prop_trt", "prog_prop_ctrl", "hr_before_after")

# merge delayed -----------------------------------------------------------

old_delayed <- readRDS("data/results/results_martin_2023-04/data/simulation_delayed_effect_ims-node2_2023-04-05_151719/results.Rds")
new_delayed <- readRDS("data/results/Additional_2023-04-24/additional_results_delayed_ims-node2_2023-04-22_144212.Rds")

merged_delayed <- merge_additional_results(
  old_delayed,
  new_delayed,
  design_names = parameter_cols,
  descriptive_regex = "^descriptive"
  )

saveRDS(merged_delayed, file = paste0(save_folder, "delayed.Rds"))

# merge crossing ----------------------------------------------------------

old_crossing <- readRDS("data/results/results_martin_2023-04/data/simulation_crossing_hazards_ims-node3_2023-04-05_151749/results.Rds")
new_crossing <- readRDS("data/results/Additional_2023-04-24/additional_results_crossing_ims-node3_2023-04-23_104958.Rds")

merged_crossing <- merge_additional_results(
  old_crossing,
  new_crossing,
  design_names = parameter_cols,
  descriptive_regex = "^descriptive"
  )

saveRDS(merged_crossing, file = paste0(save_folder, "crossing.Rds"))

# merge progression -------------------------------------------------------

old_progression <- readRDS("data/results/results_martin_2023-04/data/simulation_progression_ims-node1_2023-04-05_151339/results.Rds")
new_progression <- readRDS("data/results/Additional_2023-04-24/additional_results_progression_ims-node1_2023-04-24_022352.Rds")

merged_progression <- merge_additional_results(
  old_progression,
  new_progression,
  design_names = parameter_cols,
  descriptive_regex = "^descriptive"
  )

saveRDS(merged_progression, file = paste0(save_folder, "progression.Rds"))
