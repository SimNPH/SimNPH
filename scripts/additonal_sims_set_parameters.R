library(tidyverse)
library(SimNPH)

# months to days
m2d <- \(t) 365.25*t/12

# helper function ---------------------------------------------------------

update_parameters_with_seed <- function(old_result_filename, add_summary_stats = identity){
  # extract old nodename and read old results
  old_nodename <- str_extract(old_result_filename, "[^_]*(?=_\\d{4}-\\d{2}-\\d{2})")
  old_results <- readRDS(old_result_filename)

  # select only design columns, discard results and attributes
  design_cols <- attr(old_results, "design_names")$design

  old_results <- old_results |>
    select(all_of(design_cols), old_seed = SEED)

  attributes(old_results) <- attributes(old_results)[c("names", "row.names", "class")]

  # add summary statistics
  new_parameters <- add_summary_stats(old_results)
  # and remove the old columns, if columns were duplicated
  new_parameters <- new_parameters[-which(rev(duplicated(rev(names(new_parameters)))))]

  # add old_node_name as a column (simplest way when saving as csv)
  new_parameters$node_name <- old_nodename

  new_parameters
}


# do the actual work ------------------------------------------------------

new_parameters_delayed <- update_parameters_with_seed(
  "data/results/results_martin_2023-04/data/simulation_delayed_effect_ims-node2_2023-04-05_151719/results.Rds",
  partial(true_summary_statistics_delayed_effect, cutoff_stats=m2d(360))
)

new_parameters_crossing <- update_parameters_with_seed(
  "data/results/results_martin_2023-04/data/simulation_crossing_hazards_ims-node3_2023-04-05_151749/results.Rds",
  partial(true_summary_statistics_crossing_hazards, cutoff_stats=m2d(360))
)

new_parameters_subgroup <- update_parameters_with_seed(
  "data/results/results_martin_2023-04/data/simulation_subgroup_ims-node5_2023-04-05_151813/results.Rds",
  partial(true_summary_statistics_subgroup, cutoff_stats=m2d(360))
)

new_parameters_progression <- update_parameters_with_seed(
  "data/results/results_martin_2023-04/data/simulation_progression_ims-node1_2023-04-05_151339/results.Rds",
  partial(true_summary_statistics_progression, cutoff_stats=m2d(360))
)

# save --------------------------------------------------------------------
filename_delayed <- paste0("data/parameters/additional_delayed_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.table(new_parameters_delayed, file=filename_delayed, quote=FALSE, sep=", ", dec=".", row.names = FALSE, col.names = TRUE)

filename_crossing <- paste0("data/parameters/additional_crossing_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.table(new_parameters_crossing, file=filename_crossing, quote=FALSE, sep=", ", dec=".", row.names = FALSE, col.names = TRUE)

filename_subgroup <- paste0("data/parameters/additional_subgroup_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.table(new_parameters_subgroup, file=filename_subgroup, quote=FALSE, sep=", ", dec=".", row.names = FALSE, col.names = TRUE)

filename_os <- paste0("data/parameters/additional_progression_os_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.table(new_parameters_progression, file=filename_os, quote=FALSE, sep=", ", dec=".", row.names = FALSE, col.names = TRUE)

