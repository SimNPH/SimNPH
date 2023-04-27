## Load Data
data_folder <- "./"
results_progression <- readRDS(paste0(data_folder,"results_progression_2023-04-26.Rds")) |>
  results_pivot_longer() |>
  mutate(rejection_0.025=rejection_0.05) |> ## Fake it until you make it!
  add_meta_data(metadata,col_labels) |>
  mutate(across(starts_with(c("median_survival",
                              "recruitment",
                              "rmst",
                              "crossing",
                              "descriptive.study_time",
                              "study_time")),~.x/30.4375),
         average_study_time = coalesce(study_time,descriptive.study_time)) |> ## Transform days to months
  mutate(bias = ifelse(method %in% c("median_surv",
                                     "diff_med_weibull",
                                     "rmst_diff_6m",
                                     "rmst_diff_12m"),
                       bias_median_survival(bias,median_survival_trt,median_survival_ctrl),
                       bias)) |> ## Compute relative bias
  arrange(method) |> ## we don't need this anymore
  mutate(pfs_median_survival_ctrl = round(r2m(hazard_ctrl))) |>
  filter(method != "peto_peto") # same as fh_1_0

## Test plots
results_crossing |>
  plot_sim_results(methods = c('logrank','cox'),
                   parameter_x = c('hr_before','crossing'),
                   parameter_y = 'rejection')
