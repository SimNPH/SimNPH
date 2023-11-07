## code to prepare `combination_tests_delayed` dataset goes here

library(SimNPH)
library(parallel)

cl <- makeCluster(8)
clusterEvalQ(cl, {
  library(SimDesign)
  library(SimNPH)
  library(parallel)
})


condition <- expand.grid(
  delay=m2d(seq(0, 10, by=2)),  # delay of 0, 1, ..., 10 months
  hazard_ctrl=m2r(c(12,18,24)), # median survival control of 12,18,24 months
  hr=c(0.9, 0.8, 0.7),          # hazard ratio after onset of 0.9, 0.8, 0.7
  random_withdrawal=m2r(120),   # median time to random withdrawal 10 years
  n_pat=c(300, 600, 1200),      # 300, 600, 1200 patients
  followup=m2d(24),  # followup 2 years
  recruitment=m2d(c(0,3,6)) # recruitment time instantly, 3 months, 6 months
) |>
  within({
    n_trt  = n_pat/2
    n_ctrl = n_pat/2
    hazard_trt = hazard_ctrl * hr
    final_events = 0.75 * n_pat
  }) |>
  true_summary_statistics_delayed_effect(cutoff_stats = 15)


my_generator <- function(condition, fixed_objects=NULL){
  generate_delayed_effect(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_events(condition$final_events)
}

combination_tests_delayed <- runSimulation(
  design = condition,
  replications = 1500,
  generate = my_generator,
  analyse = list(
    logrank  = analyse_logrank(alternative = "one.sided"),
    maxcombo = analyse_maxcombo(alternative = "one.sided"),
    mwlrt    = analyse_modelstly_weighted(t_star = m2d(24))
  ),
  summarise = create_summarise_function(
    logrank = summarise_test(0.025),
    mwlrt = summarise_test(0.025),
    maxcombo = summarise_test(0.025)
  ),
  cl = cl,
  save = FALSE
)


usethis::use_data(combination_tests_delayed, overwrite = TRUE)
