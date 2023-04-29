devtools::load_all()

m2d <- \(t) 365.25*t/12

design_original <- readRDS("../final_results/merged/progression.Rds")

design_new <- design_original[, c("recruitment", "n_pat", "interim_events", "final_events", "n_ctrl",
                     "n_trt", "hazard_ctrl", "prog_prop_trt", "prog_prop_ctrl", "hr_before_after", "censoring_prop", "effect_size_ph", "hazard_after_prog", "prog_rate_ctrl", "prog_rate_trt",
                     "hazard_trt", "random_withdrawal"
)]


design_new <- design_new |>
  true_summary_statistics_progression(
    milestones   = m2d(c("6m"=6, "12m"=12)),
    cutoff_stats = m2d(c("6m"=6, "12m"=12, "360m"=360))
  )

plot(design_new[, c("gAHR_schemper_cox_12m", "gAHR_schemper_wcox_12m", "gAHR_rauch_12m")])
plot(design_new[, c("gAHR_schemper_cox_6m", "gAHR_schemper_wcox_6m", "gAHR_rauch_6m")])
plot(design_new[, c("AHR_schemper_cox_12m", "AHR_schemper_wcox_12m", "AHR_rauch_12m")])
plot(design_new[, c("AHR_schemper_cox_6m", "AHR_schemper_wcox_6m", "AHR_rauch_6m")])

saveRDS(design_new, "data/other_ahrs_progression.Rds")
