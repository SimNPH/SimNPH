library(SimNPH)

# Helper functions --------------------------------------------------------

# months to days
m2d <- \(t) 365.25*t/12

# Options -----------------------------------------------------------------

options <- expand.grid(
  recruitment = m2d(c(18, 30)),
  n_pat = c(300, 500, 1000, 1500)
) |>
  within({
    n_trt <- n_pat / 2
    n_ctrl <- n_pat / 2
    interim_events <- n_pat / 4
    final_events <- n_pat / 2
  })


# Assumptions -------------------------------------------------------------

assumptions <- expand.grid(
  hazard_ctrl = nph::m2r(c(36, 12, 6)),
  prog_prop_trt  = c(0.1, 0.2),
  prog_prop_ctrl = c(0.1, 0.2),
  hr_before_after = c(0.8, 0.5),
  censoring_prop = c(0, 0.1, 0.3)
)

assumptions$hazard_trt <- assumptions$hazard_ctrl
assumptions$hazard_after_prog <- assumptions$hazard_ctrl / assumptions$hr_before_after

# Merging Options and Assumptions -----------------------------------------

design <- merge(
  options,
  assumptions,
  by=NULL
)

# calculate hazards from PH effect size -----------------------------------

# calculate hazards (discuss options, which hazard in the treatment arm should
# be callibrated? hazard had no progression occured, hazard for pfs, ...)

# calculate progression rate from progression proportion  -----------------

design <- design |>
  within({
    followup <- log(6) / hazard_ctrl
  })

design <- SimNPH:::progression_rate_from_progression_prop(design)


# calculate random withdrawal ---------------------------------------------

# design$random_withdrawal <- 0

# calculate real statistics -----------------------------------------------

design <- design |>
  true_summary_statistics_progression()
