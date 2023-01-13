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
    final_events <- ceiling(n_pat * 0.75)
    interim_events <- ceiling(final_events * 0.5)
  })


# Assumptions -------------------------------------------------------------

assumptions <- expand.grid(
  hazard_ctrl = nph::m2r(c(36, 12, 6)),
  censoring_prop = c(0, 0.1, 0.3),
  effect_size_ph = c(0, 0.5, 0.8, 0.9),
  prevalence = seq(0.1, 0.5, by=0.2),
  hr_subgroup_relative = seq(0.7, 0.9, by=0.1)
)

# Merging Options and Assumptions -----------------------------------------

design <- merge(
  options,
  assumptions,
  by=NULL
)


# Callibrating Implicitly defined Parameters ------------------------------

# hazard rate after onset of treatment effect
# calculated such that the median survival is the same as under proportional
# hazards with the given effect size. Fails if onset of treatment effect is
# after the median survival time
design <- design |>
  hazard_subgroup_from_PH_effect_size()

# rate of random censoring
design <- design |>
  cen_rate_from_cen_prop_subgroup()

# Calculating True Summary Statistics -------------------------------------

design <- design |>
  true_summary_statistics_subgroup(milestones = m2d(c(ms_surv_8=8, ms_surv_12=12)))


# Saving ------------------------------------------------------------------

filename <- paste0("data/parameters/subgroup_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.table(design, file=filename, quote=FALSE, sep=", ", dec=".", row.names = FALSE, col.names = TRUE)

