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
  censoring_prop = c(0, 0.1, 0.3),
  delay = m2d(1:9),
  effect_size_ph = c(0, 0.5, 0.8, 0.9)
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
  within({
    followup <- log(6) / hazard_ctrl
  }) |>
  hr_after_onset_from_PH_effect_size()

# rate of random censoring
design <- design |>
  cen_rate_from_cen_prop_delayed_effect()


# Excluding Scenarios which did not give Reasonable Parameter Valu --------

# Excluding Scenarios for which a hazard in the treatment arm could not be calculated.
# This happens when the median of the survival functions is before onset of treatment effect.
design <- design |>
  subset(!is.na(hazard_trt))

# Saving ------------------------------------------------------------------

filename <- paste("data/parameters/delayed_effect_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.table(design, file=filename, quote=FALSE, sep=", ", dec=".", row.names = FALSE, col.names = TRUE)

