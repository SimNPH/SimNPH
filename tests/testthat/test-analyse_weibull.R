test_that("analyse weibull works", {

  capture_output({
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(3) |>
      tail(1)
  })

  dat <- generate_delayed_effect(condition)

  res <- analyse_weibull()(NA, dat)

  expect_named(res, c(
    "p",
    "med_trt_est", "med_trt_lower", "med_trt_upper",
    "med_ctrl_est", "med_ctrl_lower", "med_ctrl_upper",
    "diff_med_est", "diff_med_lower", "diff_med_upper",
    "CI_level",
    "N_pat", "N_evt"
    ))

})
