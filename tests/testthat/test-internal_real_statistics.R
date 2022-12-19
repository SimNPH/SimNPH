test_that("internal function for real summary statistics outputs correct values", {

  capture_output({
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      tail(1)
  })

  t_max <- max(
    log(10000) / condition$hazard_ctrl,
    log(10000) / condition$hazard_trt
  )

  treatment <- nph::pchaz(
    c(0, condition$delay, t_max),
    c(condition$hazard_ctrl, condition$hazard_trt)
  )

  control <- nph::pchaz(
    c(0, t_max),
    c(condition$hazard_ctrl)
  )

  control_low_rate <- nph::pchaz(
    c(0, t_max),
    c(condition$hazard_ctrl/1000)
  )

  result_1 <- internal_real_statistics_pchaz(treatment, control)
  result_2 <- internal_real_statistics_pchaz(treatment, treatment)
  result_3 <- internal_real_statistics_pchaz(control, control)
  result_4 <- internal_real_statistics_pchaz(control, control_low_rate)

  result_1_b <- fast_real_statistics_pchaz(treatment$Tint[-3], treatment$lambda, control$Tint[-2], control$lambda, milestones=c(msurv_5=5, msurv_10=10, msurv_15=15))
  result_2_b <- fast_real_statistics_pchaz(treatment$Tint[-3], treatment$lambda, treatment$Tint[-3], treatment$lambda, milestones=c(5, 10, 15))
  result_3_b <- fast_real_statistics_pchaz(control$Tint[-2], control$lambda, control$Tint[-2], control$lambda)

  # correct class, names, etc.
  expect_s3_class(result_1, "data.frame")
  expect_named(result_1, c("rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR"))

  expect_s3_class(result_1_b, "data.frame")
  expect_named(result_1_b, c("rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR",
                             "msurv_5_trt", "msurv_10_trt", "msurv_15_trt", "msurv_5_ctrl", "msurv_10_ctrl", "msurv_15_ctrl"))
  expect_named(result_2_b, c("rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR",
                             "milestone_surv_5_trt", "milestone_surv_10_trt", "milestone_surv_15_trt", "milestone_surv_5_ctrl",
                             "milestone_surv_10_ctrl", "milestone_surv_15_ctrl"))

  expect_equal(result_4$median_survival_ctrl, Inf)

  # AHR has to be one, if treatment and control have the same distribution
  expect_equal(result_2$gAHR, 1)
  expect_equal(result_2$AHR, 1)
  expect_equal(result_3$gAHR, 1)
  expect_equal(result_3$AHR, 1)

  expect_equal(result_2_b$gAHR, 1)
  expect_equal(result_2_b$AHR, 1)
  expect_equal(result_3_b$gAHR, 1)
  expect_equal(result_3_b$AHR, 1)

  result_4

})
