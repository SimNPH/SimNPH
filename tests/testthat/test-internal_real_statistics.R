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

  control_different_t <- nph::pchaz(
    c(0, 100),
    c(condition$hazard_ctrl)
  )

  control_different_t_2 <- control
  control_different_t_2$t[50] <- -1

  treatment_10 <- nph::pchaz(
    c(0, condition$delay, t_max)*10,
    c(condition$hazard_ctrl, condition$hazard_trt)/10
  )

  control_10 <- nph::pchaz(
    c(0, t_max)*10,
    c(condition$hazard_ctrl)/10
  )

  result_1 <- internal_real_statistics_pchaz_discrete(treatment, control, milestones=c("5m"=5, "10m"=10, "15m"=15), cutoff=c("cutoff_5"=5, "cutoff_10"=10))
  result_2 <- internal_real_statistics_pchaz_discrete(treatment, treatment, milestones=c(5, 10, 15), cutoff=c(5, 10))
  result_3 <- internal_real_statistics_pchaz_discrete(control, control, cutoff=15)
  result_4 <- internal_real_statistics_pchaz_discrete(control, control_low_rate)

  result_1_b <- fast_real_statistics_pchaz(treatment$Tint[-3], treatment$lambda, control$Tint[-2], control$lambda, milestones=c("5m"=5, "10m"=10, "15m"=15), cutoff=c("cutoff_5"=5, "cutoff_10"=10))
  result_2_b <- fast_real_statistics_pchaz(treatment$Tint[-3], treatment$lambda, treatment$Tint[-3], treatment$lambda, milestones=c(5, 10, 15), cutoff=c(5, 10))
  result_3_b <- fast_real_statistics_pchaz(control$Tint[-2], control$lambda, control$Tint[-2], control$lambda, cutoff=5)

  result_1_c <- internal_real_statistics_pchaz_discrete(treatment_10, control_10, milestones=c("5m"=5, "10m"=10, "15m"=15)*10, cutoff=c("5"=5, "10"=10)*10)
  result_1_c[, c("median_survival_trt","median_survival_ctrl","rmst_trt_5","rmst_ctrl_5","rmst_trt_10","rmst_ctrl_10")] <-
    result_1_c[, c("median_survival_trt","median_survival_ctrl","rmst_trt_5","rmst_ctrl_5","rmst_trt_10","rmst_ctrl_10")] / 10

  # correct class, names, etc.
  expect_s3_class(result_1, "data.frame")
  expect_named(result_3, c("rmst_trt_15", "median_survival_trt", "rmst_ctrl_15", "median_survival_ctrl", "gAHR_15", "AHR_15"), ignore.order = TRUE)

  expect_s3_class(result_1_b, "data.frame")
  expect_named(result_1_b, c("median_survival_trt", "median_survival_ctrl", "rmst_trt_cutoff_5",
                             "rmst_ctrl_cutoff_5", "gAHR_cutoff_5", "AHR_cutoff_5", "AHRoc_cutoff_5",
                             "AHRoc_robust_cutoff_5", "rmst_trt_cutoff_10", "rmst_ctrl_cutoff_10",
                             "gAHR_cutoff_10", "AHR_cutoff_10", "AHRoc_cutoff_10", "AHRoc_robust_cutoff_10",
                             "milestone_survival_trt_5m", "milestone_survival_ctrl_5m", "milestone_survival_trt_10m",
                             "milestone_survival_ctrl_10m", "milestone_survival_trt_15m",
                             "milestone_survival_ctrl_15m"))
  expect_named(result_2_b, c("median_survival_trt", "median_survival_ctrl", "rmst_trt_5",
                             "rmst_ctrl_5", "gAHR_5", "AHR_5", "AHRoc_5", "AHRoc_robust_5",
                             "rmst_trt_10", "rmst_ctrl_10", "gAHR_10", "AHR_10", "AHRoc_10",
                             "AHRoc_robust_10", "milestone_survival_trt_5", "milestone_survival_ctrl_5",
                             "milestone_survival_trt_10", "milestone_survival_ctrl_10", "milestone_survival_trt_15",
                             "milestone_survival_ctrl_15"))

  expect_equal(result_4$median_survival_ctrl, Inf)

  # AHR has to be one, if treatment and control have the same distribution
  expect_equal(result_2$gAHR_5, 1)
  expect_equal(result_2$AHR_5, 1)
  expect_equal(result_3$gAHR_15, 1)
  expect_equal(result_3$AHR_15, 1)

  expect_equal(result_2_b$gAHR_5, 1)
  expect_equal(result_2_b$AHR_5, 1)
  expect_equal(result_3_b$gAHR_5, 1)
  expect_equal(result_3_b$AHR_5, 1)

  expect_error(internal_real_statistics_pchaz_discrete(treatment, control_different_t))
  expect_error(internal_real_statistics_pchaz_discrete(treatment, control_different_t_2))

  # expect(all(abs((result_1_c / result_1_b) - 1) < 0.15), "relative error between discrete approximation with 1/10 interval and exact calculation is less than 15%")
})

test_that("average hazard ratios work", {

  res <- fast_real_statistics_pchaz(
    Tint_trt=0,  lambda_trt=0.5,
    Tint_ctrl=0, lambda_ctrl=1,
    N_trt=1, N_ctrl=1, cutoff=c(1, 5, 10, 100), milestones=c(1, 5, 10, 100)
  )

  max_abs_err_gahr <- max(abs(res[, c("gAHR_1", "gAHR_5", "gAHR_10", "gAHR_100")] - 0.5))
  max_abs_err_ahr  <- max(abs(res[, c("AHR_1", "AHR_5", "AHR_10", "AHR_100")] - 0.5))

  # tolerance for absolute error from integrate funciton
  expect_lt(max_abs_err_gahr, .Machine$double.eps^0.25)
  expect_lt(max_abs_err_ahr, .Machine$double.eps^0.25)

  res <- fast_real_statistics_pchaz(
    Tint_trt=0,  lambda_trt=1,
    Tint_ctrl=0, lambda_ctrl=1,
    N_trt=1, N_ctrl=1, cutoff=c(1, 5, 10, 100), milestones=c(1, 5, 10, 100)
  )

  max_abs_err_gahr <- max(abs(res[, c("gAHR_1", "gAHR_5", "gAHR_10", "gAHR_100")] - 1))
  max_abs_err_ahr  <- max(abs(res[, c("AHR_1", "AHR_5", "AHR_10", "AHR_100")] - 1))

  # tolerance for absolute error from integrate funciton
  expect_lt(max_abs_err_gahr, .Machine$double.eps^0.25)
  expect_lt(max_abs_err_ahr, .Machine$double.eps^0.25)

})
