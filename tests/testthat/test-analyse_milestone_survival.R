test_that("analyse milestone survival works", {
  capture.output({
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  })

  dat <- generate_delayed_effect(condition)

  expect_type(analyse_milestone_survival(), "closure")
  expect_error(analyse_milestone_survival(what="something else"))
  expect_error(analyse_milestone_survival(package="something else"))

  times <- c(2,5)

  result1 <- analyse_milestone_survival(times=times)(condition, dat)
  result2 <- analyse_milestone_survival(times=times, what="diff")(condition, dat)
  result3 <- analyse_milestone_survival(times=times, what="quot", package="survival")(condition, dat)
  result4 <- analyse_milestone_survival(times=times, what="diff", package="survival")(condition, dat)

  expect_type(result1, "list")
  expect_s3_class(result1, NA)
  expect_named(result1, c("p", "milestone_surv_ratio", "milestone_surv_ratio_lower", "milestone_surv_ratio_upper", "times", "N_pat", "N_evt"))

  expect_type(result2, "list")
  expect_s3_class(result1, NA)
  expect_named(result2, c("p", "milestone_surv_diff", "milestone_surv_diff_lower", "milestone_surv_diff_upper", "times", "N_pat", "N_evt"))

  expect_type(result3, "list")
  expect_s3_class(result3, NA)
  expect_named(result3, c("milestone_surv_ratio", "times", "N_pat", "N_evt"))

  expect_type(result4, "list")
  expect_s3_class(result4, NA)
  expect_named(result4, c("milestone_surv_diff", "times", "N_pat", "N_evt"))
})
