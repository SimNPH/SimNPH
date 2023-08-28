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

  # for this seed one group does not have events before t=2
  # this lead to an error earlier
  dat2 <- withr::with_seed(14, generate_delayed_effect(condition))

  expect_type(analyse_milestone_survival(), "closure")
  expect_error(analyse_milestone_survival(what="something else"))
  expect_error(analyse_milestone_survival(package="something else"))

  times <- m2d(c(2,5))

  result1 <- analyse_milestone_survival(times=times)(condition, dat)
  result2 <- analyse_milestone_survival(times=times, what="diff")(condition, dat)

  expect_type(result1, "list")
  expect_s3_class(result1, NA)
  expect_named(result1, c("p", "alternative", "milestone_surv_ratio", "milestone_surv_ratio_lower", "milestone_surv_ratio_upper", "CI_level", "times", "N_pat", "N_evt"))

  expect_type(result2, "list")
  expect_s3_class(result1, NA)
  expect_named(result2, c("p", "alternative", "milestone_surv_diff", "milestone_surv_diff_lower", "milestone_surv_diff_upper", "CI_level", "times", "N_pat", "N_evt"))

})
