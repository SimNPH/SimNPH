test_that("analyse coxph works", {
  capture.output({
    condition <- merge(
        assumptions_delayed_effect(),
        design_fixed_followup(),
        by=NULL
      ) |>
      head(1)
  })
  dat <- generate_delayed_effect(condition)
  results <- analyse_coxph()(condition, dat)
  results2 <- analyse_coxph(alternative = "one.sided")(condition, dat)

  expect_error({
    analyse_coxph(alternative="no_allowed_alternative")
  })

  expect_type(results, "list")
  expect_s3_class(results, NA)
  expect_named(results, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt"))

  expect_type(results2, "list")
  expect_s3_class(results2, NA)
  expect_named(results2, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt"))
})
