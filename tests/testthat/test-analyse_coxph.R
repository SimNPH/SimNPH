test_that("analyse coxph works", {
  capture.output({
    condition <- desing_skeleton_delayed_effect() |>
      head(1)
  })
  dat <- generate_delayed_effect(condition)
  results <- analyse_coxph(condition, dat)

  expect_type(results, "list")
  expect_s3_class(results, NA)

  expect_named(results, c("p", "coef", "hr", "hr_lower", "hr_upper", "N_pat", "N_evt"))
})
