test_that("weighted cox regression works", {
  capture.output(
  condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
    head(3) |>
    tail(1)
  )

  dat <- generate_delayed_effect(condition) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_time(condition$followup)


  analyse_s  <- analyse_coxph_weighted(type="S")
  analyse_sg <- analyse_coxph_weighted(type="SG")

  res_s  <- analyse_s (condition, dat)
  res_sg <- analyse_sg(condition, dat)

  expect_error(analyse_coxph_weighted(type="something else"))

  expect_type(res_s, "list")
  expect_s3_class(res_s, NA)
  expect_named(res_s, c("p", "coef", "hr", "hr_lower", "hr_upper", "N_pat", "N_evt"))

  expect_type(res_sg, "list")
  expect_s3_class(res_sg, NA)
  expect_named(res_sg, c("p", "coef", "hr", "hr_lower", "hr_upper", "N_pat", "N_evt"))
})
