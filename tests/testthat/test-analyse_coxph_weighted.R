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

  withr::with_seed(
    7,
    dat <- generate_delayed_effect(condition) |>
      recruitment_uniform(condition$recruitment) |>
      random_censoring_exp(condition$random_withdrawal) |>
      admin_censoring_time(condition$followup)
  )

  analyse_s  <- analyse_coxph_weighted(type="S")
  analyse_sg <- analyse_coxph_weighted(type="SG")
  analyse_g  <- analyse_coxph_weighted(type="G")

  analyse_s_12m  <- analyse_coxph_weighted(type="S", max_time=m2d(12))
  analyse_sg_12m <- analyse_coxph_weighted(type="SG", max_time=m2d(12))
  analyse_g_12m  <- analyse_coxph_weighted(type="G", max_time=m2d(12))

  res_s  <- analyse_s (condition, dat)
  res_sg <- analyse_sg(condition, dat)
  res_g  <- analyse_g (condition, dat)

  res_s_50  <- analyse_s_12m (condition, dat)
  res_sg_50 <- analyse_sg_12m(condition, dat)
  res_g_50  <- analyse_g_12m (condition, dat)

  expect_error(analyse_coxph_weighted(type="something else"))

  expect_type(res_s, "list")
  expect_s3_class(res_s, NA)
  expect_named(res_s, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt", "followup"))

  expect_type(res_sg, "list")
  expect_s3_class(res_sg, NA)
  expect_named(res_sg, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt", "followup"))

  expect_type(res_g, "list")
  expect_s3_class(res_g, NA)
  expect_named(res_g, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt", "followup"))

  expect_type(res_s_50, "list")
  expect_s3_class(res_s_50, NA)
  expect_named(res_s_50, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt", "followup"))

  expect_type(res_sg_50, "list")
  expect_s3_class(res_sg_50, NA)
  expect_named(res_sg_50, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt", "followup"))

  expect_type(res_g_50, "list")
  expect_s3_class(res_g_50, NA)
  expect_named(res_g_50, c("p", "alternative", "coef", "hr", "hr_lower", "hr_upper", "CI_level", "N_pat", "N_evt", "followup"))

})
