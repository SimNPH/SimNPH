test_that("analyse ahr works", {
  capture.output({
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  })

  dat <- generate_delayed_effect(condition)

  expect_type(analyse_ahr(), "closure")
  expect_error(analyse_ahr(type="something else"))

  results1 <- analyse_ahr()(condition, dat)
  results2 <- analyse_ahr(type="gAHR")(condition, dat)
  results3 <- analyse_ahr(max_time=50, type="AHR") (condition, dat)
  results4 <- analyse_ahr(max_time=50, type="gAHR", alternative="one.sided")(condition, dat)

  expect_type(results1, "list")
  expect_s3_class(results1, NA)
  expect_named(results1, c("p", "alternative", "AHR", "AHR_lower", "AHR_upper", "CI_level", "N_pat", "N_evt"))

  expect_type(results2, "list")
  expect_s3_class(results2, NA)
  expect_named(results2, c("p", "alternative", "gAHR", "gAHR_lower", "gAHR_upper", "CI_level", "N_pat", "N_evt"))

  expect_type(results3, "list")
  expect_s3_class(results3, NA)
  expect_named(results3, c("p", "alternative", "AHR", "AHR_lower", "AHR_upper", "CI_level", "N_pat", "N_evt"))

  expect_type(results4, "list")
  expect_s3_class(results4, NA)
  expect_named(results4, c("p", "alternative", "gAHR", "gAHR_lower", "gAHR_upper", "CI_level", "N_pat", "N_evt"))

})
