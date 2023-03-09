test_that("analyse rmst diff works", {
  capture.output({
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  })

  dat <- generate_delayed_effect(condition)
  results <- analyse_rmst_diff()(condition, dat)
  results2 <- analyse_rmst_diff(10)(condition, dat)

  expect_type(analyse_rmst_diff(), "closure")

  expect_type(results, "list")
  expect_s3_class(results, NA)
  expect_named(results, c("p", "alternative", "rmst_diff", "rmst_diff_lower", "rmst_diff_upper", "CI_level", "N_pat", "N_evt"))

  expect_type(results2, "list")
  expect_s3_class(results2, NA)
  expect_named(results2, c("p", "alternative", "rmst_diff", "rmst_diff_lower", "rmst_diff_upper", "CI_level", "N_pat", "N_evt"))
})
