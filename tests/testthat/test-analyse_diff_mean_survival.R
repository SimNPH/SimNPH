test_that("analyse diff means survival works", {
  capture.output({
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  })

  dat <- generate_delayed_effect(condition)

  expect_type(analyse_rmst_diff(), "closure")

  results <- analyse_diff_mean_survival()(condition, dat)
  results2 <- analyse_diff_mean_survival(0.7)(condition, dat)

  expect_type(results, "list")
  expect_s3_class(results, NA)
  expect_named(results, c("p", "diff_Q", "diff_Q_lower", "diff_Q_upper", "quantile", "N_pat", "N_evt"))
  expect_equal(results$quantile, 0.5)

  expect_type(results2, "list")
  expect_s3_class(results2, NA)
  expect_named(results2, c("p", "diff_Q", "diff_Q_lower", "diff_Q_upper", "quantile", "N_pat", "N_evt"))
  expect_equal(results2$quantile, 0.7)
})
