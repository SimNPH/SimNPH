test_that("analyse_modelstly_weighted outputs plausible data.frame for delayed effect", {
  capture_output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )
  dat <- generate_delayed_effect(condition)
  res <- analyse_modelstly_weighted(20)(condition, dat)

  expect_true(hasName(res, "p"), label="result has the column p")
  expect_lte(res$p, 1, label="p value le 1")
  expect_gte(res$p, 0, label="p value ge 0")

  res2 <- analyse_modelstly_weighted(max(dat$t)+1)(condition, dat)

  expect(!is.nan(res2$p), "p should not be nan, even if t* is larger than max(t)")

  expect_error(analyse_modelstly_weighted(c(10, 20)))
})
