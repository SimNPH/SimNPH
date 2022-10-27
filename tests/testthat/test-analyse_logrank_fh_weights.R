test_that("analyse_logrank_fh_weights outputs plausible data.frame for delayed effect", {
  capture_output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )
  dat <- generate_delayed_effect(condition)

  res <- analyse_logrank_fh_weights(rho=1, gamma=0)(condition, dat)
  expect_true(hasName(res, "p"), label="result has the column p")
  expect_lte(res$p, 1, label="p value le 1")
  expect_gte(res$p, 0, label="p value ge 0")
})
