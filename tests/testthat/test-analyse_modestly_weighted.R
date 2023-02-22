test_that("analyse_modelstly_weighted outputs plausible data.frame for delayed effect", {
  capture_output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )
  dat <- withr::with_seed(1, generate_delayed_effect(condition))
  dat2 <- dat
  dat2$trt <- 1-dat$trt

  res  <- analyse_modelstly_weighted(20)(condition, dat)
  res2 <- analyse_modelstly_weighted(20)(condition, dat2)

  expect_lt( res$p, 0.05)
  expect_gt(res2$p, 0.05)

  expect_true(hasName(res, "p"), label="result has the column p")
  expect_lte(res$p, 1, label="p value le 1")
  expect_gte(res$p, 0, label="p value ge 0")

  res3 <- analyse_modelstly_weighted(max(dat$t)+1)(condition, dat)

  expect(!is.nan(res3$p), "p should not be nan, even if t* is larger than max(t)")

  expect_error(analyse_modelstly_weighted(c(10, 20)))

})
