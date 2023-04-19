test_that("analyse_gehan_wilcoxon outputs plausible data.frame for delayed effect", {
  capture_output(
    condition <- merge(
      data.frame(
        delay=5,
        hazard_ctrl=0.2,
        hazard_trt =0.05,
        random_withdrawal=0
      ),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )

  dat <- withr::with_seed(1, generate_delayed_effect(condition))
  res <- analyse_gehan_wilcoxon()(condition, dat)
  res2 <- analyse_gehan_wilcoxon(alternative="one.sided")(condition, dat)

  dat2 <- dat
  dat2$trt <- 1-dat2$trt
  res3 <- analyse_gehan_wilcoxon(alternative="one.sided")(condition, dat2)

  expect_true(hasName(res, "p"), label="result has the column p")
  expect_lte(res$p, 1, label="p value le 1")
  expect_gte(res$p, 0, label="p value ge 0")

  expect_equal(res$p, res2$p*2)
  expect_equal(res3$p, 1)
})
