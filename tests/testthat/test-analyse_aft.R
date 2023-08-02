test_that("analyse_aft works", {
  capture.output({
    condition <- merge(
        assumptions_delayed_effect(),
        design_fixed_followup(),
        by=NULL
      ) |>
      head(1)
  })
  dat <- generate_delayed_effect(condition)
  results1 <- analyse_aft()(condition, dat)
  results2 <- analyse_aft(alternative = "one.sided")(condition, dat)
  results3 <- analyse_aft(dist="lognormal")(condition, dat)


  expect_type(results1, "list")
  expect_s3_class(results1, NA)

  expect_named(results1, c("p", "alternative", "coef", "lower", "upper", "CI_level", "N_pat", "N_evt"))

  expect_lte(results2$p, results1$p)
})
