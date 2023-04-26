test_that("analyse_piecewise_exponential outputs plausible data.frame for delayed effect", {
  capture_output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )

  dat <- generate_delayed_effect(condition)
  my_analyse <- analyse_piecewise_exponential(cuts=c(30))
  res <- my_analyse(condition, dat)

  dat2 <- data.frame(
    t   = c(1, 1, 1, 10, 10, 10),
    trt = c(1, 1, 0,  0,  0,  1),
    evt = TRUE
  )

  my_analyse2 <- analyse_piecewise_exponential(cuts=c(2,9,11))
  res2 <- my_analyse2(condition, dat2)

  my_analyse3 <- analyse_piecewise_exponential(cuts=c(2, 9, 11), testing_only = TRUE)
  res3 <- my_analyse3(condition, dat2)

  my_analyse4 <- analyse_piecewise_exponential(cuts=c(30), testing_only = TRUE)
  res4 <- my_analyse4(condition, dat2)

  my_analyse5 <- analyse_piecewise_exponential(cuts=c(30))
  res5 <- my_analyse5(condition, dat2)

  expect_named(res, c("p", "hr_s", "p_s", "hr_s_lower", "hr_s_upper", "N_pat", "N_evt", "interval_table"), ignore.order = TRUE)
  expect(all( (res$p_s >= 0) %in% c(TRUE, NA)), "all p values >= 0 or NA")
  expect(all( (res$p_s <= 1) %in% c(TRUE, NA)), "all p values <= 1 or NA")

  expect(
    all(
      is.na(c(
        res2$hr_s["trt:intervalI2"], res2$hr_s["trt:intervalI4"],
        res2$p_s["trt:intervalI2"], res2$p_s["trt:intervalI4"],
        res2$hr_s_upper["trt:intervalI2"], res2$hr_s_upper["trt:intervalI4"],
        res2$hr_s_lower["trt:intervalI2"], res2$hr_s_lower["trt:intervalI4"]
    ))
  ),
    "expect values for intervals without events to be NA"
  )

  expect_named(res3, "p")
})

test_that("analyse_piecewise_exponential works if there's only events in one interval", {
  testdata <- withr::with_seed(1, {
    data.frame(
      t = rexp(100, 0.1),
      evt = TRUE,
      trt = rep(c(0,1), each=50)
    )
  })

  test_p_val <- analyse_piecewise_exponential(cuts=c(50, 100), testing_only = TRUE)(NA, testdata)$p
  expect_gte(test_p_val, 0.5)
})
