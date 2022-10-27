test_that("assumptions_delayed_effect outputs correct tibble", {
  capture_output(
    expect_invisible(
      assumptions_delayed_effect(),
      label = "assumptions_delayed_effect returns invisibly"
    )
  )

  expect_output(
    assumptions_delayed_effect(),
    regexp = "^expand\\.grid.*",
    label = "assumptions_delayed_effect prints something with createDesign"
  )

  capture_output(
    test_design <- assumptions_delayed_effect()
  )

  expect_true(
    all(hasName(
      test_design,
      c("delay", "hazard_ctrl", "hazard_trt", "random_withdrawal")
    )),
    label = "output of assumptions_delayed_effect has the right columns"
  )

  expect_true(
    test_design[, c("delay", "hazard_ctrl", "hazard_trt", "random_withdrawal")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of assumptions_delayed_effect have the right datatype"
  )

})

test_that("test that generate_delayed_effect outputs correct tibble", {
  capture_output(
    scenario <- merge(assumptions_delayed_effect(), design_fixed_followup(), by=NULL)[2, ]
  )
  one_simulation <- generate_delayed_effect(scenario)

  expect_equal(
    nrow(one_simulation),
    scenario$n_trt + scenario$n_ctrl,
    label = "nrow equals treatment + control"
  )

  expect_true(
    all(hasName(
      one_simulation,
      c("t", "trt", "evt")
    )),
    label = "simulated dataset has the right columns"
  )

  expect_equal(
    sapply(one_simulation[, c("t", "trt", "evt")], class),
    c(t="numeric", trt="numeric", evt="logical"),
    label = "columns of simulated dataset have the right datatypes"
  )

})

test_that("generate delayed effect fails on delay < 0", {
  capture_output(
    scenario <- merge(assumptions_delayed_effect(), design_fixed_followup(), by=NULL)[1, ]
  )

  scenario$delay <- -1

  expect_error(
    generate_delayed_effect(scenario)
  )
})

test_that("generate delayed effect uses t_max, if given in fixed_objects", {
  capture_output(
    scenario <- merge(assumptions_delayed_effect(), design_fixed_followup(), by=NULL)[1, ]
  )

  one_simulation <- generate_delayed_effect(scenario, fixed_objects = list(t_max=10))

  expect(all(one_simulation$t <= 10), "all times lte t_max")
})

test_that("test that generate_delayed_effect outputs correct tibble with delay=0", {
  capture_output(
    scenario <- merge(assumptions_delayed_effect(), design_fixed_followup(), by=NULL)[1, ]
  )

  one_simulation <- generate_delayed_effect(scenario)

  expect_equal(
    nrow(one_simulation),
    scenario$n_trt + scenario$n_ctrl,
    label = "nrow equals treatment + control"
  )

  expect_true(
    all(hasName(
      one_simulation,
      c("t", "trt", "evt")
    )),
    label = "simulated dataset has the right columns"
  )

  expect_equal(
    sapply(one_simulation[, c("t", "trt", "evt")], class),
    c(t="numeric", trt="numeric", evt="logical"),
    label = "columns of simulated dataset have the right datatypes"
  )

})

test_that("test that true_summary_statistics_delayed_effect works", {
  test_design <- createDesign(
    n_trt=50,
    n_ctrl=50,
    delay=c(0,7),
    hazard_ctrl=0.2,
    hazard_trt=c(0.2, 0.02),
    cutoff_stats=c(7, 15)
  )

  test_design <- test_design |>
    true_summary_statistics_delayed_effect(cutoff_stats = test_design$cutoff_stats)


  expect_named(test_design, c("n_trt", "n_ctrl", "delay", "hazard_ctrl", "hazard_trt", "cutoff_stats", "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR"))

  expect(all(test_design$AHR[test_design$hazard_ctrl == test_design$hazard_trt] == 1), "all average hazard ratios should be 1 for equal hazards")
  expect(all(test_design$AHR[test_design$delay >= test_design$cutoff_stats] == 1), "all average hazard ratios should be 1 if effect starts after cutoff")
  expect(all(test_design$AHR[(test_design$delay < test_design$cutoff_stats) & (test_design$hazard_ctrl > test_design$hazard_trt)] < 1), "all average hazard ratios should be less than 1 if there's an effect before cutoff")

  expect(all(test_design$gAHR[test_design$hazard_ctrl == test_design$hazard_trt] == 1), "all geometric average hazard ratios should be 1 for equal hazards")
  expect(all(test_design$gAHR[test_design$delay >= test_design$cutoff_stats] == 1), "all geometric average hazard ratios should be 1 if effect starts after cutoff")
  expect(all(test_design$gAHR[(test_design$delay < test_design$cutoff_stats) & (test_design$hazard_ctrl > test_design$hazard_trt)] < 1), "all geometric average hazard ratios should be less than 1 if there's an effect before cutoff")
})

test_that("test that true_summary_statistics_delayed_effect fails on delay < 0", {
  capture_output(
    scenario <- merge(assumptions_delayed_effect(), design_fixed_followup(), by=NULL)[1, ]
  )

  scenario$delay <- -1

  expect_error(
    true_summary_statistics_delayed_effect(scenario)
  )
})

test_that("test that true_summary_statistics_delayed_effect works, if t_max is given", {
  capture_output(
    scenario <- merge(assumptions_delayed_effect(), design_fixed_followup(), by=NULL)[1, ]
  )

  expect_s3_class(
    true_summary_statistics_delayed_effect(scenario, fixed_objects=list(t_max=10)), "data.frame"
  )
})
