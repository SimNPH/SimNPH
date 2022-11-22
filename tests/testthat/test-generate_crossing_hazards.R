test_that("assumptions_crossing_hazards outputs correct tibble", {
  capture_output(
    expect_invisible(
      assumptions_crossing_hazards(),
      label = "assumptions_delayed_effect returns invisibly"
    )
  )

  expect_output(
    assumptions_crossing_hazards(),
    regexp = "^expand\\.grid.*",
    label = "assumptions_crossing_hazards prints something with expand.grid"
  )

  capture_output(
    test_design <- assumptions_crossing_hazards()
  )

  expect_true(
    all(hasName(
      test_design,
      c("crossing", "hazard_ctrl", "hazard_trt_before", "hazard_trt_after", "random_withdrawal")
    )),
    label = "output of assumptions_crossing_hazards has the right columns"
  )

  expect_true(
    test_design[, c("crossing", "hazard_ctrl", "hazard_trt_before", "hazard_trt_after", "random_withdrawal")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of assumptions_delayed_effect have the right datatype"
  )

})

test_that("test that generate_crossing_hazards outputs correct tibble", {
  capture_output(
    scenario <- merge(assumptions_crossing_hazards(), design_fixed_followup(), by=NULL)[2, ]
  )
  one_simulation <- generate_crossing_hazards(scenario)

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

test_that("generate_crossing_hazards fails on crossing < 0", {
  capture_output(
    scenario <- merge(assumptions_crossing_hazards(), design_fixed_followup(), by=NULL)[1, ]
  )

  scenario$crossing <- -1

  expect_error(
    generate_crossing_hazards(scenario)
  )
})

test_that("generate_crossing_hazards uses t_max, if given in fixed_objects", {
  capture_output(
    scenario <- merge(assumptions_crossing_hazards(), design_fixed_followup(), by=NULL)[1, ]
  )

  one_simulation <- generate_crossing_hazards(scenario, fixed_objects = list(t_max=10))

  expect(all(one_simulation$t <= 10), "all times lte t_max")
})

test_that("test that generate_crossing_hazards outputs correct tibble with delay=0", {
  capture_output(
    scenario <- merge(assumptions_crossing_hazards(), design_fixed_followup(), by=NULL)[1, ]
  )

  one_simulation <- generate_crossing_hazards(scenario)

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

test_that("test that true_summary_statistics_crossing_hazards works", {
  test_design <- createDesign(
    n_trt=50,
    n_ctrl=50,
    crossing=c(0,7),
    hazard_ctrl=0.2,
    hazard_trt_before=c(0.5, 0.2),
    hazard_trt_after =c(0.2, 0.02),
    cutoff_stats=c(7, 15)
  )

  test_design <- test_design |>
    true_summary_statistics_crossing_hazards(cutoff_stats = test_design$cutoff_stats)


  expect_named(test_design, c("n_trt", "n_ctrl", "crossing", "hazard_ctrl", "hazard_trt_before", "hazard_trt_after", "cutoff_stats", "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR"))

  expect(
    with(
      test_design,
      all(
        AHR[(hazard_ctrl == hazard_trt_before) & (hazard_ctrl == hazard_trt_after)] == 1
      )
    ),
    "all average hazard ratios should be 1 for equal hazards"
  )

  expect(
    with(
      test_design,
      all(
        AHR[(crossing >= cutoff_stats) & (hazard_trt_before == hazard_ctrl)] == 1
      )
    ),
    "all average hazard ratios should be 1 if hazards are only differnt after cutoff"
  )

  expect(
    with(
      test_design,
      all(
        gAHR[(hazard_ctrl == hazard_trt_before) & (hazard_ctrl == hazard_trt_after)] == 1
      )
    ),
    "all geometric average hazard ratios should be 1 for equal hazards"
  )

  expect(
    with(
      test_design,
      all(
        gAHR[(crossing >= cutoff_stats) & (hazard_trt_before == hazard_ctrl)] == 1
      )
    ),
    "all geometric average hazard ratios should be 1 if hazards are only differnt after cutoff"
  )


})

test_that("test that true_summary_statistics_crossing_hazards fails on crossing < 0", {
  capture_output(
    scenario <- merge(assumptions_crossing_hazards(), design_fixed_followup(), by=NULL)[1, ]
  )

  scenario$crossing <- -1

  expect_error(
    true_summary_statistics_crossing_hazards(scenario)
  )
})

test_that("test that true_summary_statistics_crossing_hazards works, if t_max is given", {
  capture_output(
    scenario <- merge(assumptions_crossing_hazards(), design_fixed_followup(), by=NULL)[1, ]
  )

  expect_s3_class(
    true_summary_statistics_crossing_hazards(scenario, fixed_objects=list(t_max=10)), "data.frame"
  )
})
