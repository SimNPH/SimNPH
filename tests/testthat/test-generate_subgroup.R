test_that("assumptions_subgroup outputs correct tibble", {
  capture_output(
    expect_invisible(
      assumptions_subgroup(),
      label = "assumptions_subgroup returns invisibly"
    )
  )

  expect_output(
    assumptions_subgroup(),
    regexp = "^expand.grid.*",
    label = "assumptions_subgroup prints something with createDesign"
  )

  capture_output(
    test_design <- assumptions_subgroup()
  )

  expect_true(
    all(hasName(
      test_design,
      c("hazard_ctrl", "hazard_trt", "hazard_subgroup", "prevalence", "random_withdrawal")
    )),
    label = "output of assumptions_subgroup has the right columns"
  )

  expect_true(
    test_design[, c("hazard_ctrl", "hazard_trt", "hazard_subgroup", "prevalence", "random_withdrawal")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of assumptions_subgroup have the right datatype"
  )

})

test_that("test that generate_subgroup outputs correct tibble", {
  capture_output(
    scenario <- merge(
        assumptions_subgroup(),
        design_fixed_followup(),
        by=NULL
      )[2, ]
  )
  one_simulation <- generate_subgroup(scenario)

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

test_that("generate_subgroup fails on proportion not between 0 and 1", {
  capture_output(
    scenario <- scenario <- merge(
      assumptions_subgroup(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )

  scenario$prevalence <- -1

  expect_error(
    generate_subgroup(scenario)
  )

  scenario$prevalence <- 2

  expect_error(
    generate_subgroup(scenario)
  )
})


test_that("generate_subgroup uses t_max, if given in fixed_objects", {
  capture_output(
    scenario <- merge(assumptions_subgroup(), design_fixed_followup(), by=NULL)[1, ]
  )

  one_simulation <- generate_subgroup(scenario, fixed_objects = list(t_max=10))

  expect(all(one_simulation$t <= 10), "all times lte t_max")
})




test_that("test that true_summary_statistics_subgroup works", {
  test_design <- createDesign(
    n_trt=50,
    n_ctrl=50,
    prevalence=c(0, 0.5, 1),
    hazard_ctrl=0.2,
    hazard_trt=c(0.2, 0.02),
    hazard_subgroup=1e-4,
    random_withdrawal=0.01,
    cutoff_stats=c(7, 15)
  )

  test_design <- test_design |>
    true_summary_statistics_subgroup(cutoff_stats = test_design$cutoff_stats)


  expect_named(test_design, c("n_trt", "n_ctrl", "prevalence", "hazard_ctrl", "hazard_trt", "hazard_subgroup", "random_withdrawal", "cutoff_stats", "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR"))

  expect(all(test_design$gAHR[(test_design$prevalence == 0) & (test_design$hazard_ctrl == test_design$hazard_trt)] == 1), "all gAHR should be 1 for equal hazards and prevalence == 0")
  expect(all(test_design$ AHR[(test_design$prevalence == 0) & (test_design$hazard_ctrl == test_design$hazard_trt)] == 1), "all AHR should be 1 for equal hazards and prevalence == 0")

})

test_that("test that true_summary_statistics_subgroup fails prevalence not in [0,1]", {
  capture_output(
    scenario <- merge(assumptions_subgroup(), design_fixed_followup(), by=NULL)[1, ]
  )

  scenario$prevalence <- -1

  expect_error(
    true_summary_statistics_subgroup(scenario)
  )

  scenario$prevalence <- 2

  expect_error(
    true_summary_statistics_subgroup(scenario)
  )
})

test_that("test that true_summary_statistics_subgroup works, if t_max is given", {
  capture_output(
    scenario <- merge(assumptions_subgroup(), design_fixed_followup(), by=NULL)[1, ]
  )

  expect_s3_class(
    true_summary_statistics_subgroup(scenario, fixed_objects=list(t_max=10)), "data.frame"
  )
})

