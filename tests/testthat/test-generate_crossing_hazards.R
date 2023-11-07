test_that("assumptions_crossing_hazards outputs correct tibble", {
  capture_output(
    expect_invisible(
      assumptions_crossing_hazards(),
      label = "assumptions_delayed_effect returns invisibly"
    )
  )

  expect_output(
    assumptions_crossing_hazards(print=TRUE),
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
    hazard_trt_after =c(0.2, 0.02)
  )

  test_design1 <- test_design |>
    true_summary_statistics_crossing_hazards(cutoff_stats = c(7, 15))

  test_design2 <- test_design |>
    true_summary_statistics_crossing_hazards()

  test_design4 <- test_design |>
    true_summary_statistics_crossing_hazards(cutoff_stats = c(7,15), milestones=c("milestone1"=10))

  expect_named(test_design1, c("n_trt", "n_ctrl", "crossing", "hazard_ctrl", "hazard_trt_before",
                               "hazard_trt_after", "median_survival_trt", "median_survival_ctrl",
                               "rmst_trt_7", "rmst_ctrl_7", "gAHR_7", "AHR_7", "AHRoc_7", "AHRoc_robust_7",
                               "rmst_trt_15", "rmst_ctrl_15", "gAHR_15", "AHR_15", "AHRoc_15",
                               "AHRoc_robust_15"))

  expect_named(test_design2, c("n_trt", "n_ctrl", "crossing", "hazard_ctrl", "hazard_trt_before",
                               "hazard_trt_after", "median_survival_trt", "median_survival_ctrl"))

  expect_named(test_design4, c("n_trt", "n_ctrl", "crossing", "hazard_ctrl", "hazard_trt_before",
                               "hazard_trt_after", "median_survival_trt", "median_survival_ctrl",
                               "rmst_trt_7", "rmst_ctrl_7", "gAHR_7", "AHR_7", "AHRoc_7", "AHRoc_robust_7",
                               "rmst_trt_15", "rmst_ctrl_15", "gAHR_15", "AHR_15", "AHRoc_15",
                               "AHRoc_robust_15", "milestone_survival_trt_milestone1", "milestone_survival_ctrl_milestone1"
  ))

  expect(
    with(
      test_design1,
      all(
        AHR_7[(hazard_ctrl == hazard_trt_before) & (hazard_ctrl == hazard_trt_after)] == 1
      )
    ),
    "all average hazard ratios should be 1 for equal hazards"
  )

  expect(
    with(
      test_design1,
      all(
        AHR_7[(crossing >= 7) & (hazard_trt_before == hazard_ctrl)] == 1
      )
    ),
    "all average hazard ratios should be 1 if hazards are only differnt after cutoff"
  )

  expect(
    with(
      test_design1,
      all(
        gAHR_7[(hazard_ctrl == hazard_trt_before) & (hazard_ctrl == hazard_trt_after)] == 1
      )
    ),
    "all geometric average hazard ratios should be 1 for equal hazards"
  )

  expect(
    with(
      test_design1,
      all(
        gAHR_7[(crossing >= 7) & (hazard_trt_before == hazard_ctrl)] == 1
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


test_that("test that hr_after_onset_from_PH_effect_size works", {
  capture.output(
    my_design <- merge(
      assumptions_crossing_hazards(),
      design_fixed_followup(),
      by=NULL
    )
  )

  my_design$hazard_trt_after <- NA
  my_design$hazard_ctrl <- 0.1

  my_design$followup <- NULL
  my_design$final_events <- (my_design$n_trt + my_design$n_ctrl) * 0.75

  suppressWarnings(
    expect_warning(
      my_design_B <- hr_after_crossing_from_PH_effect_size(my_design, 0.9)
    )
  )

})


test_that("cen_rate_from_cen_prop_crossing_hazards works", {
  design <- expand.grid(
    crossing=seq(0, 10, by=5),
    hazard_ctrl=0.2,
    hazard_trt_after=c(0.02, NA),
    hazard_trt_before=c(0.5, 0.2),
    censoring_prop=c(0.1, 0.25, 0.01, 0),
    followup=100,
    n_trt=50,
    n_ctrl=50
  )

  result <- cen_rate_from_cen_prop_crossing_hazards(design)

  expect(all(is.na(design$hazard_trt)==is.na(result$random_withdrawal)), "NA iff treatment hazard is NA")
  expect(all(result$random_withdrawal>=0, na.rm = TRUE), "all rates >= 0")
  expect(all(result$random_withdrawal[design$censoring_prop == 0]==0, na.rm = TRUE), "rate 0 if proportion 0")
})

