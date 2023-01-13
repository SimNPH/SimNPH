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
    cutoff_stats=c(7, 15),
    followup = 18
  )

  test_design1 <- test_design |>
    true_summary_statistics_crossing_hazards(cutoff_stats = test_design$cutoff_stats)

  test_design2 <- test_design |>
    true_summary_statistics_crossing_hazards()

  test_design3  <- test_design |>
    subset(select=c(-followup)) |>
    true_summary_statistics_crossing_hazards()

  test_design4 <- test_design |>
    true_summary_statistics_crossing_hazards(cutoff_stats = test_design$cutoff_stats, milestones=c("milestone1"=10))

  expect_named(test_design1, c("n_trt", "n_ctrl", "crossing", "hazard_ctrl", "hazard_trt_before", "hazard_trt_after", "cutoff_stats", "followup", "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))

  expect_named(test_design4, c("n_trt", "n_ctrl", "crossing", "hazard_ctrl", "hazard_trt_before", "hazard_trt_after", "cutoff_stats", "followup", "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "milestone1_trt", "milestone1_ctrl", "cutoff_used"))

  expect(
    with(
      test_design1,
      all(
        AHR[(hazard_ctrl == hazard_trt_before) & (hazard_ctrl == hazard_trt_after)] == 1
      )
    ),
    "all average hazard ratios should be 1 for equal hazards"
  )

  expect(
    with(
      test_design1,
      all(
        AHR[(crossing >= cutoff_stats) & (hazard_trt_before == hazard_ctrl)] == 1
      )
    ),
    "all average hazard ratios should be 1 if hazards are only differnt after cutoff"
  )

  expect(
    with(
      test_design1,
      all(
        gAHR[(hazard_ctrl == hazard_trt_before) & (hazard_ctrl == hazard_trt_after)] == 1
      )
    ),
    "all geometric average hazard ratios should be 1 for equal hazards"
  )

  expect(
    with(
      test_design1,
      all(
        gAHR[(crossing >= cutoff_stats) & (hazard_trt_before == hazard_ctrl)] == 1
      )
    ),
    "all geometric average hazard ratios should be 1 if hazards are only differnt after cutoff"
  )

  expect(all(test_design2$cutoff_used == test_design2$followup), "cutoff used defaults to followup")

  expect(all(!is.na(test_design3$cutoff_used)), "cutoff is calculated if cutoff_stats and followup are missing")
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


test_that("test that hr_after_onset_from_gAHR works", {
  capture.output(
    my_design <- merge(
      assumptions_crossing_hazards(),
      design_fixed_followup(),
      by=NULL
    )
  )

  my_design$hazard_trt_after <- NA
  my_design$hazard_ctrl <- 0.1

  my_design_A <- hr_after_crossing_from_gAHR(my_design, 0.95)
  my_design_B <- hr_after_crossing_from_gAHR(my_design, 0.95, cutoff = 150)

  my_design$followup <- NULL

  expect_error(hr_after_crossing_from_gAHR(my_design, 0.95))
  expect(all(!is.na(my_design_A$hazard_trt_after)), "hazard_trt_after not missing")
  expect(all(!is.na(my_design_B$hazard_trt_after)), "hazard_trt_after not missing")
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

