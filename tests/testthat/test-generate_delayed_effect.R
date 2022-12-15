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
    label = "assumptions_delayed_effect prints something with expand.grid"
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
    cutoff_stats=c(7, 15),
    followup = 18
  )

  test_design1 <- test_design |>
    true_summary_statistics_delayed_effect(cutoff_stats = test_design$cutoff_stats)

  test_design2 <- test_design |>
    true_summary_statistics_delayed_effect()

  test_design3  <- test_design |>
    subset(select=c(-followup)) |>
    true_summary_statistics_delayed_effect()

  expect_named(test_design1, c("n_trt", "n_ctrl", "delay", "hazard_ctrl", "hazard_trt", "cutoff_stats", "followup", "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))

  expect(all(test_design1$AHR[test_design1$hazard_ctrl == test_design1$hazard_trt] == 1), "all average hazard ratios should be 1 for equal hazards")
  expect(all(test_design1$AHR[test_design1$delay >= test_design1$cutoff_stats] == 1), "all average hazard ratios should be 1 if effect starts after cutoff")

  expect(all(test_design1$AHR[(test_design1$delay < test_design1$cutoff_stats) & (test_design1$hazard_ctrl > test_design1$hazard_trt)] < 1), "all average hazard ratios should be less than 1 if there's an effect before cutoff")

  expect(all(test_design1$gAHR[test_design1$hazard_ctrl == test_design1$hazard_trt] == 1), "all geometric average hazard ratios should be 1 for equal hazards")

  expect(all(test_design1$gAHR[test_design1$delay >= test_design1$cutoff_stats] == 1), "all geometric average hazard ratios should be 1 if effect starts after cutoff")

  expect(all(test_design1$gAHR[(test_design1$delay < test_design1$cutoff_stats) & (test_design1$hazard_ctrl > test_design1$hazard_trt)] < 1), "all geometric average hazard ratios should be less than 1 if there's an effect before cutoff")

  expect(all(test_design2$cutoff_used == test_design2$followup), "cutoff used defaults to followup")

  expect(all(!is.na(test_design3$cutoff_used)), "cutoff is calculated if cutoff_stats and followup are missing")
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

test_that("test that hr_after_onset_from_gAHR works", {
  capture.output(
    my_design <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    )
  )

  my_design$hazard_trt <- NA
  my_design$hazard_ctrl <- 0.1

  my_design_A <- hr_after_onset_from_gAHR(my_design, 0.95)
  my_design_B <- hr_after_onset_from_gAHR(my_design, 0.95, cutoff = 150)

  my_design$followup <- NULL

  expect_error(hr_after_onset_from_gAHR(my_design, 0.95))
  expect(all(!is.na(my_design_A$hazard_trt)), "hazard_trt not missing")
  expect(all(!is.na(my_design_B$hazard_trt)), "hazard_trt not missing")
})

test_that("test that hr_after_onset_from_PH_effect_size works", {
  capture.output(
    my_design <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    )
  )
  my_design$hazard_trt <- NA
  my_design$hazard_ctrl <- 0.1

  suppressWarnings(
    expect_warning(
      my_design_B <- hr_after_onset_from_PH_effect_size(my_design, 0.9)
    )
  )


  my_design_E <- hr_after_onset_from_PH_effect_size(my_design, 0)
  expect_equal(my_design_E$hazard_trt, my_design$hazard_ctrl)

  expect_error(
    my_design_C <- hr_after_onset_from_PH_effect_size(my_design, followup = 200)
  )

  my_design$effect_size_ph <- 0.9
  suppressWarnings(
    expect_warning(
      my_design_D <- hr_after_onset_from_PH_effect_size(my_design, followup=200)
    )
  )

  my_design$followup <- NULL
  expect_error(hr_after_onset_from_PH_effect_size(my_design, 0.9))

  expect_equal(is.na(my_design_B$hazard_trt), c(F, F, F, F, T, T))
})

test_that("cen_rate_from_cen_prop_delayed_effect works", {
  design <- expand.grid(
   delay=seq(0, 10, by=5),
   hazard_ctrl=0.2,
   hazard_trt=c(0.02, NA),
   censoring_prop=c(0.1, 0.25, 0.01, 0),
   followup=100,
   n_trt=50,
   n_ctrl=50
  )

  result <- cen_rate_from_cen_prop_delayed_effect(design)

  expect(all(is.na(design$hazard_trt)==is.na(result$random_withdrawal)), "NA iff treatment hazard is NA")
  expect(all(result$random_withdrawal>=0, na.rm = TRUE), "all rates >= 0")
  expect(all(result$random_withdrawal[design$censoring_prop == 0]==0, na.rm = TRUE), "rate 0 if proportion 0")
})
