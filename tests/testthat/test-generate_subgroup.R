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
      c("t", "trt", "evt", "subgroup")
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


test_that("test that true_summary_statistics_subgroup works", {
  test_design <- createDesign(
    n_trt=50,
    n_ctrl=50,
    prevalence=c(0, 0.5, 1),
    hazard_ctrl=0.2,
    hazard_trt=c(0.2, 0.02),
    hazard_subgroup=1e-4,
    random_withdrawal=0.01
  )

  test_design1 <- test_design |>
    true_summary_statistics_subgroup(cutoff_stats = c(7, 15), milestones = c(10))

  expect_named(test_design1, c("n_trt", "n_ctrl", "prevalence", "hazard_ctrl", "hazard_trt",
                               "hazard_subgroup", "random_withdrawal",
                               "median_survival_trt", "median_survival_ctrl", "rmst_trt_7",
                               "rmst_ctrl_7", "gAHR_7", "AHR_7", "AHRoc_7", "AHRoc_robust_7", "rmst_trt_15", "rmst_ctrl_15",
                               "gAHR_15", "AHR_15", "AHRoc_15", "AHRoc_robust_15", "milestone_survival_trt_10", "milestone_survival_ctrl_10"
  ))

  expect(all(test_design1$gAHR_7[(test_design1$prevalence == 0) & (test_design1$hazard_ctrl == test_design1$hazard_trt)] == 1), "all gAHR should be 1 for equal hazards and prevalence == 0")
  expect(all(test_design1$ AHR_15[(test_design1$prevalence == 0) & (test_design1$hazard_ctrl == test_design1$hazard_trt)] == 1), "all AHR should be 1 for equal hazards and prevalence == 0")

  expect(all(test_design1$gAHR_7[(test_design1$prevalence == 0) & (test_design1$hazard_ctrl == test_design1$hazard_trt)] == 1), "all gAHR should be 1 for equal hazards and prevalence == 0")
  expect(all(test_design1$ AHR_15[(test_design1$prevalence == 0) & (test_design1$hazard_ctrl == test_design1$hazard_trt)] == 1), "all AHR should be 1 for equal hazards and prevalence == 0")

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

test_that("hazard_subgroup_from_PH_effect_size works", {

  capture_output(
    my_design <- merge(
      assumptions_subgroup(),
      design_fixed_followup(),
      by=NULL
    )
  )

  my_design_2 <- my_design

  my_design$hazard_trt <- NA
  my_design$hazard_subgroup <- NA
  my_design$hr_subgroup_relative <- 0.9
  my_design$final_events <- ceiling((my_design$n_ctrl + my_design$n_trt) * 0.75)

  my_design_3 <- my_design
  my_design_3$effect_size_ph <- c(0, 0.7, 0.8, 0.9)

  result_1 <- hazard_subgroup_from_PH_effect_size(my_design, target_power_ph=0.9)
  result_3 <- hazard_subgroup_from_PH_effect_size(my_design_3)

  expect_named(result_1, c("hazard_ctrl", "hazard_trt", "hazard_subgroup", "prevalence",
                            "random_withdrawal", "n_trt", "n_ctrl", "followup", "recruitment",
                            "hr_subgroup_relative", "final_events", "target_median_trt"))

  expect_lt(max(abs((result_1$hazard_subgroup / result_1$hazard_trt) - result_1$hr_subgroup_relative), na.rm = TRUE), 1e-15)
  expect_true(all(!is.na(result_1$hazard_subgroup)))
  expect_true(all(!is.na(result_1$hazard_trt)))

  expect_named(result_3, c("hazard_ctrl", "hazard_trt", "hazard_subgroup", "prevalence",
                           "random_withdrawal", "n_trt", "n_ctrl", "followup", "recruitment",
                           "hr_subgroup_relative", "final_events", "effect_size_ph", "target_median_trt"))

  expect_lt(max(abs((result_3$hazard_subgroup / result_3$hazard_trt) - result_3$hr_subgroup_relative), na.rm = TRUE), 1e-15)
  expect_true(all(!is.na(result_3$hazard_subgroup)))
  expect_true(all(!is.na(result_3$hazard_trt)))


  expect_error(hazard_subgroup_from_PH_effect_size(my_design_2, target_power_ph=0.9))
  expect_error(hazard_subgroup_from_PH_effect_size(my_design))


})

test_that("cen_rate_from_cen_prop_subgroup works", {
  design <- expand.grid(
    hazard_ctrl=0.2,                   # hazard under control and before treatment effect
    hazard_trt=0.02,                   # hazard after onset of treatment effect
    hazard_subgroup=0.01,              # hazard in the subgroup in treatment
    prevalence = c(0.2, 0.5),           # subgroup prevalence
    censoring_prop=c(0.1, 0.25, 0.01), # 10%, 25%, 1% random censoring
    followup=100,                      # followup of 100 days
    n_trt=50,                          # 50 patients treatment
    n_ctrl=50                          # 50 patients control
  )

  design$censoring_prop[1] <- 0
  design$hazard_trt[2] <- NA_real_

  result <- cen_rate_from_cen_prop_subgroup(design)

  expect_named(result, c("hazard_ctrl", "hazard_trt", "hazard_subgroup", "prevalence",
                         "censoring_prop", "followup", "n_trt", "n_ctrl", "random_withdrawal"
  ))

  expect_equal(result$random_withdrawal[1], 0)
  expect_true(is.na(result$random_withdrawal[2]))
  expect_true(all(!is.na(tail(result$random_withdrawal, -2))))

})
