test_that("assumptions_progression outputs correct tibble", {
  capture_output(
    expect_invisible(
      assumptions_progression(),
      label = "assumptions_progression returns invisibly"
    )
  )

  expect_output(
    assumptions_progression(),
    regexp = "^expand\\.grid.*",
    label = "assumptions_progression prints something with expand.grid"
  )

  capture_output(
    test_design <- assumptions_progression()
  )

  expect_true(
    all(hasName(
      test_design,
      c("hazard_ctrl", "hazard_trt", "hazard_after_prog", "prog_rate_ctrl", "prog_rate_trt", "random_withdrawal")
    )),
    label = "output of assumptions_delayed_effect has the right columns"
  )

  expect_true(
    test_design[, c("hazard_ctrl", "hazard_trt", "hazard_after_prog", "prog_rate_ctrl", "prog_rate_trt", "random_withdrawal")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of assumptions_delayed_effect have the right datatype"
  )

})


test_that("test that generate_progression outputs correct tibble", {
  capture_output(
    scenario <- merge(assumptions_progression(), design_fixed_followup(), by=NULL)[2, ]
  )
  one_simulation <- generate_progression(scenario)

  expect_equal(
    nrow(one_simulation),
    scenario$n_trt + scenario$n_ctrl,
    label = "nrow equals treatment + control"
  )

  expect_true(
    all(hasName(
      one_simulation,
      c("t", "trt", "evt", "t_ice", "ice")
    )),
    label = "simulated dataset has the right columns"
  )

  expect_equal(
    sapply(one_simulation[, c("t", "trt", "evt", "t_ice", "ice")], class),
    c(t="numeric", trt="integer", evt="logical", t_ice="numeric", ice="logical"),
    label = "columns of simulated dataset have the right datatypes"
  )

})


test_that("true summary statistics progression works", {
  capture_output(
    design <- merge(assumptions_progression(), design_fixed_followup(), by=NULL)
  )

  design_2 <- design
  design_2$followup <- NULL

  summaries_os     <- true_summary_statistics_progression(design, what="os")
  summaries_pfs    <- true_summary_statistics_progression(design, what="pfs")
  summaries_os_2   <- true_summary_statistics_progression(design, what="os", fixed_objects = list(t_max=10000))
  summaries_pfs_2  <- true_summary_statistics_progression(design, what="pfs", fixed_objects = list(t_max=10000))
  summaries_os_3   <- true_summary_statistics_progression(design_2, what="os")
  summaries_pfs_3  <- true_summary_statistics_progression(design_2, what="pfs")

  expect_error(true_summary_statistics_progression(design, what="something else"))

  expect_named(summaries_pfs, c(names(design), "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))
  expect_named(summaries_os , c(names(design), "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))

  expect_named(summaries_pfs_2, c(names(design), "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))
  expect_named(summaries_os_2 , c(names(design), "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))

  expect_named(summaries_pfs_3, c(names(design_2), "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))
  expect_named(summaries_os_3 , c(names(design_2), "rmst_trt", "median_survival_trt", "rmst_ctrl", "median_survival_ctrl", "gAHR", "AHR", "cutoff_used"))
})

test_that("censoring rate from censoring proportion for disease progression works", {
  design <- expand.grid(
    hazard_ctrl         = 0.001518187, # hazard under control (med. survi. 15m)
    hazard_trt          = 0.001265156, # hazard under treatment (med. surv. 18m)
    hazard_after_prog   = 0.007590934, # hazard after progression (med. surv. 3m)
    prog_rate_ctrl      = 0.001897734, # hazard rate for disease progression under control (med. time to progression 12m)
    prog_rate_trt       = c(0.001897734, 0.001423300, 0.001265156), # hazard rate for disease progression unter treatment (med. time to progression 12m, 16m, 18m)
    censoring_prop      = 0.1,         # rate of random withdrawal
    followup            = 100,         # follow up time
    n_trt               = 50,          # patients in treatment arm
    n_ctrl              = 50           # patients in control arm
  )

  design_2 <- design
  design_2$censoring_prop <- 0

  res  <- cen_rate_from_cen_prop_progression(design)
  res2 <- cen_rate_from_cen_prop_progression(design_2)

  expect(all(!is.na(res$random_withdrawal)), "some values for random_withdrawal are missing")
  expect_equal(res2$random_withdrawal, c(0,0,0))
})

test_that("progression rate from progression prop works", {
  capture_output(
    my_design <- merge(
      assumptions_progression(),
      design_fixed_followup(),
      by=NULL
    )
  )
  my_design$prog_rate_ctrl <- NULL
  my_design$prog_rate_trt <- NULL
  my_design$prog_prop_trt <- 0.2
  my_design$prog_prop_ctrl <- 0.3

  res  <- progression_rate_from_progression_prop(my_design)

  expect_named(res, c(names(my_design), c("prog_rate_trt", "prog_rate_ctrl")))
  expect_equal(order(res$prog_rate_trt), order(res$prog_prop_trt))
  expect_equal(order(res$prog_rate_ctrl), order(res$prog_prop_ctrl))
})
