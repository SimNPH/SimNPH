test_that("desing_skeleton_delayed_effect outputs correct tibble", {
  capture_output(
    expect_invisible(
      desing_skeleton_delayed_effect(),
      label = "desing_skeleton_delayed_effect returns invisibly"
    )
  )

  expect_output(
    desing_skeleton_delayed_effect(),
    regexp = "^createDesign.*",
    label = "desing_skeleton_delayed_effect prints something with createDesign"
  )

  capture_output(
    test_design <- desing_skeleton_delayed_effect()
  )

  expect_true(
    all(hasName(
      test_design,
      c("n_trt", "n_ctrl", "delay", "hazard_ctrl", "hazard_trt", "t_max")
    )),
    label = "output of desing_skeleton_delayed_effect has the right columns"
  )

  expect_true(
    test_design[, c("n_trt", "n_ctrl", "delay", "hazard_ctrl", "hazard_trt", "t_max")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of desing_skeleton_delayed_effect have the right datatype"
  )

})

test_that("test that generate_delayed_effect outputs correct tibble", {
  capture_output(
    scenario <- desing_skeleton_delayed_effect()[2, ]
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

  expect_lte(
    max(one_simulation$t, na.rm = TRUE),
    scenario$t_max,
    label = "no times after t_max"
  )

  expect_equal(
    sapply(one_simulation[, c("t", "trt", "evt")], class),
    c(t="numeric", trt="numeric", evt="logical"),
    label = "columns of simulated dataset have the right datatypes"
  )

})

test_that("generate delayed effect fails on delay < 0", {
  capture_output(
    scenario <- desing_skeleton_delayed_effect() |>
      head(1)
  )

  scenario$delay <- -1

  expect_error(
    generate_delayed_effect(scenario)
  )
})

test_that("test that generate_delayed_effect outputs correct tibble with delay=0", {
  capture_output(
    scenario <- desing_skeleton_delayed_effect() |>
      head(1)
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

  expect_lte(
    max(one_simulation$t, na.rm = TRUE),
    scenario$t_max,
    label = "no times after t_max"
  )

  expect_equal(
    sapply(one_simulation[, c("t", "trt", "evt")], class),
    c(t="numeric", trt="numeric", evt="logical"),
    label = "columns of simulated dataset have the right datatypes"
  )

})
