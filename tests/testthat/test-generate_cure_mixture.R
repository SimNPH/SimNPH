test_that("desing_skeleton_cure_mixture outputs correct tibble", {
  capture_output(
    expect_invisible(
      desing_skeleton_cure_mixture(),
      label = "desing_skeleton_cure_mixture returns invisibly"
    )
  )

  expect_output(
    desing_skeleton_cure_mixture(),
    regexp = "^createDesign.*",
    label = "desing_skeleton_cure_mixture prints something with createDesign"
  )

  capture_output(
    test_design <- desing_skeleton_cure_mixture()
  )

  expect_true(
    all(hasName(
      test_design,
      c("n_trt", "n_ctrl", "hazard_ctrl", "hazard_trt", "hazard_cured", "cured_prop", "t_max")
    )),
    label = "output of desing_skeleton_cure_mixture has the right columns"
  )

  expect_true(
    test_design[, c("n_trt", "n_ctrl", "hazard_ctrl", "hazard_trt", "hazard_cured", "cured_prop", "t_max")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of desing_skeleton_cure_mixture have the right datatype"
  )

})

test_that("test that generate_cure_mixture outputs correct tibble", {
  capture_output(
    scenario <- desing_skeleton_cure_mixture()[2, ]
  )
  one_simulation <- generate_cure_mixture(scenario)

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

test_that("generate_cure_mixture fails on proportion not between 0 and 1", {
  capture_output(
    scenario <- desing_skeleton_cure_mixture() |>
      head(1)
  )

  scenario$cured_prop <- -1

  expect_error(
    generate_cure_mixture(scenario)
  )

  scenario$cured_prop <- 2

  expect_error(
    generate_cure_mixture(scenario)
  )
})
