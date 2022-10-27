test_that("design_fixed_followup works", {
  capture_output(
    expect_invisible(
      design_fixed_followup(),
      label = "design_fixed_followup returns invisibly"
    )
  )

  expect_output(
    design_fixed_followup(),
    regexp = "^expand\\.grid.*",
    label = "design_fixed_followup prints something with createDesign"
  )

  capture_output(
    test_design <- design_fixed_followup()
  )

  expect_true(
    all(hasName(
      test_design,
      c("n_trt", "n_ctrl", "followup", "recruitment")
    )),
    label = "output of design_fixed_followup has the right columns"
  )

  expect_true(
    test_design[, c("n_trt", "n_ctrl", "followup", "recruitment")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of design_fixed_followup have the right datatype"
  )
})


test_that("design_group_sequential works", {
  capture_output(
    expect_invisible(
      design_group_sequential(),
      label = "design_group_sequential returns invisibly"
    )
  )

  expect_output(
    design_group_sequential(),
    regexp = "^expand\\.grid.*",
    label = "design_group_sequential prints something with createDesign"
  )

  capture_output(
    test_design <- design_group_sequential()
  )

  expect_true(
    all(hasName(
      test_design,
      c("n_trt", "n_ctrl", "followup", "recruitment", "interim_events", "final_events")
    )),
    label = "output of design_group_sequential has the right columns"
  )

  expect_true(
    test_design[, c("n_trt", "n_ctrl", "followup", "recruitment", "interim_events", "final_events")] |>
      sapply(is.numeric) |>
      all(),
    label = "columns of output of design_group_sequential have the right datatype"
  )
})
