test_that("descriptive statistics work (summary of many lines)", {

  capture_output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      tail(4) |>
      head(1)
  )

  summarise_all <- create_summarise_function(
    describe=summarise_describe()
  )

  # TODO: remove caputre warnings when sessioninfo is fixed
  tmp_output <- capture_warnings(
    capture_messages(
      capture_output(
        # runs simulations
        sim_results <- runSimulation(
          design=condition,
          replications=100,
          generate=generate_delayed_effect,
          analyse=list(
            describe=analyse_describe()
          ),
          summarise = summarise_all,
          save=FALSE
        )
      )
    )
  )

  my_names <- c("describe.n_pat", "describe.max_followup", "describe.evt", "describe.evt_ctrl",
                "describe.evt_trt", "describe.study_time", "describe.sd_n_pat", "describe.sd_max_followup",
                "describe.sd_evt", "describe.sd_evt_ctrl", "describe.sd_evt_trt",
                "describe.sd_study_time")

  expect(all(my_names %in% names(sim_results)), "some expected names are missing")
})


test_that("descriptive statistics work (one line)", {

  capture_output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )

  withr::with_seed(123,{
    dat <- generate_delayed_effect(condition) |>
      within({
        rec_time <- 0
      }) |>
      admin_censoring_time(2000)
  })

  desc <- analyse_describe()(condition, dat)
  expect_equal(desc$study_time, 2000)
})


test_that("descriptive statistics work (subgroup)", {
  capture_output(
    condition <- merge(
      assumptions_subgroup(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )

  withr::with_seed(123,{
    dat <- generate_subgroup(condition)
  })

  desc <- analyse_describe()(condition, dat)
  expect_named(
    desc,
    c("n_pat", "max_followup", "evt", "evt_ctrl", "evt_trt", "study_time",
      "subgroup", "subgroup_ctrl", "subgroup_trt"),
    ignore.order = TRUE
  )
})


test_that("descriptive statistics work (intercurrent event)", {
  capture_output(
    condition <- merge(
      assumptions_progression(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )

  withr::with_seed(123,{
    dat <- generate_progression(condition)
  })

  desc <- analyse_describe()(condition, dat)
  expect_named(
    desc,
    c("n_pat", "max_followup", "evt", "evt_ctrl", "evt_trt", "study_time",
      "ice", "ice_ctrl", "ice_trt"),
    ignore.order = TRUE
  )
})

