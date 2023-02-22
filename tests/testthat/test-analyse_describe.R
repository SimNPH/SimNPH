test_that("descriptive statistics work", {

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
        summarise = summarise_all
      )
    )
  )

  my_names <- c("describe.n_pat", "describe.max_followup", "describe.evt", "describe.evt_ctrl",
                "describe.evt_trt", "describe.study_time", "describe.sd_n_pat", "describe.sd_max_followup",
                "describe.sd_evt", "describe.sd_evt_ctrl", "describe.sd_evt_trt",
                "describe.sd_study_time")

  expect(all(my_names %in% names(sim_results)), "some expected names are missing")
})
