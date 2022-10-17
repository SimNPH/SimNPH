test_that("group sequential tests work", {
  N_sim <- 10

  capture.output({
    condition <- desing_skeleton_delayed_effect() |>
      head(1)
  })

  condition$followup <- 200
  condition$recruitment <- 50
  condition$interim_events <- 25

  my_generator <- function(condition, fixed_objects=NULL){
    generate_delayed_effect(condition, fixed_objects) |>
      recruitment_uniform(condition$recruitment)
  }

  analyse_logrank_sequential <- analyse_group_sequential(
    followup = c(condition$interim_events, condition$followup),
    followup_type = c("event", "time"),
    alpha = c(0.025, 0.05),
    analyse_functions = analyse_logrank
  )

  dat <- withr::with_seed(0, my_generator(condition))

  result <- analyse_logrank_sequential(condition, dat)

  expect_named(result, c("rejected_at_stage", "N_pat", "N_evt", "followup", "results_stages"))
  expect_type(result$rejected_at_stage, "integer")
  expect_type(result$N_pat, "integer")
  expect_type(result$N_evt, "integer")
  expect_type(result$followup, "double")
  expect_type(result$results_stages, "list")

  expect(all(result$rejected_at_stage %in% c(1, 2, Inf)), "not all rejected_at_stage between 1 and the number of interims or Inf")
  expect_lte(result$N_pat, condition$n_trt + condition$n_ctrl)
  expect_gte(result$N_pat, 0)
  expect_lte(result$N_evt, condition$n_trt + condition$n_ctrl)
  expect_gte(result$N_evt, 0)
  expect_lte(result$followup, condition$followup)
  expect_gte(result$followup, 0)

  expect_length(result$results_stages[[1]][[1]], 2)
  expect_named(result$results_stages[[1]][[1]][[1]], c("p", "N_pat", "N_evt"))
  expect_named(result$results_stages[[1]][[1]][[2]], c("p", "N_pat", "N_evt"))
})
