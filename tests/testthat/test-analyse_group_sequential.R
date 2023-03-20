test_that("group sequential tests work", {

  capture.output({
    condition <-  merge(
        assumptions_delayed_effect(),
        design_fixed_followup(),
        by=NULL
      ) |>
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
    analyse_functions = analyse_logrank()
  )

  dat <- withr::with_seed(0, my_generator(condition))

  result <- analyse_logrank_sequential(condition, dat)

  expect_named(result, c("rejected_at_stage", "N_pat", "N_evt", "study_time", "max_followup", "results_stages"))
  expect_type(result$rejected_at_stage, "integer")
  expect_type(result$N_pat, "integer")
  expect_type(result$N_evt, "integer")
  expect_type(result$study_time, "double")
  expect_type(result$max_followup, "double")
  expect_type(result$results_stages, "list")

  expect(all(result$rejected_at_stage %in% c(1, 2, Inf)), "not all rejected_at_stage between 1 and the number of interims or Inf")
  expect_lte(result$N_pat, condition$n_trt + condition$n_ctrl)
  expect_gte(result$N_pat, 0)
  expect_lte(result$N_evt, condition$n_trt + condition$n_ctrl)
  expect_gte(result$N_evt, 0)
  expect_lte(result$max_followup, condition$followup)
  expect_gte(result$max_followup, 0)

  expect_length(result$results_stages[[1]][[1]], 2)
  expect_named(result$results_stages[[1]][[1]][[1]], c("p", "alternative", "N_pat", "N_evt"))
  expect_named(result$results_stages[[1]][[1]][[2]], c("p", "alternative", "N_pat", "N_evt"))
})

test_that("summarising results from group sequential function works", {
  capture.output({
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
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
    analyse_functions = analyse_logrank()
  )

  dat <- withr::with_seed(0, replicate(2, my_generator(condition), simplify=FALSE))

  results <- lapply(dat, \(dat){analyse_logrank_sequential(condition, dat)})

  results <- results |>
    purrr::map(tibble::as_tibble) |>
    dplyr::bind_rows()

  aggregate_results <- summarise_group_sequential()(condition, results)

  expect_s3_class(aggregate_results, "data.frame")
  expect_named(aggregate_results, c("rejection", "n_pat", "n_evt", "study_time", "followup", "sd_nevt", "sd_study_time", "sd_followup", "sd_npat", "N_missing_rejection", "N_missing_npat", "N_missing_nevt", "N_missing_study_time", "N_missing_followup", "N"), ignore.order = TRUE)
  expect_equal(nrow(aggregate_results), 1)
})

test_that("summarise group sequential deals with missings correctly", {
  tmp_results <- tibble::tribble(
    ~rejected_at_stage,   ~N_pat,   ~N_evt, ~max_followup, ~study_time,
                     1,       10,        5,        50,              60,
                     2,       20,       10,       100,             110,
              NA_real_, NA_real_, NA_real_,  NA_real_,        NA_real_,
                   Inf,       30,       15,       150,        NA_real_,
  )

  my_summarise <- summarise_group_sequential()
  results <- my_summarise(NA, tmp_results)

  expect_equal(
    results,
    data.frame(
      rejection = 2/3,
      n_pat = 20,
      n_evt = 10,
      study_time = 170/2,
      followup = 100,
      sd_npat = 10,
      sd_nevt = 5,
      sd_study_time = sd(c(60, 110)),
      sd_followup = 50,
      N = 4,
      N_missing_rejection = 1,
      N_missing_npat = 1,
      N_missing_nevt = 1,
      N_missing_study_time = 2,
      N_missing_followup = 1
    )
  )
})














