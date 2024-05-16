test_that("adding recruitment times works", {

  capture.output({
    data <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(2) |>
      tail(1) |>
      generate_delayed_effect()
  })

  data$additional_column <- NA

  data1 <- data |>
    recruitment_uniform(50, 10)

  expect(all(data1$rec_time <= 50), "max recruitment time not smaller than given in argument")

  expect(all(data1$rec_time >= 10), "min recruitment time not larger than given in argument")

  expect_named(data1, c("t", "trt", "evt", "additional_column", "rec_time"))

  expect_equal(round(data1$rec_time), data1$rec_time)

  withr::with_seed(123, {
    data2 <- data |>
      recruitment_uniform(50, 10, discrete = FALSE)
  })

  expect(all(data2$rec_time <= 50), "max recruitment time not smaller than given in argument")

  expect(all(data2$rec_time >= 10), "min recruitment time not larger than given in argument")

  expect_named(data2, c("t", "trt", "evt", "additional_column", "rec_time"))

  expect_true(any(data2$rec_time != round(data2$rec_time)))
})

test_that("administrative censoring after fixed time works", {
  test_data <- tibble::tribble(
    ~t, ~evt, ~rec_time,
    10, TRUE,        11,
    40, TRUE,        31,
    50, TRUE,         5,
    60, TRUE,        40,
    100, TRUE,         1,
    110, TRUE,       105,
  )

  test_data$additional_column <- NA

  followup <- 100

  # keeping the non recruited patients
  test_res <- admin_censoring_time(test_data, followup, keep_non_recruited = TRUE)

  expect_equal(test_res$t, c(10, 40, 50, 60, 99, NA), info="testing failure times in administrative censoring")
  expect_equal(test_res$evt, c(TRUE, TRUE, TRUE, TRUE, FALSE, NA), info="testing event indicators in administrative censoring")
  expect_equal(test_res$rec_time, test_data$rec_time, info="check that rec_times is not changed")
  expect_named(test_res, c("t", "evt", "rec_time", "additional_column"))
  expect(!is.null(attr(test_res, "followup")), "check followup attribute exists")
  expect(attr(test_res, "followup") == followup, "check followup attribute has the correct value")

  # without keeping the non recruited patients
  test_res <- admin_censoring_time(test_data, followup)

  expect_equal(test_res$t, c(10, 40, 50, 60, 99), info="testing failure times in administrative censoring")
  expect_equal(test_res$evt, c(TRUE, TRUE, TRUE, TRUE, FALSE), info="testing event indicators in administrative censoring")
  expect_named(test_res, c("t", "evt", "rec_time", "additional_column"))
  expect(!is.null(attr(test_res, "followup")), "check followup attribute exists")
  expect(attr(test_res, "followup") == followup, "check followup attribute has the correct value")


  test_data_ice <- test_data
  test_data_ice$t_ice <- test_data_ice$t
  test_data_ice$ice <- test_data$evt

  test_res <- admin_censoring_time(test_data_ice, followup, keep_non_recruited = TRUE)

  expect_equal(test_res$t_ice, c(10, 40, 50, 60, 99, NA), info="testing failure times in administrative censoring")
  expect_equal(test_res$ice, c(TRUE, TRUE, TRUE, TRUE, FALSE, NA), info="testing event indicators in administrative censoring")
  expect_equal(test_res$rec_time, test_data_ice$rec_time, info="check that rec_times is not changed")
  expect_named(test_res, c("t", "evt", "rec_time", "additional_column", "t_ice", "ice"))
  expect(!is.null(attr(test_res, "followup")), "check followup attribute exists")
  expect(attr(test_res, "followup") == followup, "check followup attribute has the correct value")
})

test_that("administrative censoring after fixed number of events", {
  test_data <- tibble::tribble(
    ~t, ~evt, ~rec_time,
    10, TRUE,        11,
    40, TRUE,        30,
    50, TRUE,         5,
    60, TRUE,        41,
    100, TRUE,         1,
    110, TRUE,       105,
  )

  test_data$additional_column <- NA

  # keeping the non-recruited patients
  test_res <- admin_censoring_events(test_data, 3, keep_non_recruited = TRUE)

  expect_equal(test_res$t, c(10, 40, 50, 29, 69, NA), info="testing failure times in administrative censoring")
  expect_equal(test_res$evt, c(TRUE, TRUE, TRUE, FALSE, FALSE, NA), info="testing event indicators in administrative censoring")
  expect_equal(test_res$rec_time, test_data$rec_time, info="check that rec_times is not changed")
  expect_named(test_res, c("t", "evt", "rec_time", "additional_column"))
  expect(!is.null(attr(test_res, "followup")), "check followup attribute exists")
  expect(attr(test_res, "followup") == 70, "check followup attribute has the correct value")

  # without keeping the non recruited patients
  test_res2 <- admin_censoring_events(test_data, 3)

  expect_equal(test_res2$t, c(10, 40, 50, 29, 69), info="testing failure times in administrative censoring")
  expect_equal(test_res2$evt, c(TRUE, TRUE, TRUE, FALSE, FALSE), info="testing event indicators in administrative censoring")
  expect_named(test_res2, c("t", "evt", "rec_time", "additional_column"))
  expect(!is.null(attr(test_res2, "followup")), "check followup attribute exists")
  expect(attr(test_res2, "followup") == 70, "check followup attribute has the correct value")

  # with more events than are present in the dataset
  test_res_3 <- admin_censoring_events(test_data, nrow(test_data)+1)
  expect(all(test_res_3$evt), "If less events than targeted, no event should be excluded")

  # with warning if there are not enough events
  expect_warning(admin_censoring_events(test_data, nrow(test_data)+1, on_incomplete = "warn"))

  # with error if there are not enough events
  expect_error(admin_censoring_events(test_data, nrow(test_data)+1, on_incomplete = "stop"))
})
