test_that("random censoring works", {
  test_data <- tibble::tribble(
    ~t, ~evt,
    10, TRUE,
    40, TRUE,
    50, TRUE,
    60, TRUE,
    100, TRUE,
    110, TRUE,
  )

  test_data$additional_column <- NA

  rate <- 0.01
  test_res <- random_censoring_exp(test_data, rate)

  expect(all(test_res$t <= test_data$t), "not all censored times less equal uncensored times")
  expect_named(test_res, c("t", "evt", "additional_column"))
})
