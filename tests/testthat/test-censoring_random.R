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

  test_data2 <- test_data
  test_data2$ice <- test_data2$evt
  test_data2$t_ice <- test_data2$t

  rate <- 0.01
  test_res <- random_censoring_exp(test_data, rate)

  expect(all(test_res$t <= test_data$t), "not all censored times less equal uncensored times")
  expect_named(test_res, c("t", "evt", "additional_column"))

  expect_error(random_censoring_exp(test_data, -1))

  test_res2 <- random_censoring_exp(test_data, 0)

  expect_equal(test_data, test_res2)

  test_res3 <- random_censoring_exp(test_data2, rate)
  expect(all(test_res3$t_ice <= test_data2$t_ice), "not all censored times less equal uncensored times")
  expect_named(test_res3, c("t", "evt", "additional_column", "ice", "t_ice"))
})
