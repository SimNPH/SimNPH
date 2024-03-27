test_that("if labs_from_labels works", {
  test <- mtcars

  attr(test$wt, "label") <- "weight"

  gg <- ggplot2::ggplot(test, ggplot2::aes(x=wt, y=mpg)) +
    ggplot2::geom_point()

  gg2 <- labs_from_labels(gg)

  expect_equal(gg$labels, list(x="wt", y="mpg"))
  expect_equal(gg2$labels, list(x="weight", y="mpg"))
})

test_that("results pivot longer works", {
  data_wide <- combination_tests_delayed

  data_long <- data_wide |>
    results_pivot_longer()

  data_long_2 <- data_wide |>
    results_pivot_longer(exclude_from_methods = "logrank")

  n_methods <- length(unique(data_long$method))
  cols_per_method <- (ncol(data_long_2) - ncol(data_long))

  # check dimensions of output
  expect_equal(n_methods * nrow(data_wide), nrow(data_long))
  expect_equal((n_methods-1) * nrow(data_wide), nrow(data_long_2))

  expect_equal(
    ncol(data_wide) - cols_per_method * n_methods,
    ncol(data_long) - cols_per_method - 1
  )
})

test_that("combined_plot works", {
  results_long <- combination_tests_delayed |>
    results_pivot_longer()

  expect_error(
    combined_plot(
      results_long,
      c("logrank", "mwlrt", "maxcombo"),
      c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
      "rejection_0.025",
      grid_level=0
    )
  )

  expect_error(
    combined_plot(
      results_long,
      c("logrank", "mwlrt", "maxcombo"),
      c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
      "rejection_0.025",
      split_var=6
    )
  )

  expect_error(
    combined_plot(
      results_long,
      c("logrank", "mwlrt", "maxcombo"),
      c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
      "rejection_0.025",
      split_var=-1
    )
  )

})
