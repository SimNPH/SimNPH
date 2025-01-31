test_that("if labs_from_labels works", {
  test <- mtcars

  attr(test$wt, "label") <- "weight"

  gg <- ggplot2::ggplot(test, ggplot2::aes(x=wt, y=mpg)) +
    ggplot2::geom_point()

  gg2 <- labs_from_labels(gg)

  if ("get_labs" %in% getNamespaceExports("ggplot2")) {
    # ggplot2 3.6.0 also extracts label attribute by default
    expect_equal(ggplot2::get_labs(gg)[c("x", "y")], list(x = "weight", y = "mpg"))
    expect_equal(ggplot2::get_labs(gg2)[c("x", "y")], list(x = "weight", y = "mpg"))
  } else {
    expect_equal(gg$labels, list(x="wt", y="mpg"))
    expect_equal(gg2$labels, list(x="weight", y="mpg"))
  }
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

  expect_warning(
    combined_plot(
      results_long,
      c("logrank", "mwlrt", "maxcombo"),
      c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
      "rejection_0.025",
      scale_stairs=0.7
    )
  )

  gg <- combined_plot(
    results_long,
    c("logrank", "mwlrt", "maxcombo"),
    c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
    "rejection_0.025",
    grid_level=2
  )

  labels <- if ("get_labs" %in% getNamespaceExports("ggplot2")) {
    ggplot2::get_labs(gg[[1]])
  } else {
    gg[[1]]$labels
  }

  expect_equal(labels$y, "rejection_0.025")
  expect_equal(labels$colour, "method")

  expect_equal(gg[[2]][[1]]$labels$y, "hr")
  expect_equal(gg[[2]][[2]]$labels$y, "n_pat_design")
  expect_equal(gg[[2]][[3]]$labels$y, "delay")
  expect_equal(gg[[2]][[4]]$labels$y, "hazard_ctrl")
  expect_equal(gg[[2]][[5]]$labels$y, "recruitment")

  expect_s3_class(gg, "ggplot")
  expect_s3_class(gg, "patchwork")

  expect_s3_class(gg[[1]], "ggplot")
  expect_s3_class(gg[[2]], "patchwork")
  expect_s3_class(gg[[2]], "ggplot")

  expect_s3_class(gg[[2]][[1]], "ggplot")
  expect_s3_class(gg[[2]][[2]], "ggplot")
  expect_s3_class(gg[[2]][[3]], "ggplot")
  expect_s3_class(gg[[2]][[4]], "ggplot")
  expect_s3_class(gg[[2]][[5]], "ggplot")

  expect_equal(sort(unique(gg[[1]]$data$x_split)), 1:3)

  gg2 <- combined_plot(
    results_long,
    c("logrank", "mwlrt", "maxcombo"),
    c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
    "rejection_0.025",
    grid_level=2,
    split_var = 0
  )

  expect_equal(sort(unique(gg2[[1]]$data$x_split)), 1)

  my_colours <- c(
    logrank="black",
    mwlrt="blue",
    maxcombo="green"
  )

  my_shapes <- c(
    logrank=1,
    mwlrt=2,
    maxcombo=2
  )

  gg3 <- combined_plot(
    results_long,
    c("logrank", "mwlrt", "maxcombo"),
    c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
    "rejection_0.025",
    grid_level=2,
    use_colours = my_colours,
    use_shapes = my_shapes
  )

  g3 <- ggplot2::ggplot_build(gg3[[1]])

  expect_equal(unique(g3$data[[1]][["shape"]]), 1:2)
  expect_equal(unique(g3$data[[1]][["colour"]]), c("black", "green", "blue"))
})
