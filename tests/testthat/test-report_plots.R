test_that("if labs_from_labels works", {
  test <- mtcars

  attr(test$wt, "label") <- "weight"

  gg <- ggplot2::ggplot(test, ggplot2::aes(x=wt, y=mpg)) +
    ggplot2::geom_point()

  gg2 <- labs_from_labels(gg)

  expect_equal(gg$labels, list(x="wt", y="mpg"))
  expect_equal(gg2$labels, list(x="weight", y="mpg"))
})
