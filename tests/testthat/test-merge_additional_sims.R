test_that("upsert_merge works", {

  a <- data.frame(x=5:2, y=5:2, a=5:2)
  b <- data.frame(x=1:4, y=1:4+10, b=1:4*10)
  c <- upsert_merge(a, b, by="x")

  expected <- data.frame(
    x = c(5L, 4L, 3L, 2L, 1L),
    y = c(5L, 14L, 13L, 12L, 11L),
    a = c(5L, 4L, 3L, 2L, NA),
    b = c(NA, 40, 30, 20, 10)
  )

  expect_equal(c, expected)
})

test_that("merge results works", {
  a <- data.frame(x=5:2, y=5:2, a=5:2)
  b <- data.frame(x=1:4, y=1:4+10, b=1:4*10)
  c <- merge_additional_results(a, b, design_names = "x")

  expected <- data.frame(
    x = c(5L, 4L, 3L, 2L, 1L),
    y = c(5L, 14L, 13L, 12L, 11L),
    a = c(5L, 4L, 3L, 2L, NA),
    b = c(NA, 40, 30, 20, 10)
  )
  attr(expected, "design_names")  <- list()

  expect_equal(c, expected)
})

test_that("merge results works without design names", {
  a <- data.frame(x=5:2, y=5:2, a=5:2)
  b <- data.frame(x=1:4, y=1:4+10, b=1:4*10)

  attr(a, "design_names") <- list(design="x")
  attr(b, "design_names") <- list(design="x")

  c <- merge_additional_results(a, b)

  expected <- data.frame(
    x = c(5L, 4L, 3L, 2L, 1L),
    y = c(5L, 14L, 13L, 12L, 11L),
    a = c(5L, 4L, 3L, 2L, NA),
    b = c(NA, 40, 30, 20, 10)
  )

  attr(expected, "design_names") <- list(design="x")

  expect_equal(c, expected)
})

test_that("merge gives warnings if no exact replication", {
  a <- data.frame(x=5:2, y=5:2, a=5:2, d=5:2)
  b <- data.frame(x=1:4, y=1:4+10, b=1:4*10, d=1:4+100)

  expect_warning(merge_additional_results(a, b, design_names = "x", descriptive_regex = "^d"))
})
