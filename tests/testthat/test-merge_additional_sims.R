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


test_that("merge gives warning on different design names", {
  old <- new <- combination_tests_delayed
  names_old <- attr(combination_tests_delayed, "design_names")$design
  attr(old, "design_names")$design <- names_old[1:17]
  attr(new, "design_names")$design <- names_old[c(1:15, 18:19)]

  expect_warning(merge_additional_results(old, new))
})

test_that("renaming results columns works", {
  new <- combination_tests_delayed |>
    rename_results_column(c(
      logrank.rejection_0.025 = "lr_rej025",
      delay = "onset"
    ))

  # colnames should be changed
  expect_equal(names(new)[1], "onset")
  expect_equal(names(new)[20], "lr_rej025")

  # attr should be changed
  expect_equal(attr(new, "design_names")$design[1], "onset")
  expect_equal(attr(new, "design_names")$sim[1], "lr_rej025")

  # everything else should stay the same
  old <- unname(combination_tests_delayed)
  attributes(old) <- NULL

  new <- unname(new)
  attributes(new) <- NULL

  expect_equal(old, new)
})

test_that("renaming results columns by pattern works", {
  new <- combination_tests_delayed |>
    rename_results_column_pattern(
      "mwlrt",
      "modestly_weighted_logrank_test"
    )

  # colnames should be changed
  expect_equal(names(new)[29:37], c("modestly_weighted_logrank_test.rejection_0.025", "modestly_weighted_logrank_test.N_missing_0.025",
                                    "modestly_weighted_logrank_test.N", "modestly_weighted_logrank_test.mean_n_pat",
                                    "modestly_weighted_logrank_test.sd_n_pat", "modestly_weighted_logrank_test.mean_n_evt",
                                    "modestly_weighted_logrank_test.sd_n_evt", "modestly_weighted_logrank_test.N_missing_n_pat",
                                    "modestly_weighted_logrank_test.N_missing_n_evt"))

  # attr should be changed
  expect_equal(attr(new, "design_names")$sim[10:18], c("modestly_weighted_logrank_test.rejection_0.025", "modestly_weighted_logrank_test.N_missing_0.025",
                                                       "modestly_weighted_logrank_test.N", "modestly_weighted_logrank_test.mean_n_pat",
                                                       "modestly_weighted_logrank_test.sd_n_pat", "modestly_weighted_logrank_test.mean_n_evt",
                                                       "modestly_weighted_logrank_test.sd_n_evt", "modestly_weighted_logrank_test.N_missing_n_pat",
                                                       "modestly_weighted_logrank_test.N_missing_n_evt"))

  # everything else should stay the same
  old <- unname(combination_tests_delayed)
  attributes(old) <- NULL

  new <- unname(new)
  attributes(new) <- NULL

  expect_equal(old, new)
})
