test_that("wrap all in tryCatch works", {
  list_of_functions <- list(
    a=\(){stop("test")},
    b=\(){1}
  )

  expect_error(lapply(list_of_functions, \(f){f()} ))

  list_of_functions_2 <- wrap_all_in_trycatch(list_of_functions)
  list_of_functions_3 <- wrap_all_in_trycatch(list_of_functions, error=\(e){-1})

  expect_warning(res2 <- lapply(list_of_functions_2, \(f){f()} ))
  res3 <- lapply(list_of_functions_3, \(f){f()} )

  expect_equal(res2, list(a=NA, b=1))
  expect_equal(res3, list(a=-1, b=1))
})

test_that("wrap all in preserve seed works", {
  funs1 <- list(
    b = \(){runif(1)}
  )

  funs2 <- list(
    a = \(){runif(1)},
    b = \(){runif(1)}
  ) |> wrap_all_in_preserve_seed()

  withr::with_seed(1, {
    x1 <- lapply(funs1, \(f){f()})
  })

  withr::with_seed(1, {
    x2 <- lapply(funs2, \(f){f()})
  })

  expect_equal(x1$b, x2$b)

})
