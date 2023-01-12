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
