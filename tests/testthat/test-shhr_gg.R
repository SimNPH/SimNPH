test_that("shhr ggplot works", {

  B <- nph::pchaz(c(0, 10, 100), c(0.1, 0.05))
  A <- nph::pchaz(c(0, 100), c(0.1))
  withr::with_package("ggplot2", {
    withr::with_package("patchwork", {
      my_obj <- shhr_gg(A, B)
    })
  })

  expect_s3_class(my_obj, "patchwork")

  withr::with_package("ggplot2", {
    withr::with_package("patchwork", {
      my_obj <- shhr_gg(A, B, lab_time="Months", trafo_time=d2m)
    })
  })

  expect_s3_class(my_obj, "patchwork")

  expect_equal(my_obj[[1]]$labels$y, "Survival")
  expect_equal(my_obj[[2]]$labels$y, "Hazard")
  expect_equal(my_obj[[3]]$labels$y, "Hazard ratio")


  withr::with_package("ggplot2", {
    withr::with_package("patchwork", {
      my_obj2 <- shhr_gg(A, B, lab_time="Months", trafo_time=d2m, as_list=TRUE)
    })
  })

  expect_equal(my_obj2[[1]]$labels$y, "Survival")
  expect_equal(my_obj2[[2]]$labels$y, "Hazard")
  expect_equal(my_obj2[[3]]$labels$y, "Hazard ratio")

  expect_s3_class(my_obj2, NA)

  my_obj3 <- shhr_gg(A, B)
})

test_that("shhr gives warning if packages are missing", {
  B <- nph::pchaz(c(0, 10, 100), c(0.1, 0.05))
  A <- nph::pchaz(c(0, 100), c(0.1))

  # overwrite requireNamespace
  old_fn <- base::requireNamespace
  myrequireNamespace <- function(...) FALSE
  unlockBinding("requireNamespace", as.environment("package:base"))
  assign("requireNamespace",myrequireNamespace, "package:base")

  expect_message(my_obj <- shhr_gg(A, B))
  expect_null(my_obj)

  # restore normal search path
  assign("requireNamespace",old_fn, "package:base")
  lockBinding("requireNamespace", as.environment("package:base"))
  rm(old_fn, myrequireNamespace)
})
