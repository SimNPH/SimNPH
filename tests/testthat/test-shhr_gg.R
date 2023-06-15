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
})
