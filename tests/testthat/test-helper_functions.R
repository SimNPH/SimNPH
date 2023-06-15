test_that("trycatch_nphparams works", {
 expect_warning(
   res <- trycatch_nphparams(stop("a")),
   "a"
  )

  expect_equal(res, list(
    tab=list(
      p_unadj=NA_real_,
      lwr_unadj=NA_real_,
      upr_unadj=NA_real_,
      Estimate=NA_real_
    )
  ))
})

test_that("utility functions work", {
  expect_equal(r2m((12*log(2))/365.25), 1)
  expect_equal(m2r((12*log(2)/365.25)), 1)
  expect_equal(m2d(12), 365.25)
  expect_equal(d2m(365.25), 12)
})
