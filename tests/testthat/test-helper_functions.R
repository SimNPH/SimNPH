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
