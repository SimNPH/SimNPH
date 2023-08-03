test_that("internal function for progression works", {
  expect_error(subpop_hazVfun_simnph(c(0,10), 0.01, 0.1, 0.02, timezero=FALSE))

  # result from SimNPH functions
  result_1 <- subpop_hazVfun_simnph(
    Tint=m2d(c(0,48)),
    lambda1=m2r(36),
    lambda2=m2r(6),
    lambdaProg=m2r(24),
    timezero = TRUE
  )

  # result from reference implementation
  result_nph <- nph:::subpop_hazVFfun(
    Tint=m2d(c(0,48)),
    lambda1=m2r(36),
    lambda2=m2r(6),
    lambdaProg=m2r(24),
    timezero = TRUE
  )

  # compare names
  expect_named(result_1, names(result_nph))

  # compare values copied from input to result
  expect_identical(result_1$t, result_nph$t)
  expect_identical(result_1$t, result_nph$t)
  expect_identical(result_1$t, result_nph$t)
  expect_identical(result_1$t, result_nph$t)
  expect_identical(result_1$t, result_nph$t)
  expect_identical(result_1$t, result_nph$t)

  # compare calcluated values
  # those can be different due to different implementations
  # but should be numerically close
  expect_lt(max(abs(result_1$haz - result_nph$haz), na.rm=TRUE), 0.001)
  expect_lt(max(abs(result_1$cumhaz - result_nph$cumhaz), na.rm=TRUE), 0.01)
  expect_lt(max(abs(result_1$S - result_nph$S), na.rm=TRUE), 0.001)
  expect_lt(max(abs(result_1$F - result_nph$F), na.rm=TRUE), 0.001)

})
