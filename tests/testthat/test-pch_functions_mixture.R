test_that("checking p for mixutres works", {
  expect_error(prepare_p(c(-1,1)))

  p1 <- prepare_p(c(2,8))
  expect_equal(p1, c(0.2,0.8))
})

test_that("checking lists in mixtures works",{
  expect_true(check_lists(c(1,2), list(1,2)), "check list should return true")
  expect_invisible(check_lists(c(1,2), list(1,2)), "check list should return invisibly")

  expect_error(check_lists(c(0,1), list(1)))
})

test_that("mixture_haz_fun works", {
  # hazard is f/S
  # using functions that are not densities or survival functions but easy to
  # check
  fun_1 <- mixture_haz_fun(
    c(2, 8),
    pdfs = list(
      \(x){rep(5, length(x))},
      \(x){rep(1.25, length(x))}
    ),
    survs = list(
      \(x){rep(5, length(x))},
      \(x){rep(1.25, length(x))}
    )
  )

  expect_equal(fun_1(c(0)), 1)
  expect_equal(fun_1(c(0,10)), c(1,1))
})

test_that("mixture_cumhaz_fun works", {
  # cumhaz is -log(S)
  # using functions that are not survival functions but easy to
  # check
  fun_1 <- mixture_cumhaz_fun(
    c(2, 8),
    survs = list(
      \(x){rep(5/2, length(x))},
      \(x){rep(1.25/2, length(x))}
    )
  )

  expect_equal(fun_1(10), 0)
  expect_equal(fun_1(c(0,10)), c(0,0))
})

test_that("mixture_cdf_fun works", {
  # using functions that are not cdf functions but easy to
  # check
  fun_1 <- mixture_cdf_fun(
    c(2, 8),
    cdfs = list(
      \(x){rep(5/2, length(x))},
      \(x){rep(1.25/2, length(x))}
    )
  )

  expect_equal(fun_1(10), 1)
  expect_equal(fun_1(c(0,10)), c(1,1))
})


test_that("mixture_pdf_fun works", {
  # using functions that are not densities functions but easy to
  # check
  fun_1 <- mixture_pdf_fun(
    c(2, 8),
    pdfs = list(
      \(x){rep(5/2, length(x))},
      \(x){rep(1.25/2, length(x))}
    )
  )

  expect_equal(fun_1(10), 1)
  expect_equal(fun_1(c(0,10)), c(1,1))
})

test_that("mixture_surv_fun works", {
  # using functions that are not survival functions but easy to
  # check
  fun_1 <- mixture_surv_fun(
    c(2, 8),
    survs = list(
      \(x){rep(5/2, length(x))},
      \(x){rep(1.25/2, length(x))}
    )
  )

  expect_equal(fun_1(10), 1)
  expect_equal(fun_1(c(0,10)), c(1,1))
})

test_that("mixture_quant_fun works", {
  fun_1 <- mixture_quant_fun(c(1,3), cdfs=list(\(x){2*x},\(x){x*4/6}), quants=list(\(x){x/2},\(x){x*6/4}))
  # cdf with a jump at 0.5
  fun_2 <- mixture_quant_fun(
    1,
    cdfs=  list(\(x){0.5*x*(x<=0.5) + (0.5 + 0.5*x)*(x>0.5)}),
    quants=list(\(y){2*y*(y<=0.25) + (0.5 + 2*(y-0.75))*(y>0.75) + 0.5 *(y>0.25)*(y<=0.75)})
  )

  expect_equal(fun_1(0.5), 0.5)
  expect_equal(fun_1(c(0.25, 0.75)), c(0.25, 0.75))
  expect_equal(fun_2(0.5), 0.5)
})

test_that("mixture_rng_fun works", {
  n <- 10
  p <- c(0.2, 0.8)
  rngs <- list(
    \(n){ rnorm(n, mean=-1)},
    \(n){ rnorm(n, mean= 1)}
  )

  withr::with_seed(1, {
    numbers <- rmultinom(1, n, p)
    result_1 <- sample(c(rngs[[1]](numbers[1]), rngs[[2]](numbers[2])))
  })

  withr::with_seed(1,{
    result_2 <- mixture_rng_fun(p, rngs)(n)
  })

  expect_equal(result_1, result_2)
})
