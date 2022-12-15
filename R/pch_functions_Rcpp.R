check_tint_lambda <- function(Tint, lambda){
  stopifnot(!is.unsorted(Tint))
  stopifnot(all(lambda >= 0))
  stopifnot(length(lambda) == length(Tint))
  invisible(TRUE)
}



#' Fast C implementation of hazard, cumulative hazard, ... for piecewise constant hazards
#'
#' @param Tint left points of the time intervals the hazards are defined on
#' @param lambda hazards in the time intervals
#'
#' @return
#' A function with one parameter, a vector of times where the function should be evaluated.
#' @export
#'
#' @describeIn fast_haz_fun fast hazard function
#'
#' @details the last time interval extends to +Inf
#'
#' @examples
#' haz <- fast_haz_fun(c(0, 10, 20), c(10, 15, 7))
#' haz(seq(0, 30, by=0.1))
fast_haz_fun    <- function(Tint, lambda){
  check_tint_lambda(Tint, lambda)
  function(v){
    .Call("_SimNPH_hazFunCpp", PACKAGE="SimNPH", Tint, lambda, v)
  }
}

#' @export
#'
#' @describeIn fast_haz_fun fast cumulative hazard function
#'
#' @examples
#' cumhaz <- fast_cumhaz_fun(c(0, 10, 20), c(10, 15, 7))
#' cumhaz(seq(0, 30, by=0.1))
fast_cumhaz_fun <- function(Tint, lambda){
  check_tint_lambda(Tint, lambda)
  function(v){
    .Call("_SimNPH_cumhazFunCpp", PACKAGE="SimNPH", Tint, lambda, v)
  }
}

#' @export
#'
#' @describeIn fast_haz_fun fast cumulative density function function
#'
#' @examples
#' cdf <- fast_cdf_fun(c(0, 10, 20), c(10, 15, 7))
#' cdf(seq(0, 30, by=0.1))
fast_cdf_fun    <- function(Tint, lambda){
  check_tint_lambda(Tint, lambda)
  function(v){
    .Call("_SimNPH_cdfFunCpp", PACKAGE="SimNPH", Tint, lambda, v)
  }
}

#' @export
#'
#' @describeIn fast_haz_fun fast probability density function
#'
#' @examples
#' pdf <- fast_pdf_fun(c(0, 10, 20), c(10, 15, 7))
#' pdf(seq(0, 30, by=0.1))
fast_pdf_fun    <- function(Tint, lambda){
  check_tint_lambda(Tint, lambda)
  function(v){
    .Call("_SimNPH_pdfFunCpp", PACKAGE="SimNPH", Tint, lambda, v)
  }
}

#' @export
#'
#' @describeIn fast_haz_fun fast survival function
#'
#' @examples
#' surv <- fast_surv_fun(c(0, 10, 20), c(10, 15, 7))
#' surv(seq(0, 30, by=0.1))
fast_surv_fun   <- function(Tint, lambda){
  check_tint_lambda(Tint, lambda)
  function(v){
    .Call("_SimNPH_survFunCpp", PACKAGE="SimNPH", Tint, lambda, v)
  }
}

#' @export
#'
#' @describeIn fast_haz_fun fast quantile function
#'
#' @examples
#' quant <- fast_quant_fun(c(0, 10, 20), c(10, 15, 7))
#' quant(seq(0, 1, by=0.01))
fast_quant_fun   <- function(Tint, lambda){
  check_tint_lambda(Tint, lambda)
  function(v){
    .Call("_SimNPH_quantFunCpp", PACKAGE="SimNPH", Tint, lambda, v)
  }
}

#' @export
#'
#' @param discrete=TRUE should surivial times be rounded down to the next whole day
#'
#' @return for fast_rng_fun: a function with the argument that draws `n` samples
#'   from the survival distribution
#'
#' @describeIn fast_haz_fun fast random numbers from survival distribution
#'
#' @examples
#' rng <- fast_rng_fun(c(0, 10, 20), c(10, 15, 7))
#' rng(100)
fast_rng_fun <- function(Tint, lambda, discrete=TRUE){
  check_tint_lambda(Tint, lambda)
  if(discrete){
    function(n){
      floor(.Call("_SimNPH_quantFunCpp", PACKAGE="SimNPH", Tint, lambda, runif(n)) + 1)
    }
  } else {
    function(n){
      .Call("_SimNPH_quantFunCpp", PACKAGE="SimNPH", Tint, lambda, runif(n))
    }
  }
}
