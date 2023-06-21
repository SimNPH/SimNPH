prepare_p <- function(p){
  stopifnot(all(p >= 0))
  p/sum(p)
}

check_lists <- function(p, ...){
  lists <- list(...)
  stopifnot(all(sapply(lists, length)==length(p)))
  invisible(TRUE)
}

#' Fast implementation of hazard, cumulative hazard, ... for mixtures of subpopulations
#'
#' @param p vector of probabilities of the mixture
#' @param pdfs list of probability density functions of the mixture components
#' @param survs list of survuval functions of the mixture components
#'
#' @return
#' A function with one parameter, a vector of times where the function should be evaluated.
#' @export
#'
#' @describeIn mixture_haz_fun hazard function of mixture
#'
#' @details the last time interval extends to +Inf
#'
#' @examples
#' haz <- mixture_haz_fun(
#'   p = c(0.3, 0.7),
#'   pdfs = list(
#'     fast_pdf_fun(0, 0.1),
#'     fast_pdf_fun(c(0,5), c(0.1, 0.12))
#'   ),
#'   survs = list(
#'     fast_surv_fun(0, 0.1),
#'     fast_surv_fun(c(0,5), c(0.1, 0.12))
#'   )
#' )
#' plot(haz(seq(0, 30, by=0.15)), ylim=c(0, 0.2), type="l")
#' abline(h=0)
mixture_haz_fun    <- function(p, pdfs, survs){
  p <- prepare_p(p)
  check_lists(p, pdfs, survs)

  function(v){
    pdf <- sapply(pdfs, \(pdf) pdf(v)) %*% p
    surv <- sapply(survs, \(surv) surv(v)) %*% p
    pdf/surv
  }
}

#' @describeIn mixture_haz_fun cumulative hazard function of mixture
#'
#' @export
#'
#' @examples
#' cumhaz <- mixture_cumhaz_fun(
#'   p = c(0.3, 0.7),
#'   survs = list(
#'     fast_surv_fun(0, 0.1),
#'     fast_surv_fun(c(0,5), c(0.1, 0.12))
#'   )
#' )
#' plot(cumhaz(seq(0, 30, by=0.15)), type="l")
mixture_cumhaz_fun <- function(p, survs){
  p <- prepare_p(p)
  function(v){
    surv <- sapply(survs, \(surv) surv(v)) %*% p
    -log(surv)
  }
}

#' @describeIn mixture_haz_fun cumulative density function of mixture
#'
#' @export
#'
#' @examples
#' cdf <- mixture_cdf_fun(
#'   p = c(0.3, 0.7),
#'   cdfs = list(
#'     fast_cdf_fun(0, 0.1),
#'     fast_cdf_fun(c(0,5), c(0.1, 0.12))
#'   )
#' )
#' plot(cdf(seq(0, 30, by=0.15)), type="l")
mixture_cdf_fun    <- function(p, cdfs){
  p <- prepare_p(p)
  function(v){
    sapply(cdfs, \(cdf) cdf(v)) %*% p
  }
}

#' @describeIn mixture_haz_fun probability density function of mixture
#'
#' @export
#'
#' @examples
#' pdf <- mixture_pdf_fun(
#'   p = c(0.3, 0.7),
#'   pdfs = list(
#'     fast_pdf_fun(0, 0.1),
#'     fast_pdf_fun(c(0,5), c(0.1, 0.12))
#'   )
#' )
#' plot(pdf(seq(0, 30, by=0.15)), type="l")
mixture_pdf_fun    <- function(p, pdfs){
  p <- prepare_p(p)
  function(v){
    sapply(pdfs, \(pdf) pdf(v)) %*% p
  }
}

#' @describeIn mixture_haz_fun survival function of mixture
#'
#' @export
#'
#' @examples
#' surv <- mixture_pdf_fun(
#'   p = c(0.3, 0.7),
#'   survs = list(
#'     fast_surv_fun(0, 0.1),
#'     fast_surv_fun(c(0,5), c(0.1, 0.12))
#'   )
#' )
#' plot(surv(seq(0, 30, by=0.15)), type="l")
mixture_surv_fun   <- function(p, survs){
  p <- prepare_p(p)
  function(v){
    sapply(survs, \(surv) surv(v)) %*% p
  }
}

#' @describeIn mixture_haz_fun quantile function of mixture
#'
#' @export
#'
#' @param cdfs list of cumulative density functions of the mixture components
#' @param quants list of quantile functions of the mixture components
#'
#' @details mixture_quant_fun relies on numeric root finding and is therefore
#'   not as fast as fast_quant_fun.
#'
#' @examples
#' quant <- mixture_quant_fun(
#'   p = c(0.3, 0.7),
#'   cdfs = list(
#'     fast_cdf_fun(0, 0.1),
#'     fast_cdf_fun(c(0,5), c(0.1, 0.12))
#'   )
#' )
#' plot(surv(seq(0, 30, by=0.15)), type="l")
mixture_quant_fun   <- function(p, cdfs, quants){
  p <- prepare_p(p)
  check_lists(p, cdfs, quants)

  mix_cdf <- mixture_cdf_fun(p, cdfs)

  target_fun <- function(x, v){
    mix_cdf(x) - v
  }

  function(v){
    sapply(v, \(y){
      lims <- range(sapply(quants, \(q){q(y)}))
      uniroot(
        target_fun,
        v=y,
        interval=lims,
        extendInt="upX",
        tol=2*.Machine$double.eps
      )$root
    })
  }
}

#' @describeIn mixture_haz_fun quantile function of mixture
#'
#' @export
#'
#' @param rngs random number generating functions of the components
#'
#' @details in contrast to fast_rng_fun mixture_rng fun does not rely on the
#'   quantile function but samples the counts from the respective mixtures from
#'   a multinomial distribution with parameter `p` and then samples from the
#'   components and shuffles the result. This is faster, since the quantile
#'   function relies on numeric root finding for each point.
#'
#' @examples
#' rng <- mixture_rng_fun(
#'   p = c(0.3, 0.7),
#'   rngs = list(
#'     fast_rng_fun(0, 0.1),
#'     fast_rng_fun(c(0,5), c(0.1, 0.12))
#'   )
#' )
#' hist(rng(100))
mixture_rng_fun <- function(p, rngs){
  p <- prepare_p(p)
  m <- length(p)
  function(n){
    # how many draws from wich distribution
    i <- rmultinom(1, n, p)
    # draw from rngs, shuffle
    mapply(\(r,i){r(i)}, rngs, i) |>
      unlist() |>
      sample()
  }
}
