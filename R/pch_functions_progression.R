#' Fast implementation of cumulative density function, survival function, ... for scenarios with progression
#'
#' @param hazard_before hazard for death before progression
#' @param prog_rate hazard rate for progression
#' @param hazard_after hazard for death after progression
#'
#' @return
#' A function with one parameter, a vector of times/probabilities where the function should be evaluated.
#'
#' @details
#' Calculations are done by viewing the disease process as a three state
#' (non-progressed disease, progressed disease, death) continuous time markov
#' chain. Calculations can then easily be done using the matrix exponential
#' function and Q-matrices.
#'
#' @export
#'
#' @describeIn progression_cdf_fun cumulative density function for progression scenario
#'
#' @examples
#' cdf <- progression_cdf_fun(
#'   hazard_before = m2r(48),
#'   prog_rate = m2r(18),
#'   hazard_after = m2r(6)
#' )
#' t <- 0:1000
#' plot(t, cdf(t), type="l")
progression_cdf_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
                               0, -hazard_after,  hazard_after,
                               0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    sapply(v, \(v_){
      as.numeric(pi %*% Matrix::expm(v_*Q) %*% c(0,0,1))
    })
  }
}

#' @export
#'
#' @describeIn progression_cdf_fun survival function for progression scenario
#'
#' @examples
#' surv <- progression_surv_fun(
#'   hazard_before = m2r(48),
#'   prog_rate = m2r(18),
#'   hazard_after = m2r(6)
#' )
#' t <- 0:1000
#' plot(t, surv(t), type="l")
progression_surv_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
    0, -hazard_after,  hazard_after,
    0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    sapply(v, \(v_){
      1-as.numeric(pi %*% Matrix::expm(v_*Q) %*% c(0,0,1))
    })
  }
}

#' @export
#'
#' @describeIn progression_cdf_fun probability density function for progression scenario
#'
#' @examples
#' pdf <- progression_pdf_fun(
#'   hazard_before = m2r(48),
#'   prog_rate = m2r(18),
#'   hazard_after = m2r(6)
#' )
#' t <- 0:1000
#' plot(t, pdf(t), type="l")
progression_pdf_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
    0, -hazard_after,  hazard_after,
    0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    sapply(v, \(v_){
      as.numeric(pi %*% Matrix::expm(v_*Q) %*% Q %*% c(0,0,1))
    })
  }
}

#' @export
#'
#' @describeIn progression_cdf_fun hazard function for progression scenario
#'
#' @examples
#' haz <- progression_haz_fun(
#'   hazard_before = m2r(48),
#'   prog_rate = m2r(18),
#'   hazard_after = m2r(6)
#' )
#' t <- 0:1000
#' plot(t, haz(t), type="l")
progression_haz_fun <- function(hazard_before, prog_rate, hazard_after){
  pi <- c(1,0,0)
  Q <- matrix(c(
    -(hazard_before + prog_rate),     prog_rate, hazard_before,
    0, -hazard_after,  hazard_after,
    0,             0,             0
  ), 3,3, byrow=TRUE)

  function(v){
    S <- sapply(v, \(v_){
      1-as.numeric(pi %*% Matrix::expm(v_*Q) %*% c(0,0,1))
    })

    f <- sapply(v, \(v_){
      as.numeric(pi %*% Matrix::expm(v_*Q) %*% Q %*% c(0,0,1))
    })

    S/f
  }
}

#' @export
#'
#' @describeIn progression_cdf_fun quantile function for progression scenario
#'
#' @examples
#' quant <- progression_quant_fun(
#'   hazard_before = m2r(48),
#'   prog_rate = m2r(18),
#'   hazard_after = m2r(6)
#' )
#' p <- seq(0,0.99, by=.01)
#' plot(p, quant(p), type="l")
progression_quant_fun <- function(hazard_before, prog_rate, hazard_after){
  F <- progression_cdf_fun(hazard_before, prog_rate, hazard_after)
  target_fun <- function(q, p){
    F(q) - p
  }

  function(v){
    sapply(v, \(v_){
      uniroot(target_fun, interval=c(0,1), p=v_, extendInt = "upX")$root
    })
  }
}
