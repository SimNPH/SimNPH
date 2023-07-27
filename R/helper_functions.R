# run nphparams with trycatch, used in various analyse functions
trycatch_nphparams <- function(call){
  tryCatch(
    call,
    error = function(e){
      warning(e$message)
      list(
        tab=list(
          p_unadj=NA_real_,
          lwr_unadj=NA_real_,
          upr_unadj=NA_real_,
          Estimate=NA_real_
        )
      )
    })
}

# hazard ratio required, inverted Schönfeld sample size formula
#
# Nevt: number of events
# alpha: target alpha
# beta: target beta
# p: allocation ratio
hr_required_schoenfeld <- function(Nevt, alpha=0.05, beta=0.2, p=0.5){
  if(beta == 1) return(1)
  exp( (qnorm(beta) + qnorm(alpha) ) / sqrt(p*(1-p)*Nevt) )
}

# censoring until t_max from cumhaz
#
# n_trt: patients in the treatment arm
# n_ctrl: patients in the control arm
# censoring_prop: target proportion of censored patients up until t_max
# cumhaz_ctrl: value of the cumulative hazard function for the control arm at t_max
# cumhaz_trt: value of the cumulative hazard function for the treatment arm at t_max
# t_max: time at which the proportion is evaluated
censoring_prop_from_cumhaz <- function(n_trt, n_ctrl, censoring_prop, cumhaz_ctrl, cumhaz_trt, t_max){
  alloc_ratio <- n_trt / (n_trt + n_ctrl)

  a <- 1-censoring_prop
  b <- alloc_ratio*cumhaz_trt + (1-alloc_ratio)*cumhaz_ctrl - (cumhaz_ctrl+cumhaz_trt)*censoring_prop
  c <- -cumhaz_ctrl*cumhaz_trt*censoring_prop

  (-b + sqrt(b^2 - 4*a*c))/(2*a*t_max)
}

#' Functions to Convert Between Days and Months and Medians and Rates
#'
#' Some functions to convert between days and months and rates and medians.
#'
#' @describeIn r2m daily rate to median in months
#'
#' @param lambda hazard rate
#'
#' @return median survival time in months (`r2m`)
#'
#' @examples
#' r2m(0.002)
#' @export
r2m <- function(lambda){
  12*log(2)/(365.25*lambda)
}

#' @describeIn r2m median to months to daily rate
#'
#' @param med median in months
#'
#' @return hazard rate per day (`m2r`)
#' @export
#'
#' @examples
#' m2r(12)
m2r <- function(med){
  (12*log(2))/(med*365.25)
}

#' @describeIn r2m months to days
#'
#' @param mon time in months
#'
#' @return time in days (`m2d`)
#' @export
#'
#' @examples
#' m2d(1)
m2d <- function(mon){
  365.25*mon/12
}

#' @describeIn r2m days to months
#'
#' @param day time in days
#'
#' @return time in months (`d2m`)
#' @export
#'
#' @examples
#' d2m(31)
d2m <- function(day){
  12*day/365.25
}
