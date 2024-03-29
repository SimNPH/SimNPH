#' Apply Random Exponentially Distributed Censoring
#'
#' @param dat the dataset to apply the random censoring to
#' @param rate time of end of enrollment
#' @param discrete should the censoring times be rounded to whole days?
#'
#' @return
#' Returns a Function with one argument `dat` that modifies a dataset generated
#' by the generate functions by censoring the times and setting the event
#' indicator to `FALSE` for censored observations.
#' @export
#'
#' @examples
#'
# generate a sample with delayed effect
#' one_simulation <- merge(
#'   assumptions_delayed_effect(),
#'   design_fixed_followup(),
#'   by=NULL
#' ) |>
#'   head(1) |>
#'   generate_delayed_effect()
#'
#' # apply censoring to dataset
#' censored_sim <- random_censoring_exp(one_simulation, 0.01)
#'
#' # plot
#' # uncensored (blue) observations are the same for original and modified
#' # dataset
#' # censored (red) observations are smaller than the uncensored ones
#' plot(
#'   one_simulation$t,
#'   censored_sim$t,
#'   col=ifelse(censored_sim$evt, "blue", "red"),
#'   xlab = "uncensored times",
#'   ylab = "censored times"
#' )
#' abline(0,1)
random_censoring_exp <- function(dat, rate, discrete=TRUE){
  censor_fixed_time_internal <- function(dat, time_var, evt_var, cen_time){
    if(all(c(time_var, evt_var) %in% names(dat))){
      dat[[evt_var]][ dat[[time_var]] > cen_time ] <- FALSE
      dat[[time_var]] <- pmin(dat[[time_var]], cen_time)
    }
    dat
  }

  if(rate > 0){
    if(discrete){
      censoring_time <- floor(rexp(nrow(dat), rate = rate)) + 1
    } else {
      censoring_time <- rexp(nrow(dat), rate = rate)
    }
    dat <- dat |>
      censor_fixed_time_internal("t", "evt", censoring_time) |>
      censor_fixed_time_internal("t_ice", "ice", censoring_time)
  } else if (rate < 0){
    stop(gettext("Rate of random censoring has to be >= 0"))
  }
  dat
}



