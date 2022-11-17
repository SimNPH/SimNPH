#' Analyse the Dataset using the difference in RMST
#'
#' @param condition condition of the simulation
#' @param dat generated datasets
#' @param fixed_objects other constants
#'
#' @return a data frame with the columns
#' * `p` p value of the test, see Details
#' * `rmst_diff` estimated differnce in RMST
#' * `rmst_diff_lower` adjusted lower bound of the confidence interval for differnce in RMST
#' * `rmst_diff_upper` adjusted upper bound of the confidence interval for differnce in RMST
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @details
#' The implementation from the nph package is used, see the documentation there
#' for details.
#'
#' @seealso
#' \link[nph:nphparams]
#'
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_rmst_diff(condition, dat)
analyse_rmst_diff <- function(condition, dat, fixed_objects = NULL){
  model <- nph::nphparams(dat$t, dat$evt, dat$trt, param_type="RMST")

  list(
    p = model$tab$p_adjusted,
    rmst_diff = model$tab$Estimate,
    rmst_diff_lower = model$tab$lwr_adjusted,
    rmst_diff_upper = model$tab$upr_adjusted,
    N_pat=nrow(dat),
    N_evt=sum(dat$evt)
  )
}

