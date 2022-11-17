#' Analyse the Dataset using the difference in RMST
#'
#' @param max_time time for which the RMST is calculated
#'
#' @return Returns an analysis function, that can be used in runSimulations
#'
#' @export
#'
#' @details
#' The implementation from the nph package is used, see the documentation there
#' for details.
#'
#' The data.frame returned by the created function includes the follwing
#' columns:
#'
#' * `p` p value of the test, see Details
#' * `rmst_diff` estimated differnce in RMST
#' * `rmst_diff_lower` unadjusted lower bound of the confidence interval for differnce in RMST
#' * `rmst_diff_upper` unadjusted upper bound of the confidence interval for differnce in RMST
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @seealso
#' [nph:nphparams]
#'
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_rmst_diff()(condition, dat)
analyse_rmst_diff <- function(max_time=NA){
  function(condition, dat, fixed_objects = NULL){
    model <- nph::nphparams(dat$t, dat$evt, dat$trt, param_type="RMST", param_par=max_time)

    list(
      p = model$tab$p_unadj,
      rmst_diff = model$tab$Estimate,
      rmst_diff_lower = model$tab$lwr_unadj,
      rmst_diff_upper = model$tab$upr_unadj,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )
  }
}
