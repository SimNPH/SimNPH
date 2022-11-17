#' Analyse the dataset using differnce in mena survival
#'
#' @param quant=0.5 quantile for which the difference should be calculated, defaults to the median
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
#' * `diff_Q` estimated differnce in quantile of the suvivla functions
#' * `diff_Q_lower` adjusted lower bound of the confidence interval for the differnce in quantile of the suvivla functions
#' * `diff_Q_upper` adjusted upper bound of the confidence interval for the differnce in quantile of the suvivla functions
#' * `quantile` quantile used for extimation
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
#' analyse_diff_mean_survival()(condition, dat)
analyse_diff_mean_survival <- function(quant=0.5){
  function(condition, dat, fixed_objects = NULL){
    model <- nph::nphparams(dat$t, dat$evt, dat$trt, param_type="Q", param_par=quant)

    list(
      p = model$tab$p_adjusted,
      diff_Q = model$tab$Estimate,
      diff_Q_lower = model$tab$lwr_adjusted,
      diff_Q_upper = model$tab$upr_adjusted,
      quantile = quant,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )
  }
}

