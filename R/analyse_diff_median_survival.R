#' Analyse the dataset using differnce in median survival
#'
#' @param quant quantile for which the difference should be calculated, defaults to the median
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
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
#' * `diff_Q_lower` unadjusted lower bound of the confidence interval for the differnce in quantile of the suvivla functions
#' * `diff_Q_upper` unadjusted upper bound of the confidence interval for the differnce in quantile of the suvivla functions
#' * `quantile` quantile used for extimation
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @seealso
#' [nph::nphparams]
#'
#' @examples
#' condition <- merge(
#'   assumptions_delayed_effect(),
#'   design_fixed_followup(),
#'   by = NULL
#' ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_diff_median_survival()(condition, dat)
analyse_diff_median_survival <- function(quant = 0.5, level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  alt_ <- switch(alternative,
    two.sided = "two.sided",
    one.sided = "greater",
    stop(gettext("'alternative' has to be either 'two.sided' or 'one.sided'."))
  )

  function(condition, dat, fixed_objects = NULL) {
    model <- trycatch_nphparams(nph::nphparams(
      dat$t, dat$evt, dat$trt,
      param_type = "Q",
      param_par = quant,
      lvl = level,
      alternative_test = alt_
    ))

    list(
      p = model$tab$p_unadj,
      alternative = alternative,
      diff_Q = model$tab$Estimate,
      diff_Q_lower = model$tab$lwr_unadj,
      diff_Q_upper = model$tab$upr_unadj,
      CI_level = level,
      quantile = quant,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )
  }
}
