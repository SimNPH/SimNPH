#' Analyse the Dataset using the difference in RMST
#'
#' @param max_time time for which the RMST is calculated
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return Returns an analysis function, that can be used in runSimulations
#'
#' @export
#'
#' @details
#' The implementation from the nph package is used, see the documentation there
#' for details.
#'
#' `alternative` can be "two.sided" for a two sided test of equality of the
#' summary statistic or "one.sided" for a one sided test testing H0: treatment
#' has equal or shorter survival than control vs. H1 treatment has longer
#' survival than control.
#'
#' The data.frame returned by the created function includes the follwing
#' columns:
#'
#' * `p` p value of the test, see Details
#' * `alternative` the alternative used
#' * `rmst_diff` estimated differnce in RMST
#' * `rmst_diff_lower` unadjusted lower bound of the confidence interval for differnce in RMST
#' * `rmst_diff_upper` unadjusted upper bound of the confidence interval for differnce in RMST
#' * `CI_level` the CI level used
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
#' analyse_rmst_diff()(condition, dat)
analyse_rmst_diff <- function(max_time = NA, level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  alt_ <- switch(alternative,
    two.sided = "two.sided",
    one.sided = "greater",
    stop(gettext("'alternative' has to be either 'two.sided' or 'one.sided'."))
  )


  function(condition, dat, fixed_objects = NULL) {
    model <- trycatch_nphparams(nph::nphparams(
      dat$t, dat$evt, dat$trt,
      param_type = "RMST",
      param_par = max_time,
      lvl = level,
      param_alternative = alt_,
      alternative_test = alt_
    ))

    list(
      p = model$tab$p_unadj,
      alternative = alternative,
      rmst_diff = model$tab$Estimate,
      rmst_diff_lower = model$tab$lwr_unadj,
      rmst_diff_upper = model$tab$upr_unadj,
      CI_level = level,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )
  }
}
