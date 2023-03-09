#' Analyse the dataset using extimators for the the average hazard ratio
#'
#' @param max_time time for which the RMST is calculated
#' @param type "AHR" for average hazard ratio "gAHR" for geometric average hazard ratio
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
#' * `AHR`/`gAHR` estimated (geometric) average hazard ratio
#' * `AHR_lower`/`gAHR_lower` unadjusted lower bound of the confidence interval for the (geometric) average hazard ratio
#' * `AHR_upper`/`gAHR_upper` unadjusted upper bound of the confidence interval for the (geometric) average hazard ratio
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
#' analyse_ahr()(condition, dat)
#' analyse_ahr(type = "gAHR")(condition, dat)
#' analyse_ahr(max_time = 50, type = "AHR")(condition, dat)
#' analyse_ahr(max_time = 50, type = "gAHR")(condition, dat)
analyse_ahr <- function(max_time = NA, type = "AHR", level = 0.95, alternative = "two.sided") {
  alt_ <- switch(alternative,
    two.sided = "two.sided",
    one.sided = "less",
    stop(gettext("'alternative' has to be either 'two.sided' or 'one.sided'."))
  )

  switch(type,
    AHR = {
      function(condition, dat, fixed_objects = NULL) {
        model <- trycatch_nphparams(nph::nphparams(
          dat$t, dat$evt, dat$trt,
          param_type = "avgHR",
          param_par = max_time,
          lvl = level,
          param_alternative = alt_,
          alternative_test = alternative
        ))

        list(
          p = model$tab$p_unadj,
          alternative = alternative,
          AHR = model$tab$Estimate,
          AHR_lower = model$tab$lwr_unadj,
          AHR_upper = model$tab$upr_unadj,
          CI_level = level,
          N_pat = nrow(dat),
          N_evt = sum(dat$evt)
        )
      }
    },
    gAHR = {
      function(condition, dat, fixed_objects = NULL) {
        model <- trycatch_nphparams(nph::nphparams(
          dat$t, dat$evt, dat$trt,
          param_type = "HR",
          param_par = max_time,
          lvl = level,
          param_alternative = alt_,
          alternative_test = alternative
        ))

        list(
          p = model$tab$p_unadj,
          alternative = alternative,
          gAHR = model$tab$Estimate,
          gAHR_lower = model$tab$lwr_unadj,
          gAHR_upper = model$tab$upr_unadj,
          CI_level = level,
          N_pat = nrow(dat),
          N_evt = sum(dat$evt)
        )
      }
    },
    {
      stop(gettext("Invalid type, use AHR for the average hazard ratio or gAHR for the geometric average hazard ratio."))
    }
  )
}
