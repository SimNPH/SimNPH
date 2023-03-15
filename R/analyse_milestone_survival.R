#' Analyse the Dataset using difference or quotient of milestone survival
#'
#' @param times followup times at which the the survival should be compared
#' @param what "quot" for quotient and "diff" for differnce of surival probabilities
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
#' * `milestone_surv_ratio` / `milestone_surv_diff` ratio or differnce of survival probabilities
#' * `times` followup times at which the the survival are compared
#' * `N_pat` number of patients
#' * `N_evt` number of events
#' * `p` p value for the H0 that the ratios are 1 or the differnce is 0 respectively
#' * `alternative` the alternative used
#' * `milestone_surv_ratio_lower` / `milestone_surv_diff_lower` upper/lower CI for the estimate
#' * `milestone_surv_ratio_upper` / `milestone_surv_diff_upper` upper/lower CI for the estimate
#' * `CI_level` the CI level used
#'
#' @seealso
#' [nph::nphparams]
#'
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_milestone_survival(3:5)(condition, dat)
#' analyse_milestone_survival(3:5, what="diff")(condition, dat)
analyse_milestone_survival <- function(times, what="quot", level=0.95, alternative="two.sided"){
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  alt_ <- switch(alternative,
                 two.sided = "two.sided",
                 one.sided = "greater",
                 stop(gettext("'alternative' has to be either 'two.sided' or 'one.sided'."))
  )


  if (what == "quot") {
    function(condition, dat, fixed_objects = NULL){
      model <- trycatch_nphparams(nph::nphparams(
        dat$t, dat$evt, dat$trt,
        param_type="logS",
        param_par=times,
        lvl=level,
        alternative_test = alt_,
        param_alternative = alt_
        ))

      list(
        p = model$tab$p_unadj,
        alternative = alternative,
        milestone_surv_ratio = model$tab$Estimate,
        milestone_surv_ratio_lower = model$tab$lwr_unadj ,
        milestone_surv_ratio_upper = model$tab$upr_unadj,
        CI_level=level,
        times = times,
        N_pat=nrow(dat),
        N_evt=sum(dat$evt)
      )
    }
  } else if(what == "diff") {
    function(condition, dat, fixed_objects = NULL){
      model <- trycatch_nphparams(nph::nphparams(
        dat$t, dat$evt, dat$trt,
        param_type="S",
        param_par=times,
        lvl=level,
        alternative_test = alt_,
        param_alternative = alt_
        ))

      list(
        p = model$tab$p_unadj,
        alternative = alternative,
        milestone_surv_diff = model$tab$Estimate,
        milestone_surv_diff_lower = model$tab$lwr_unadj ,
        milestone_surv_diff_upper = model$tab$upr_unadj,
        CI_level=level,
        times = times,
        N_pat=nrow(dat),
        N_evt=sum(dat$evt)
      )
    }
  } else {
    stop(gettext("Invalid selection, choose 'quot' or 'diff' for argument 'what'."))
  }

}
