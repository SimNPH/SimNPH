#' Analyse the Dataset using difference or quotient of milestone survival
#'
#' @param times followup times at which the the survival should be compared
#' @param what="quot" "quot" for quotient and "diff" for differnce of surival probabilities
#' @param package="nph" package to be used for computations: "nph" or "survival"
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
#' * `milestone_surv_ratio` / `milestone_surv_diff` ratio or differnce of survival probabilities
#' * `times` followup times at which the the survival are compared
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' If the package nph is used the following columns are also returned:
#' * `p` p value for the H0 that the ratios are 1 or the differnce is 0 respectively
#' * `milestone_surv_ratio_lower` / `milestone_surv_diff_lower` upper/lower CI for the estimate
#' * `milestone_surv_ratio_upper` / `milestone_surv_diff_upper` upper/lower CI for the estimate
#'
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
#' analyse_milestone_survival(3:5)(condition, dat)
#' analyse_milestone_survival(3:5, what="diff")(condition, dat)
#' analyse_milestone_survival(3:5, what="quot", package="survival")(condition, dat)
#' analyse_milestone_survival(3:5, what="diff", package="survival")(condition, dat)
analyse_milestone_survival <- function(times, what="quot", package="nph"){

  if (package == "nph" & what == "quot") {
    function(condition, dat, fixed_objects = NULL){
      model <- nph::nphparams(dat$t, dat$evt, dat$trt, param_type="logS", param_par=times)

      list(
        p = model$tab$p_unadj,
        milestone_surv_ratio = model$tab$Estimate,
        milestone_surv_ratio_lower = model$tab$lwr_unadj ,
        milestone_surv_ratio_upper = model$tab$upr_unadj,
        times = times,
        N_pat=nrow(dat),
        N_evt=sum(dat$evt)
      )
    }
  } else if(package == "nph" & what == "diff") {
    function(condition, dat, fixed_objects = NULL){
      model <- nph::nphparams(dat$t, dat$evt, dat$trt, param_type="S", param_par=times)

      list(
        p = model$tab$p_unadj,
        milestone_surv_diff = model$tab$Estimate,
        milestone_surv_diff_lower = model$tab$lwr_unadj ,
        milestone_surv_diff_upper = model$tab$upr_unadj,
        times = times,
        N_pat=nrow(dat),
        N_evt=sum(dat$evt)
      )
    }
  } else if(package == "survival" & what == "quot") {
    function(condition, dat, fixed_objects = NULL){
      model <- summary(survival::survfit(survival::Surv(t, evt)~trt, dat), times=times)
      ratios <- model$surv[model$strata=="trt=1"] / model$surv[model$strata=="trt=0"]

      list(
        milestone_surv_ratio = ratios,
        times = times,
        N_pat=nrow(dat),
        N_evt=sum(dat$evt)
      )
    }
  } else if(package == "survival" & what == "diff") {
    function(condition, dat, fixed_objects = NULL){
      model <- summary(survival::survfit(survival::Surv(t, evt)~trt, dat), times=times)
      diffs <- model$surv[model$strata=="trt=1"] - model$surv[model$strata=="trt=0"]

      list(
        milestone_surv_diff = diffs,
        times = times,
        N_pat=nrow(dat),
        N_evt=sum(dat$evt)
      )
    }
  } else {
    stop(gettext("Invalid selection, choose 'quot' or 'diff' for argument 'what' and 'nph' or 'survival' for argument package"))
  }

}
