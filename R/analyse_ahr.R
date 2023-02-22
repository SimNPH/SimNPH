#' Analyse the dataset using extimators for the the average hazard ratio
#'
#' @param max_time time for which the RMST is calculated
#' @param type "AHR" for average hazard ratio "gAHR" for geometric average hazard ratio
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
#' * `AHR`/`gAHR` estimated (geometric) average hazard ratio
#' * `AHR_lower`/`gAHR_lower` unadjusted lower bound of the confidence interval for the (geometric) average hazard ratio
#' * `AHR_upper`/`gAHR_upper` unadjusted upper bound of the confidence interval for the (geometric) average hazard ratio
#' * `N_pat` number of patients
#' * `N_evt` number of events
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
#' analyse_ahr()(condition, dat)
#' analyse_ahr(type="gAHR")(condition, dat)
#' analyse_ahr(max_time=50, type="AHR")(condition, dat)
#' analyse_ahr(max_time=50, type="gAHR")(condition, dat)
analyse_ahr <- function(max_time=NA, type="AHR"){

  switch(
    type,
    AHR = {
      function(condition, dat, fixed_objects = NULL){
        model <- trycatch_nphparams(nph::nphparams(dat$t, dat$evt, dat$trt, param_type="avgHR", param_par=max_time))

        list(
          p = model$tab$p_unadj,
          AHR = model$tab$Estimate,
          AHR_lower = model$tab$lwr_unadj,
          AHR_upper = model$tab$upr_unadj,
          N_pat=nrow(dat),
          N_evt=sum(dat$evt)
        )
      }
    },
    gAHR = {
      function(condition, dat, fixed_objects = NULL){
        model <- trycatch_nphparams(nph::nphparams(dat$t, dat$evt, dat$trt, param_type="HR", param_par=max_time))

        list(
          p = model$tab$p_unadj,
          gAHR = model$tab$Estimate,
          gAHR_lower = model$tab$lwr_unadj,
          gAHR_upper = model$tab$upr_unadj,
          N_pat=nrow(dat),
          N_evt=sum(dat$evt)
        )
      }
    },
    {
      stop(gettext("Invalid type, use AHR for the average hazard ratio or gAHR for the geometric average hazard ratio."))
    }
  )

}
