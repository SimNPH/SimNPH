#' Create Analyse function for the modestly weighted logrank test
#'
#' @param t_star parameter t* of the modestly weighted logrank test
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @examples
#' condition <- merge(
#'    assumptions_delayed_effect(),
#'    design_fixed_followup(),
#'    by=NULL
#'  ) |>
#'  head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_modelstly_weighted(20)(condition, dat)
analyse_modelstly_weighted <- function(t_star){
  function(condition, dat, fixed_objects = NULL){
    # estimated survival functions
    model_km <- survival::survfit(survival::Surv(t, evt)~1, dat)
    s_t_star <- summary(model_km, times=t_star)$surv
    s_t <- summary(model_km, times=sort(unique(dat$t)))$surv

    # event time weights
    w_t <- 1/pmax(s_t_star, s_t)

    result_tmp <- list(
      p = nph::logrank.test(dat$t, dat$evt, dat$trt, event_time_weights = w_t)$test$p,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )

    result_tmp
  }
}
