#' Analyse Dataset with Weibull Regression
#'
#' @return an analysis function that returns a data.frame with the columns
#' * `est_med_trt` estimated median survival in the treatment arm
#' * `est_med_ctrl` estimated median survival in the control arm
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(3) |>
#'   tail(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_weibull()(condition, dat)
analyse_weibull <- function() {
  function(condition, dat, fixed_objects = NULL){
    model_trt <- survreg(Surv(t, evt) ~ 1, dat, dist="weibull", subset = (trt==1))
    model_ctrl <- survreg(Surv(t, evt) ~ 1, dat, dist="weibull", subset = (trt==0))

    med_trt <- exp(model_trt$coefficients["(Intercept)"]) * (log(2)^(model_trt$scale)) - 1
    med_ctrl <- exp(model_ctrl$coefficients["(Intercept)"]) * (log(2)^(model_ctrl$scale)) - 1

    list(
      est_med_trt = med_trt,
      est_med_ctrl = med_ctrl,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )
  }
}
