#' Analyse Dataset with the Cox Protportional Hazards Model
#'
#' @param condition condition of the simulation
#' @param dat generated datasets
#' @param fixed_objects other constants
#'
#' @return a data frame with the columns
#' * `p` p value of the score test
#' * `coef` coefficient for `trt`
#' * `hr` hazard ratio for `trt`
#' * `hr_lower` lower 95% confidence intervall boundary for the hazard ratio for `trt`
#' * `hr_upper`lower 95% confidence intervall boundary for the hazard ratio for `trt`
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
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_coxph(condition, dat)
analyse_coxph <- function(condition, dat, fixed_objects = NULL){
  model <- survival::coxph(survival::Surv(t, evt) ~ trt, dat, robust=FALSE)
  model_summary <- summary(model)

  list(
    p = 1-pchisq(model$score, 1),
    coef = coefficients(model)["trt"],
    hr   = exp(coefficients(model)["trt"]),
    hr_lower = model_summary$conf.int[, "lower .95"],
    hr_upper = model_summary$conf.int[, "upper .95"],
    N_pat=nrow(dat),
    N_evt=sum(dat$evt)
  )
}

