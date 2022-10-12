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
#' * `confint` 95% confidence intervall for the hazard ratio for `trt`
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @examples
#' condition <- desing_skeleton_delayed_effect() |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_coxph(condition, dat)
analyse_coxph <- function(condition, dat, fixed_objects = NULL){
  model <- survival::coxph(survival::Surv(t, evt) ~ trt, dat, robust=FALSE)

  list(
    p = 1-pchisq(model$score, 1),
    coef = coefficients(model)["trt"],
    hr   = exp(coefficients(model)["trt"]),
    confint = summary(model)$conf.int[, c("lower .95", "upper .95")],
    N_pat=nrow(dat),
    N_evt=sum(dat$evt)
  )
}

