#' Analyse Dataset with the Cox Protportional Hazards Model
#'
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return an analyse function that returns a list with the elements
#' * `p` p value of the score test (two.sided) or the Wald test (one.sided)
#' * `alternative` the alternative used
#' * `coef` coefficient for `trt`
#' * `hr` hazard ratio for `trt`
#' * `hr_lower` lower 95% confidence intervall boundary for the hazard ratio for `trt`
#' * `hr_upper`lower 95% confidence intervall boundary for the hazard ratio for `trt`
#' * `CI_level` the CI level used
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @details
#'
#' `alternative` can be "two.sided" for a two sided test of equality of the
#' summary statistic or "one.sided" for a one sided test testing H0: treatment
#' has equal or shorter survival than control vs. H1 treatment has longer
#' survival than control.
#'
#' @examples
#' condition <- merge(
#'   assumptions_delayed_effect(),
#'   design_fixed_followup(),
#'   by = NULL
#' ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_coxph()(condition, dat)
analyse_coxph <- function(level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  function(condition, dat, fixed_objects = NULL) {
    model <- survival::coxph(survival::Surv(t, evt) ~ trt, dat, robust = FALSE)
    model_summary <- summary(model, conf.int = level)

    p_value <- switch(alternative,
      two.sided = {
        1 - pchisq(model$score, 1)
      },
      one.sided = {
        1 - pnorm(model$wald.test)
      }
    )

    list(
      p = p_value,
      alternative = alternative,
      coef = coefficients(model)["trt"],
      hr = exp(coefficients(model)["trt"]),
      hr_lower = model_summary$conf.int[, 3],
      hr_upper = model_summary$conf.int[, 4],
      CI_level = level,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )
  }
}
