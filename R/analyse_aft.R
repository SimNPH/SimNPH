#' Analyse Dataset with accelarated failure time models
#'
#' @param level confidence level for CI computation
#' @param dist passed to survival::survreg
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return an analyse function that returns a list with the elements
#' * `p` p value of the score test (two.sided) or the Wald test (one.sided)
#' * `alternative` the alternative used
#' * `coef` coefficient for `trt`
#' * `lower` lower 95% confidence intervall boundary for the coefficient
#' * `upper`lower 95% confidence intervall boundary for the coefficient
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
#' analyse_aft()(condition, dat)
#' analyse_aft(dist="lognormal")(condition, dat)
analyse_aft <- function(level = 0.95, dist="weibull", alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  function(condition, dat, fixed_objects = NULL) {
    model <- survival::survreg(survival::Surv(t, evt) ~ trt, dat, dist=dist, robust = FALSE)
    model_summary <- summary(model, conf.int = level)
    model_confint <- confint(model, level = level)

    p_value <- switch(alternative,
      two.sided = {
        1 - pchisq((-model_summary$table["trt", "z"])^2, 1)
      },
      one.sided = {
        1 - pnorm(model_summary$table["trt", "z"])
      }
    )

    list(
      p = p_value,
      alternative = alternative,
      coef = coefficients(model)["trt"],
      lower = model_confint["trt", 1],
      upper = model_confint["trt", 2],
      CI_level = level,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )
  }
}
