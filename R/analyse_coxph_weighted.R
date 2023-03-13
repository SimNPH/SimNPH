#' Analyse Dataset with weighted Cox regression
#'
#' @param type type of weights, see Details
#' @param max_time cutoff time for estimation, defaults to last event time
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return an analyse function that returns a data.frame with the columns
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
#' Type can be "SG", "S" or "G". For "SG" the weights are set to
#' `S(t)*G(t)^(-1)`, for "S" `S(t)` is used, for "G", `G(t)^(-1)` is used. Here
#' `S(t)` denotes the Kaplan Meier estimator for the event Times and `G(t)`
#' denotes the Kaplan Meier estimator for the censoring times.
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
#' analyse_coxph_weighted()(condition, dat)
analyse_coxph_weighted <- function(type = "SG", max_time = NA_real_, level = 0.95, alternative = "two.sided") {
  if (!(type %in% c("SG", "S", "G"))) {
    stop('in analyse_coxph_weighted: invalid "type", currently implemented: "SG" and "G".')
  }

  stopifnot(alternative %in% c("two.sided", "one.sided"))

  function(condition, dat, fixed_objects = NULL) {
    # TODO: check this
    # censor at cutoff time
    if (!is.na(max_time)) {
      dat$evt[dat$t > max_time] <- FALSE
      dat$t <- pmin(dat$t, max_time)
    }

    # set t_max to last non-censored observation to avoid infinite weights
    t_max <- max(dat$t[dat$evt])
    dat <- dat[dat$t <= t_max, ]

    switch(type,
      S = {
        S <- survfit(Surv(t, evt) ~ 1, data = dat)
        weights <- S$surv[match(dat$t, S$time)]
      },
      SG = {
        S <- survfit(Surv(t, evt) ~ 1, data = dat)
        weights <- S$surv[match(dat$t, S$time)]
        W <- survfit(Surv(t, !evt) ~ 1, data = dat)
        weights <- weights * ((W$surv[match(dat$t, W$time)])^(-1))
      },
      G = {
        W <- survfit(Surv(t, !evt) ~ 1, data = dat)
        weights <- ((W$surv[match(dat$t, W$time)])^(-1))
      }
    )

    model <- survival::coxph(survival::Surv(t, evt) ~ trt, dat, robust = TRUE, weights = weights, subset = (weights != 0))
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
      N_evt = sum(dat$evt),
      followup = t_max
    )
  }
}
