#' Analyse Dataset with weighted Cox regression
#'
#' @param type="SG" type of weights, see Details
#' @param max_time=NA_real_ cutoff time for estimation, defaults to last event time
#'
#' @return an analyse function that returns a data.frame with the columns
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
#' @details
#'
#' Type can be "SG", "S" or "G". For "SG" the weights are set to
#' `S(t)*G(t)^(-1)`, for "S" `S(t)` is used, for "G", `G(t)^(-1)` is used. Here
#' `S(t)` denotes the Kaplan Meier estimator for the event Times and `G(t)`
#' denotes the Kaplan Meier estimator for the censoring times.
#'
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_coxph_weighted()(condition, dat)
analyse_coxph_weighted <- function(type="SG", max_time=NA_real_){
  if(!(type %in% c("SG", "S", "G"))){
    stop('in analyse_coxph_weighted: invalid "type", currently implemented: "SG" and "G".')
  }

  function(condition, dat, fixed_objects = NULL){
    # TODO: check this
    # censor at cutoff time
    if(!is.na(max_time)){
      dat$evt[dat$t > max_time] <- FALSE
      dat$t <- pmin(dat$t, max_time)
    }

    # set t_max to last non-censored observation to avoid infinite weights
    t_max <- max(dat$t[dat$evt])
    dat <- dat[dat$t <= t_max, ]

    switch(
      type,
      S = {
        S <- survfit(Surv(t, evt)~1, data=dat)
        weights <- S$surv[match(dat$t, S$time)]
      },
      SG = {
        S <- survfit(Surv(t, evt)~1, data=dat)
        weights <- S$surv[match(dat$t, S$time)]
        W <- survfit(Surv(t, !evt)~1, data=dat)
        weights <- weights * ((W$surv[match(dat$t, W$time)])^(-1))
      },
      G={
        W <- survfit(Surv(t, !evt)~1, data=dat)
        weights <- ((W$surv[match(dat$t, W$time)])^(-1))
      }
    )

    model <- survival::coxph(survival::Surv(t, evt) ~ trt, dat, robust=TRUE, weights=weights, subset=(weights!=0))
    model_summary <- summary(model)

    list(
      p = 1-pchisq(model$score, 1),
      coef = coefficients(model)["trt"],
      hr   = exp(coefficients(model)["trt"]),
      hr_lower = model_summary$conf.int[, "lower .95"],
      hr_upper = model_summary$conf.int[, "upper .95"],
      N_pat=nrow(dat),
      N_evt=sum(dat$evt),
      followup=t_max
    )
  }
}
