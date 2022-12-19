#' Create Analyse function for piecewise exponential model
#'
#' @param cuts interval boundaries for the piecewise exponential model
#'
#' @return an analyse function that can be used in runSimulation
#'
#' @details If there's any time interval no patients ever enter, NA is returned
#'   for all time intervals. This behavior will likely change in future package
#'   versions.
#'
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
#' analyse_piecewise_exponential(cuts=c(90, 360))(condition, dat)
analyse_piecewise_exponential <- function(cuts){

  n_intervals <- length(cuts) + 1
  interval_labels <- paste0("I", 1:n_intervals)
  coef_labels <- paste0("trt:interval", interval_labels)

  function(condition, dat, fixed_objects = NULL){
    dat2 <- survSplit(Surv(t, evt)~trt, data=dat, cut=cuts)
    dat2$interval <- factor(dat2$tstart, levels=c(0, cuts), labels=interval_labels)
    dat2$t_interval = dat2$t-dat2$tstart

    # exclude intervals with no events
    intervals_with_events <- tapply(dat2$evt, dat2$interval, sum, default=0) |>
      Filter(f=\(x){x>0}) |>
      names()

    dat3 <- dat2 |>
      subset(interval %in% intervals_with_events)

    model <- glm(evt ~ trt*interval-trt-1+offset(log(t_interval)), data=dat3, family=poisson())

    summary <- summary(model)

    # preparing output
    # coef_labels also included the ones with no events
    # since those are not in the model summary or confint, they are NA in the
    # output

    tmp_confint <- confint(model)

    lower <- tmp_confint[, 1]
    lower <- lower[coef_labels]
    names(lower) <- coef_labels

    upper <- tmp_confint[, 2]
    upper <- upper[coef_labels]
    names(upper) <- coef_labels

    hr_s <- summary$coefficients[, "Estimate"]
    hr_s <- hr_s[coef_labels]
    names(hr_s) <- coef_labels
    hr_s <- exp(hr_s)

    p_s <- summary$coefficients[, "Pr(>|z|)"]
    p_s <- p_s[coef_labels]
    names(p_s) <- coef_labels


    list(
      hr_s = tibble::as_tibble(t(hr_s)),
      p_s = tibble::as_tibble(t(p_s)),
      hr_s_lower = tibble::as_tibble(t(lower)),
      hr_s_upper = tibble::as_tibble(t(upper)),
      N_pat = nrow(dat),
      N_evt = sum(dat$evt),
      interval_table = table(dat2$interval, dat2$evt)
    )
  }


}
