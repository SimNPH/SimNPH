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

    # if no patients enter any interval, return NA for all intervals
    # TODO: calculate values for intervals where it's possible
    if(length(unique(dat2$tstart)) < n_intervals){

      na_tibble <- rep(NA_real_, n_intervals) |>
        setNames(coef_labels) |>
        t() |>
        tibble::as_tibble()

      list(
        hr_s = na_tibble,
        p_s =  na_tibble,
        hr_s_lower = na_tibble,
        hr_s_upper = na_tibble,
        N_pat = nrow(dat),
        N_evt = sum(dat$evt),
        interval_table = table(dat2$tstart, dat2$evt)
      )
    } else {
      dat2$interval <- factor(dat2$tstart, labels=interval_labels)
      dat2$t_interval = dat2$t-dat2$tstart

      model <- glm(evt ~ trt*interval-trt-1+offset(log(t_interval)), data=dat2, family=poisson())

      summary <- summary(model)

      # preparing output

      lower <- confint(model)[, 1]
      lower <- lower[coef_labels]
      names(lower) <- coef_labels

      upper <- confint(model)[, 2]
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
        interval_table = list(list(table(dat2$tstart, dat2$evt)))
      )
    }
  }


}
