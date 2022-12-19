#' Create Analyse function for piecewise exponential model
#'
#' @param cuts interval boundaries for the piecewise exponential model
#' @param testing_only=FALSE if set to true omits all statistics in the
#'   intervals and just returns the p value of the global test.
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
analyse_piecewise_exponential <- function(cuts, testing_only=FALSE){

  n_intervals <- length(cuts) + 1
  interval_labels <- paste0("I", 1:n_intervals)
  coef_labels <- paste0("trt:interval", interval_labels)

  if(testing_only){

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

      if( length(unique(dat3$interval)) > 1 ) {
        model <- glm(evt ~ trt*interval-trt-1+offset(log(t_interval)), data=dat3, family=poisson())
        null_model <- glm(evt ~ interval-1+offset(log(t_interval)), data=dat3, family=poisson())

      } else {
        model <- glm(evt ~ trt-1+offset(log(t_interval)), data=dat3, family=poisson())
        null_model <- glm(evt ~ -1+offset(log(t_interval)), data=dat3, family=poisson())
      }

      overall_p_value <- anova(null_model, model, test="Chisq")[["Pr(>Chi)"]][2]

      list(
        p = overall_p_value
      )
    }

  } else {

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

      if( length(unique(dat3$interval)) > 1 ) {
        model <- glm(evt ~ trt*interval-trt-1+offset(log(t_interval)), data=dat3, family=poisson())
        null_model <- glm(evt ~ interval-1+offset(log(t_interval)), data=dat3, family=poisson())
        summary <- summary(model)

        tmp_confint <- tryCatch(
          confint(model),
          error=function(e){
            matrix(NA_real_, 0, 2)
          })

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
      } else {
        only_interval_label <-  paste0("trt:interval", unique(dat3$interval))

        model <- glm(evt ~ trt-1+offset(log(t_interval)), data=dat3, family=poisson())
        null_model <- glm(evt ~ -1+offset(log(t_interval)), data=dat3, family=poisson())
        summary <- summary(model)

        tmp_confint <- tryCatch(
          confint(model),
          error=function(e){
            numeric(0)
          })

        lower <- tmp_confint[1]
        names(lower) <- only_interval_label
        lower <- lower[coef_labels]
        names(lower) <- coef_labels

        upper <- tmp_confint[2]
        names(upper) <- only_interval_label
        upper <- upper[coef_labels]
        names(upper) <- coef_labels

        hr_s <- summary$coefficients[, "Estimate"]
        names(hr_s) <- only_interval_label
        hr_s <- hr_s[coef_labels]
        names(hr_s) <- coef_labels
        hr_s <- exp(hr_s)

        p_s <- summary$coefficients[, "Pr(>|z|)"]
        names(p_s) <- only_interval_label
        p_s <- p_s[coef_labels]
        names(p_s) <- coef_labels
      }

      overall_p_value <- anova(null_model, model, test="Chisq")[["Pr(>Chi)"]][2]

      # preparing output
      # coef_labels also included the ones with no events
      # since those are not in the model summary or confint, they are NA in the
      # output



      list(
        p = overall_p_value,
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


}
