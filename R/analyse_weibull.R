#' Analyse Dataset with Weibull Regression
#'
#' @param level confidence level for CI computation
#'
#' @return an analysis function that returns a data.frame
#'
#' @export
#'
#' @details the columns in the return are the two-sided p-value for the test of
#'   equal medians. The estimated medians in the treatment and control group and
#'   the estimated difference in median survival with confidence intervals.
#'
#'   The estimates and tests are comstructed by fitting seperate Weibull
#'   regression models in the treatment and control groups and then estimating
#'   the medians and respective variances with the delta-method.
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
analyse_weibull <- function(level=0.95) {
  function(condition, dat, fixed_objects = NULL){
    # Calculations provided by Andrew Hooker

    reg_1 <- survreg(Surv(t, evt) ~ 1, dat, dist="weibull", subset = (trt==1))
    reg_0 <- survreg(Surv(t, evt) ~ 1, dat, dist="weibull", subset = (trt==0))


    params_surv_0 <- c(
      "intercept"=reg_0$coefficients[[1]],
      "log_scale"=reg_0$scale
    )
    cov_0 <- reg_0$var
    rownames(cov_0) <- names(params_surv_0)
    colnames(cov_0) <- names(params_surv_0)

    params_surv_1 <- c(
      "intercept"=reg_1$coefficients[[1]],
      "log_scale"=reg_1$scale
    )
    cov_1 <- reg_1$var
    rownames(cov_1) <- names(params_surv_1)
    colnames(cov_1) <- names(params_surv_1)

    med_tte_0 <- car::deltaMethod(
      params_surv_0,
      "exp(intercept)*log(2)^(log_scale)",
      cov_0,
      level=level
    )

    med_tte_1 <- car::deltaMethod(
      params_surv_1,
      "exp(intercept)*log(2)^(log_scale)",
      cov_1,
      level=level
    )

    params_surv_both <- c(params_surv_0,params_surv_1)
    names(params_surv_both) <- c("intercept_0","log_scale_0","intercept_1","log_scale_1")

    cov_both <- diag(4)
    cov_both[1:2,1:2] <- cov_0
    cov_both[3:4,3:4] <- cov_1
    rownames(cov_both) <- names(params_surv_both)
    colnames(cov_both) <- names(params_surv_both)

    diff_med_tte <- car::deltaMethod(
      params_surv_both,
      "exp(intercept_1)*log(2)^(log_scale_1) - exp(intercept_0)*log(2)^(log_scale_0)",
      cov_both,
      level=level
    )

    z_val <- (diff_med_tte$Estimate)/diff_med_tte$SE
    p_val_1sided <- 1-pnorm(z_val)
    p_val_2sided <- 1-pchisq(z_val^2, 1)

    list(
      p = p_val_2sided,
      med_trt_est = med_tte_0[, "Estimate"],
      med_trt_lower = med_tte_0[, 3],
      med_trt_upper = med_tte_0[, 4],
      med_ctrl_est = med_tte_1[, "Estimate"],
      med_ctrl_lower = med_tte_1[, 3],
      med_ctrl_upper = med_tte_1[, 4],
      diff_med_est = diff_med_tte[, "Estimate"],
      diff_med_lower = diff_med_tte[, 3],
      diff_med_upper = diff_med_tte[, 4],
      CI_level = level,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )
  }
}
