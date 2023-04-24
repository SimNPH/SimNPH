#' Create Analyse function for Gehan Wilcoxon test
#'
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return an analyse function that can be used in runSimulation
#'
#' @details
#'
#' `alternative` can be "two.sided" for a two sided test of equality of the
#' summary statistic or "one.sided" for a one sided test testing H0: treatment
#' has equal or shorter survival than control vs. H1 treatment has longer
#' survival than control.
#'
#' @export
#'
#' @examples
#' condition <- merge(
#'   assumptions_delayed_effect(),
#'   design_fixed_followup(),
#'   by = NULL
#' ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_gehan_wilcoxon()(condition, dat)
analyse_gehan_wilcoxon <- function(alternative="two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  function(condition, dat, fixed_objects = NULL) {
    # according to survival:::print.survdiff
    model <- survival::survdiff(survival::Surv(t, evt) ~ trt, dat, rho = 1)
    df <- sum(model$exp > 0) - 1

    p_value <- switch(alternative,
                      two.sided = {
                        pchisq(model$chisq, df, lower.tail = FALSE)
                      },
                      one.sided = {
                        ifelse(
                          (model$obs - model$exp)[2] < 0,
                          pchisq(model$chisq, df, lower.tail = FALSE)/2,
                          1
                        )
                      }
    )

    result_tmp <- list(
      p = p_value,
      alternative=alternative,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )

    result_tmp
  }
}
