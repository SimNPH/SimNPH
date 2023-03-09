#' Create Analyse function for Gehan Wilcoxon test
#'
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return an analyse function that can be used in runSimulation
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
analyse_gehan_wilcoxon <- function(alternative) {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  function(condition, dat, fixed_objects = NULL) {
    # according to survival:::print.survdiff
    model <- survival::survdiff(survival::Surv(t, evt) ~ trt, dat, rho = 1)
    df <- sum(model$exp > 0) - 1

    p_value <- switch(alternative,
                      two.sided = {
                        1 - pchisq(model$score, 1)
                      },
                      one.sided = {
                        1 - pnorm(model$wald.test)
                      }
    )

    result_tmp <- list(
      p = pchisq(model$chisq, df, lower.tail = FALSE),
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )

    result_tmp
  }
}
