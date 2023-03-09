#' Analyse Dataset with the Fleming Harrington weighted Logrank Test
#'
#' @param rho rho for the rho-gamma family of weights
#' @param gamma gamma for the rho-gamma family of weights
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return a function with the arguments condition, dat and fixed_objects that
#'   returns a dataframe with the p-value of the weighted logrank test in the
#'   column p. See ?SimDesign::Analyse for details on the arguments condition,
#'   dat, fixed_arguments.
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
#' # create two functions with different weights
#' analyse_01 <- analyse_logrank_fh_weights(rho = 0, gamma = 1)
#' analyse_10 <- analyse_logrank_fh_weights(rho = 1, gamma = 0)
#' # run the tests created before
#' analyse_01(condition, dat)
#' analyse_10(condition, dat)
analyse_logrank_fh_weights <- function(rho, gamma, alternative="two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  alt_ <- switch(alternative,
    two.sided = "two.sided",
    one.sided = "greater",
    stop(gettext("'alternative' has to be either 'two.sided' or 'one.sided'."))
  )

  function(condition, dat, fixed_objects = NULL) {
    list(
      p = nph::logrank.test(dat$t, dat$evt, dat$trt, rho = rho, gamma = gamma, alternative = alt_)$test$p,
      alternative = alternative,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )
  }
}
