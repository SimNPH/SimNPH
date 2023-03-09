#' Analyse Dataset with the Logrank Test
#'
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return an analysis function that returns a data.frame with the columns
#' * `p` p-value of the logrank test
#' * `alternative` the alternative used
#' * `N_pat` number of patients
#' * `N_evt` number of events
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
#' analyse_logrank()(condition, dat)
analyse_logrank <- function(alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  alt_ <- switch(alternative,
    two.sided = "two.sided",
    one.sided = "greater",
    stop(gettext("'alternative' has to be either 'two.sided' or 'one.sided'."))
  )

  function(condition, dat, fixed_objects = NULL) {
    list(
      p = nph::logrank.test(dat$t, dat$evt, dat$trt, alternative = alt_)$test$p,
      alternative = alternative,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )
  }
}
