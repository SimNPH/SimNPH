#' Analyse Dataset with the Maxcombo Test
#'
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
#'
#' @return an analyse function that returns a data.frame with the combined
#'   p-value of the max combo test in the column p
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
#' analyse_maxcombo()(condition, dat)
analyse_maxcombo <- function(alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  alt_ <- switch(alternative,
    two.sided = "two.sided",
    one.sided = "greater",
    stop(gettext("'alternative' has to be either 'two.sided' or 'one.sided'."))
  )

  function(condition, dat, fixed_objects = NULL) {
    list(
      p = nph::logrank.maxtest(dat$t, dat$evt, dat$trt, alternative = alt_)$pmult,
      N_pat = nrow(dat),
      N_evt = sum(dat$evt)
    )
  }
}
