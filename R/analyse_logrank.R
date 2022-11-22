#' Analyse Dataset with the Logrank Test
#'
#' @return an analysis function that returns a data.frame with the columns
#' * `p` p-value of the logrank test
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_logrank()(condition, dat)
analyse_logrank <- function() {
  function(condition, dat, fixed_objects = NULL){
    list(
      p=nph::logrank.test(dat$t, dat$evt, dat$trt)$test$p,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )
  }
}
