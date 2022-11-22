#' Analyse Dataset with an ANCOVA Type test on the RMST
#'
#' @return an analysis function that returns a data.frame with the p-value of
#'   the test in the column p
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
#' analyse_rmst()(condition, dat)
analyse_rmst <- function(){
  function(condition, dat, fixed_objects = NULL){
    data.frame(
      p=survRM2::rmst2(dat$t, dat$evt, dat$trt)$unadjusted.result[1, "p"]
    )
  }
}
