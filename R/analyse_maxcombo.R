#' Analyse Dataset with the Maxcombo Test
#'
#' @return an analyse function that returns a data.frame with the combined
#'   p-value of the max combo test in the column p
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
#' analyse_maxcombo()(condition, dat)
analyse_maxcombo <- function(){
  function(condition, dat, fixed_objects = NULL){
    list(
      p=nph::logrank.maxtest(dat$t, dat$evt, dat$trt)$pmult,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )
  }
}
