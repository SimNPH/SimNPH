#' Create Analyse function for the moderately weighted logrank test
#'
#' @return an analyse function that can be used in runSimulation
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
#' analyse_moderately_weighted()(condition, dat)
analyse_moderately_weighted <- function(){
  function(condition, dat, fixed_objects = NULL){
    result_tmp <- list(
      p = 1,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )

    result_tmp
  }
}
