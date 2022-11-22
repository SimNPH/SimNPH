#' Create Analyse function for Gehan Wilcoxon test
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
#' analyse_gehan_wilcoxon()(condition, dat)
analyse_gehan_wilcoxon <- function(){
  function(condition, dat, fixed_objects = NULL){
    result_tmp <- list(
      p = 1,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )

    result_tmp
  }
}
