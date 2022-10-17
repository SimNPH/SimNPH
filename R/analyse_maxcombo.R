#' Analyse Dataset with the Maxcombo Test
#'
#' @param condition condition of the simulation
#' @param dat generated datasets
#' @param fixed_objects other constants
#'
#' @return a dataframe with the combined p-value of the max combo test in the
#'   column p
#' @export
#'
#' @examples
#' condition <- desing_skeleton_delayed_effect() |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_maxcombo(condition, dat)
analyse_maxcombo <- function(condition, dat, fixed_objects = NULL){
  list(
    p=nph::logrank.maxtest(dat$t, dat$evt, dat$trt)$pmult,
    N_pat=nrow(dat),
    N_evt=sum(dat$evt)
  )
}
