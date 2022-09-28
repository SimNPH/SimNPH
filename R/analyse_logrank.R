#' Analyse Dataset with the Logrank Test
#'
#' @param condition condition of the simulation
#' @param dat generated datasets
#' @param fixed_objects other constants
#'
#' @return a dataframe with the p-value of the logrank test in the column p
#' @export
#'
#' @examples
#' condition <- desing_skeleton_delayed_effect() |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_logrank(condition, dat)
analyse_logrank <- function(condition, dat, fixed_objects = NULL){
  list(
    p=nph::logrank.test(dat$t, dat$evt, dat$trt)$test$p,
    N_pat=nrow(dat),
    N_evt=sum(dat$evt)
  )
}
