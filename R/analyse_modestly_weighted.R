#' Create Analyse function for the modestly weighted logrank test
#'
#' @param t_star parameter t* of the modestly weighted logrank test
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
#' analyse_modelstly_weighted(20)(condition, dat)
analyse_modelstly_weighted <- function(t_star){

  if(length(t_star) != 1){
    stop("length(t_star) != 1, please provide a scalar cutoff time.")
  }

  function(condition, dat, fixed_objects = NULL){
    z_score <- nphRCT::wlrt(Surv(t, evt)~trt, dat, t_star=t_star, method="mw")$z
    p_value <- (pnorm(-z_score, lower.tail = FALSE))

    result_tmp <- list(
      p = p_value,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )

    result_tmp
  }
}
