#' Analyse Dataset with the Fleming Harrington weighted Logrank Test
#'
#' @param rho rho for the rho-gamma family of weights
#' @param gamma gamma for the rho-gamma family of weights
#'
#' @return a function with the arguments condition, dat and fixed_objects that
#'   returns a dataframe with the p-value of the weighted logrank test in the
#'   column p. See ?SimDesign::Analyse for details on the arguments condition,
#'   dat, fixed_arguments.
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
#' # create two functions with different weights
#' analyse_01 <- analyse_logrank_fh_weights(rho=0, gamma=1)
#' analyse_10 <- analyse_logrank_fh_weights(rho=1, gamma=0)
#' # run the tests created before
#' analyse_01(condition, dat)
#' analyse_10(condition, dat)
analyse_logrank_fh_weights <- function(rho, gamma){
  function(condition, dat, fixed_objects = NULL){
    data.frame(
      p=nph::logrank.test(dat$t, dat$evt, dat$trt, rho=rho, gamma=gamma)$test$p
    )
  }}
