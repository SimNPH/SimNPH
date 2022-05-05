#' Generate Dataset with mixture of cured patients
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condidtion has to contain the following columns:
#'
#'   * n_trt number of paitents in treatment arm
#'   * n_ctrl number of patients in control arm
#'   * hazard_ctrl hazard in the control arm
#'   * hazard_trt hazard in the treatment arm for not cured patients
#'   * hazard_cured hazard in for cured patients
#'   * cured_prop proportion of cured patients
#'   * t_max the maximal (cutoff) time
#'
#' @return
#' For generate_cure_mixture: A dataset with the columns t (time) and trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_cure_mixture simulates a dataset with a mixture of cured
#'   patients
#'
#' @examples
#' one_simulation <- desing_skeleton_cure_mixture() |>
#'   head(1) |>
#'   generate_cure_mixture()
#' head(one_simulation)
#' tail(one_simulation)
generate_cure_mixture <- function(condition, fixed_objects=NULL){
  if (condition$cured_prop < 0 || condition$cured_prop > 1) {
    stop(gettext("Cured proportion has to be between 0 and 1"))
  }

  data_trt <- data.frame(
    t = nph::rSurv_fun(
      condition$n_trt,
      nph::pop_pchaz(
        Tint = c(0, condition$t_max),
        lambdaMat1 = matrix(c(
          condition$hazard_trt,
          condition$hazard_cured
        ), ncol=1),
        lambdaMat2 = matrix(c(0, 0), ncol=1),
        lambdaProgMat = matrix(c(0, 0), ncol=1),
        p=c(1-condition$cured_prop, condition$cured_prop)
      )
    ),
    trt = 1,
    evt = TRUE
  )

  data_ctrl <- data.frame(
    t = nph::rSurv_fun(
      condition$n_ctrl,
      nph::pchaz(
        c(0, condition$t_max),
        c(condition$hazard_ctrl)
      )
    ),
    trt = 0,
    evt = TRUE
  )

  rbind(data_trt, data_ctrl)
}

#' Create design skeleton for use with generate_cure_mixture
#'
#' @return For desing_skeleton_cure_mixture: a design tibble with default values invisibly
#'
#' @details desing_skeleton_cure_mixture prints the code to generate a default
#'   design tibble for use with generate_cure_mixture and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_cure_mixture generate default design tibble
#'
#' @examples
#' Design <- desing_skeleton_cure_mixture()
#' Design
desing_skeleton_cure_mixture <- function(){
  skel <- "createDesign(
  n_trt=50,                         # 100 patients in the treatment arm
  n_ctrl=50,                        # 100 patients in the control arm
  hazard_ctrl  = 0.003795467,       # hazard under control (med. survi. 6m)
  hazard_trt   = 0.001265156,       # hazard under treatment (med. surv. 18m)
  hazard_cured = 9.488668e-05,      # hazard cured (med. surv. 20y)
  cured_prop   = seq(0, 1, by=0.2), # proportion of cured patients
  t_max=5*365.25                    # cutoff time
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}
