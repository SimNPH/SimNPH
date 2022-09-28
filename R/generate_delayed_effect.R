#' Generate Dataset with delayed effect
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condidtion has to contain the following columns:
#'
#'   * n_trt number of paitents in treatment arm
#'   * n_ctrl number of patients in control arm
#'   * delay time until onset of effect
#'   * hazard_ctrl hazard in the control arm = hazard before onset of treatment
#'     effect
#'   * hazard_trt hazard in the treatment arm afert onset of treatment effect
#'
#' If fixed_objects is given and contains an element `t_max`, then this is used
#' as the cutoff for the simulation used internally. If t_max is not given in
#' this way the 1-(1/10000) quantile of the survival distribution in the control
#' or treatment arm is used (which ever is larger).
#'
#' @return
#' For generate_delayed_effect: A dataset with the columns t (time) and trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_delayed_effect simulates a dataset with delayed
#'   treatment effect
#'
#' @examples
#' one_simulation <- desing_skeleton_delayed_effect() |>
#'   head(1) |>
#'   generate_delayed_effect()
#' head(one_simulation)
#' tail(one_simulation)
generate_delayed_effect <- function(condition, fixed_objects=NULL){

  # if t_max is not given in fixed_objects
  if(is.null(fixed_objects) || (!hasName(fixed_objects, t_max))){
    # set t_max to 1-1/10000 quantile of control or treatment survival function
    # whichever is later
    t_max <- max(
      log(10000) / condition$hazard_ctrl,
      log(10000) / condition$hazard_trt
    )
  }

  # simulate treatment group
  if (condition$delay < 0){
  # if delay is smaller than 0 stop with error
    stop(gettext("Delay has to be >= 0"))
  } else if (condition$delay == 0){
  # if delay is 0 leave out period bevore treatment effect
  # (times have to be strictly monotonous for rSurv_fun)
    data_trt <- data.frame(
      t = nph::rSurv_fun(
        condition$n_trt,
        nph::pchaz(
          c(0, t_max),
          c(condition$hazard_trt)
        )
      ),
      trt = 1,
      evt = TRUE
    )
  } else {
  # if delay is positive simulate in the time intervals bevore and after
  # treatment effect
    data_trt <- data.frame(
      t = nph::rSurv_fun(
        condition$n_trt,
        nph::pchaz(
          c(0, condition$delay, t_max),
          c(condition$hazard_ctrl, condition$hazard_trt)
        )
      ),
      trt = 1,
      evt = TRUE
    )
  }

  # simulate control group with constant hazard from 0 to t_max
  data_ctrl <- data.frame(
    t = nph::rSurv_fun(
      condition$n_ctrl,
      nph::pchaz(
        c(0, t_max),
        c(condition$hazard_ctrl)
      )
    ),
    trt = 0,
    evt = TRUE
  )

  rbind(data_trt, data_ctrl)
}

#' Create design skeleton for use with generate_delayed_effect
#'
#' @return For desing_skeleton_delayed_effect: a design tibble with default values invisibly
#'
#' @details desing_skeleton_delayed_effect prints the code to generate a default
#'   design tibble for use with generate_delayed_effect and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_delayed_effect generate default design tibble
#'
#' @examples
#' Design <- desing_skeleton_delayed_effect()
#' Design
desing_skeleton_delayed_effect <- function(){
  skel <- "createDesign(
  n_trt=50,               # 100 patients in the treatment arm
  n_ctrl=50,              # 100 patients in the control arm
  delay=seq(0, 10, by=2), # delay of 0, 1, ..., 10 days
  hazard_ctrl=0.2,        # hazard under control and before treatment effect
  hazard_trt=0.02         # hazard after onset of treatment effect
)
"

  cat(skel)
  invisible(
    skel |>
      str2expression() |>
      eval()
  )
}

#' Calculate hr after onset of treatment effect from gAHR
#'
#' @param condition condition data.frame
#' @param target_gAHR target geometric average hazard ratio
#'
#' @return For prepare_design_delayed_effect_gAHR: the condition
#'   data.frame passed as argument with the additional column hazard_trt.
#' @export
#'
#' @describeIn generate_delayed_effect  Calculate hr after onset of treatment effect from gAHR
#'
#' @examples
#' my_condition <- desing_skeleton_delayed_effect()
#' my_condition$hr_trt <- NA
#' my_condition <- prepare_design_delayed_effect_gAHR(my_condition, 0.8)
#' my_condition
prepare_design_delayed_effect_gAHR <- function(condition, target_gAHR){

}
