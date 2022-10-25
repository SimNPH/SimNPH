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
#' one_simulation <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1) |>
#'   generate_delayed_effect()
#' head(one_simulation)
#' tail(one_simulation)
generate_delayed_effect <- function(condition, fixed_objects=NULL){

  # if t_max is not given in fixed_objects
  if(is.null(fixed_objects) || (!hasName(fixed_objects, "t_max"))){
    # set t_max to 1-1/10000 quantile of control or treatment survival function
    # whichever is later
    t_max <- max(
      log(10000) / condition$hazard_ctrl,
      log(10000) / condition$hazard_trt
    )
  } else {
    t_max <- fixed_objects$t_max
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

#' Create an empty assumtions data.frame for generate_delayed_effect
#'
#' @return For assumptions_delayed_effect: a design tibble with default values invisibly
#'
#' @details assumptions_delayed_effect prints the code to generate a default
#'   design tibble for use with generate_delayed_effect and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_delayed_effect generate default design tibble
#'
#' @examples
#' Design <- assumptions_delayed_effect()
#' Design
assumptions_delayed_effect <- function(){
  skel <- "expand.grid(
  delay=seq(0, 10, by=2), # delay of 0, 1, ..., 10 days
  hazard_ctrl=0.2,        # hazard under control and before treatment effect
  hazard_trt=0.02,        # hazard after onset of treatment effect
  random_withdrawal=0.01  # rate of random withdrawal
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
#' @param design design data.frame
#' @param target_gAHR target geometric average hazard ratio
#'
#' @return For hr_after_onset_from_gAHR: the design data.frame passed as
#'   argument with the additional column hazard_trt.
#' @export
#'
#' @describeIn generate_delayed_effect  Calculate hr after onset of treatment effect from gAHR
#'
#' @examples
#' my_design <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design$hr_trt <- NA
#' my_design <- hr_after_onset_from_gAHR(my_design, 0.8)
#' my_design
hr_after_onset_from_gAHR <- function(design, target_gAHR){

}

#' Calculate true summary statistics for scenarios with delayed treatment effect
#'
#' @param Design Design data.frame for delayed effect
#' @param cutoff_stats Cutoff time for rmst and average hazard ratios
#' @param fixed_objects=NULL fixed objects not used for now
#'
#' @return For true_summary_statistics_delayed_effect: the design data.frame
#'   passed as argument with the additional columns:
#' * `rmst_trt` rmst in the treatment group
#' * `median_surv_trt` median survival in the treatment group
#' * `rmst_ctrl` rmst in the control group
#' * `median_surv_ctrl` median survial in the control group
#' * `gAHR` geometric average hazard ratio
#' * `AHR` average hazard ratio
#'
#' @export
#'
#' @describeIn generate_delayed_effect  calculate true summary statistics for delayed effect
#'
#' @examples
#' my_design <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design$follwup <- 15
#' my_design <- true_summary_statistics_delayed_effect(my_design, cutoff_stats=my_design$followup)
#' my_design
true_summary_statistics_delayed_effect <- function(Design, cutoff_stats=10, fixed_objects=NULL){

  true_summary_statistics_delayed_effect_rowwise <- function(condition, cutoff_stats){

    # if t_max is not given in fixed_objects
    if(is.null(fixed_objects) || (!hasName(fixed_objects, "t_max"))){
      # set t_max to 1-1/10000 quantile of control or treatment survival function
      # whichever is later
      t_max <- max(
        log(10000) / condition$hazard_ctrl,
        log(10000) / condition$hazard_trt
      )
    } else {
      t_max <- fixed_objects$t_max
    }

    # create functions for treatment group
    if (condition$delay < 0){
      # if delay is smaller than 0 stop with error
      stop(gettext("Delay has to be >= 0"))
    } else if (condition$delay == 0){
      # if delay is 0 leave out period bevore treatment effect
      # (times have to be strictly monotonous for rSurv_fun)
      data_generating_model_trt <- nph::pchaz(
        c(0, t_max),
        c(condition$hazard_trt)
      )

    } else {
      # if delay is positive create piecewise constant hazards and respective
      # functions in the time intervals bevore and after treatment effect
      data_generating_model_trt <- nph::pchaz(
        c(0, condition$delay, t_max),
        c(condition$hazard_ctrl, condition$hazard_trt)
      )

    }

    # create functions for control group with constant hazard from 0 to t_max
    data_generating_model_ctrl <- nph::pchaz(
      c(0, t_max),
      c(condition$hazard_ctrl)
    )

    res <- cbind(
      condition,
      internal_real_statistics_pchaz(
        data_generating_model_trt,
        data_generating_model_ctrl,
        N_trt=condition$n_trt,
        N_ctrl=condition$n_ctrl,
        cutoff = cutoff_stats
      )
    )

    row.names(res) <- NULL
    res
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    mapply(FUN=true_summary_statistics_delayed_effect_rowwise, cutoff_stats = cutoff_stats, SIMPLIFY = FALSE)

  Design <- do.call(rbind, Design)

  Design
}
