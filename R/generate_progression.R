#' Create an empty assumtions data.frame for generate_progression
#'
#' @return For generate_progression: a design tibble with default values invisibly
#'
#' @details assumptions_progression prints the code to generate a default
#'   design tibble for use with generate_progression and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_progression generate default assumptions tibble
#'
#' @examples
#' Design <- assumptions_progression()
#' Design
assumptions_progression <- function(){
  skel <- "expand.grid(
  hazard_ctrl         = 0.001518187, # hazard under control (med. survi. 15m)
  hazard_trt          = 0.001265156, # hazard under treatment (med. surv. 18m)
  hazard_after_prog   = 0.007590934, # hazard after progression (med. surv. 3m)
  prog_rate_ctrl      = 0.001897734, # hazard rate for disease progression under control (med. time to progression 12m)
  prog_rate_trt       = c(0.001897734, 0.001423300, 0.001265156), # hazard rate for disease progression unter treatment (med. time to progression 12m, 16m, 18m)
  random_withdrawal   = 0.01         # rate of random withdrawal
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}

#' Generate Dataset with changing hazards after disease progression
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
#'   * hazard_after_prog hazard after disease progression
#'   * prog_rate_ctrl hazard rate for disease progression unter control
#'   * prog_rate_trt hazard rate for disease progression unter treatment
#'
#' @return
#' For generate_progression: A dataset with the columns t (time) and trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations),
#' t_ice (time of intercurrent event), ice (intercurrent event)
#'
#' @export
#' @describeIn generate_progression simulates a dataset with changing hazards after disease progression
#'
#' @examples
#' one_simulation <- merge(
#'     assumptions_progression(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   tail(1) |>
#'   generate_progression()
#' head(one_simulation)
#' tail(one_simulation)

# TODO: update t_max calculation such that lower hazard after intercurrent event (cure) gives reasonable t_max
generate_progression <- function(condition, fixed_objects=NULL){
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

  t_evt_ctrl <- nph::rSurv_fun(
      condition$n_ctrl,
      nph::pchaz(
        c(0, t_max),
        c(condition$hazard_ctrl)
      )
    )

  t_evt_trt <- nph::rSurv_fun(
    condition$n_trt,
    nph::pchaz(
      c(0, t_max),
      c(condition$hazard_trt)
    )
  )

  t_prog_ctrl <- nph::rSurv_fun(
    condition$n_ctrl,
    nph::pchaz(
      c(0, t_max),
      c(condition$prog_rate_ctrl)
    )
  )

  t_prog_trt <- nph::rSurv_fun(
    condition$n_trt,
    nph::pchaz(
      c(0, t_max),
      c(condition$prog_rate_trt)
    )
  )

  t_evt_after_prog_ctrl <- nph::rSurv_fun(
    condition$n_ctrl,
    nph::pchaz(
      c(0, t_max),
      c(condition$hazard_after_prog)
    )
  )

  t_evt_after_prog_ctrl <- t_prog_ctrl + t_evt_after_prog_ctrl

  t_evt_after_prog_trt <- nph::rSurv_fun(
    condition$n_ctrl,
    nph::pchaz(
      c(0, t_max),
      c(condition$hazard_after_prog)
    )
  )

  t_evt_after_prog_trt <- t_prog_trt + t_evt_after_prog_trt

  data_trt <- data.frame(
    t = ifelse(t_prog_trt < t_evt_trt, t_evt_after_prog_trt, t_evt_trt),
    trt = 1L,
    evt = TRUE,
    t_ice = ifelse(t_prog_trt < t_evt_trt, t_prog_trt, Inf),
    ice   = (t_prog_trt < t_evt_trt)
  )

  data_ctrl <- data.frame(
    t = ifelse(t_prog_ctrl < t_evt_ctrl, t_evt_after_prog_ctrl, t_evt_ctrl),
    trt = 0L,
    evt = TRUE,
    t_ice = ifelse(t_prog_ctrl < t_evt_ctrl, t_prog_ctrl, Inf),
    ice   = (t_prog_ctrl < t_evt_ctrl)
  )

  rbind(data_trt, data_ctrl)
}

#' Calculate true summary statistics for scenarios with differential treatment effect in subgroup
#'
#' @param Design Design data.frame for subgroup
#' @param what="os" True summary statistics for which estimand
#' @param cutoff_stats Cutoff time for rmst and average hazard ratios
#' @param fixed_objects=NULL fixed objects not used for now
#'
#' @return For true_summary_statistics_subgroup: the design data.frame
#'   passed as argument with the additional columns:
#' * `rmst_trt` rmst in the treatment group
#' * `median_surv_trt` median survival in the treatment group
#' * `rmst_ctrl` rmst in the control group
#' * `median_surv_ctrl` median survial in the control group
#' * `gAHR` geometric average hazard ratio
#' * `AHR` average hazard ratio
#'
#' @details
#'
#' `what` can be `"os"` for overall survival and `"pfs"` for progression free
#' survival.
#'
#' @export
#'
#' @describeIn generate_subgroup  calculate true summary statistics for subgroup
#'
#' @examples
#' my_design <- merge(
#'     assumptions_progression(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design$follwup <- 15
#' my_design_os  <- true_summary_statistics_subgroup(my_design, "os",  cutoff_stats=my_design$followup)
#' my_design_pfs <- true_summary_statistics_subgroup(my_design, "pfs", cutoff_stats=my_design$followup)
#' my_design_os
#' my_design_pfs
true_summary_statistics_subgroup <- function(Design, what="os", cutoff_stats=10, fixed_objects=NULL){

  true_summary_statistics_progression_rowwise_pfs <- function(condition, cutoff_stats){
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

    # create functions for control group with constant hazard from 0 to t_max
    # minimum of two exponential (constant hazards) distributions is again
    # exponential with sum of the rates
    data_generating_model_ctrl <- nph::pchaz(
      c(0, t_max),
      c(condition$hazard_ctrl + condition$prog_rate_ctrl)
    )

    # minimum of two exponential (constant hazards) distributions is again
    # exponential with sum of the rates
    # create functions for control group with constant hazard from 0 to t_max
    data_generating_model_trt <- nph::pchaz(
      c(0, t_max),
      c(condition$hazard_trt + condition$prog_rate_trt)
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

  true_summary_statistics_progression_rowwise_os <- function(condition, cutoff_stats){
    data_generating_model_ctrl <- nph::pchaz(
      c(0, t_max),
      lambda1 = condition$hazard_ctrl,
      lambda2 = condition$hazard_after_prog,
      lambdaProg = condition$prog_rate_ctrl
    )

    data_generating_model_trt <- nph::pchaz(
      c(0, t_max),
      lambda1 = condition$hazard_trt,
      lambda2 = condition$hazard_after_prog,
      lambdaProg = condition$prog_rate_trt
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

  true_summary_statistics_progression_rowwise <- switch(what,
    "os"  = true_summary_statistics_progression_rowwise_os,
    "pfs" = true_summary_statistics_progression_rowwise_pfs,
    {stop(paste0(gettext("Invalid value for"), " what: ", what, " ", gettext('use "os" for overall survival or "pfs" for progression free survival.')))}
  )

  Design <- Design |>
    split(1:nrow(Design)) |>
    mapply(FUN=true_summary_statistics_progression_rowwise, cutoff_stats = cutoff_stats, SIMPLIFY = FALSE)

  Design <- do.call(rbind, Design)

  Design
}



# TODO: function to callibrate progression rate to percentage who experience progression
