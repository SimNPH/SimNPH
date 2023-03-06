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
  hazard_ctrl       = 0.001518187, # hazard under control (med. survi. 15m)
  hazard_trt        = 0.001265156, # hazard under treatment (med. surv. 18m)
  hazard_after_prog = 0.007590934, # hazard after progression (med. surv. 3m)
  prog_rate_ctrl    = 0.001897734, # hazard rate for disease progression under control
  prog_rate_trt     = c(0.001897734, 0.001423300, 0.001265156), # haz. rate for progression, trt.
  random_withdrawal = 0.01         # rate of random withdrawal
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
generate_progression <- function(condition, fixed_objects=NULL){

  t_evt_ctrl <- fast_rng_fun(
      c(0),
      c(condition$hazard_ctrl)
    )(condition$n_ctrl)

  t_evt_trt <- fast_rng_fun(
      c(0),
      c(condition$hazard_trt)
    )(condition$n_trt)

  t_prog_ctrl <- fast_rng_fun(
      c(0),
      c(condition$prog_rate_ctrl)
    )(condition$n_ctrl)

  t_prog_trt <- fast_rng_fun(
      c(0),
      c(condition$prog_rate_trt)
    )(condition$n_trt)

  t_evt_after_prog_ctrl <- fast_rng_fun(
      c(0),
      c(condition$hazard_after_prog)
    )(condition$n_ctrl)

  t_evt_after_prog_ctrl <- t_prog_ctrl + t_evt_after_prog_ctrl

  t_evt_after_prog_trt <- fast_rng_fun(
      c(0),
      c(condition$hazard_after_prog)
    )(condition$n_ctrl)

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

#' @param Design Design data.frame for subgroup
#' @param what True summary statistics for which estimand
#' @param cutoff_stats (optionally named) cutoff time, see details
#' @param milestones (optionally named) vector of times at which milestone survival should be calculated
#' @param fixed_objects additional settings, see details
#'
#' @return For true_summary_statistics_subgroup: the design data.frame
#'   passed as argument with the additional columns
#'
#' @details
#'
#' `what` can be `"os"` for overall survival and `"pfs"` for progression free
#' survival.
#'
#' The if `fixed_objects` contains `t_max` then this value is used as the
#' maximum time to calculate function like survival, hazard, ... of the data
#' generating models. If this is not given `t_max` is choosen as the minimum of
#' the `1-(1/10000)` quantile of all survival distributions in the model.
#'
#' `cutoff_stats` are the times used to calculate the statistics like average
#' hazard ratios and RMST, that are only calculated up to a certain point.
#'
#' @export
#'
#' @describeIn generate_progression calculate true summary statistics for scenarios with disease progression
#'
#' @examples
#' my_design <- merge(
#'     assumptions_progression(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design_os  <- true_summary_statistics_subgroup(my_design, "os")
#' my_design_pfs <- true_summary_statistics_subgroup(my_design, "pfs")
#' my_design_os
#' my_design_pfs
true_summary_statistics_progression <- function(Design, what="os", cutoff_stats=NULL, fixed_objects=NULL, milestones=NULL){

  true_summary_statistics_progression_rowwise_pfs <- function(condition, cutoff_stats, milestones){

    real_stats <- fast_real_statistics_pchaz(
      Tint_trt =  0, lambda_trt  = condition$hazard_trt  + condition$prog_rate_trt,
      Tint_ctrl = 0, lambda_ctrl = condition$hazard_ctrl + condition$prog_rate_ctrl,
      cutoff = cutoff_stats, N_trt = condition$n_trt, N_ctrl = condition$n_ctrl, milestones=milestones
    )

    res <- cbind(
      condition,
      real_stats
    )

    row.names(res) <- NULL
    res
  }

  true_summary_statistics_progression_rowwise_os <- function(condition, cutoff_stats, milestones){

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

    data_generating_model_ctrl <- subpop_hazVfun_simnph(
      c(0, t_max),
      lambda1 = condition$hazard_ctrl,
      lambda2 = condition$hazard_after_prog,
      lambdaProg = condition$prog_rate_ctrl,
      timezero = TRUE
    )

    data_generating_model_trt <- subpop_hazVfun_simnph(
      c(0, t_max),
      lambda1 = condition$hazard_trt,
      lambda2 = condition$hazard_after_prog,
      lambdaProg = condition$prog_rate_trt,
      timezero = TRUE
    )

    res <- cbind(
      condition,
      internal_real_statistics_pchaz_discrete(
        data_generating_model_trt,
        data_generating_model_ctrl,
        N_trt=condition$n_trt,
        N_ctrl=condition$n_ctrl,
        cutoff = cutoff_stats,
        milestones = milestones
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
    purrr::map(true_summary_statistics_progression_rowwise, cutoff_stats = cutoff_stats, milestones=milestones, .progress = TRUE) |>
    do.call(what=rbind)

  Design
}



#' Calculate progression rate from proportion of patients who progress
#'
#' @param design design data.frame
#'
#' @describeIn generate_progression Calculate progression rate from proportion of patients who progress
#'
#' @return For progression_rate_from_progression_prop: the design data.frame passed as
#'   argument with the additional columns prog_rate_trt, prog_rate_ctrl
#'
#' @details For progression_rate_from_progression_prop, the design data.frame,
#'   has to contain the columns `prog_prop_trt` and `prog_prop_ctrl` with the
#'   proportions of patients, who progress in the respective arms.
#'
#' @export
#'
#' @examples
#' my_design <- merge(
#'     assumptions_progression(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design$prog_rate_ctrl <- NA_real_
#' my_design$prog_rate_trt <- NA_real_
#' my_design$prog_prop_trt <- 0.2
#' my_design$prog_prop_ctrl <- 0.3
#' my_design <- progression_rate_from_progression_prop(my_design)
#' my_design
progression_rate_from_progression_prop <- function(design){

  rowwise_fun <- function(condition){
    # set t_max to 1-1/1000 quantile of control or treatment survival function
    # whichever is later
    t_max <- max(
      log(1000) / condition$hazard_ctrl,
      log(1000) / condition$hazard_trt
    )

    cumhaz_trt <- fast_cumhaz_fun(
      c(                   0),
      c(condition$hazard_trt)
    )

    cumhaz_ctrl <- fast_cumhaz_fun(
      c(                    0),
      c(condition$hazard_ctrl)
    )

    target_fun_trt <- Vectorize(\(r){
      cumhaz_prog_trt <- fast_cumhaz_fun(0, r)
      prob_prog_trt  <- cumhaz_prog_trt(t_max)/(cumhaz_prog_trt(t_max) + cumhaz_trt(t_max))
      prob_prog_trt-condition$prog_prop_trt
    })

    target_fun_ctrl <- Vectorize(\(r){
      cumhaz_prog_ctrl <- fast_cumhaz_fun(0, r)
      prob_prog_ctrl  <- cumhaz_prog_ctrl(t_max)/(cumhaz_prog_ctrl(t_max) + cumhaz_ctrl(t_max))
      prob_prog_ctrl-condition$prog_prop_ctrl
    })

    condition$prog_rate_trt  <- uniroot(target_fun_trt,  interval=c(0, 1e-6), extendInt = "upX", tol=.Machine$double.eps)$root
    condition$prog_rate_ctrl <- uniroot(target_fun_ctrl, interval=c(0, 1e-6), extendInt = "upX", tol=.Machine$double.eps)$root

    condition
  }

  result <- design |>
    split(1:nrow(design)) |>
    purrr::map(rowwise_fun, .progress = TRUE) |>
    do.call(what=rbind)

  result
}

#' @describeIn generate_progression calculate censoring rate from censoring proportion
#'
#' @return for cen_rate_from_cen_prop_progression: design data.frame with the
#'   additional column random_withdrawal
#' @export
#'
#' @details cen_rate_from_cen_prop_progression takes the proportion of
#'   censored patients from the column `censoring_prop`. This column describes
#'   the proportion of patients who are censored randomly before experiencing an
#'   event, without regard to administrative censoring.
#'
#' @examples
#' design <- expand.grid(
#' hazard_ctrl         = 0.001518187, # hazard under control (med. survi. 15m)
#' hazard_trt          = 0.001265156, # hazard under treatment (med. surv. 18m)
#' hazard_after_prog   = 0.007590934, # hazard after progression (med. surv. 3m)
#' prog_rate_ctrl      = 0.001897734, # hazard rate for disease progression under control (med. time to progression 12m)
#' prog_rate_trt       = c(0.001897734, 0.001423300, 0.001265156), # hazard rate for disease progression unter treatment (med. time to progression 12m, 16m, 18m)
#' censoring_prop      = 0.1,         # rate of random withdrawal
#' followup            = 100,         # follow up time
#' n_trt               = 50,          # patients in treatment arm
#' n_ctrl              = 50           # patients in control arm
#' )
#' cen_rate_from_cen_prop_progression(design)
cen_rate_from_cen_prop_progression <- function(design){

  rowwise_fun <- function(condition){
    if(condition$censoring_prop == 0){
      condition$random_withdrawal <- 0.
      return(condition)
    }

    # set t_max to 1-1/1000 quantile of control or treatment survival function
    # whichever is later
    t_max <- max(
      log(1000) / condition$hazard_ctrl,
      log(1000) / condition$hazard_trt
    )

    a <- condition$n_trt / (condition$n_trt + condition$n_ctrl)
    b <- 1-a

    data_generating_model_ctrl <- subpop_hazVfun_simnph(
      c(0, t_max),
      lambda1 = condition$hazard_ctrl,
      lambda2 = condition$hazard_after_prog,
      lambdaProg = condition$prog_rate_ctrl,
      timezero = TRUE
    )

    data_generating_model_trt <- subpop_hazVfun_simnph(
      c(0, t_max),
      lambda1 = condition$hazard_trt,
      lambda2 = condition$hazard_after_prog,
      lambdaProg = condition$prog_rate_trt,
      timezero = TRUE
    )

    cumhaz_trt_tmax  <- tail(data_generating_model_trt$cumhaz, 1)
    cumhaz_ctrl_tmax <- tail(data_generating_model_trt$cumhaz, 1)

    target_fun <- Vectorize(\(r){
      cumhaz_censoring <- fast_cumhaz_fun(0, r)
      prob_cen_ctrl <- cumhaz_censoring(t_max)/(cumhaz_censoring(t_max) + cumhaz_ctrl_tmax)
      prob_cen_trt  <- cumhaz_censoring(t_max)/(cumhaz_censoring(t_max) + cumhaz_trt_tmax)
      prob_cen <- a*prob_cen_trt + b*prob_cen_ctrl
      prob_cen-condition$censoring_prop
    })

    condition$random_withdrawal <- uniroot(target_fun, interval=c(0, 1e-6), extendInt = "upX", tol=.Machine$double.eps)$root

    condition
  }

  result <- design |>
    split(1:nrow(design)) |>
    purrr::map(rowwise_fun, .progress = TRUE) |>
    do.call(what=rbind)

  result

}




#' Calculate hr after onset of treatment effect
#'
#' @param design design data.frame
#' @param target_power_ph target power under proportional hazards
#' @param final_events target events for inversion of Schönfeld Formula, defaults to `condition$final_events`
#' @param target_alpha target alpha level for the power calculation
#'
#' @return For hazard_before_progression_from_PH_effect_size: the design
#'   data.frame passed as argument with the additional column hazard_trt.
#' @export
#'
#' @describeIn generate_progression Calculate hazard in the treatment arm before progression from PH effect size
#'
#' @details `hazard_before_progression_from_PH_effect_size` calculates the
#'   hazard ratio after onset of treatment effect as follows: First calculate
#'   the hazard in the control arm that would give the same median survival
#'   under an exponential model. Then calculate the median survival in the
#'   treatment arm that would give the desired power of the logrank test under
#'   exponential models in control and treatment arm. Then callibrate the hazard
#'   before progression in the treatment arm to give the same median survival
#'   time.
#'
#'   This is a heuristic and to some extent arbitrary approach to calculate
#'   hazard ratios that correspond to reasonable and realistic scenarios.
#'
#' @examples
#' \dontrun{
#' my_design <- hazard_before_progression_from_PH_effect_size(my_design, target_power_ph=0.9)
#' }
hazard_before_progression_from_PH_effect_size <- function(design, target_power_ph=NA_real_, final_events=NA_real_, target_alpha=0.05){

  get_hr_after <- function(condition, target_power_ph=NA_real_, final_events=NA_real_){

    if(condition$effect_size_ph == 0){
      condition$hazard_trt <- condition$hazard_ctrl
      return(condition)
    }

    if(is.na(final_events)){
      if(hasName(condition, "final_events")){
        final_events <- condition$final_events
      } else {
        stop("final_events not given and not present in condition")
      }
    }

    if(is.na(target_power_ph)){
      if(hasName(condition, "effect_size_ph")){
        target_power_ph <- condition$effect_size_ph
      } else {
        stop(gettext("target_ph_power not given and effect_size_ph not present in design"))
      }
    }

    # set t_max to 1/5000 quantile of control arm
    t_max <- log(500) / condition$hazard_ctrl

    model_control <- subpop_hazVfun_simnph(
      c(0, t_max),
      lambda1 = condition$hazard_ctrl,
      lambda2 = condition$hazard_after_prog,
      lambdaProg = condition$prog_rate_ctrl,
      timezero = TRUE
    )

    # setting median to max(t) if less than half die
    # only for uniroot later, not for general use!
    median_progression <- function(mod){
      med <- mod$t[mod$S <= 0.5][1] + 1
      if(is.na(med)){
        med <- max(mod$t)
      }
      med
    }

    ph_hr <- hr_required_schoenfeld(
      final_events,
      alpha=target_alpha,
      beta=(1-target_power_ph),
      p=(condition$n_ctrl/(condition$n_ctrl + condition$n_trt))
    )
    median_ctrl <- median_progression(model_control)

    hazard_ctrl_ph <- uniroot(
      \(h){
        fast_quant_fun(0, h)(0.5) - median_ctrl
      },
      interval=c(1e-8, 0.0001),
      extendInt="downX"
    )$root

    median_trt_ph  <- fast_quant_fun(0, hazard_ctrl_ph * ph_hr)(0.5)

    target_fun_hazard_trt <- function(hazard_after){
      sapply(hazard_after, \(h){
        mod_trt <- SimNPH:::subpop_hazVfun_simnph(
          c(0, t_max),
          lambda1 = h,
          lambda2 = condition$hazard_after_prog,
          lambdaProg = condition$prog_rate_trt,
          timezero = TRUE
        )
        median_trt <- median_progression(mod_trt)
        median_trt_ph - median_trt
      })
    }

    condition$hazard_trt <- uniroot(target_fun_hazard_trt, interval=c(1e-8, 0.0001), extendInt = "upX")$root
    condition
  }

  result <- design |>
    split(1:nrow(design)) |>
    purrr::map(get_hr_after, target_power_ph=target_power_ph, final_events=final_events, .progress=TRUE) |>
    do.call(what=rbind)

  result

}
