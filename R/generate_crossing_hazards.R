#' Generate Dataset with crossing hazards
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condidtion has to contain the following columns:
#'
#'   * n_trt number of paitents in treatment arm
#'   * n_ctrl number of patients in control arm
#'   * crossing time of crossing of the hazards
#'   * hazard_ctrl hazard in the control arm = hazard before onset of treatment
#'     effect
#'   * hazard_trt_before hazard in the treatment arm before onset of treatment effect
#'   * hazard_trt_after hazard in the treatment arm afert onset of treatment effect
#'
#' If fixed_objects is given and contains an element `t_max`, then this is used
#' as the cutoff for the simulation used internally. If t_max is not given in
#' this way the 1-(1/10000) quantile of the survival distribution in the control
#' or treatment arm is used (which ever is larger).
#'
#' @return
#' For generate_crossing_hazards: A dataset with the columns t (time) and trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_crossing_hazards simulates a dataset with crossing
#'   hazards
#'
#' @examples
#' one_simulation <- merge(
#'     assumptions_crossing_hazards(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1) |>
#'   generate_crossing_hazards()
#' head(one_simulation)
#' tail(one_simulation)
generate_crossing_hazards <- function(condition, fixed_objects=NULL){

  # simulate treatment group
  if (condition$crossing < 0){
    # if crossing is smaller than 0 stop with error
    stop(gettext("Time of crossing has to be >= 0"))
  } else if (condition$crossing == 0){
    # if crossing is 0 leave out period bevore treatment effect
    data_trt <- data.frame(
      t = miniPCH::rpch_fun(
        c(0),
        c(condition$hazard_trt_after),
        discrete = TRUE
      )(condition$n_trt),
      trt = 1,
      evt = TRUE
    )
  } else {
    # if crossing is positive simulate in the time intervals bevore and after
    # treatment effect
    data_trt <- data.frame(
      t = miniPCH::rpch_fun(
        c(0, condition$crossing),
        c(condition$hazard_trt_before, condition$hazard_trt_after),
        discrete = TRUE
      )(condition$n_trt),
      trt = 1,
      evt = TRUE
    )
  }

  # simulate control group with constant hazard from 0
  data_ctrl <- data.frame(
    t = miniPCH::rpch_fun(
      c(0),
      c(condition$hazard_ctrl),
      discrete = TRUE
    )(condition$n_trt),
    trt = 0,
    evt = TRUE
  )

  rbind(data_trt, data_ctrl)
}

#' Create an empty assumtions data.frame for generate_crossing_hazards
#'
#' @return For assumptions_crossing_hazards: a design tibble with default values invisibly
#'
#' @details assumptions_crossing_hazards prints the code to generate a default
#'   design tibble for use with generate_crossing_hazards and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_crossing_hazards generate default design tibble
#'
#' @examples
#' Design <- assumptions_crossing_hazards()
#' Design
assumptions_crossing_hazards <- function(){
  skel <- "expand.grid(
  crossing=m2d(seq(0, 10, by=2)), # crossing after of 0, 1, ..., 10 months
  hazard_ctrl=m2r(24),            # median survival control of 24 months
  hazard_trt_before=m2r(18),      # median survival before crossing 18 months
  hazard_trt_after=m2r(36),       # median survival after crossing 36 months
  random_withdrawal=m2r(120)      # median time to random withdrawal 10 years
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}


#' Calculate hr after crossing the hazard functions
#'
#' @param design design data.frame
#' @param target_power_ph target power under proportional hazards
#' @param final_events target events for inversion of Schönfeld Formula, defaults to `condition$final_events`
#' @param target_alpha target one-sided alpha level for the power calculation
#'
#' @return For hr_after_crossing_from_PH_effect_size: the design data.frame passed as
#'   argument with the additional column hazard_trt.
#' @export
#'
#' @describeIn generate_crossing_hazards Calculate hr after crossing of the hazards from PH effect size
#'
#' @details `hr_after_crossing_from_PH_effect_size` calculates the hazard ratio
#'   after crossing of hazards as follows: First, the hazard ratio needed
#'   to archive the desired power under proportional hazards is calculated by
#'   inverting Schönfeld's sample size formula. Second the median survival times
#'   for both arm under this hazard ratio and proportional hazards are
#'   calculated. Finally the hazard rate of the treatment arm after crossing of
#'   hazards is set such that the median survival time is the same as the one
#'   calculated under proportional hazards.
#'
#'   This is a heuristic and to some extent arbitrary approach to calculate
#'   hazard ratios that correspond to reasonable and realistic scenarios.
#'
#' @examples
#' my_design <- merge(
#'     assumptions_crossing_hazards(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#'
#' my_design$final_events <- ceiling((my_design$n_trt + my_design$n_ctrl)*0.75)
#' my_design$hazard_trt <- NA
#' my_design <- hr_after_crossing_from_PH_effect_size(my_design, target_power_ph=0.9)
#' my_design
hr_after_crossing_from_PH_effect_size <- function(design, target_power_ph=NA_real_, final_events=NA_real_, target_alpha=0.025){

  get_hr_after <- function(condition, target_power_ph=NA_real_, final_events=NA_real_){

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

    scale <- 1/condition$hazard_ctrl
    median_ctrl <- miniPCH::qpch_fun(0, scale*condition$hazard_ctrl)(0.5)

    if(target_power_ph == 0){
      median_trt <- median_ctrl
      ph_hr <- 1
    } else {
      ph_hr <- hr_required_schoenfeld(
        final_events,
        alpha=target_alpha,
        beta=(1-target_power_ph),
        p=(condition$n_ctrl/(condition$n_ctrl + condition$n_trt))
      )

      median_trt <- miniPCH::qpch_fun(0, scale*condition$hazard_ctrl*ph_hr)(0.5)
    }

    median_trt_before <- miniPCH::qpch_fun(0, scale*condition$hazard_trt_before)(0.5)

    if(scale*median_ctrl <= condition$crossing ||
       scale*median_trt_before <= condition$crossing
    ){
      warning("Median survival reached before crossing of the hazards curves, calculation not possible")
      condition$hazard_trt_after  <- NA_real_
      condition$target_median_trt <- median_trt * scale
      condition$target_hr         <- ph_hr
      return(condition)
    }

    if(condition$crossing != 0){
      target_fun_hazard_after <- function(hazard_after){
        sapply(hazard_after, \(h){
          median_trt -
            miniPCH::qpch_fun(
              c(0, condition$crossing/scale),
              c(condition$hazard_trt_before*scale, h)
            )(0.5)
        })
      }
    } else {
      target_fun_hazard_after <- function(hazard_after){
        sapply(hazard_after, \(h){
          median_trt -
            miniPCH::qpch_fun(
              c(0),
              c(h)
            )(0.5)
        })
      }
    }

    # setting the lower interval bound to 0 and f.lower to -Inf
    # and extendInt="upX" guarantees, that the root is searched on all positives
    # but also that the target function is never evaluated at non-positive values
    my_root <- uniroot(
      target_fun_hazard_after,
      interval=c(0, 2),
      f.lower = -Inf,
      extendInt = "upX",
      tol=2*.Machine$double.eps
    )

    condition$target_median_trt <- median_trt * scale
    condition$hazard_trt_after  <- my_root$root / scale
    condition$target_hr         <- ph_hr
    condition
  }

  result <- design |>
    split(1:nrow(design)) |>
    lapply(get_hr_after, target_power_ph=target_power_ph, final_events=final_events) |>
    do.call(what=rbind)

  result

}

#' @describeIn generate_crossing_hazards calculate censoring rate from censoring proportion
#'
#' @return for cen_rate_from_cen_prop_crossing_hazards: design data.frame with the
#'   additional column random_withdrawal
#' @export
#'
#' @details cen_rate_from_cen_prop_crossing_hazards takes the proportion of
#'   censored patients from the column `censoring_prop`. This column describes
#'   the proportion of patients who are censored randomly before experiencing an
#'   event, without regard to administrative censoring.
#'
#' @examples
#' design <- data.frame(
#'   crossing = c(2, 4, 6),
#'   hazard_ctrl = c(0.05, 0.05, 0.05),
#'   hazard_trt_before = c(0.025, 0.025, 0.025),
#'   hazard_trt_after = c(0.1, 0.1, 0.1),
#'   censoring_prop = c(0.1, 0.3, 0.2),
#'   n_trt = c(50, 50, 50),
#'   n_ctrl = c(50, 50, 50),
#'   followup = c(200, 200, 200),
#'   recruitment = c(50, 50, 50)
#' )
#' cen_rate_from_cen_prop_crossing_hazards(design)
cen_rate_from_cen_prop_crossing_hazards <- function(design){

  rowwise_fun <- function(condition){
    if(is.na(condition$hazard_trt_after)){
      return(NA_real_)
    }

    if(condition$censoring_prop == 0){
      condition$random_withdrawal <- 0.
      return(condition)
    }

    # set t_max to 1-1/10000 quantile of control or treatment survival function
    # whichever is later
    t_max <- max(
      log(10000) / condition$hazard_ctrl,
      log(10000) / condition$hazard_trt
    )

    if(condition$crossing != 0){
      cumhaz_trt <- miniPCH::chpch_fun(
        c(                          0,         condition$crossing),
        c(condition$hazard_trt_before, condition$hazard_trt_after)
      )(t_max)
    } else {
      cumhaz_trt <- miniPCH::chpch_fun(
        c(condition$crossing),
        c(condition$hazard_trt_after)
      )(t_max)
    }

    cumhaz_ctrl <- miniPCH::chpch_fun(
      c(                    0),
      c(condition$hazard_ctrl)
    )(t_max)

    condition$random_withdrawal <- censoring_prop_from_cumhaz(
      n_trt          = condition$n_trt,
      n_ctrl         = condition$n_ctrl,
      censoring_prop = condition$censoring_prop,
      cumhaz_ctrl    = cumhaz_ctrl,
      cumhaz_trt     = cumhaz_trt,
      t_max          = t_max
    )

    condition
  }

  result <- design |>
    split(1:nrow(design)) |>
    lapply(rowwise_fun) |>
    do.call(what=rbind)

  result

}

#' Calculate true summary statistics for scenarios with crossing hazards
#'
#' @param Design Design data.frame for crossing hazards
#' @param cutoff_stats (optionally named) cutoff time, see details
#' @param milestones (optionally named) vector of times at which milestone survival should be calculated
#' @param fixed_objects additional settings, see details
#'
#' @return For true_summary_statistics_crossing_hazards: the design data.frame
#'   passed as argument with additional columns,
#'
#' @export
#'
#' @details
#'
#' `cutoff_stats` are the times used to calculate the statistics like average
#' hazard ratios and RMST, that are only calculated up to a certain point.
#'
#' @describeIn generate_crossing_hazards  calculate true summary statistics for crossing hazards
#'
#' @examples
#' my_design <- merge(
#'     assumptions_crossing_hazards(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design$follwup <- 15
#' my_design <- true_summary_statistics_crossing_hazards(my_design)
#' my_design
true_summary_statistics_crossing_hazards <- function(Design, cutoff_stats=NULL, milestones=NULL, fixed_objects=NULL){

  true_summary_statistics_crossing_hazards_rowwise <- function(condition, cutoff_stats, milestones){
    # create functions for treatment group
    if (condition$crossing < 0){
      # if crossing is smaller than 0 stop with error
      stop(gettext("Time of crossing has to be >= 0"))
    } else if (condition$crossing == 0){
      # if crossing is 0 leave out period bevore treatment effect
      real_stats <- fast_real_statistics_pchaz(
        Tint_trt = 0, lambda_trt  = condition$hazard_trt_after,
        Tint_ctrl= 0, lambda_ctrl = condition$hazard_ctrl,
        cutoff = cutoff_stats, N_trt = condition$n_trt, N_ctrl = condition$n_ctrl, milestones = milestones
      )

    } else {
      # if crossing is positive create piecewise constant hazards and respective
      # functions in the time intervals bevore and after treatment effect
      real_stats <- fast_real_statistics_pchaz(
        Tint_trt = c(0, condition$crossing), lambda_trt  = c(condition$hazard_trt_before, condition$hazard_trt_after),
        Tint_ctrl= 0, lambda_ctrl = condition$hazard_ctrl,
        cutoff = cutoff_stats, N_trt = condition$n_trt, N_ctrl = condition$n_ctrl, milestones = milestones
      )
    }


    res <- cbind(
      condition,
      real_stats
    )

    row.names(res) <- NULL
    res
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    lapply(true_summary_statistics_crossing_hazards_rowwise, cutoff_stats = cutoff_stats, milestones=milestones)

  Design <- do.call(rbind, Design)

  Design
}
