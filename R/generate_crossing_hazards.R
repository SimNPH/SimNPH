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
      t = fast_rng_fun(
        c(0),
        c(condition$hazard_trt_after)
      )(condition$n_trt),
      trt = 1,
      evt = TRUE
    )
  } else {
    # if crossing is positive simulate in the time intervals bevore and after
    # treatment effect
    data_trt <- data.frame(
      t = fast_rng_fun(
        c(0, condition$crossing),
        c(condition$hazard_trt_before, condition$hazard_trt_after)
      )(condition$n_trt),
      trt = 1,
      evt = TRUE
    )
  }

  # simulate control group with constant hazard from 0
  data_ctrl <- data.frame(
    t = fast_rng_fun(
        c(0),
        c(condition$hazard_ctrl)
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
  crossing=seq(0, 10, by=2),  # crossing of hazards after 0, 1, ..., 10 days
  hazard_ctrl=0.05,        # hazard under control and before treatment effect
  hazard_trt_before=0.025, # hazard before crossing of the hazard curves
  hazard_trt_after=0.1,    # hazard after crossing of the hazard curves
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

#' Calculate hr after crossing of hazards from gAHT
#'
#' @param design design data.frame
#' @param target_gAHR target geometric average hazard ratio
#' @param cutoff=NA_real_ time until which the gAHR should be calculated, defaults to `condition$followup`
#'
#' @return For hr_after_crossing_from_gAHR: the design data.frame passed as
#'   argument with the additional column hazard_trt.
#' @export
#'
#' @describeIn generate_crossing_hazards  Calculate hr after crossing of hazards from gAHR
#'
#' @examples
#' my_design <- merge(
#'     assumptions_crossing_hazards(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design$hazard_trt <- NA
#' my_design <- hr_after_crossing_from_gAHR(my_design, 0.8, 200)
#' my_design
hr_after_crossing_from_gAHR <- function(design, target_gAHR, cutoff=NA_real_){

  fast_gAHR <- function(hazard_trt_after, condition, cutoff, target_gAHR=1, N_trt=1, N_ctrl=1){
    h1 <- fast_haz_fun(c(0, condition$crossing), c(condition$hazard_trt_before, hazard_trt_after))
    h0 <- fast_haz_fun(c(0), c(condition$hazard_ctrl))

    f1 <- fast_pdf_fun(c(0, condition$crossing), c(condition$hazard_trt_before, hazard_trt_after))
    f0 <- fast_pdf_fun(c(0), c(condition$hazard_ctrl))

    f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*f1(t) + N_ctrl*f0(t))}
    w  <- \(t){1} # Cox

    gAHR <- exp(integrate(\(t){log(h1(t) / h0(t)) * f(t) * w(t)}, 0, cutoff)$value)

    gAHR-target_gAHR
  }

  get_hr_after <- function(condition, cutoff=cutoff){
    if(is.na(cutoff)){
      if(hasName(condition, "followup")){
        cutoff <- condition$followup
      } else {
        stop(gettext("cutoff not given and followup not present in design"))
      }
    }

    condition$hazard_trt_after <- uniroot(fast_gAHR, interval = c(1e-8, 1), condition=condition, cutoff=cutoff, target_gAHR=target_gAHR)$root
    condition
  }

  result <- design |>
    split(1:nrow(design)) |>
    lapply(get_hr_after, cutoff=cutoff) |>
    do.call(what=rbind)

  result
}


#' Calculate hr after crossing the hazard functions
#'
#' @param design design data.frame
#' @param target_power_ph target power under proportional hazards
#' @param final_events=NA_real_ target events for inversion of Schönfeld Formula, defaults to `condition$final_events`
#' @param target_alpha=0.05 target alpha level for the power calculation
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
hr_after_crossing_from_PH_effect_size <- function(design, target_power_ph=NA_real_, final_events=NA_real_, target_alpha=0.05){
  # hazard ratio required, inverted Schönfeld sample size formula
  hr_required_schoenfeld <- function(Nevt, alpha=0.05, beta=0.2, p=0.5){
    exp( (qnorm(beta) + qnorm(alpha)) / sqrt(p*(1-p)*Nevt) )
  }

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

    ph_hr <- hr_required_schoenfeld(final_events, alpha=target_alpha, beta=(1-target_power_ph), p=(condition$n_ctrl/(condition$n_ctrl + condition$n_trt)))

    median_trt <- fast_quant_fun(0, condition$hazard_ctrl * ph_hr)(0.5)
    median_ctrl <- fast_quant_fun(0, condition$hazard_ctrl        )(0.5)
    median_trt_before <- fast_quant_fun(0, condition$hazard_trt_before)(0.5)

    if(target_power_ph == 0){
      median_trt <- median_ctrl
    }

    if(median_trt <= condition$crossing ||
       median_ctrl <= condition$crossing ||
       median_trt_before <= condition$crossing
       ){
      warning("Median survival reached before crossing of the hazards curves, calculation not possible")
      condition$hazard_trt_after <- NA_real_
      return(condition)
    }

    target_fun_hazard_after <- function(hazard_after){
      sapply(hazard_after, \(h){
        median_trt - fast_quant_fun(c(0, condition$crossing), c(condition$hazard_trt_before, h))(0.5)
      })
    }

    condition$hazard_trt_after  <- uniroot(target_fun_hazard_after, interval=c(1e-8, condition$hazard_ctrl))$root
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
#' design <- expand.grid(
#'   crossing=seq(0, 10, by=5),         # crossing at 0, 1, ..., 10 days
#'   hazard_ctrl=0.2,                   # hazard under control and before treatment effect
#'   hazard_trt=0.02,                   # hazard after onset of treatment effect
#'   censoring_prop=c(0.1, 0.25, 0.01), # 10%, 25%, 1% random censoring
#'   followup=100,                      # followup of 100 days
#'   n_trt=50,                          # 50 patients treatment
#'   n_ctrl=50                          # 50 patients control
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

    t_max <- condition$followup * 10

    a <- condition$n_trt / (condition$n_trt + condition$n_ctrl)
    b <- 1-a

    cumhaz_trt <- fast_cumhaz_fun(
      c(                          0,         condition$crossing),
      c(condition$hazard_trt_before, condition$hazard_trt_after)
    )

    cumhaz_ctrl <- fast_cumhaz_fun(
      c(                    0),
      c(condition$hazard_ctrl)
    )

    target_fun <- Vectorize(\(r){
      cumhaz_censoring <- fast_cumhaz_fun(0, r)
      prob_cen_ctrl <- cumhaz_censoring(t_max)/(cumhaz_censoring(t_max) + cumhaz_ctrl(t_max))
      prob_cen_trt  <- cumhaz_censoring(t_max)/(cumhaz_censoring(t_max) + cumhaz_trt(t_max))
      prob_cen <- a*prob_cen_trt + b*prob_cen_ctrl
      prob_cen-condition$censoring_prop
    })

    condition$random_withdrawal <- uniroot(target_fun, interval=c(0, 1e-6), extendInt = "upX", tol=.Machine$double.eps)$root

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
#' @param cutoff_stats=NA_real_ cutoff time, see details
#' @param milestones=NULL (optionally named) vector of times at which milestone survival should be calculated
#' @param fixed_objects=NULL additional settings, see details
#'
#' @return For true_summary_statistics_crossing_hazards: the design data.frame
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
#' @details
#'
#' The if `fixed_objects` contains `t_max` then this value is used as the
#' maximum time to calculate function like survival, hazard, ... of the data
#' generating models. If this is not given `t_max` is choosen as the minimum of
#' the `1-(1/10000)` quantile of all survival distributions in the model.
#'
#' `cutoff_stats` is the time used to calculate the statistics like average
#' hazard ratios and RMST, that are only calculated up to a certain point. It
#' defaults to `NA_real_` in which case the variable `followup` from the Design
#' dataset is used. If `followup` is also not set it uses `t_max`.
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
true_summary_statistics_crossing_hazards <- function(Design, cutoff_stats=NA_real_, milestones=NULL, fixed_objects=NULL){

  true_summary_statistics_crossing_hazards_rowwise <- function(condition, cutoff_stats, milestones){

    # if t_max is not given in fixed_objects
    if(is.null(fixed_objects) || (!hasName(fixed_objects, "t_max"))){
      # set t_max to 1-1/10000 quantile of control or treatment survival function
      # whichever is later
      t_max <- max(
        log(10000) / condition$hazard_ctrl,
        log(10000) / condition$hazard_trt_after,
        log(10000) / condition$hazard_trt_before
      )
    } else {
      t_max <- fixed_objects$t_max
    }

    if(is.na(cutoff_stats)){
      if(hasName(condition, "followup")){
        cutoff_stats <- condition$followup
      } else {
        cutoff_stats <- t_max
      }
    }


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
      real_stats,
      cutoff_used=cutoff_stats
    )

    row.names(res) <- NULL
    res
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    mapply(FUN=true_summary_statistics_crossing_hazards_rowwise, cutoff_stats = cutoff_stats, MoreArgs = list(milestones=milestones), SIMPLIFY = FALSE)

  Design <- do.call(rbind, Design)

  Design
}
