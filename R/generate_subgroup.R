#' Generate Dataset with different treatment effect in subgroup
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
#'   * hazard_subgroup hazard in the subgroup in the treatment arm
#'   * prevalence proportion of cured patients
#'
#' @return
#' For generate_subgroup: A dataset with the columns t (time) and trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_subgroup simulates a dataset with a mixture of cured
#'   patients
#'
#' @examples
#' one_simulation <- merge(
#'     assumptions_subgroup(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1) |>
#'   generate_subgroup()
#' head(one_simulation)
#' tail(one_simulation)
generate_subgroup <- function(condition, fixed_objects=NULL){
  if (condition$prevalence < 0 || condition$prevalence > 1) {
    stop(gettext("Subgroup prevalence has to be between 0 and 1"))
  }

  counts <- rmultinom(1, condition$n_trt, prob=c(condition$prevalence, 1-condition$prevalence))

  data_subgroup <- data.frame(
    t = fast_rng_fun(0, condition$hazard_subgroup)(counts[1]),
    trt = 1,
    evt = TRUE,
    subgroup = 1
  )

  data_trt <- data.frame(
    t = fast_rng_fun(0, condition$hazard_trt)(counts[2]),
    trt = 1,
    evt = TRUE,
    subgroup = 0
  )

  data_ctrl <- data.frame(
    t = fast_rng_fun(0, condition$hazard_ctrl)(condition$n_ctrl),
    trt = 0,
    evt = TRUE,
    subgroup = rbinom(condition$n_ctrl, 1, condition$prevalence)
  )

  rbind(data_trt, data_subgroup, data_ctrl)
}

#' Create an empty assumtions data.frame for generate_subgroup
#'
#' @return For assumptions_subgroup: a design tibble with default values invisibly
#'
#' @details assumptions_subgroup prints the code to generate a default
#'   design tibble for use with generate_subgroup and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_subgroup generate default assumptions tibble
#'
#' @examples
#' Design <- assumptions_subgroup()
#' Design
assumptions_subgroup <- function(){
  skel <- "expand.grid(
  hazard_ctrl       = 0.003795467,       # hazard under control (med. survi. 6m)
  hazard_trt        = 0.001265156,       # hazard under treatment (med. surv. 18m)
  hazard_subgroup   = 9.488668e-05,      # hazard for subgroup under treatment (med. surv. 20y)
  prevalence        = seq(0, 1, by=0.2), # proportion of patients belonging to subgroup
  random_withdrawal = 0.01               # rate of random withdrawal
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}



#' Calculate true summary statistics for scenarios with differential treatment effect in subgroup
#'
#' @param Design Design data.frame for subgroup
#' @param cutoff_stats (optionally named) cutoff times, see details
#' @param milestones (optionally named) vector of times at which milestone survival should be calculated
#' @param fixed_objects additional settings, see details
#'
#' @return For true_summary_statistics_subgroup: the design data.frame
#'   passed as argument with the additional columns
#'
#' @export
#'
#' @details
#'
#' `cutoff_stats` are the times used to calculate the statistics like average
#' hazard ratios and RMST, that are only calculated up to a certain point.
#'
#' @describeIn generate_subgroup  calculate true summary statistics for subgroup
#'
#' @examples
#' my_design <- merge(
#'     assumptions_subgroup(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design <- true_summary_statistics_subgroup(my_design)
#' my_design
true_summary_statistics_subgroup <- function(Design, cutoff_stats=NULL, milestones=NULL, fixed_objects=NULL){

  true_summary_statistics_subgroup_rowwise <- function(condition, cutoff_stats, milestones){

    if (condition$prevalence < 0 || condition$prevalence > 1) {
      stop(gettext("Subgroup prevalence has to be between 0 and 1"))
    }

    haz_trt   <-   mixture_haz_fun(
      c(1-condition$prevalence, condition$prevalence),
      pdfs = list(
        fast_pdf_fun(0, condition$hazard_trt),
        fast_pdf_fun(0, condition$hazard_subgroup)
      ),
      survs = list(
        fast_surv_fun(0, condition$hazard_trt),
        fast_surv_fun(0, condition$hazard_subgroup)
      )
    )

    pdf_trt   <-   mixture_pdf_fun(
      c(1-condition$prevalence, condition$prevalence),
      list(
        fast_pdf_fun(0, condition$hazard_trt),
        fast_pdf_fun(0, condition$hazard_subgroup)
      )
    )

    surv_trt  <-  mixture_surv_fun(
      c(1-condition$prevalence, condition$prevalence),
      list(
        fast_surv_fun(0, condition$hazard_trt),
        fast_surv_fun(0, condition$hazard_subgroup)
      )
    )

    quant_trt <- mixture_quant_fun(
      c(1-condition$prevalence, condition$prevalence),
      cdfs=list(
        fast_cdf_fun(0, condition$hazard_trt),
        fast_cdf_fun(0, condition$hazard_subgroup)
      ),
      quants=list(
        fast_quant_fun(0, condition$hazard_trt),
        fast_quant_fun(0, condition$hazard_subgroup)
      )
    )

    haz_ctrl   <-   fast_haz_fun(0, condition$hazard_ctrl)
    pdf_ctrl   <-   fast_pdf_fun(0, condition$hazard_ctrl)
    surv_ctrl  <-  fast_surv_fun(0, condition$hazard_ctrl)
    quant_ctrl <- fast_quant_fun(0, condition$hazard_ctrl)

    real_stats <- fast_real_statistics(
      haz_trt,  pdf_trt,  surv_trt, quant_trt,
      haz_ctrl, pdf_ctrl, surv_ctrl, quant_ctrl,
      N_trt=condition$n_trt, N_ctrl=condition$n_ctrl, cutoff=cutoff_stats, milestones=milestones
    )

    res <- cbind(
      condition,
      real_stats
    )

    row.names(res) <- NULL
    res
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    lapply(true_summary_statistics_subgroup_rowwise, cutoff_stats = cutoff_stats, milestones=milestones)

  Design <- do.call(rbind, Design)

  Design
}


#' Calculate hazards in treatment arm in subgroup and compliment
#'
#' @param design design data.frame
#' @param target_power_ph target power under proportional hazards
#' @param final_events target events for inversion of Schönfeld Formula, defaults to `condition$final_events`
#' @param target_alpha target alpha level for the power calculation
#'
#' @return For hazard_subgroup_from_PH_effect_size: the design data.frame passed as
#'   argument with the additional columns hazard_trt and hazard_subgroup.
#' @export
#'
#' @describeIn generate_subgroup Calculate hazards in treatement arm
#'
#' @details `hazard_subgroup_from_PH_effect_size` calculates the hazard rate in
#'   the subgroup and the compliment of the subgroup in the treatment arm as
#'   follows: First, the hazard ratio needed to archive the desired power under
#'   proportional hazards is calculated by inverting Schönfeld's sample size
#'   formula. Second the median survival times for both arms under this hazard
#'   ratio and proportional hazards are calculated. Finally the hazard rate of
#'   the treatment arm in the subgroup and its complement are set such that the
#'   median survival time is the same as the one calculated under proportional
#'   hazards.
#'
#'   This is a heuristic and to some extent arbitrary approach to calculate
#'   hazard ratios that correspond to reasonable and realistic scenarios.
#'
#' @examples
#' my_design <- merge(
#'     assumptions_subgroup(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design$hazard_trt <- NA
#' my_design$hazard_subgroup <- NA
#' my_design$hr_subgroup_relative <- 0.9
#' my_design <- hazard_subgroup_from_PH_effect_size(my_design, target_power_ph=0.9)
#' my_design
hazard_subgroup_from_PH_effect_size <- function(design, target_power_ph=NA_real_, final_events=NA_real_, target_alpha=0.025){

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

    if(target_power_ph == 0){
      condition$hazard_trt <- condition$hazard_ctrl
    }

    ph_hr <- hr_required_schoenfeld(
      final_events,
      alpha=target_alpha,
      beta=(1-target_power_ph),
      p=(condition$n_ctrl/(condition$n_ctrl + condition$n_trt))
    )

    scale <- 1/condition$hazard_ctrl

    median_trt  <- fast_quant_fun(0, scale * condition$hazard_ctrl * ph_hr)(0.5)
    median_ctrl <- fast_quant_fun(0, scale * condition$hazard_ctrl        )(0.5)

    if(target_power_ph == 0){
      median_trt <- median_ctrl
    }

    target_fun_hazards_subgroups <- function(hazard_compliment){
      sapply(hazard_compliment, \(h){
        my_quant_fun <- mixture_quant_fun(
          c(condition$prevalence, 1-condition$prevalence),
          cdfs = list(
            fast_cdf_fun(0, h*condition$hr_subgroup_relative),
            fast_cdf_fun(0, h)
          ),
          quants = list(
            fast_quant_fun(0, h*condition$hr_subgroup_relative),
            fast_quant_fun(0, h)
          )
          )
        median_trt - my_quant_fun(0.5)
      })
    }

    my_root <- uniroot(
      target_fun_hazards_subgroups,
      interval=c(0, 2),
      f.lower = -Inf,
      extendInt = "upX",
      tol=2*.Machine$double.eps
    )

    condition$hazard_trt <- my_root$root / scale
    condition$target_median_trt <- median_trt * scale
    condition$hazard_subgroup <- condition$hazard_trt * condition$hr_subgroup_relative
    condition
  }

  result <- design |>
    split(1:nrow(design)) |>
    purrr::map(get_hr_after, target_power_ph=target_power_ph, final_events=final_events, .progress=TRUE) |>
    do.call(what=rbind)

  result
}

#' @describeIn generate_subgroup calculate censoring rate from censoring proportion
#'
#' @return for cen_rate_from_cen_prop_subgroup: design data.frame with the
#'   additional column random_withdrawal
#' @export
#'
#' @details cen_rate_from_cen_prop_subgroup takes the proportion of
#'   censored patients from the column `censoring_prop`. This column describes
#'   the proportion of patients who are censored randomly before experiencing an
#'   event, without regard to administrative censoring.
#'
#' @examples
#' design <- expand.grid(
#'   hazard_ctrl=0.2,                   # hazard under control and before treatment effect
#'   hazard_trt=0.02,                   # hazard after onset of treatment effect
#'   hazard_subgroup=0.01,              # hazard in the subgroup in treatment
#'   prevalence = c(0.2, 0.5),           # subgroup prevalence
#'   censoring_prop=c(0.1, 0.25, 0.01), # 10%, 25%, 1% random censoring
#'   followup=100,                      # followup of 100 days
#'   n_trt=50,                          # 50 patients treatment
#'   n_ctrl=50                          # 50 patients control
#' )
#' cen_rate_from_cen_prop_subgroup(design)
cen_rate_from_cen_prop_subgroup <- function(design){

  rowwise_fun <- function(condition){
    if(is.na(condition$hazard_trt)){
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

    cumhaz_trt <- mixture_cumhaz_fun(
      c(condition$prevalence, 1-condition$prevalence),
      survs = list(
        fast_surv_fun(0, condition$hazard_subgroup),
        fast_surv_fun(0, condition$hazard_trt)
      )
    )(t_max)

    cumhaz_ctrl <- fast_cumhaz_fun(
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
    purrr::map(rowwise_fun, .progress = TRUE) |>
    do.call(what=rbind)

  result

}
