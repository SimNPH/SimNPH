# function to calculate the real summary statistics for piecewise constant hazards
#
# arguments:
#  * haz_trt     hazard function of the treatment arm
#  * pdf_trt     probability density function of the treatment arm
#  * surv_trt    survival function of the treatment arm
#  * quant_trt   quantile function of the treatment arm
#  * haz_ctrl    hazard function of the control arm
#  * pdf_ctrl    probability density function of the control arm
#  * surv_ctrl   survival function of the control arm
#  * quant_ctrl  quantile function of the control arm
#  * N_trt       number of patients in the treatment arm
#  * N_ctrl      number of patients in the control arm
#  * cutoff      cutoff used to calculate rmst and average hazard ratios
#  * milestones  milestones at which to compute the milestone survival
fast_real_statistics <- function(haz_trt,
                                 pdf_trt,
                                 surv_trt,
                                 quant_trt,
                                 haz_ctrl,
                                 pdf_ctrl,
                                 surv_ctrl,
                                 quant_ctrl,
                                 N_trt = 1,
                                 N_ctrl = 1,
                                 cutoff = NULL,
                                 milestones = NULL) {
  res <- data.frame(median_survival_trt = quant_trt(0.5),
                    median_survival_ctrl = quant_ctrl(0.5))

  # if cutoff time for rmst and gahr are given
  if (!is.null(cutoff)) {
    if (is.null(names(cutoff))) {
      names(cutoff) <- as.character(cutoff)
    } else {
      names(cutoff)[names(cutoff) %in% c("", NA_character_)] <-
        as.character(cutoff[names(cutoff) %in% c("", NA_character_)])
    }

    # define helper functions for the calculation of the average hazard ratios
    h  <- \(t) {
      haz_trt(t) + haz_ctrl(t)
    }
    f  <-
      \(t) {
        (1 / (N_trt + N_ctrl)) * (N_trt * pdf_trt(t) + N_ctrl * pdf_ctrl(t))
      }

    myint<-function(f,lower,upper,steps=1000) {
      delta<-(upper-lower)/steps
      x<-seq(lower+delta/2,upper-delta/2,delta)
      sum(f(x)*delta)
    }
    rmst_ahr <- purrr::imap(cutoff, function(cutoff, label) {
      #AHRoc
      true_avg_HR_fun0<-function(x) surv_ctrl(x)*pdf_trt(x)
      true_avg_HR_fun1<-function(x) surv_trt(x)*pdf_ctrl(x)
      #integrate can run into numeric problems with NPH objects, myint is more robust
      Int0<-myint(true_avg_HR_fun0,lower=0,upper=cutoff)
      Int1<-myint(true_avg_HR_fun1,lower=0,upper=cutoff)
      AHRoc_myint<-Int0/Int1

      Int0A<-integrate(true_avg_HR_fun0,lower=0,upper=cutoff)
      Int1A<-integrate(true_avg_HR_fun1,lower=0,upper=cutoff)
      AHRoc_integrate<-Int0A$value/Int1A$value



      tmp <- data.frame(
        rmst_trt = integrate(surv_trt, 0, cutoff)$value,
        rmst_ctrl = integrate(surv_ctrl, 0, cutoff)$value,
        gAHR = exp(
          integrate(\(t) {
            log(haz_trt(t) / haz_ctrl(t)) * f(t)
          }, 0, cutoff, stop.on.error = FALSE)$value /
            integrate(\(t) {
              f(t)
            }, 0, cutoff, stop.on.error = FALSE)$value
        ),
        AHR = integrate(\(t) {
          (haz_trt(t) / h(t)) * f(t)
        }, 0, cutoff, stop.on.error = FALSE)$value /
          integrate(\(t) {
            (haz_ctrl(t) / h(t)) * f(t)
          }, 0, cutoff, stop.on.error = FALSE)$value,
        #### Add correctly weighted AHR here
        AHRoc = AHRoc_integrate,
        AHRoc_robust = AHRoc_myint
      )

      names(tmp) <- paste0(names(tmp), "_", label)
      tmp
    }) |>
      unname() |>
      do.call(what = cbind)

    res <- cbind(res, rmst_ahr)
  }

  # if times for milestone survival are given calculate milestone survival at given times
  if (!is.null(milestones)) {
    if (is.null(names(milestones))) {
      names(milestones) <- as.character(milestones)
    } else {
      names(milestones)[names(milestones) %in% c("", NA_character_)] <-
        as.character(milestones[names(milestones) %in% c("", NA_character_)])
    }

    milestones <- purrr:::imap(milestones, function(v, label) {
      tmp <- data.frame(
        milestone_survival_trt  = surv_trt(v),
        milestone_survival_ctrl = surv_ctrl(v)
      )

      names(tmp) <- paste0(names(tmp), "_", label)
      tmp
    }) |>
      unname() |>
      do.call(what = cbind)

    res <- cbind(res,
                 milestones)
  }

  res
}

fast_real_statistics_pchaz <- function(
    Tint_trt,  lambda_trt,
    Tint_ctrl, lambda_ctrl,
    N_trt=1, N_ctrl=1, cutoff=NULL, milestones=NULL
){

  fast_real_statistics(
    haz_trt    =   fast_haz_fun( Tint_trt,  lambda_trt),
    pdf_trt    =   fast_pdf_fun( Tint_trt,  lambda_trt),
    surv_trt   =  fast_surv_fun( Tint_trt,  lambda_trt),
    quant_trt  = fast_quant_fun( Tint_trt,  lambda_trt),
    haz_ctrl   =   fast_haz_fun(Tint_ctrl, lambda_ctrl),
    pdf_ctrl   =   fast_pdf_fun(Tint_ctrl, lambda_ctrl),
    surv_ctrl  =  fast_surv_fun(Tint_ctrl, lambda_ctrl),
    quant_ctrl = fast_quant_fun(Tint_ctrl, lambda_ctrl),
    N_trt=N_trt, N_ctrl=N_ctrl, cutoff=cutoff, milestones=milestones
  )
}

#' Calculate true summary statistics for scenarios with delayed treatment effect
#'
#' @param Design Design data.frame for delayed effect
#' @param cutoff_stats (optionally named) cutoff times, see details
#' @param milestones (optionally named) vector of times at which milestone survival should be calculated
#' @param fixed_objects additional settings, see details
#'
#' @return For true_summary_statistics_delayed_effect: the design data.frame
#'   passed as argument with additional columns
#'
#' @export
#'
#' @details
#'
#' `cutoff_stats` are the times used to calculate the statistics like average
#' hazard ratios and RMST, that are only calculated up to a certain point.
#'
#' @describeIn generate_delayed_effect calculate true summary statistics for delayed effect
#'
#' @examples
#' my_design <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   )
#' my_design <- true_summary_statistics_delayed_effect(my_design)
#' my_design
true_summary_statistics_delayed_effect <- function(Design, cutoff_stats=NULL, milestones=NULL, fixed_objects=NULL){

  true_summary_statistics_delayed_effect_rowwise <- function(condition, cutoff_stats, milestones){

    # create functions for treatment group
    if (condition$delay < 0){
      # if delay is smaller than 0 stop with error
      stop(gettext("Delay has to be >= 0"))
    } else if (condition$delay == 0){
      # if delay is 0 leave out period bevore treatment effect
      real_stats <- fast_real_statistics_pchaz(
        Tint_trt = 0,  lambda_trt = condition$hazard_trt,
        Tint_ctrl = 0, lambda_ctrl = condition$hazard_ctrl,
        cutoff = cutoff_stats, N_trt = condition$n_trt, N_ctrl = condition$n_ctrl, milestones=milestones
      )

    } else {
      # if delay is positive create piecewise constant hazards and respective
      # functions in the time intervals bevore and after treatment effect
      real_stats <- fast_real_statistics_pchaz(
        Tint_trt = c(0, condition$delay),  lambda_trt = c(condition$hazard_ctrl, condition$hazard_trt),
        Tint_ctrl = 0, lambda_ctrl = condition$hazard_ctrl,
        cutoff = cutoff_stats, N_trt = condition$n_trt, N_ctrl = condition$n_ctrl, milestones=milestones
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
    lapply(true_summary_statistics_delayed_effect_rowwise, milestones=milestones, cutoff_stats=cutoff_stats)

  Design <- do.call(rbind, Design)

  Design
}

