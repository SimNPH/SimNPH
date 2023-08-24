internal_real_statistics_pchaz_discrete <- function(data_gen_model_trt, data_gen_model_ctrl, N_trt=1, N_ctrl=1, cutoff=NULL, milestones=NULL){
  stopifnot(
    length(data_gen_model_trt$t) == length(data_gen_model_ctrl$t)
  )

  stopifnot(
    all(data_gen_model_trt$t == data_gen_model_ctrl$t)
  )

  res <- data.frame(
    median_survival_trt = ifelse(
      any(data_gen_model_trt$F >= 0.5),
      data_gen_model_trt$t[data_gen_model_trt$F >= 0.5][1],
      Inf
    ),
    median_survival_ctrl= ifelse(
      any(data_gen_model_ctrl$F >= 0.5),
      data_gen_model_ctrl$t[data_gen_model_ctrl$F >= 0.5][1],
      Inf
    )
  )

  if(!is.null(cutoff)){

    if(is.null(names(cutoff))){
      names(cutoff) <- as.character(cutoff)
    } else {
      names(cutoff)[names(cutoff) %in% c("", NA_character_)] <- as.character(cutoff[names(cutoff) %in% c("", NA_character_)])
    }

    h <- data_gen_model_trt$haz + data_gen_model_ctrl$haz
    f1 <- c(diff(data_gen_model_trt$F), 0)
    f2 <- c(diff(data_gen_model_ctrl$F), 0)
    f <- (1/(N_trt+N_ctrl))*(N_trt*f1 + N_ctrl*f2)

    #AHRoc
    true_avg_HR_fun0 <- data_gen_model_ctrl$S * c(diff(data_gen_model_trt$F), 0)
    true_avg_HR_fun1 <- data_gen_model_trt$S * c(diff(data_gen_model_ctrl$F), 0)


    rmst_ahr <- purrr::imap(cutoff, function(cutoff, label){
      ind <- (data_gen_model_trt$t <= cutoff)

      Int0 <- sum(true_avg_HR_fun0[ind])
      Int1 <- sum(true_avg_HR_fun1[ind])

      log_quot <- suppressWarnings({
        log(data_gen_model_trt$haz[ind] / data_gen_model_ctrl$haz[ind])
      })

      if(any(is.nan(log_quot))){
        warning("NaN values in log hazard ratio, check cutoff value.")
      }

      tmp <- data.frame(
        rmst_trt            = sum(data_gen_model_trt$S[ind]),
        rmst_ctrl           = sum(data_gen_model_ctrl$S[ind]),
        gAHR                = exp(sum( (log_quot * f[ind]) , na.rm=TRUE)),
        AHR                 = sum(((data_gen_model_trt$haz[ind]/h[ind]) * f[ind]), na.rm=TRUE) /
          sum(((data_gen_model_ctrl$haz[ind]/h[ind]) * f[ind]), na.rm=TRUE),
        AHRoc = Int0/Int1,
        AHRoc_robust = Int0/Int1
      )

      names(tmp) <- paste0(names(tmp), "_", label)
      tmp
    }) |>
      unname() |>
      do.call(what=cbind)

    res <- cbind(res, rmst_ahr)
  }

  if(!is.null(milestones)){

    if(is.null(names(milestones))){
      names(milestones) <- as.character(milestones)
    } else {
      names(milestones)[names(milestones) %in% c("", NA_character_)] <- as.character(milestones[names(milestones) %in% c("", NA_character_)])
    }

    milestones <- purrr::imap(milestones, function(v, label){
      tmp <- data.frame(
        milestone_survival_trt  = tail(data_gen_model_trt$S[data_gen_model_trt$t < v], 1),
        milestone_survival_ctrl = tail(data_gen_model_ctrl$S[data_gen_model_ctrl$t < v], 1)
      )

      names(tmp) <- paste0(names(tmp), "_", label)
      tmp
    }) |>
      unname() |>
      do.call(what=cbind)

    res <- cbind(
      res,
      milestones
    )
  }

  res
}



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
fast_real_statistics <- function(
    haz_trt,  pdf_trt,  surv_trt, quant_trt,
    haz_ctrl, pdf_ctrl, surv_ctrl, quant_ctrl,
    N_trt=1, N_ctrl=1, cutoff=NULL, milestones=NULL
){

  res <- data.frame(
    median_survival_trt = quant_trt(0.5),
    median_survival_ctrl= quant_ctrl(0.5)
  )

  # if cutoff time for rmst and gahr are given
  if(!is.null(cutoff)){

    if(is.null(names(cutoff))){
      names(cutoff) <- as.character(cutoff)
    } else {
      names(cutoff)[names(cutoff) %in% c("", NA_character_)] <- as.character(cutoff[names(cutoff) %in% c("", NA_character_)])
    }

    # define helper functions for the calculation of the average hazard ratios
    h  <- \(t){haz_trt(t)+haz_ctrl(t)}
    f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*pdf_trt(t) + N_ctrl*pdf_ctrl(t))}

    # rectangle rule (more robust in some settings where integrate fails)
    myint<-function(f,lower,upper,steps=1000) {
      delta<-(upper-lower)/steps
      x<-seq(lower+delta/2,upper-delta/2,delta)
      sum(f(x)*delta)
    }

    rmst_ahr <- purrr::imap(cutoff, function(cutoff, label){
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
          integrate(\(t){log(haz_trt(t) / haz_ctrl(t)) * f(t) }, 0, cutoff, stop.on.error = FALSE)$value /
            integrate(\(t){f(t) }, 0, cutoff, stop.on.error = FALSE)$value
          ),
        AHR = integrate(\(t){(haz_trt(t)/h(t)) * f(t) }, 0, cutoff, stop.on.error = FALSE)$value /
          integrate(\(t){(haz_ctrl(t)/h(t)) * f(t)}, 0, cutoff, stop.on.error = FALSE)$value,
        #### Add correctly weighted AHR here
        AHRoc = AHRoc_integrate,
        AHRoc_robust = AHRoc_myint
        )

      names(tmp) <- paste0(names(tmp), "_", label)
      tmp
    }) |>
      unname() |>
      do.call(what=cbind)

    res <- cbind(res, rmst_ahr)
  }

  # if times for milestone survival are given calculate milestone survival at given times
  if(!is.null(milestones)){

    if(is.null(names(milestones))){
      names(milestones) <- as.character(milestones)
    } else {
      names(milestones)[names(milestones) %in% c("", NA_character_)] <- as.character(milestones[names(milestones) %in% c("", NA_character_)])
    }

    milestones <- purrr::imap(milestones, function(v, label){
      tmp <- data.frame(
        milestone_survival_trt  = surv_trt(v),
        milestone_survival_ctrl = surv_ctrl(v)
      )

      names(tmp) <- paste0(names(tmp), "_", label)
      tmp
    }) |>
      unname() |>
      do.call(what=cbind)

    res <- cbind(
      res,
      milestones
    )
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

