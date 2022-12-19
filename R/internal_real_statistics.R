internal_real_statistics_pchaz <- function(data_gen_model_trt, data_gen_model_ctrl, N_trt=1, N_ctrl=1, cutoff=1000){

  real_stats_one_group <- function(group, cutoff, milestones){
    res <- data.frame(
      rmst      = integrate(group$funs$surv_f, 0, cutoff)$value,
      median_survival = if(group$funs$surv_f(max(group$t, cutoff)) < 0.5) {
        uniroot(\(x) group$funs$surv_f(x) - 0.5, lower=0, upper=12, extendInt = "downX")$root
      } else {
        Inf
      }
    )

    res
  }

  real_stats_trt  <- real_stats_one_group(data_gen_model_trt, cutoff = cutoff, milestones = milestones)
  real_stats_ctrl <- real_stats_one_group(data_gen_model_ctrl, cutoff = cutoff, milestones = milestones)

  # helper functions for average hazard ratios
  # notation corresponds to Schemper 2009
  h1 <- data_gen_model_trt$funs$haz_f
  h2 <- data_gen_model_ctrl$funs$haz_f
  h  <- \(t){h1(t)+h2(t)}
  f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*data_gen_model_trt$funs$pdf_f(t) + N_ctrl*data_gen_model_ctrl$funs$pdf_f(t))}
  w  <- \(t){1} # Cox
  # w  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*group_1$funs$surv_f(t) + N_ctrl*group_2$funs$surv_f(t))} # WCox

  gAHR <- exp(integrate(\(t){log(h1(t) / h2(t)) * f(t) * w(t)}, 0, cutoff)$value)
  AHR  <- integrate(\(t){(h1(t)/h(t)) * f(t) * w(t)}, 0, cutoff)$value / integrate(\(t){(h2(t)/h(t)) * f(t) * w(t)}, 0, cutoff)$value


  names(real_stats_trt) <- paste0(names(real_stats_trt), "_trt")
  names(real_stats_ctrl) <- paste0(names(real_stats_ctrl), "_ctrl")

  cbind(real_stats_trt, real_stats_ctrl, gAHR=gAHR, AHR=AHR)
}


fast_real_statistics <- function(
    haz_trt,  pdf_trt,  surv_trt, quant_trt,
    haz_ctrl, pdf_ctrl, surv_ctrl, quant_ctrl,
    N_trt=1, N_ctrl=1, cutoff=1000, milestones=NULL
){

  h  <- \(t){haz_trt(t)+haz_ctrl(t)}
  f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*pdf_trt(t) + N_ctrl*pdf_ctrl(t))}

  res <- data.frame(
    rmst_trt            = integrate(surv_trt, 0, cutoff)$value,
    median_survival_trt = quant_trt(0.5),
    rmst_ctrl           = integrate(surv_ctrl, 0, cutoff)$value,
    median_survival_ctrl= quant_ctrl(0.5),
    gAHR                = exp(integrate(\(t){log(haz_trt(t) / haz_ctrl(t)) * f(t) }, 0, cutoff)$value),
    AHR                 = integrate(\(t){(haz_trt(t)/h(t)) * f(t) }, 0, cutoff)$value /
      integrate(\(t){(haz_ctrl(t)/h(t)) * f(t)}, 0, cutoff)$value
  )

  if(!is.null(milestones)){
    if(is.null(names(milestones))){
      my_colnames <- paste0("milestone_surv_", milestones)
    } else {
      my_colnames <- names(milestones)
    }

    milestones_trt <- milestones |>
      surv_trt() |>
      setNames(paste0(my_colnames, "_trt")) |>
      t() |>
      as.data.frame()

    milestones_ctrl <- milestones |>
      surv_ctrl() |>
      setNames(paste0(my_colnames, "_ctrl")) |>
      t() |>
      as.data.frame()

    res <- cbind(
      res,
      milestones_trt,
      milestones_ctrl
    )

  }

  res
}

fast_real_statistics_pchaz <- function(
    Tint_trt,  lambda_trt,
    Tint_ctrl, lambda_ctrl,
    N_trt=1, N_ctrl=1, cutoff=1000, milestones=NULL
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

