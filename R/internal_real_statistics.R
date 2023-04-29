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

    rmst_ahr <- purrr::imap(cutoff, function(cutoff, label){

      tmp <- data.frame(
        rmst_trt  = sum(data_gen_model_trt$S[data_gen_model_trt$t <= cutoff]),
        rmst_ctrl = sum(data_gen_model_ctrl$S[data_gen_model_ctrl$t <= cutoff]),

        gAHR_schemper_cox  = gahr_schemper_cox_discrete(data_gen_model_trt, data_gen_model_ctrl, N_trt, N_ctrl, cutoff),
        gAHR_schemper_wcox = gahr_schemper_wcox_discrete(data_gen_model_trt, data_gen_model_ctrl, N_trt, N_ctrl, cutoff),
        AHR_schemper_cox  = ahr_schemper_cox_discrete(data_gen_model_trt, data_gen_model_ctrl, N_trt, N_ctrl, cutoff),
        AHR_schemper_wcox = ahr_schemper_wcox_discrete(data_gen_model_trt, data_gen_model_ctrl, N_trt, N_ctrl, cutoff),
        AHR_rauch         = ahrs_rauch_discrete(data_gen_model_trt, data_gen_model_ctrl, cutoff, gamma=1),
        gAHR_rauch         = ahrs_rauch_discrete(data_gen_model_trt, data_gen_model_ctrl, cutoff, gamma=0.5)
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

    milestones <- purrr:::imap(milestones, function(v, label){
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

    rmst_ahr <- purrr::imap(cutoff, function(cutoff, label){
      tmp <- data.frame(
        rmst_trt = integrate(surv_trt, 0, cutoff)$value,
        rmst_ctrl = integrate(surv_ctrl, 0, cutoff)$value,
        gAHR_schemper_cox  = gahr_schemper_cox(haz_trt, pdf_trt, haz_ctrl, pdf_ctrl, N_trt, N_ctrl, cutoff),
        gAHR_schemper_wcox = gahr_schemper_wcox(haz_trt, pdf_trt, surv_trt, haz_ctrl, pdf_ctrl, surv_ctrl, N_trt, N_ctrl, cutoff),
        AHR_schemper_cox  = ahr_schemper_cox(haz_trt, pdf_trt, haz_ctrl, pdf_ctrl, N_trt, N_ctrl, cutoff),
        AHR_schemper_wcox = ahr_schemper_wcox(haz_trt, pdf_trt, surv_trt, haz_ctrl, pdf_ctrl, surv_ctrl, N_trt, N_ctrl, cutoff),
        AHR_rauch         = ahrs_rauch(haz_trt, pdf_trt, surv_trt, haz_ctrl, pdf_ctrl, surv_ctrl, cutoff, gamma=1),
        gAHR_rauch         = ahrs_rauch(haz_trt, pdf_trt, surv_trt, haz_ctrl, pdf_ctrl, surv_ctrl, cutoff, gamma=0.5)
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

    milestones <- purrr:::imap(milestones, function(v, label){
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

# helper functions for different AHRs -------------------------------------

gahr_schemper_cox <- function(
    haz_trt,  pdf_trt,
    haz_ctrl, pdf_ctrl,
    N_trt, N_ctrl,
    cutoff
){
  h  <- \(t){haz_trt(t)+haz_ctrl(t)}
  f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*pdf_trt(t) + N_ctrl*pdf_ctrl(t))}

  ahr <- exp(
    integrate(\(t){log(haz_trt(t) / haz_ctrl(t)) * f(t) }, 0, cutoff, stop.on.error = FALSE)$value /
      integrate(\(t){f(t) }, 0, cutoff, stop.on.error = FALSE)$value
  )

  ahr
}

ahr_schemper_cox <- function(
    haz_trt,  pdf_trt,
    haz_ctrl, pdf_ctrl,
    N_trt, N_ctrl,
    cutoff
){
  h  <- \(t){haz_trt(t)+haz_ctrl(t)}
  f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*pdf_trt(t) + N_ctrl*pdf_ctrl(t))}

  ahr <- integrate(\(t){(haz_trt(t)/h(t)) * f(t) }, 0, cutoff, stop.on.error = FALSE)$value /
    integrate(\(t){(haz_ctrl(t)/h(t)) * f(t)}, 0, cutoff, stop.on.error = FALSE)$value

  ahr
}

gahr_schemper_wcox <- function(
    haz_trt,  pdf_trt,  surv_trt,
    haz_ctrl, pdf_ctrl, surv_ctrl,
    N_trt, N_ctrl,
    cutoff
){
  h  <- \(t){haz_trt(t)+haz_ctrl(t)}
  f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*pdf_trt(t) + N_ctrl*pdf_ctrl(t))}
  w  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*surv_trt(t) + N_ctrl*surv_ctrl(t))}

  ahr <- exp(
    integrate(\(t){log(haz_trt(t) / haz_ctrl(t)) * f(t) * w(t)}, 0, cutoff, stop.on.error = FALSE)$value /
      integrate(\(t){f(t) * w(t)}, 0, cutoff, stop.on.error = FALSE)$value
  )

  ahr
}

ahr_schemper_wcox <- function(
    haz_trt,  pdf_trt,  surv_trt,
    haz_ctrl, pdf_ctrl, surv_ctrl,
    N_trt, N_ctrl,
    cutoff
){
  h  <- \(t){haz_trt(t)+haz_ctrl(t)}
  f  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*pdf_trt(t) + N_ctrl*pdf_ctrl(t))}
  w  <- \(t){(1/(N_trt+N_ctrl))*(N_trt*surv_trt(t) + N_ctrl*surv_ctrl(t))}

  ahr <- integrate(\(t){(haz_trt(t)/h(t)) * f(t) * w(t)}, 0, cutoff, stop.on.error = FALSE)$value /
    integrate(\(t){(haz_ctrl(t)/h(t)) * f(t) * w(t)}, 0, cutoff, stop.on.error = FALSE)$value

  ahr
}

# gamma=0.5: gAHR
# gamma=1     AHR
ahrs_rauch <- function(
    haz_trt,  pdf_trt,  surv_trt,
    haz_ctrl, pdf_ctrl, surv_ctrl,
    cutoff, gamma
){
  g <- function(t){
    gamma * surv_trt(t)^(gamma-1) * surv_ctrl(t)^gamma * (-pdf_trt(t)) +
      gamma * surv_trt(t)^gamma * surv_ctrl(t)^(gamma-1) * (-pdf_ctrl(t))
  }

  numerator <- function(t){
    (haz_trt(t) / (haz_trt(t)+haz_ctrl(t))) * g(t)
  }

  denominator <- function(t){
    (haz_ctrl(t) / (haz_trt(t)+haz_ctrl(t))) * g(t)
  }

  ahr <- integrate(numerator,   0, cutoff)$value /
    integrate(denominator, 0, cutoff)$value

  ahr
}


# helper functions for different AHRs, discrete ---------------------------


gahr_schemper_cox_discrete <- function(
    data_gen_model_trt, data_gen_model_ctrl, N_trt=1, N_ctrl=1, cutoff=NULL
){
  h <- data_gen_model_trt$haz + data_gen_model_ctrl$haz
  f1 <- c(diff(data_gen_model_trt$F), 0)
  f2 <- c(diff(data_gen_model_ctrl$F), 0)
  f <- (1/(N_trt+N_ctrl))*(N_trt*f1 + N_ctrl*f2)

  ahr <- exp(sum(
    (log(data_gen_model_trt$haz /
           data_gen_model_ctrl$haz) * f)[data_gen_model_trt$t <= cutoff] ,
    na.rm=TRUE))

  ahr
}

ahr_schemper_cox_discrete <- function(
    data_gen_model_trt, data_gen_model_ctrl, N_trt=1, N_ctrl=1, cutoff=NULL
){
  h <- data_gen_model_trt$haz + data_gen_model_ctrl$haz
  f1 <- c(diff(data_gen_model_trt$F), 0)
  f2 <- c(diff(data_gen_model_ctrl$F), 0)
  f <- (1/(N_trt+N_ctrl))*(N_trt*f1 + N_ctrl*f2)

  ahr <- sum(((data_gen_model_trt$haz/h) * f)[data_gen_model_trt$t <= cutoff], na.rm=TRUE) /
    sum(((data_gen_model_ctrl$haz/h) * f)[data_gen_model_ctrl$t <= cutoff], na.rm=TRUE)

  ahr
}

gahr_schemper_wcox_discrete <- function(
    data_gen_model_trt, data_gen_model_ctrl, N_trt=1, N_ctrl=1, cutoff=NULL
){
  h <- data_gen_model_trt$haz + data_gen_model_ctrl$haz
  f1 <- c(diff(data_gen_model_trt$F), 0)
  f2 <- c(diff(data_gen_model_ctrl$F), 0)
  f <- (1/(N_trt+N_ctrl))*(N_trt*f1 + N_ctrl*f2)
  w <- (1/(N_trt+N_ctrl))*(N_trt*data_gen_model_trt$S + N_ctrl*data_gen_model_ctrl$S)

  ahr <- exp(sum(
    (log(data_gen_model_trt$haz /
           data_gen_model_ctrl$haz) * f * w)[data_gen_model_trt$t <= cutoff] ,
    na.rm=TRUE))

  ahr
}

ahr_schemper_wcox_discrete <- function(
    data_gen_model_trt, data_gen_model_ctrl, N_trt=1, N_ctrl=1, cutoff=NULL
){
  h <- data_gen_model_trt$haz + data_gen_model_ctrl$haz
  f1 <- c(diff(data_gen_model_trt$F), 0)
  f2 <- c(diff(data_gen_model_ctrl$F), 0)
  f <- (1/(N_trt+N_ctrl))*(N_trt*f1 + N_ctrl*f2)
  w <- (1/(N_trt+N_ctrl))*(N_trt*data_gen_model_trt$S + N_ctrl*data_gen_model_ctrl$S)

  ahr <- sum(((data_gen_model_trt$haz/h) * f * w)[data_gen_model_trt$t <= cutoff], na.rm=TRUE) /
    sum(((data_gen_model_ctrl$haz/h) * f * w)[data_gen_model_ctrl$t <= cutoff], na.rm=TRUE)

  ahr
}

# gamma=0.5: gAHR
# gamma=1     AHR
ahrs_rauch_discrete <- function(
    data_gen_model_trt, data_gen_model_ctrl, cutoff=NULL, gamma
){

  f_trt <- c(diff(data_gen_model_trt$F), 0)
  f_ctrl <- c(diff(data_gen_model_ctrl$F), 0)

  g <- (
    gamma * data_gen_model_trt$S^(gamma-1) * data_gen_model_ctrl$S^gamma * (-f_trt) +
      gamma * data_gen_model_trt$S^gamma * data_gen_model_ctrl$S^(gamma-1) * (-f_ctrl)
  )

  numerator <- (
    (data_gen_model_trt$haz / (data_gen_model_trt$haz + data_gen_model_ctrl$haz)) * g
  )

  denominator <- (
    (data_gen_model_ctrl$haz / (data_gen_model_trt$haz + data_gen_model_ctrl$haz)) * g
  )

  ahr <- sum(numerator[data_gen_model_ctrl$t <= cutoff], na.rm=TRUE) /
    sum(denominator[data_gen_model_ctrl$t <= cutoff], na.rm=TRUE)

  ahr
}

