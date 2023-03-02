# run nphparams with trycatch, used in various analyse functions
trycatch_nphparams <- function(call){
  tryCatch(
    call,
    error = function(e){
      warning(e$message)
      list(
        tab=list(
          p_unadj=NA_real_,
          lwr_unadj=NA_real_,
          upr_unadj=NA_real_,
          Estimate=NA_real_
        )
      )
    })
}

# hazard ratio required, inverted Schönfeld sample size formula
hr_required_schoenfeld <- function(Nevt, alpha=0.05, beta=0.2, p=0.5){
  exp( (qnorm(beta) + qnorm(alpha) ) / sqrt(p*(1-p)*Nevt) )
}
