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
