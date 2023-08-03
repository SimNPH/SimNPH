# Tint: time intervals
# lambda1: hazard for death before progression
# lambda2: hazard for death after progression
# lambdaProg: hazard for progression
subpop_hazVfun_simnph <- function (Tint, lambda1, lambda2, lambdaProg,
                                   timezero = FALSE){
  stopifnot(timezero)

  t_max <- max(Tint)
  Tint_ <- Tint[-length(Tint)]

  t <- seq(0, t_max-1, 1)
  t_ <- seq(1, t_max, 1)
  nt <- length(t)

  f_tod             <-    fast_pdf_fun(Tint_,    lambda1)(t_)
  S_tod             <-   fast_surv_fun(Tint_,    lambda1)(t_)
  f_tod_nach_switch <-    fast_pdf_fun(Tint_,    lambda2)(t_)
  S_switch          <-   fast_surv_fun(Tint_, lambdaProg)(t)
  f_switch          <-    fast_pdf_fun(Tint_, lambdaProg)(t)

  f_switch_vor_tod <- f_switch * S_tod
  f_tod_vor_switch <- f_tod * S_switch

  f_switch_und_tod <- convolve(f_switch_vor_tod, rev(f_tod_nach_switch), type="open")[1:nt]
  f_mix <- f_switch_und_tod + f_tod_vor_switch

  F_mix <- cumsum(c(0,f_mix[-nt]))
  S_mix <- 1-F_mix
  cummixhaz <- -log(S_mix)
  mixhaz <- f_mix/S_mix

  out <- list(haz = mixhaz, cumhaz = cummixhaz, S = S_mix,
              F = F_mix, t = t, Tint = Tint, lambda1 = lambda1, lambda2 = lambda2,
              lambdaProg = lambdaProg, timezero = TRUE)
  class(out) <- c("mixpch", "subpop")
  out

}
