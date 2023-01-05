subpop_pzhaz_simnph_internal <- function(Tint, lambda1, lambda2, lambdaProg, int_control){

  cumhaz_tod <- fast_cumhaz_fun(Tint, lambda1)
  cumhaz_tod_nach_switch <- fast_cumhaz_fun(Tint, lambda2)
  pdf_switch <- fast_pdf_fun(Tint, lambdaProg)
  haz_tod <- fast_haz_fun(Tint, lambda1)
  haz_tod_nach_switch <- fast_haz_fun(Tint, lambda2)

  Smix_ <- function(v){
    int_fun <- function(s, t){
      t1 = cumhaz_tod(pmin(s, t))
      t2 = (cumhaz_tod_nach_switch(t - s)) * (s < t)
      exp(-(t1 + t2)) * pdf_switch(s)
    }

    mapply(
      function(x){
        integrate(
          f=int_fun,
          lower=0,
          upper=Inf,
          t=x,
          rel.tol = int_control$rel.tol, abs.tol = int_control$abs.tol
        )$value
      },
      v
    )
  }

  Smix <- memoise::memoise(Smix_)

  hmix <- function(v){
    int_fun <- function(s,t){
      t1 = cumhaz_tod(pmin(s,t))
      t2 = cumhaz_tod_nach_switch(t-s) * (s < t)

      exp(-(t2+t1)) * pdf_switch(s) * haz_tod(t) * (s > t) +
        haz_tod_nach_switch(t-s) * (s < t)
    }

    mapply(
      function(x){
        invS <- (1/Smix(x))
        Ss = integrate(
          f=int_fun,
          lower=0,
          upper=Inf,
          t=x,
          rel.tol = int_control$rel.tol, abs.tol = int_control$abs.tol
        )$value

        invS * Ss
      },
      v
    )
  }

  cummixhaz <- function(v){
    -log(Smix(v))
  }

  Fmix <- function(v){
    1 - Smix(v)
  }

  out <- list(haz_f = hmix, cumhaz_f = cummixhaz, surv_f = Smix,
              cdf_f = Fmix)
  out
}



subpop_pchaz_simnph <- function (Tint, lambda1, lambda2, lambdaProg,
          int_control = list(rel.tol = .Machine$double.eps^0.4, abs.tol = 1e-09), t)
{
  call <- match.call()
  if (any(diff(Tint) <= 0)) {
    stop("Tint should be non-decreasing")
  }

  if (!all(c("rel.tol", "abs.tol") %in% names(int_control))) {
    int_control = list(rel.tol = .Machine$double.eps^0.4,
                       abs.tol = 1e-09)
  }
  if (all(lambdaProg == Inf)) {
    funs = internal_pchaz(Tint, lambda2)
  }
  else if (all(lambdaProg == 0)) {
    funs = internal_pchaz(Tint, lambda1)
  }
  else {
    funs = subpop_pzhaz_simnph_internal(Tint, lambda1, lambda2,
                                 lambdaProg, int_control)
  }
  haz = funs$haz_f(t)
  cumhaz = funs$cumhaz_f(t)
  S = funs$surv_f(t)
  F = funs$cdf_f(t)
  out <- list(haz = haz, cumhaz = cumhaz, S = S, F = F, t = t,
              Tint = Tint, lambda1 = lambda1, lambda2 = lambda2, lambdaProg = lambdaProg,
              timezero = FALSE, call = call, funs = funs)
  class(out) <- c("mixpch", "subpop")
  out
}

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
