# This is an adapted version of the pop_pchaz function from the nph package
# this function also returns the funs object that is included in pchaz but not
# in pop_pchaz in the nph package.
#
# nph is GPL3, since we rely heavily on the nph package we should licence this
# package under GPL3 anyways, so copying the code should not be a problem.
#
# Tobias is in contact with the author of nph to find out:
#   1. if the code here is correct
#   2. if the changes can be adapted in the nph package instead of this
#      workaround
alternative_pop_pchaz <- function (Tint, lambdaMat1, lambdaMat2, lambdaProgMat, p, timezero = FALSE,
                          int_control = list(rel.tol = .Machine$double.eps^0.4, abs.tol = 1e-09),
                          discrete_approximation = FALSE)
{
  call <- match.call()
  if (!identical(dim(as.matrix(lambdaMat1)), dim(as.matrix(lambdaMat2)),
                 dim(as.matrix(lambdaProgMat)))) {
    stop("Dimensions of all lambda matrices should be the same")
  }
  if (any(diff(Tint) <= 0)) {
    stop("Tint should be non-decreasing")
  }
  if (any(c(lambdaMat1, lambdaMat2, lambdaProgMat) < 0)) {
    stop("Only positive values are allowed for the hazard rates")
  }
  if ((length(Tint) - 1) != ncol(as.matrix(lambdaMat1))) {
    stop("The length of Tint should be equal to the number of column of the lambda matrices + 1")
  }
  if (length(p) != nrow(as.matrix(lambdaMat1))) {
    stop("The length of p should be equal to the number of column of the lambda matrices")
  }
  if (abs(sum(p) - 1) > .Machine$double.eps^0.5) {
    stop("The prevalences in p should sum up to 1")
  }
  if (length(timezero) == 1) {
    timezero = rep(timezero, length(p))
  }
  if (length(timezero) != length(p)) {
    stop("timezero should be of length 1 or the same length as p")
  }
  if (discrete_approximation) {
    out = nph:::pop_hazVFfun(Tint, lambdaMat1, lambdaMat2, lambdaProgMat,
                       p, timezero)
    return(out)
  }
  if (!all(c("rel.tol", "abs.tol") %in% names(int_control))) {
    int_control = list(rel.tol = .Machine$double.eps^0.4,
                       abs.tol = 1e-09)
  }

  t <- seq(0, max(Tint) - 1, 1)
  nsubpop <- length(p)
  Fmix <- rep(0, length(t))

  # this is a list of lists
  # for each subpopulation a the list contains the list of all functions
  all_funs <- lapply(1:nsubpop, \(i){
    lambdaProg = lambdaProgMat[i, ]
    lambda1 = lambdaMat1[i, ]
    lambda2 = lambdaMat2[i, ]
    if (all(lambdaProg == Inf)) {
      funs = nph:::internal_pchaz(Tint, lambda2)
    }
    else if (all(lambdaProg == 0)) {
      funs = nph:::internal_pchaz(Tint, lambda1)
    }
    else {
      funs = nph:::internal_subpop_pchaz(Tint, lambda1, lambda2,
                                         lambdaProg, timezero[i], int_control)
    }

    funs
  })

  # combine the functions of the subpopulations
  # to get the functions of the whole population

  # returns a function of t
  Fmix_f <- function(t) {
    vapply(all_funs, \(funs_i){funs_i$cdf_f(t)}, numeric(length(t))) %*% p |>
      as.numeric()
  }

  Smix_f <- \(t) 1-Fmix_f(t)
  cummixhaz_f <- \(t) -log(Smix_f(t))

  # returns a function of t
  fmix_f <- function(t){
    vapply(all_funs, \(funs_i){funs_i$haz_f(t)*funs_i$surv_f(t)}, numeric(length(t))) %*% p |>
      as.numeric()
  }

  mixhaz_f <- \(t) fmix_f(t)/Smix_f(t)

  Fmix <- Fmix_f(t)
  Smix <- 1-Fmix
  cummixhaz <- -log(Smix)
  mixhaz <- c(diff(cummixhaz), NA_real_)

  funs <- list(
    haz_f = mixhaz_f,
    cumhaz_f = cummixhaz_f,
    surv_f = Smix_f,
    cdf_f = Fmix_f,
    pdf_f = fmix_f
  )

  out <- list(haz = mixhaz, cumhaz = cummixhaz, S = Smix,
              F = Fmix, t = t, Tint = Tint, lambdaMat1 = lambdaMat1,
              lambdaMat2 = lambdaMat2, lambdaProgMat = lambdaProgMat,
              p = p, timezero = timezero, call = call, funs=funs)

  out
}


