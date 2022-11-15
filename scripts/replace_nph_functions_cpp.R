library(nph)
library(Rcpp)

sourceCpp("C/pch_functions.cpp")

new_haz_fun    <- function(Tint, lambda){
  function(v){
    hazFunCpp(Tint, lambda, v)
  }
}

new_cumhaz_fun <- function(Tint, lambda){
  function(v){
    cumhazFunCpp(Tint, lambda, v)
  }
}

new_cdf_fun    <- function(Tint, lambda){
  function(v){
    cdfFunCpp(Tint, lambda, v)
  }
}

new_pdf_fun    <- function(Tint, lambda){
  function(v){
    pdfFunCpp(Tint, lambda, v)
  }
}

new_surv_fun   <- function(Tint, lambda){
  function(v){
    survFunCpp(Tint, lambda, v)
  }
}



assignInNamespace("haz_fun"   , new_haz_fun   , ns="nph")
assignInNamespace("cumhaz_fun", new_cumhaz_fun, ns="nph")
assignInNamespace("cdf_fun"   , new_cdf_fun   , ns="nph")
assignInNamespace("pdf_fun"   , new_pdf_fun   , ns="nph")
assignInNamespace("surv_fun"  , new_surv_fun  , ns="nph")

