library(nph)
library(SimNPH)

assignInNamespace("haz_fun"   , SimNPH::fast_haz_fun   , ns="nph")
assignInNamespace("cumhaz_fun", SimNPH::fast_cumhaz_fun, ns="nph")
assignInNamespace("cdf_fun"   , SimNPH::fast_cdf_fun   , ns="nph")
assignInNamespace("pdf_fun"   , SimNPH::fast_pdf_fun   , ns="nph")
assignInNamespace("surv_fun"  , SimNPH::fast_surv_fun  , ns="nph")

