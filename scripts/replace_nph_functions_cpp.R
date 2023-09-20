library(nph)
library(SimNPH)

assignInNamespace("haz_fun"   , miniPCH::hpch_fun,  ns="nph")
assignInNamespace("cumhaz_fun", miniPCH::chpch_fun, ns="nph")
assignInNamespace("cdf_fun"   , miniPCH::dpch_fun,  ns="nph")
assignInNamespace("pdf_fun"   , miniPCH::dpch_fun,  ns="nph")
assignInNamespace("surv_fun"  , miniPCH::spch_fun,  ns="nph")

