library(SimNPH)
library(nph)

# Options -----------------------------------------------------------------

options <- expand.grid(
  n_pat = c(300, 500, 1000, 1500)
) |>
  within({
    n_trt <- n_pat / 2
    n_ctrl <- n_pat / 2
    final_events <- ceiling(n_pat * 0.75)
  })


# Assumptions -------------------------------------------------------------

assumptions <- expand.grid(
  hazard_ctrl = nph::m2r(c(36, 12, 6)),
  effect_size_ph = c(0, 0.5, 0.8, 0.9)
)


# Merge Options and Assumptions -------------------------------------------

design <- merge(
  options,
  assumptions,
  by=NULL
)

# Schoenfeld Formula ------------------------------------------------------

# hazard ratio required, inverted Schönfeld sample size formula
hr_required_schoenfeld <- function(Nevt, alpha=0.05, beta=0.2, p=0.5){
  exp( (qnorm(beta) + qnorm(alpha)) / sqrt(p*(1-p)*Nevt) )
}


# calculate hr and medians ------------------------------------------------

target_alpha <- 0.05

get_medians <- function(condition){

  if(condition$effect_size_ph == 0){
    ph_hr <- 1
  } else {
    ph_hr <- hr_required_schoenfeld(
      condition$final_events,
      alpha=target_alpha,
      beta=(1-condition$effect_size_ph),
      p=(condition$n_ctrl/(condition$n_ctrl + condition$n_trt)) # allocation ratio
    )
  }

  condition$hr_ph <- ph_hr
  condition$median_trt  <- log(2) / (condition$hazard_ctrl * ph_hr)
  condition$median_ctrl <- log(2) / (condition$hazard_ctrl)

  condition
}

design <- design |>
  split(1:nrow(design)) |>
  lapply(get_medians) |>
  do.call(what=rbind)


# save --------------------------------------------------------------------

filename <- paste0("data/parameters/joint_models_", format(Sys.Date(), "%Y-%m-%d"), ".csv")
write.table(design, file=filename, quote=FALSE, sep=", ", dec=".", row.names = FALSE, col.names = TRUE)

