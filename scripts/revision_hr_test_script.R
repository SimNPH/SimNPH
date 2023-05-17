library(SimNPH)
library(SimDesign)
library(parallel)

if(packageVersion("SimNPH") != "0.3.2"){
  stop("Please run the simulations with the correct vesion of the SimNPH package for reproducability.")
}
## Here we would need to implement the alternative AHR
source('scripts/revision_ahr_functions.R')

# Helper functions --------------------------------------------------------

# months to days
m2d <- \(t) 365.25*t/12

# Options -----------------------------------------------------------------

options <- expand.grid(
  recruitment = m2d(18),#m2d(c(18, 30)),
  n_pat = 1000#c(300, 500, 1000, 1500)
) |>
  within({
    n_trt <- n_pat / 2
    n_ctrl <- n_pat / 2
    final_events <- ceiling(n_pat * 0.75)
    interim_events <- ceiling(final_events * 0.5)
  })


# Assumptions -------------------------------------------------------------

assumptions <- expand.grid(
  hazard_ctrl = nph::m2r(12),#nph::m2r(c(36, 12, 6)),
  censoring_prop = 0.1,#c(0, 0.1, 0.3),
  delay = m2d(seq(0,9, by=4)),
  effect_size_ph = c(0, 0.8)
)


# Merging Options and Assumptions -----------------------------------------

design <- merge(
  options,
  assumptions,
  by=NULL
)


# Calibrating Implicitly defined Parameters ------------------------------

# hazard rate after onset of treatment effect
# calculated such that the median survival is the same as under proportional
# hazards with the given effect size. Fails if onset of treatment effect is
# after the median survival time
design <- design |>
  hr_after_onset_from_PH_effect_size()

# rate of random censoring
design <- design |>
  cen_rate_from_cen_prop_delayed_effect()


# Excluding Scenarios which did not give Reasonable Parameter Value--------

# Excluding Scenarios for which a hazard in the treatment arm could not be calculated.
# This happens when the median of the survival functions is before onset of treatment effect.
design <- design |>
  subset(!is.na(hazard_trt))

# Calculating True Summary Statistics -------------------------------------
rqm <- \(lambda,q=.5) 12*log(1/q)/(365.25*lambda)


design <- design |> #rowwise() |>
  true_summary_statistics_delayed_effect(
    milestones   = m2d(c("6m"=6, "12m"=12)),
    cutoff_stats = m2d(c("6m"=6, "12m"=12, "360m" = 360))
  )

# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  ahr_6m   = analyse_ahr(type="AHR", max_time=m2d(6), alternative = "one.sided"),
  ahr_12m  = analyse_ahr(type="AHR", max_time=m2d(12), alternative = "one.sided"),
  gahr_6m  = analyse_ahr(type="gAHR", max_time=m2d(6), alternative = "one.sided"),
  gahr_12m = analyse_ahr(type="gAHR", max_time=m2d(12), alternative = "one.sided"),
  cox = analyse_coxph(alternative = "one.sided"),
  # tests
  descriptive = analyse_describe()
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  # estimation
  cox_gahr = summarise_estimator(est=hr, real=gAHR_360m, lower=hr_lower, upper=hr_upper, null=1),
  cox_ahr = summarise_estimator(est=hr, real=AHR_360m, lower=hr_lower, upper=hr_upper, null=1),
  ## Turn this on to obtain results for ahroc
  cox_ahroc = summarise_estimator(est=hr, real=AHRoc_360m, lower=hr_lower, upper=hr_upper, null=1),
  ahr_6m  = summarise_estimator(est=AHR, real=AHR_6m, lower=AHR_lower, upper=AHR_upper, null=1),
  ahr_12m = summarise_estimator(est=AHR, real=AHR_12m, lower=AHR_lower, upper=AHR_upper, null=1),
  gahr_6m  = summarise_estimator(est=gAHR, real=gAHR_6m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  gahr_12m = summarise_estimator(est=gAHR, real=gAHR_12m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  descriptive = summarise_describe()
)



N_sim <- 50

# Helper functions --------------------------------------------------------
# months to days
m2d <- \(t) 365.25*t/12

# setup cluster -----------------------------------------------------------

n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(SimDesign)
  library(SimNPH)
  library(parallel)
})


# setup data generation ---------------------------------------------------

# define generator
my_generator <- function(condition, fixed_objects=NULL){
  generate_delayed_effect(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_events(condition$final_events)
}

alpha <- 0.025
nominal_alpha <- ldbounds::ldBounds(c(0.5,1), sides=1, alpha = 0.025)$nom.alpha

clusterExport(cl, "nominal_alpha")


# run ---------------------------------------------------------------------

save_folder <- "tmp_test_hrs"
dir.create(save_folder)

results <- runSimulation(
  design,
  replications = N_sim,
  generate = my_generator,
  analyse = my_analyse,
  summarise = my_summarise,
  seed = design$old_seed,
  cl = cl,
  parallel = TRUE,
  save_details = list(
    out_rootdir = save_folder
  ),
  extra_options = list(
    store_warning_seeds = TRUE,
    allow_na = TRUE
  )
)


