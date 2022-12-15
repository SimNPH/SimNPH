library(SimNPH)
library(SimDesign)
library(parallel)

N_sim <- 10

# setup cluster -----------------------------------------------------------

n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(SimDesign)
  library(SimNPH)
  library(parallel)
})


# setup data generation ---------------------------------------------------

# load parameters
design <- read.table("data/parameters/delayed_effect_2022-12-15.csv", sep=",", dec=".", header=TRUE)

design <- design[sample(1:nrow(design), 1), ]

# define generator
my_generator <- function(condition, fixed_objects=NULL){
  generate_delayed_effect(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_events(condition$final_events)
}

alpha <- 0.05
nominal_alpha <- ldbounds::ldBounds(c(0.5,1))$nom.alpha

clusterExport(cl, "nominal_alpha")


# define analyse functions ------------------------------------------------

# TODO: weighted cox
# TODO: piecewise exponential

my_analyse <- list(
  # estimation
  ahr  = analyse_ahr(type="AHR"),
  gahr = analyse_ahr(type="gAHR"),
  mean_surv = analyse_diff_mean_survival(),
  milestone = analyse_milestone_survival(times=c(8, 12)),
  rmst = analyse_rmst_diff(),
  cox = analyse_coxph(),
  # tests
  peto_peto = analyse_gehan_wilcoxon(),
  fh = analyse_logrank_fh_weights(rho=0, gamma=1),
  logrank = analyse_logrank(),
  max_combo = analyse_maxcombo(),
  modest_6 = analyse_modelstly_weighted(t_star=6),
  modest_8 = analyse_modelstly_weighted(t_star=8),
  # tests group sequential
  peto_peto_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_gehan_wilcoxon()
  ),
  fh_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=0, gamma=1)
  ),
  logrank_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank()
  ),
  max_combo_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_maxcombo()
  ),
  modest_gs_6 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_modelstly_weighted(t_star=6)
  ),
  modest_gs_8 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_modelstly_weighted(t_star=8)
  )
)


# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  ahr = summarise_estimator(est=AHR, real=0),
  gahr = summarise_estimator(est=gAHR, real=0),
  mean_surv = summarise_estimator(est=diff_Q, real=0),
  milestone = summarise_estimator(est=milestone_surv_ratio[1], real=0, name="milestone_1"),
  milestone = summarise_estimator(est=milestone_surv_ratio[2], real=0, name="milestone_2"),
  rmst = summarise_estimator(est=rmst_diff, real=0),
  cox = summarise_estimator(est=hr, lower = hr_lower, upper=hr_upper),
  peto_peto = summarise_test(alpha),
  fh = summarise_test(alpha),
  logrank = summarise_test(alpha),
  max_combo = summarise_test(alpha),
  modest_6 = summarise_test(alpha),
  modest_8 = summarise_test(alpha),
  peto_peto_gs = summarise_group_sequential(),
  fh_gs = summarise_group_sequential(),
  logrank_gs = summarise_group_sequential(),
  max_combo_gs = summarise_group_sequential(),
  modest_gs_6 = summarise_group_sequential(),
  modest_gs_8 = summarise_group_sequential()
)


# run ---------------------------------------------------------------------

save_folder <- paste0(paste0("data/simulation_delayed_effect_", Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S")))
dir.create(save_folder)

results <- runSimulation(
  design,
  replications = N_sim,
  generate = my_generator,
  analyse = my_analyse,
  summarise = my_summarise,
  cl = cl,
  save_details = list(
    out_rootdir = save_folder
  )
)

saveRDS(results, paste0(save_folder, "/results.Rds"))

