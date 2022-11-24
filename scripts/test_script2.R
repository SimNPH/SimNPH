# How To use this test-script
#
# The installation of the packages can prompt for input in the console if some
# packages are not yet installed.
# It's best to run the lines with install.packages and remotes::install_git
# interactively one time and then sourcing the whole file.
#
# Please don't run the rest of the file line by line, but source the whole
# script in one go since the script also measures running time of the
# simulations.
#
# The number of simulations is set in the variable N_sim and is currently 250.
# This is very low and should therefore run quickly. If the script runs without
# error for 250 simulations, please increase N_sim to 2500 and rund the script
# again.

# don't stop with an error if a package masks another function
options(conflicts.policy = NULL)

# installation of packages ------------------------------------------------

# install.packages(c("SimDesign", "remotes", "ldbounds"))
# remotes::install_git("https://github.com/SimNPH/SimNPH.git")

# preparation -------------------------------------------------------------

library(SimDesign)
library(SimNPH)
library(parallel)

N_sim <- 250

n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(SimDesign)
  library(SimNPH)
  library(parallel)
})

times <- numeric(0)

# Test 1 simple tests -----------------------------------------------------

times <- c(times, "start"=Sys.time())

Assumptions <- assumptions_delayed_effect()
Options <- design_group_sequential()
Design <- merge(Assumptions, Options, by=NULL)

my_generator <- function(condition, fixed_objects=NULL){
  generate_delayed_effect(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_time(condition$followup)
}

alpha <- 0.05
nominal_alpha <- ldbounds::ldBounds(c(0.5,1))$nom.alpha

clusterExport(cl, "nominal_alpha")

my_analyse <- list(
  # estimation
  ahr  = analyse_ahr(type="AHR"),
  gahr = analyse_ahr(type="gAHR"),
  mean_surv = analyse_diff_mean_survival(),
  milestone = analyse_milestone_survival(times=c(8, 12)),
  rmst = analyse_rmst_diff(),
  # tests
  peto_peto = analyse_gehan_wilcoxon(),
  fh = analyse_logrank_fh_weights(rho=0, gamma=1),
  logrank = analyse_logrank(),
  max_combo = analyse_maxcombo(),
  # modest = analyse_modelstly_weighted(t_star=c(6, 8)),
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
  ) #,
  # modest_gs = analyse_group_sequential(
  #   followup = c(condition$interim_events, condition$final_events),
  #   followup_type = c("event", "event"),
  #   alpha = nominal_alpha,
  #   analyse_functions = analyse_modelstly_weighted(t_star=c(6, 8))
  # )
)

my_summarise <- create_summarise_function(
  ahr = summarise_estimator(est=AHR, real=0),
  gahr = summarise_estimator(est=gAHR, real=0),
  mean_surv = summarise_estimator(est=diff_Q, real=0),
  milestone = summarise_estimator(est=milestone_surv_ratio[1], real=0, name="milestone_1"),
  milestone = summarise_estimator(est=milestone_surv_ratio[2], real=0, name="milestone_2"),
  rmst = summarise_estimator(est=rmst_diff, real=0),
  peto_peto = summarise_test(alpha),
  fh = summarise_test(alpha),
  logrank = summarise_test(alpha),
  max_combo = summarise_test(alpha),
  # modest = summarise_test(alpha),
  # peto_peto_gs = summarise_group_sequential(),
  fh_gs = summarise_group_sequential(),
  logrank_gs = summarise_group_sequential(),
  max_combo_gs = summarise_group_sequential()#,
  # modest_gs = summarise_group_sequential()
)

times <- c(times, "test_1_setup"=Sys.time())

result_1 <- runSimulation(
  Design,
  replications = N_sim,
  generate = my_generator,
  analyse = my_analyse,
  summarise = my_summarise,
  cl = cl,
  save=FALSE
)

times <- c(times, "test_1_sim"=Sys.time())

# output performance measures ---------------------------------------------

performance <- list(
  n_cores = n_cores,
  sessioninfo = sessionInfo(),
  duration = diff(times),
  N_sim = N_sim,
  N_methods = length(my_analyse),
  N_scenarios = nrow(Design)
)

filename <- paste0("script_2_", Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S"), ".Rds")
saveRDS(performance, filename)
