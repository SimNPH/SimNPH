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

install.packages(c("SimDesign", "remotes"))
remotes::install_git("https://github.com/SimNPH/SimNPH.git")

# preparation -------------------------------------------------------------

library(SimDesign)
library(SimNPH)
library(parallel)

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

N_sim <- 250
Assumptions <- assumptions_delayed_effect()
Options <- design_fixed_followup()
Design <- merge(Assumptions, Options, by=NULL)

my_generator <- function(condition, fixed_objects=NULL){
  generate_delayed_effect(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_time(condition$followup)
}

alpha <- 0.05

Summarise <- create_summarise_function(
  maxcombo = function(condition, results, fixed_objects=NULL){
    data.frame("rejection"=mean(results$p < alpha))
  },
  logrank = function(condition, results, fixed_objects=NULL){
    data.frame("rejection"=mean(results$p < alpha))
  }
)

times <- c(times, "test_1_setup"=Sys.time())

result_1 <- runSimulation(
  Design,
  replications = N_sim,
  generate = my_generator,
  analyse = list(
    logrank  = analyse_logrank,
    maxcombo = analyse_maxcombo
  ),
  summarise = Summarise,
  cl = cl,
  save=FALSE
)

times <- c(times, "test_1_sim"=Sys.time())

# Test 2 group sequential -------------------------------------------------

Options <- design_group_sequential()
Design <- merge(Assumptions, Options, by=NULL)


nominal_alpha <- ldbounds::ldBounds(c(0.5,1))$nom.alpha

clusterExport(cl, "nominal_alpha")

Analyse <-  list(
  logrank_seq  = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank
  ),
  maxcombo_seq = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_maxcombo
  )
)

Summarise <- create_summarise_function(
  maxcombo_seq = summarise_group_sequential(),
  logrank_seq = summarise_group_sequential()
)


times <- c(times, "test_2_setup"=Sys.time())

result_2 <- runSimulation(
  Design,
  replications = N_sim,
  generate = my_generator,
  analyse = Analyse,
  summarise = Summarise,
  cl = cl,
  save=FALSE
)

times <- c(times, "test_2_sim"=Sys.time())

# Test 3 estimator --------------------------------------------------------

Options <- design_fixed_followup()
Design <- merge(Assumptions, Options, by=NULL)
Design <- Design |>
  true_summary_statistics_delayed_effect()

Summarise <- create_summarise_function(
    coxph=summarise_estimator(hr, gAHR, hr_lower, hr_upper, name="gAHR"),
    coxph=summarise_estimator(hr, hazard_trt/hazard_ctrl, hr_lower, hr_upper, name="HR")
  )

Analyse <- list(
  coxph = analyse_coxph
)

times <- c(times, "test_3_setup"=Sys.time())

result_3 <- runSimulation(
  Design,
  replications = N_sim,
  generate = my_generator,
  analyse = Analyse,
  summarise = Summarise,
  cl = cl,
  save=FALSE
)

times <- c(times, "test_3_sim"=Sys.time())


# output performance measures ---------------------------------------------

performance <- list(
  n_cores = n_cores,
  sessioninfo = sessionInfo(),
  duration = diff(times),
  N_sim = N_sim
)

filename <- paste0(Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S"), ".Rds")
saveRDS(performance, filename)


# output of the simulations -----------------------------------------------
# this script is mainly to check if the simulations run on each computer
# but some simulations are done (yet with unrealistic parameters)
# if you want to inspect the results, use the commented code below.

# result_1 |>
#   subset(select=c("delay", "maxcombo.rejection", "logrank.rejection"))
#
# result_2 |>
#   subset(select=c(
#     "delay",
#     "maxcombo_seq.rejection", "logrank_seq.rejection",
#     "maxcombo_seq.n_pat", "logrank_seq.n_pat",
#     "maxcombo_seq.n_evt", "logrank_seq.n_evt",
#     "maxcombo_seq.followup", "logrank_seq.followup"
#   ))
#
# result_3 |>
#   subset(select=c(
#     "delay", "coxph.HR.bias", "coxph.gAHR.bias", "coxph.HR.mse",
#     "coxph.gAHR.mse", "coxph.HR.coverage", "coxph.gAHR.coverage"
#   ))
