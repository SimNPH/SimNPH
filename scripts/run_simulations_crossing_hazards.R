package_url <- "https://github.com/SimNPH/SimNPH/archive/refs/tags/sims_crossing.zip"
download.file(package_url, destfile = "sims_crossing.zip")
dir.create("sims_crossing")
unzip("sims_crossing.zip", exdir = "sims_crossing")
setwd("sims_crossing/SimNPH-sims_crossing/")

# installing package to local library to not interfere with other versions
# running in the same filesystem
dir.create("local_lib")
withr::with_libpaths(
  "local_lib",
  devtools::install(".", build=TRUE, quick=TRUE, dependencies=TRUE),
  action="prefix"
)

withr::with_libpaths("local_lib", library("SimNPH"))
library(SimDesign)
library(parallel)

if(packageVersion("SimNPH") != "0.1.2"){
  stop("Please run the simulations with the correct vesion of the SimNPH package for reproducability.")
}

N_sim <- 2500

# Helper functions --------------------------------------------------------

# months to days
m2d <- \(t) 365.25*t/12

# setup cluster -----------------------------------------------------------

n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores)

clusterEvalQ(cl, {
  library(SimDesign)
  withr::with_libpaths("local_lib", library("SimNPH"))
  library(parallel)
})


# setup data generation ---------------------------------------------------

# load parameters
design <- read.table("data/parameters/crossing_hazards_2022-12-21.csv", sep=",", dec=".", header=TRUE)

# define generator
my_generator <- function(condition, fixed_objects=NULL){
  generate_crossing_hazards(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_events(condition$final_events)
}

alpha <- 0.05
nominal_alpha <- ldbounds::ldBounds(c(0.5,1))$nom.alpha

clusterExport(cl, "nominal_alpha")


# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  ahr  = analyse_ahr(type="AHR"),
  gahr = analyse_ahr(type="gAHR"),
  mean_surv = analyse_diff_mean_survival(),
  milestone = analyse_milestone_survival(times=m2d(c(8, 12))),
  rmst = analyse_rmst_diff(),
  cox = analyse_coxph(),
  weighted_cox = analyse_coxph_weighted(type="G"),
  # tests
  pw_exp_3  = analyse_piecewise_exponential(cuts=m2d(seq(0, 240, by= 3)), testing_only=TRUE),
  pw_exp_12 = analyse_piecewise_exponential(cuts=m2d(seq(0, 240, by=12)), testing_only=TRUE),
  peto_peto = analyse_gehan_wilcoxon(),
  fh = analyse_logrank_fh_weights(rho=0, gamma=1),
  logrank = analyse_logrank(),
  max_combo = analyse_maxcombo(),
  modest_6 = analyse_modelstly_weighted(t_star=m2d(6)),
  modest_8 = analyse_modelstly_weighted(t_star=m2d(8)),
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
    analyse_functions = analyse_modelstly_weighted(t_star=m2d(6))
  ),
  modest_gs_8 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_modelstly_weighted(t_star=m2d(8))
  )
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  ahr = summarise_estimator(est=AHR, real=AHR),
  gahr = summarise_estimator(est=gAHR, real=gAHR),
  mean_surv = summarise_estimator(est=diff_Q, real=median_survival_trt/median_survival_ctrl),
  milestone = summarise_estimator(est=milestone_surv_ratio[1], real= ms_surv_8_ctrl/ ms_surv_8_trt, name="milestone_8" ),
  milestone = summarise_estimator(est=milestone_surv_ratio[2], real=ms_surv_12_ctrl/ms_surv_12_trt, name="milestone_12"),
  rmst = summarise_estimator(est=rmst_diff, real=rmst_trt-rmst_ctrl),
  cox = summarise_estimator(est=hr, real=gAHR, lower=hr_lower, upper=hr_upper),
  weighted_cox = summarise_estimator(est=hr, real=gAHR, lower=hr_lower, upper=hr_upper),
  pw_exp_3 = summarise_test(alpha),
  pw_exp_12 = summarise_test(alpha),
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

save_folder <- paste0(paste0("data/simulation_crossing_hazards_", Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S")))
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

