
library(SimNPH)
library(SimDesign)
library(parallel)

if(packageVersion("SimNPH") != "0.3.2"){
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
  library(SimNPH)
  library(parallel)
})


# setup data generation ---------------------------------------------------

# load parameters
design <- read.table("data/parameters/additional_progression_os_2023-04-19.csv", sep=",", dec=".", header=TRUE)

# define generator
my_generator <- function(condition, fixed_objects=NULL){
  generate_progression(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_events(condition$final_events)
}

alpha <- 0.025
nominal_alpha <- ldbounds::ldBounds(c(0.5,1), sides=1, alpha = 0.025)$nom.alpha

clusterExport(cl, "nominal_alpha")


# define analysis and summarise functions ---------------------------------
source("scripts/additional_sims_common.R")

# run ---------------------------------------------------------------------

save_folder <- "data/results/progression"
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


saveRDS(results, paste0(save_folder, "/additional_results_progression_", Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S"), ".Rds"))


