
library(SimNPH)
library(SimDesign)
library(parallel)

if(packageVersion("SimNPH") != "0.5.4"){
  stop("Please run the simulations with the correct vesion of the SimNPH package for reproducability.")
}

N_sim <- 2500

# recalculate summary statistics ------------------------------------------

# load parameters
design <- read.table("data_sim_study/parameters/revision_crossing_hazards_2023-05-22.csv", sep=",", dec=".", header=TRUE)

design <- design |>
  dplyr::select(-c(median_survival_trt:milestone_survival_ctrl_12m))

design <- purrr::map_dfr(1:nrow(design), function(n){
  data <- design[n,]
  true_summary_statistics_crossing_hazards(
    data,
    cutoff_stats = m2d(c("24m"=24, "36m"=36))
  )})

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
  generate_crossing_hazards(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    random_censoring_exp(condition$random_withdrawal) |>
    admin_censoring_events(condition$final_events)
}

alpha <- 0.025
nominal_alpha <- ldbounds::ldBounds(c(0.5,1), sides=1, alpha = 0.025)$nom.alpha

clusterExport(cl, "nominal_alpha")


# define analysis and summarise functions ---------------------------------
# this is the same for all scenarios

source("scripts/revision_2_sims_common.R")

# run ---------------------------------------------------------------------

save_folder <- paste0(paste0("data/revision_2_crossing_hazards_", Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S")))
dir.create(save_folder)

results <- runSimulation(
  design,
  replications = N_sim,
  generate = my_generator,
  analyse = my_analyse,
  summarise = my_summarise,
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

saveRDS(results, paste0(save_folder, "/results.Rds"))

