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

if(packageVersion("SimNPH") != "0.1.4"){
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
design <- read.table("data/parameters/crossing_hazards_2023-01-17.csv", sep=",", dec=".", header=TRUE)

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



# define analysis and summarise functions ---------------------------------
# this is the same for all scenarios

source("scripts/run_simulations_common.R")


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
  ),
  extra_options = list(
    store_warning_seeds = TRUE,
    allow_na = TRUE
  )
)

saveRDS(results, paste0(save_folder, "/results.Rds"))

