# Load the library
# install first by running:
# remotes::install_git("https://github.com/SimNPH/SimNPH.git")
library(SimNPH)
library(parallel)

if(packageVersion("SimNPH") != "0.2.0"){
  stop("Please run the simulations with the correct vesion of the SimNPH package for reproducability.")
}

# Script options/settings --------------------------------------------------------

N_sim <- 2500 # number of simulations per scenario

run_parallel <- TRUE # should we parallelize?
n_cores <- parallel::detectCores() - 4

save_folder <- here::here(
  "data",
  paste0("simulation_common_ext_data_",
         Sys.info()["nodename"], "_",
         strftime(Sys.time(), "%Y-%m-%d_%H%M%S")))
sim_data_dir <- here::here("..","..","parametric_simulations","sim_data")

alpha <- 0.025
nominal_alpha <- ldbounds::ldBounds(c(0.5,1), sides=1, alpha=0.025)$nom.alpha

# Helper functions --------------------------------------------------------

# to convert from months to days
m2d <- \(t) 365.25*t/12

# setup cluster -----------------------------------------------------------

cl <- NULL
if(run_parallel){
  cl <- makeCluster(n_cores)

  clusterEvalQ(cl, {
    library(SimDesign)
    library(SimNPH)
    library(parallel)
  })

  clusterExport(cl, c("sim_data_dir","save_folder","alpha","nominal_alpha"))

}

# Design data.frame -------------------------------------------------------

# In this data.frame describe the design, either read in my dataset and merge
# or edit the code here.
# The columns final_events and recruitment are needed to add administrative
# censoring.
# design_1 <- expand.grid(
#   scenario     = 1,
#   final_events = 300 * 0.75,
#   recruitment  = m2d(c(18, 30))
# )
#
# Design <- design_1
# # |> dplyr::arrange(scenario)
Design <- read.table("data/parameters/ext_data_2023-03-14.csv", sep=",", dec=".", header=TRUE)
# Design <- Design[1:2,]

# Function to read in dataset ---------------------------------------------

# This function reads the dataset and converts it to the format used internally.
# The path of the file is built from the Design data.frame.

# condition is a one line data.frame containing the line of Design that is
# currently simulated. Additionally condition has the column REPLICATIONS that
# contains the number of the replication.
generate_read <- function(condition, fixed_objects=NULL){
  # TODO: edit to fit your folder structure and filename pattern
  filename <- paste0(
    # folder name relative to the folder the simulation is started
    sim_data_dir,
    "/scenario_",
    # format the scenario number, add leading 0 to a width of 4 digits
    formatC(condition$scenario, width=NULL, flag="0"),
    # start of filename
    "/simtab",
    # format number of replication, add leading 0 to a width of 4 digits
    formatC(condition$REPLICATION, width=NULL, flag="0"),
    # end of filename / extension
    ".dat.zip"
  )

  # dat <- read.csv(filename)
  dat <- readr::read_table(filename,show_col_types=FALSE)

  dat <- dat |> dplyr::filter(FLAG!=1) |>
    dplyr::select("id"=ID,"t"=TIME,"evt"=DV,"trt"=TRT)

  # rename the columns
  # names(dat) <- c("ID"="id", "TIME"="t", "DV"="evt", "TRT"="trt")[names(dat)]

  # convert the columns to the format used by SimNPH
  dat <- within(dat, {
    trt = 1-trt
    evt = as.logical(evt)
  })

  # return the data
  dat
}


# generate, analyse, summarise functions ----------------------------------

# generate
# 1. read data
# 2. add time of recruitment
# 3. add administrative censoring
my_generator <- function(condition, fixed_objects=NULL){
  generate_read(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    admin_censoring_events(condition$final_events)
}

# analyse
# example with just logrank test and cox regression
# will later be replaced by sourcing scripts/run_simulations_common.R
# my_analyse <- list(
#   logrank     = analyse_logrank(),
#   coxph       = analyse_coxph(),
#   descriptive = analyse_describe()
# )

# summarise
# example with just logrank test, cox regression and descriptive statistics
# will later be replaced by sourcing scripts/run_simulations_common.R
# my_summarise <- create_summarise_function(
#   logrank     = summarise_test(alpha),
#   coxph       = summarise_estimator(coef, 0),
#   descriptive = summarise_describe()
# )

source("scripts/run_simulations_common.R")


# run simulations ---------------------------------------------------------

dir.create(save_folder)

results <- runSimulation(
  Design,
  replications = N_sim,
  generate = my_generator,
  analyse = my_analyse,
  summarise = my_summarise,
  cl=cl,
  save_details = list(
    out_rootdir = save_folder
  ),
  #save_seeds = T, # if MPI
  # include_replication_index is crucial
  # if this argument is missing condition will not contain REPLICATIONS column
  extra_options = list(include_replication_index=TRUE,
                       store_warning_seeds = TRUE,
                       allow_na = TRUE)
)

if(run_parallel) stopCluster(cl)


saveRDS(results, paste0(save_folder, "/results.Rds"))

results
names(results)
results$SEED
results$logrank.rejection_0.025

