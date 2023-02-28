# Load the library
# install first by running:
# remotes::install_git("https://github.com/SimNPH/SimNPH.git")
library(SimNPH)

# number of simulations
# TODO: set this to the number of datasets of replications
N_sim <- 2

# to convert from months to days
m2d <- \(t) 365.25*t/12


# Design data.frame -------------------------------------------------------

# In this data.frame describe the design, either read in my dataset and merge
# or edit the code here.
# The columns final_events and recruitment are needed to add administrative
# censoring.
Design <- expand.grid(
  scenario     = 1,
  final_events = 1500 * 0.75,
  recruitment  = m2d(c(18, 30))
)

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
    "testdata/scenario_",
    # format the scenario number, add leading 0 to a width of 4 digits
    formatC(condition$scenario, width=1, flag="0"),
    # start of filename
    "/Dataset",
    # format number of replication, add leading 0 to a width of 4 digits
    formatC(condition$REPLICATION, width=4, flag="0"),
    # end of filename / extension
    ".csv"
  )

  dat <- read.csv(filename)

  # rename the columns
  names(dat) <- c("ID"="id", "TIME"="t", "DV"="evt", "TRT"="trt")[names(dat)]

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
my_analyse <- list(
  logrank     = analyse_logrank(),
  coxph       = analyse_coxph(),
  descriptive = analyse_describe()
)

# summarise
# example with just logrank test, cox regression and descriptive statistics
# will later be replaced by sourcing scripts/run_simulations_common.R
my_summarise <- create_summarise_function(
  logrank     = summarise_test(0.05),
  coxph       = summarise_estimator(coef, 0),
  descriptive = summarise_describe()
)


# run simulations ---------------------------------------------------------

res <- runSimulation(
  Design,
  replications = N_sim,
  generate = my_generator,
  analyse = my_analyse,
  summarise = my_summarise,
  # include_replication_index is crucial
  # if this argument is missing condition will not contain REPLICATIONS column
  extra_options = list(include_replication_index=TRUE)
)

res
