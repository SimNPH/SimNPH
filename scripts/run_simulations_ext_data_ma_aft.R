# Load the library
# install first by running:
# remotes::install_git("https://github.com/SimNPH/SimNPH.git")
library(SimNPH)
library(SimDesign)
library(parallel)

if(packageVersion("SimNPH") != "0.3.2"){
  stop("Please run the simulations with the correct vesion of the SimNPH package for reproducability.")
}

source("./R/analyse_ma_aft.R")

# Script options/settings --------------------------------------------------------
N_sim <- 5 # number of simulations per scenario

run_parallel <- TRUE # should we parallelize?
n_cores <- parallel::detectCores() - 4

# assume we are in base directory of the package
save_folder <- file.path(
  #"..",
  "data",
  paste0("simulation_ext_data_ma_aft",
         Sys.info()["nodename"], "_",
         strftime(Sys.time(), "%Y-%m-%d_%H%M%S")))

sim_data_dir <- file.path("..","..","parametric_simulations","sim_data_3")

design_table <- file.path("data","parameters","ext_data_2023-04-22.csv")

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
Design <- read.table(
  design_table,
  sep=",", dec=".", header=TRUE)
Design <- Design[c(
  1
  #,5,9,13,17,21,25,29,33,37,41,45
  ),]

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
    dplyr::select("id"=ID,"t"=TIME,"evt"=DV,"trt"=TRT) |>
    dplyr::mutate(t=t/24,trt=1-trt,evt=as.logical(evt))

  if(length(unique(dat$trt))==1) {
    warning(filename," has only one group\n")
    return(NULL)
  }

  # rename the columns
  # names(dat) <- c("ID"="id", "TIME"="t", "DV"="evt", "TRT"="trt")[names(dat)]

  # # convert the columns to the format used by SimNPH
  # dat <- within(dat, {
  #   trt = 1-trt
  #   evt = as.logical(evt)
  # })

  # return the data
  dat
}


# generate, analyse, summarise functions ----------------------------------

# generate
# 1. read data
# 2. add time of recruitment
# 3. add administrative censoring
my_generator <- function(condition, fixed_objects=NULL){
  dat <- generate_read(condition, fixed_objects) |>
    recruitment_uniform(condition$recruitment) |>
    admin_censoring_events(condition$final_events)

  dat$t <- round(dat$t, 5)
  dat
}


# analyse
# example with just logrank test and cox regression
# will later be replaced by sourcing scripts/run_simulations_common.R
my_analyse <- list(
  ma_aft = analyse_ma_aft(alternative = "one.sided",n_boot=10)
  # logrank     = analyse_logrank(alternative = "one.sided"),
  # cox = analyse_coxph(alternative = "one.sided"),
  # aft_weibull = analyse_aft(dist="weibull", alternative = "one.sided"),
  # ahr_6m   = analyse_ahr(type="AHR", max_time=m2d(6), alternative = "one.sided")
  )
my_analyse <- wrap_all_in_trycatch(my_analyse)

# summarise
# example with just logrank test, cox regression and descriptive statistics
# will later be replaced by sourcing scripts/run_simulations_common.R
my_summarise <- create_summarise_function(
  ma_aft = summarise_estimator(est=taf_est, real=NA_real_,lower=taf_lower,upper=taf_upper,null=1)#,
  # ma = summarise_estimator(est=surv_6m_ratio_est, real=milestone_surv_ratio_6m,lower=surv_6m_ratio_lower,upper=surv_6m_ratio_upper,null=1,name="surv_ratio_6m"),
  # ma = summarise_estimator(est=surv_12m_ratio_est, real=milestone_surv_ratio_12m,lower=surv_12m_ratio_lower,upper=surv_12m_ratio_upper,null=1,name="surv_ratio_12m"),
  # ma = summarise_estimator(est=surv_24m_ratio_est, real=milestone_surv_ratio_24m,lower=surv_24m_ratio_lower,upper=surv_24m_ratio_upper,null=1,name="surv_ratio_24m"),
  # ma = summarise_estimator(est=AIC_diff_est, real=NA_real_,lower=AIC_diff_lower,upper=AIC_diff_upper,null=0,name="AIC_diff")
  # logrank     = summarise_test(alpha),
  # cox = summarise_estimator(est=hr, real=gAHR_12m, lower=hr_lower, upper=hr_upper, null=1),
  # # cox = summarise_estimator(est=hr, real=gAHR_12m, lower=hr_lower, upper=hr_upper, null=1),
  # aft_weibull = summarise_estimator(est=coef, real=NA_real_, lower=lower, upper=upper, null=0),
  # ahr_6m  = summarise_estimator(est=AHR, real=AHR_6m, lower=AHR_lower, upper=AHR_upper, null=1)
)

# source("scripts/run_simulations_common.R")
#
# my_summarise <- create_summarise_function(
#   # estimation
#   ahr_6m  = summarise_estimator(est=AHR, real=AHR_6m, lower=AHR_lower, upper=AHR_upper, null=1),
#   ahr_12m = summarise_estimator(est=AHR, real=AHR_12m, lower=AHR_lower, upper=AHR_upper, null=1),
#   gahr_6m  = summarise_estimator(est=gAHR, real=gAHR_6m, lower=gAHR_lower, upper=gAHR_upper, null=1),
#   gahr_12m = summarise_estimator(est=gAHR, real=gAHR_12m, lower=gAHR_lower, upper=gAHR_upper, null=1),
#   median_surv = summarise_estimator(est=diff_Q, real=median_survival_trt-median_survival_ctrl, lower=diff_Q_lower, upper=diff_Q_upper, null=0),
#   milestone_6m = summarise_estimator(est=milestone_surv_ratio, real= milestone_surv_ratio_6m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
#   milestone_12m = summarise_estimator(est=milestone_surv_ratio, real=milestone_surv_ratio_12m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
#   rmst_diff_6m  = summarise_estimator(est=rmst_diff, real=RMST_diff_6m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
#   rmst_diff_12m = summarise_estimator(est=rmst_diff, real=RMST_diff_12m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
#   cox = summarise_estimator(est=hr, real=gAHR_12m, lower=hr_lower, upper=hr_upper, null=1),
#   weighted_cox_6m  = summarise_estimator(est=hr, real=gAHR_6m, lower=hr_lower, upper=hr_upper, null=1),
#   weighted_cox_12m = summarise_estimator(est=hr, real=gAHR_12m, lower=hr_lower, upper=hr_upper, null=1),
#   aft_weibull = summarise_estimator(est=coef, real=NA_real_, lower=lower, upper=upper, null=0),
#   aft_lognormal = summarise_estimator(est=coef, real=NA_real_, lower=lower, upper=upper, null=0),
#   diff_med_weibull = summarise_estimator(est=diff_med_est, real=med_diff,lower=diff_med_lower, upper=diff_med_upper, null=0),
#   # tests
#   pw_exp_3 = summarise_test(alpha),
#   pw_exp_12 = summarise_test(alpha),
#   peto_peto = summarise_test(alpha),
#   fh_0_0 = summarise_test(alpha),
#   fh_0_1 = summarise_test(alpha),
#   fh_1_0 = summarise_test(alpha),
#   fh_1_1 = summarise_test(alpha),
#   logrank = summarise_test(alpha),
#   max_combo = summarise_test(alpha),
#   modest_6 = summarise_test(alpha),
#   modest_8 = summarise_test(alpha),
#   # tests group sequential
#   peto_peto_gs = summarise_group_sequential(),
#   fh_gs_0_0 = summarise_group_sequential(),
#   fh_gs_0_1 = summarise_group_sequential(),
#   fh_gs_1_0 = summarise_group_sequential(),
#   fh_gs_1_1 = summarise_group_sequential(),
#   logrank_gs = summarise_group_sequential(),
#   max_combo_gs = summarise_group_sequential(),
#   modest_gs_6 = summarise_group_sequential(),
#   modest_gs_8 = summarise_group_sequential(),
#   descriptive = summarise_describe()
# )

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
  extra_options = list(
    include_replication_index=TRUE,
    store_warning_seeds = TRUE,
    allow_na = TRUE
  )
)

if(run_parallel) stopCluster(cl)


saveRDS(results, paste0(save_folder, "/results.Rds"))

