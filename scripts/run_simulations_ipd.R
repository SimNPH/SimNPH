library(SimNPH)
library(SimDesign)
library(parallel)

N_sim <- 500

# Helper functions --------------------------------------------------------

# months to days
m2d <- \(t) 365.25*t/12

# setup cluster -----------------------------------------------------------

## Be careful If you do anything else on your computer this may slow things down considerably, maybe leave some cores for other work; there is also a workaround that looks at free cores;
n_cores <- parallel::detectCores()
cl <- makeCluster(n_cores-2)

clusterEvalQ(cl, {
  library(SimDesign)
  library(SimNPH)
  library(parallel)
})


# setup data generation ---------------------------------------------------

# files for IPD
ipd_files <- data.frame(
  file = c(
    "data/ipd/506_time_to_relapse.Rdata",
    "data/ipd/6481_os.Rdata",
    "data/ipd/9883_os.Rdata"
  ),
  name = c(
    "Terifluomide",
    "Nivolumab",
    "Prembrolizumap"
  ), ## We need to add "true" values
  n_pat = c(166,419,559),
  n_trt = c(109,210,281),
  final_events = c(65,328,447),
  interim_events = c(65,328,447),
  AHR_6m = 0,
  AHR_12m = 0,
  gAHR_6m = 0,
  gAHR_12m = 0,
  median_survival_ctrl = 0,
  median_survival_trt = 0,
  median_survival_diff = 0,
  milestone_survival_ctrl_6m = 0,
  milestone_survival_trt_6m = 0,
  milestone_survival_ratio_6m = 1,
  milestone_survival_ctrl_12m = 0,
  milestone_survival_trt_12m =0,
  milestone_survival_ratio_12m = 1,
  rmst_trt_6m = 0,
  rmst_ctrl_6m = 0,
  rmst_diff_6m = 0,
  rmst_trt_12m = 0,
  rmst_ctrl_12m = 0,
  rmst_diff_12m = 0
)

# assumptions (null or real data)
assumptions <- data.frame(
  null = c(TRUE, FALSE)
)

# TODO: add options (sample size, recruitment, admininstrative censoring)?

design <- merge(
  ipd_files,
  assumptions,
  by=NULL
)

# read in the actual data
ipd_data <- lapply(split(ipd_files, 1:nrow(ipd_files)), \(i){
  load(i$file)
  names(all_patients) <- c("t", "evt", "trt")
  all_patients$evt <- as.logical(all_patients$evt)
  all_patients
})
names(ipd_data) <- ipd_files$name
# distribute data to all cluster nodes
clusterExport(cl, "ipd_data")

# define generator
my_generator <- function(condition, fixed_objects=NULL){
  # select dataset according to column in condition
  tmp_data <- ipd_data[[condition$name]]
  N <- nrow(tmp_data)

  # under the null modify the dataset:
  #   sample only from control
  #   sample treatment from full dataset
  if(condition$null){
    tmp_trt <- tmp_data$trt
    tmp_data <- tmp_data[tmp_data$trt == 0, ]
    tmp_data$trt <- sample(tmp_trt, nrow(tmp_data), replace=TRUE)
  }

  # return a sample of the same size as the original dataset
  tmp_data#[sample(1:nrow(tmp_data), N, replace = TRUE), ]
}

alpha <- 0.025
nominal_alpha <- ldbounds::ldBounds(c(0.5,1), sides=1, alpha = 0.025)$nom.alpha

clusterExport(cl, "nominal_alpha")


# define analysis and summarise functions ---------------------------------
# this is the same for all scenarios

source("scripts/run_simulations_common_ipd.R")

# run ---------------------------------------------------------------------

## we first run this with one replacition and no resampling to obtain "true" statistics which correspond to the estimates of the original data
save_folder <- paste0(paste0("data/simulation_ipd_", Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S")))
dir.create(save_folder)

results <- runSimulation(
  design,
  replications = 1,
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

## we update the design (for now)

design_true <- results |>
  mutate(AHR_6m = ahr_6m.mean_est,
         AHR_12m = ahr_12m.mean_est,
         gAHR_6m = gahr_6m.mean_est,
         gAHR_12m = gahr_12m.mean_est,
         median_survival_diff = median_surv.mean_est,
         milestone_survival_ratio_6m = milestone_6m.mean_est,
         milestone_survival_ratio_12m = milestone_12m.mean_est,
         rmst_diff_6m = rmst_diff_6m.mean_est,
         rmst_diff_12m = rmst_diff_12m.mean_est
         )


design_sim <-  design_true[names(design)]
## To assess the bias we could also just compute the differences (and change analyse/summarise functions accordingly)

# define generator
my_generator_sim <- function(condition, fixed_objects=NULL){
  # select dataset according to column in condition
  tmp_data <- ipd_data[[condition$name]]
  N <- nrow(tmp_data)

  # under the null modify the dataset:
  #   sample only from control
  #   sample treatment from full dataset
  if(condition$null){
    tmp_trt <- tmp_data$trt
    tmp_data <- tmp_data[tmp_data$trt == 0, ]
    tmp_data$trt <- sample(tmp_trt, nrow(tmp_data), replace=TRUE)
  }

  # return a sample of the same size as the original dataset
  tmp_data[sample(1:nrow(tmp_data), N, replace = TRUE), ]
}

# run ---------------------------------------------------------------------


## now we run it with N_sim replications and actual resampling

results <- runSimulation(
  design_sim,
  replications = N_sim,
  generate = my_generator_sim,
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

## For testing purposes

# analyse_ahr(type="AHR", max_time=m2d(6), alternative = "one.sided")(design,ipd_data[[1]])
# dat <- ipd_data[[1]]
# nph::nphparams(dat$t,
#                dat$evt,
#                dat$trt,
#                param_type=c("Q"),
#                param_par=c(0.4))
# my_analyse[[1]](design[1,],ipd_data[[1]])
#
#
