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


# ipd, assumptions, etc. --------------------------------------------------

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


# scale ipd data to days
# weeks
ipd_data[[1]]$t <- ipd_data[[1]]$t * 7
# months
ipd_data[[2]]$t <- m2d(ipd_data[[2]]$t)
# months
ipd_data[[3]]$t <- m2d(ipd_data[[3]]$t)

# distribute data to all cluster nodes
clusterExport(cl, "ipd_data")


# data generation for "true" summary statistics ---------------------------
# this data generator only takes one sample from the ipd datasets
# to calculate the values of the summary statistics in the original data

# define generator
my_generator <- function(condition, fixed_objects=NULL){
  # select dataset according to column in condition
  tmp_data <- ipd_data[[condition$name]]
  N <- nrow(tmp_data)
  ntrt=sum(tmp_data$trt)
  nc=N-ntrt
  # under the null:  sample both groups from control
  # under the alternative: sample trt from trt and control from control
  if(condition$null){
    c_data <- tmp_data[tmp_data$trt == 0, ]
    boot_data=rbind(c_data,c_data)
    boot_data$trt=c(rep(0,nc),rep(1,nc))
  }
  else
  {
    c_data <- tmp_data[tmp_data$trt == 0, ]
    trt_data <- tmp_data[tmp_data$trt == 1, ]
    boot_c_data=c_data
    boot_trt_data=trt_data
    boot_data=rbind(boot_c_data,boot_trt_data)
  }
  # return a sample of the same size as the original dataset
  return(boot_data)
}


alpha <- 0.025
nominal_alpha <- ldbounds::ldBounds(c(0.5,1), sides=1, alpha = 0.025)$nom.alpha

clusterExport(cl, "nominal_alpha")


# define analysis and summarise functions ---------------------------------
# this is the same for all scenarios

source("scripts/run_simulations_common.R")

# run with only one replication -------------------------------------------

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


# update the design with "true" summary statistics ------------------------
design_true <- results |>
  dplyr::mutate(AHR_6m = ahr_6m.mean_est,
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


# data generation for actual bootstrap ------------------------------------

# define generator
my_generator_sim <- function(condition, fixed_objects=NULL){
  # select dataset according to column in condition
  tmp_data <- ipd_data[[condition$name]]
  N <- nrow(tmp_data)
  ntrt=sum(tmp_data$trt)
  nc=N-ntrt
  # under the null:  sample both groups from control
  # under the alternative: sample trt from trt and control from control
  if(condition$null){
    c_data <- tmp_data[tmp_data$trt == 0, ]
    boot_data=c_data[sample(nc, N, replace=TRUE),]
    boot_data$trt=c(rep(0,nc),rep(1,ntrt))
  }
  else
  {
    c_data <- tmp_data[tmp_data$trt == 0, ]
    trt_data <- tmp_data[tmp_data$trt == 1, ]
    boot_c_data=c_data[sample(nc, nc, replace=TRUE),]
    boot_trt_data=trt_data[sample(ntrt, ntrt, replace=TRUE),]
    boot_data=rbind(boot_c_data,boot_trt_data)
  }
  # return a sample of the same size as the original dataset
  return(boot_data)
}



# run bootstrap -----------------------------------------------------------
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
