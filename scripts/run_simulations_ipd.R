
library(SimNPH)
library(SimDesign)
library(parallel)

if(packageVersion("SimNPH") != "0.3.0"){
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
  )
)

assumptions <- data.frame(
  null = c(T,F)
)

design <- merge(
  ipd_files,
  assumptions,
  by=NULL
)

ipd_data <- lapply(split(ipd_files, 1:nrow(ipd_files)), \(i){
  load(i$file)
  names(all_patients) <- c("t", "evt", "trt")
  all_patients$evt <- as.logical(all_patients$evt)
  all_patients
})
names(ipd_data) <- ipd_files$name
clusterExport(cl, "ipd_data")

# define generator
my_generator <- function(condition, fixed_objects=NULL){
  tmp_data <- ipd_data[[condition$name]]
  N <- nrow(tmp_data)
  if(condition$null){
    tmp_trt <- tmp_data$trt
    tmp_data <- tmp_data[tmp_data$trt == 0, ]
    tmp_data$trt <- sample(tmp_trt, nrow(tmp_data), replace=TRUE)
  }

  tmp_data[sample(1:nrow(tmp_data), N, replace = TRUE), ]
}

alpha <- 0.025
nominal_alpha <- ldbounds::ldBounds(c(0.5,1), sides=1, alpha = 0.025)$nom.alpha

clusterExport(cl, "nominal_alpha")


# define analysis and summarise functions ---------------------------------
# this is the same for all scenarios

source("scripts/run_simulations_common.R")

# run ---------------------------------------------------------------------

save_folder <- paste0(paste0("data/simulation_ipd_", Sys.info()["nodename"], "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S")))
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

