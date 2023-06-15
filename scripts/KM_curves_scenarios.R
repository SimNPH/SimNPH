
library(SimNPH)
library(SimDesign)
library(survival)
library(tidyverse)

N_sim <- 1
N_sample <- 100000

# delayed -----------------------------------------------------------------
design <- read.table("data/parameters/delayed_effect_2023-03-13.csv", sep=",", dec=".", header=TRUE)

design <- design |>
  mutate(
    n_ctrl = N_sample/2,
    n_trt  = N_sample/2
  ) |>
  select(delay, hazard_ctrl, hazard_trt,  n_trt, n_ctrl, effect_size_ph, recruitment) |>
  unique()

results <- runSimulation(
  design,
  replications = N_sim,
  generate = generate_delayed_effect,
  analyse = function(condition, dat, fixed_objects=NULL){
    survival::survfit(survival::Surv(t, evt)~trt, dat)
  },
  extra_options = list(
    store_warning_seeds = TRUE,
    allow_na = TRUE
  )
)

results_delayed <- bind_cols(
  design,
  tibble(
    survfit = map(results, \(x) x[[1]])
  )
)

# subgroup ----------------------------------------------------------------

design <- read.table("data/parameters/subgroup_2023-03-13.csv", sep=",", dec=".", header=TRUE)

design <- design |>
  mutate(
    n_ctrl = N_sample/2,
    n_trt  = N_sample/2
  ) |>
  select(hazard_ctrl, hazard_trt, hazard_subgroup, prevalence, n_trt, n_ctrl, effect_size_ph, recruitment, hr_subgroup_relative) |>
  unique()

results <- runSimulation(
  design,
  replications = N_sim,
  generate = generate_subgroup,
  analyse = function(condition, dat, fixed_objects=NULL){
    survival::survfit(survival::Surv(t, evt)~trt, dat)
  },
  extra_options = list(
    store_warning_seeds = TRUE,
    allow_na = TRUE
  )
)

results_subgroup <- bind_cols(
  design,
  tibble(
    survfit = map(results, \(x) x[[1]])
  )
)

# crossing ----------------------------------------------------------------


design <- read.table("data/parameters/crossing_hazards_2023-03-13.csv", sep=",", dec=".", header=TRUE)

design <- design |>
  mutate(
    n_ctrl = N_sample/2,
    n_trt  = N_sample/2
  ) |>
  select(crossing, hazard_ctrl, hazard_trt_before, hazard_trt_after, n_trt, n_ctrl, effect_size_ph, recruitment, hr_before) |>
  unique()

results <- runSimulation(
  design,
  replications = N_sim,
  generate = generate_crossing_hazards,
  analyse = function(condition, dat, fixed_objects=NULL){
    survival::survfit(survival::Surv(t, evt)~trt, dat)
  },
  extra_options = list(
    store_warning_seeds = TRUE,
    allow_na = TRUE
  )
)

results_crossing <- bind_cols(
  design,
  tibble(
    survfit = map(results, \(x) x[[1]])
  )
)

# progression -------------------------------------------------------------

design <- read.table("data/parameters/progression_os_2023-03-20.csv", sep=",", dec=".", header=TRUE)

design <- design |>
  mutate(
    n_ctrl = N_sample/2,
    n_trt  = N_sample/2
  ) |>
  select(hazard_ctrl, hazard_trt, hazard_after_prog, prog_rate_ctrl, prog_rate_trt, n_trt, n_ctrl, effect_size_ph, recruitment, prog_prop_trt, prog_prop_ctrl) |>
  unique()

results <- runSimulation(
  design,
  replications = N_sim,
  generate = generate_progression,
  analyse = function(condition, dat, fixed_objects=NULL){
    survival::survfit(survival::Surv(t, evt)~trt, dat)
  },
  extra_options = list(
    store_warning_seeds = TRUE,
    allow_na = TRUE
  )
)

results_progression <- bind_cols(
  design,
  tibble(
    survfit = map(results, \(x) x[[1]])
  )
)


# save --------------------------------------------------------------------

filename <- paste0(
  "data/theoretical_vs_simulated_",
  format(Sys.Date(), "%Y-%m-%d"),
  ".Rdata"
)

save(
  results_delayed,
  results_crossing,
  results_subgroup,
  results_progression,
  file=filename
)
