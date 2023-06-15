## Load Data
library(tidyverse)
library(SimNPH)
library(ggplot2)
library(patchwork)

data_folder <- "C:/Users/floria14/Nextcloud/NPH Project/Simulations/data/merged_2023-04-26/"

source("./scripts/present_data_common.R")
results_progression <- readRDS(paste0(data_folder,"results_progression_2023-04-26.Rds")) |>
  results_pivot_longer() |>
  mutate(rejection_0.025=rejection_0.05) |> ## Fake it until you make it!
  add_meta_data(metadata,col_labels) |>
  mutate(across(starts_with(c("median_survival",
                              "recruitment",
                              "rmst",
                              "crossing",
                              "descriptive.study_time",
                              "study_time")),~.x/30.4375),
         average_study_time = coalesce(study_time,descriptive.study_time)) |> ## Transform days to months
  mutate(pfs_median_survival_ctrl = round(r2m(hazard_ctrl))) |>
  filter(method != "peto_peto") # same as fh_1_0

results_crossing <- readRDS(paste0(data_folder,"results_crossing_2023-04-26.Rds")) |>
  results_pivot_longer() |>
  add_meta_data(metadata,col_labels) |>
  mutate(across(starts_with(c("median_survival",
                              "recruitment",
                              "rmst",
                              "crossing",
                              "descriptive.study_time",
                              "study_time")),~.x/30.4375),
         average_study_time = coalesce(study_time,descriptive.study_time)) |> ## Transform days to months
  arrange(method) |> ## we don't need this anymore
  mutate(pfs_median_survival_ctrl = round(r2m(hazard_ctrl))) |>
  filter(method != "peto_peto") # same as fh_1_0


results_delay <- readRDS(paste0(data_folder,"results_delayed_2023-04-26.Rds")) |>
  results_pivot_longer() |>
  add_meta_data(metadata,col_labels) |>
  mutate(across(starts_with(c("median_survival",
                              "recruitment",
                              "rmst",
                              "delay",
                              "descriptive.study_time",
                              "study_time")),~.x/30.4375),
         average_study_time = coalesce(study_time,descriptive.study_time)) |> ## Transform days to months
  arrange(method) |> ## we don't need this anymore
  mutate(pfs_median_survival_ctrl = round(r2m(hazard_ctrl))) |>
  filter(method != "peto_peto") # same as fh_1_0

## Test plots
methods_ci_power <- c('cox','logrank',
                      'ahr_12m','ahr_6m',
                      'diff_med_weibull','median_surv',
                      'gahr_12m','gahr_6m',
                      'rmst_diff_12m','rmst_diff_6m',
                      'milestone_12m','milestone_milestone_6m',
                      "aft_lognormal","aft_weibull")
results_progression |>
  plot_sim_results(methods = methods_ci_power,
                   parameter_x = c('hr_before_after','prog_prop_ctrl','prog_prop_trt'),
                   parameter_y = 'rejection',
                   split_var = 1)

results_crossing |>
  plot_sim_results(methods = methods_ci_power,
                   parameter_x = c('hr_before','crossing'),
                   parameter_y = 'rejection')

results_delay |>
  plot_sim_results(methods = methods_ci_power,
                   parameter_x = c('delay'),
                   parameter_y = 'rejection')

# n <- length(xvars)
# lastvar <- xvars[[n]]
# data <- data |>
#   group_by(method,!!!facet_vars_x_sym,!!!facet_vars_y_sym,!!!xvars[-n]) |>
#   group_modify(~add_row(.x,.before = 1)) |>
#   mutate(!!lastvar := ifelse(is.na(!!yvar),!!lastvar + .0,!!lastvar)) |>
#   ungroup()
