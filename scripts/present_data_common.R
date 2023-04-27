library(tidyverse)
library(SimNPH)


# metadata for methods ----------------------------------------------------

metadata <- tibble::tribble(
# method variable name      direction (one sided test)                                  name for output   estimator/test/gs-test  procedure category
  ~method,                 ~direction,                                                     ~method_name,                   ~type,                        ~category,
  "ahr_6m",                   "lower",                                        "average hazard hatio 6m",             "estimator",           "Average Hazard Ratio",
  "ahr_12m",                  "lower",                                       "average hazard ratio 12m",             "estimator",           "Average Hazard Ratio",
  "gahr_6m",                  "lower",                              "geometric average hazard ratio 6m",             "estimator", "Geometric Average Hazard Ratio",
  "gahr_12m",                 "lower",                             "geometric average hazard ratio 12m",             "estimator", "Geometric Average Hazard Ratio",
  "median_surv",             "higher",                                  "difference in median survival",             "estimator", "Diff. Median Survival",
  "milestone.milestone_6",   "higher",                                    "milestone survival ratio 6m",             "estimator", "Ratio of Milestone Survival",
  "milestone.milestone_12",  "higher",                                   "milestone survival ratio 12m",             "estimator", "Ratio of Milestone Survival",
  "rmst_diff_6m",            "higher",                                           "RMST (6m) difference",             "estimator", "Diff. RMST",
  "rmst_diff_12m",           "higher",                                          "RMST (12m) difference",             "estimator", "Diff. RMST",
  "cox",                      "lower",                                                 "Cox regression",             "estimator", "Cox regression",
  "weighted_cox_6m",          "lower",                                     "weighted Cox regression 6m",             "estimator", "Cox regression",
  "weighted_cox_12m",         "lower",                                    "weighted Cox regression 12m",             "estimator", "Cox regression",
  "aft_weibull",             "higher",                                "AFT model, Weibull distribution",             "estimator", "Accelerated Failure Time Model",
  "aft_lognormal",           "higher",                             "AFT model, log-normal distribution",             "estimator", "Accelerated Failure Time Model",
  "diff_med_weibull",        "higher",             "difference in median survival based on Weibull glm",             "estimator", "Diff. Median Survival",
  "pw_exp_3",                      NA,                 "test bases on a piecewise exponential model 3m",                  "test", "Piece-wise Exponential Model",
  "pw_exp_12",                     NA,                "test bases on a piecewise exponential model 12m",                  "test", "Piece-wise Exponential Model",
  "peto_peto",                     NA,                                                 "Peto-Peto test",                  "test", "Peto-Peto Test",
  "fh_0_0",                        NA,                          "weighted logrank test, FH weights 0,1",                  "test", "Fleming-Harrington Test",
  "fh_0_1",                        NA,                          "weighted logrank test, FH weights 0,1",                  "test", "Fleming-Harrington Test",
  "fh_1_0",                        NA,                          "weighted logrank test, FH weights 0,1",                  "test", "Fleming-Harrington Test",
  "fh_1_1",                        NA,                          "weighted logrank test, FH weights 0,1",                  "test", "Fleming-Harrington Test",
  "logrank",                       NA,                                                   "logrank test",                  "test", "Log-Rank Test",
  "max_combo",                     NA,                                                 "max-combo test",                  "test", "Max-Combo Test",
  "modest_6",                      NA,                          "modestly weighted logrank test, t*=6m",                  "test", "Modestly Weighted Test",
  "modest_8",                      NA,                          "modestly weighted logrank test, t*=8m",                  "test", "Modestly Weighted Test",
  "peto_peto_gs",                  NA,                        "Peto-Peto test, group sequential design", "group sequential test", "Fleming-Harrington Test",
  "fh_gs_0_0",                     NA, "weighted logrank test, FH weights 0,1, group sequential design", "group sequential test", "Fleming-Harrington Test",
  "fh_gs_0_1",                     NA, "weighted logrank test, FH weights 0,1, group sequential design", "group sequential test", "Fleming-Harrington Test",
  "fh_gs_1_0",                     NA, "weighted logrank test, FH weights 0,1, group sequential design", "group sequential test", "Fleming-Harrington Test",
  "fh_gs_1_1",                     NA, "weighted logrank test, FH weights 0,1, group sequential design", "group sequential test", "Fleming-Harrington Test",
  "logrank_gs",                    NA,                          "logrank test, group sequential design", "group sequential test", "Log-Rank Test",
  "max_combo_gs",                  NA,                        "max-combo test, group sequential design", "group sequential test", "Max-Combo Test",
  "modest_gs_6",                   NA, "modestly weighted logrank test, t*=6m, group sequential design", "group sequential test", "Modestly Weighted Test",
  "modest_gs_8",                   NA, "modestly weighted logrank test, t*=8m, group sequential design", "group sequential test", "Modestly Weighted Test"
)

# metadata for output columns ---------------------------------------------

col_labels <- tibble::tribble(
#                   column name,                                            added column label, where's the var from?,          unit
                       ~colname,                                                    ~col_label,                 ~from,         ~unit,
                  "recruitment",                                     "duration of recruitment",           "parameter",        "days",
                 "n_pat_design",                                          "number of patients",           "parameter",      "number",
               "interim_events",                       "number of events for interim analysis", "parameter (derived)",      "number",
                 "final_events",                         "number of events for final analysis", "parameter (derived)",      "number",
                       "n_ctrl",                           "number of patients in control arm", "parameter (derived)",      "number",
                        "n_trt",                         "number of patients in treatment arm", "parameter (derived)",      "number",
                  "hazard_ctrl",                              "hazard rate in the control arm",           "parameter",      "1/days",
               "censoring_prop",                    "proportion of randomly censored patients",           "parameter",    "fraction",
                        "delay",                          "delay of onset of treatment effect",           "parameter",        "days",
               "effect_size_ph",             "effect size (power of the logranktest under PH)",           "parameter",    "fraction",
                     "crossing",                       "time of crossing of the hazard curves",           "parameter",        "days",
                    "hr_before",           "hazard ratio before crossing of the hazard curves",           "parameter",          "HR",
                "prog_prop_trt",    "proportion of patients who progress in the treatment arm",           "parameter",    "fraction",
               "prog_prop_ctrl",      "proportion of patients who progress in the control arm",           "parameter",    "fraction",
              "hr_before_after",           "hazard ratio between before and after progression",           "parameter",          "HR",
                   "prevalence",                                  "prevalence of the subgroup",           "parameter",    "fraction",
         "hr_subgroup_relative",    "ratio of the HRs between the subgroup and the complement",           "parameter",       "HR/HR",
                   "hazard_trt", "hazard in the treatment arm after onset of treatment effect", "parameter (derived)",      "1/days",
            "target_median_trt",               "targeted median survival in the treatment arm", "parameter (derived)",        "days",
                    "target_hr",                                    "targeted PH hazard ratio", "parameter (derived)",          "HR",
            "random_withdrawal",                                   "rate of random withdrawal", "parameter (derived)",      "1/days",
          "median_survival_trt",                        "median survival in the treatment arm",   "true summary stat",        "days",
         "median_survival_ctrl",                          "median survival in the control arm",   "true summary stat",        "days",
                  "rmst_trt_6m",                              "RMST (6m) in the treatment arm",   "true summary stat",        "days",
                 "rmst_ctrl_6m",                                "RMST (6m) in the control arm",   "true summary stat",        "days",
                      "gAHR_6m",                         "geometric average hazard ratio (6m)",   "true summary stat",          "HR",
                       "AHR_6m",                                   "average hazard ratio (6m)",   "true summary stat",          "HR",
                 "rmst_trt_12m",                             "RMST (12m) in the treatment arm",   "true summary stat",        "days",
                "rmst_ctrl_12m",                               "RMST (12m) in the control arm",   "true summary stat",        "days",
                     "gAHR_12m",                        "geometric average hazard ratio (12m)",   "true summary stat",          "HR",
                      "AHR_12m",                                  "average hazard ratio (12m)",   "true summary stat",          "HR",
    "milestone_survival_trt_6m",                "milestone survival in the treatment arm (6m)",   "true summary stat",    "fraction",
   "milestone_survival_ctrl_6m",                  "milestone survival in the control arm (6m)",   "true summary stat",    "fraction",
   "milestone_survival_trt_12m",               "milestone survival in the treatment arm (12m)",   "true summary stat",    "fraction",
  "milestone_survival_ctrl_12m",                 "milestone survival in the control arm (12m)",   "true summary stat",    "fraction",
            "descriptive.n_pat",                                  "average number of patients",    "descriptive stat",      "number",
     "descriptive.max_followup",                                  "average max follow up time",    "descriptive stat",      "number",
              "descriptive.evt",                                    "average number of events",    "descriptive stat",      "number",
         "descriptive.evt_ctrl",                 "average number of events in the control arm",    "descriptive stat",      "number",
          "descriptive.evt_trt",               "average number of events in the treatment arm",    "descriptive stat",      "number",
       "descriptive.study_time",                                     "average length of study",    "descriptive stat",        "days",
         "descriptive.sd_n_pat",                    "standard deviation of number of patients",    "descriptive stat",      "number",
  "descriptive.sd_max_followup",                    "standard deviation of max follow up time",    "descriptive stat",      "number",
           "descriptive.sd_evt",                      "standard deviation of number of events",    "descriptive stat",      "number",
      "descriptive.sd_evt_ctrl",   "standard deviation of number of events in the control arm",    "descriptive stat",      "number",
       "descriptive.sd_evt_trt", "standard deviation of number of events in the treatment arm",    "descriptive stat",      "number",
    "descriptive.sd_study_time",                       "standard deviation of length of study",    "descriptive stat",        "days",
                 "REPLICATIONS",                                      "number of replications",  "SimDesign metadata",      "number",
                     "SIM_TIME",                                    "scenario simulation time",  "SimDesign metadata",     "seconds",
                    "COMPLETED",                                      "simulation finished at",  "SimDesign metadata",        "date",
                         "SEED",                                                 "random seed",  "SimDesign metadata", NA_character_,
                       "method",                                              "estimator/test",               "pivot", NA_character_,
                     "mean_est",                                         "mean point estimate",           "summarise", NA_character_,
                   "median_est",                                       "median point estimate",           "summarise", NA_character_,
                       "sd_est",                        "standard deviation of point estimate",           "summarise", NA_character_,
                         "bias",                                                        "bias",           "summarise", NA_character_,
                      "sd_bias",                                  "standard deviation of bias",           "summarise", NA_character_,
                          "mse",                                          "mean squared error",           "summarise", NA_character_,
                       "sd_mse",                    "standard deviation of mean squared error",           "summarise", NA_character_,
                          "mae",                                         "mean absolute error",           "summarise", NA_character_,
                       "sd_mae",                   "standard deviation of mean absolute error",           "summarise", NA_character_,
                     "coverage",                                                 "CI coverage",           "summarise",    "fraction",
                   "null_cover",              "proportions of CIs that contain the null value",           "summarise",    "fraction",
                  "cover_lower",      "proportion of lower CI boundaries below the true value",           "summarise",    "fraction",
                  "cover_upper",      "proportion of upper CI boundaries above the true value",           "summarise",    "fraction",
                   "null_lower",      "proportion of lower CI boundaries below the null value",           "summarise",    "fraction",
                   "null_upper",      "proportion of upper CI boundaries above the null value",           "summarise",    "fraction",
                        "width",                                            "average CI width",           "summarise", NA_character_,
                     "sd_width",                              "standard deviation of CI width",           "summarise", NA_character_,
                      "mean_sd",                        "average estimated standard deviation",           "summarise", NA_character_,
                        "sd_sd",          "standard deviation of estimated standard deviation",           "summarise", NA_character_,
                   "mean_n_pat",                                  "average number of patients",           "summarise",      "number",
                     "sd_n_pat",                    "standard deviation of number of patients",           "summarise",      "number",
                   "mean_n_evt",                                    "average number of events",           "summarise",      "number",
                     "sd_n_evt",                      "standard deviation of number of events",           "summarise",      "number",
                    "N_missing",                                 "number of missing estimates",           "summarise",      "number",
                            "N",                                      "number of replications",           "summarise",      "number",
                 "N_missing_CI",                                       "number of missing CIs",           "summarise",      "number",
              "N_missing_upper",                           "number of missing upper CI limits",           "summarise",      "number",
              "N_missing_lower",                           "number of missing lower CI limits",           "summarise",      "number",
                 "N_missing_sd",      "number of missing estimates for the standard deviation",           "summarise",      "number",
              "N_missing_n_pat",      "number or replications with missing number of patients",           "summarise",      "number",
              "N_missing_n_evt",        "number of replications with missing number of events",           "summarise",      "number",
              "N_missing_0.025",                                  "number of missing p-values",           "summarise",      "number",
                    "rejection",                                              "rejection rate",           "summarise",    "fraction",
                   "study_time",                                                  "study time",           "summarise",        "days",
                     "followup",                                                "max followup",           "summarise",        "days",
                "sd_study_time",                            "standard deviation of study time",           "summarise",        "days",
                  "sd_followup",                              "standard deviation of followup",           "summarise",        "days",
          "N_missing_rejection",      "number of missing rejection (group sequential designs)",           "summarise",      "number",
               "N_missing_npat",       "missing number of patients (group sequential designs)",           "summarise",      "number",
               "N_missing_nevt",         "missing number of events (group sequential designs)",           "summarise",      "number",
         "N_missing_study_time",               "missing study time (group sequential designs)",           "summarise",      "number",
           "N_missing_followup",             "missing max followup (group sequential designs)",           "summarise",      "number",
 "ci_based_one_sided_rejection",                "rejection rate of one sided test based on CI",                "post",    "fraction",
                    "direction",                     "direction of one sided test based on CI",                "post",        "text",
                  "method_name",                                    "estimator/test procedure",                "post",        "text",
                         "type",           "estimator/fixed sample test/group sequential test",                "post",        "text",
)

# transform data ----------------------------------------------------------

## First transform to long format: Not run
# results_long <- results |>
#   results_pivot_longer()

#' Add metadata and column labels to simulation results
#'
#' This function adds metadata and column labels to a long-format results table.
#'
#' @param results_long A data frame containing the long-format results table.
#' @param metadata A data frame containing the metadata to be added to the results table.
#' @param col_labels A data frame containing the column labels to be added to the results table.
#'
#' @return A data frame with added metadata and column labels.
#' @export
#'

#' @import dplyr
#' @importFrom tidyr left_join
#' @importFrom purrr imap
#' @importFrom stringr str_c
add_meta_data <- function(results_long,metadata,col_labels){
  # sort those columsn first
  colnames_first <- col_labels |>
    filter(from == "parameter") |>
    pull(colname)

  results_long <- results_long |>
    left_join(metadata, by = "method") |> # add metadata
    mutate(
      # one sided test based on CI
      ci_based_one_sided_rejection = case_when(
        direction == "lower" ~ 1 - null_upper,
        direction == "higher" ~ 1 - null_lower,
        TRUE ~ NA_real_
      ),
      # combine columsn with different names acros estimators/tests/gs-tests
      # coalesce: take first value, if it is missing take second value
      mean_n_pat = coalesce(mean_n_pat, n_pat),
      mean_n_evt = coalesce(mean_n_evt, n_evt),
      sd_n_pat   = coalesce(sd_n_pat, sd_npat),
      sd_n_evt   = coalesce(sd_n_evt, sd_nevt),
      study_time = coalesce(study_time, descriptive.study_time),
      followup   = coalesce(followup, descriptive.max_followup),
      rejection  = coalesce(rejection, rejection_0.025),
      rejection  = coalesce(rejection, ci_based_one_sided_rejection)
    ) |>
    select( # deselect columsn
      -n_pat, -n_evt, -sd_npat, -sd_nevt, rejection_0.025
    ) |>
    select( # reorder columns
      method_name, any_of(colnames_first), everything()
    )


  # add labels to columns (shown in rstudio viewer) -------------------------

  results_long <- results_long |>
    imap(\(x,i){
      label <- col_labels |>
        filter(colname==i) |>
        pull(col_label)
      attr(x, "label") <- label
      x
    }) |>
    as_tibble()

  return(results_long)
}
