# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  milestone_6m = analyse_milestone_survival(times=m2d(6), alternative = "one.sided"),
  milestone_12m = analyse_milestone_survival(times=m2d(12), alternative = "one.sided"),
  # tests
  peto_peto = analyse_gehan_wilcoxon(alternative = "one.sided"),
  # tests group sequential
  peto_peto_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_gehan_wilcoxon(alternative = "one.sided")
  ),
  descriptive = analyse_describe()
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  # estimation
  milestone_6m = summarise_estimator(est=milestone_surv_ratio, real= milestone_survival_ctrl_6m/ milestone_survival_trt_6m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
  milestone_12m = summarise_estimator(est=milestone_surv_ratio, real=milestone_survival_ctrl_12m/milestone_survival_trt_12m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
  # tests
  peto_peto = summarise_test(alpha),
  # tests group sequential
  peto_peto_gs = summarise_group_sequential(),
  descriptive = summarise_describe()
)

