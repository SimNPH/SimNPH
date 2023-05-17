# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  ahr_6m   = analyse_ahr(type="AHR", max_time=m2d(6), alternative = "one.sided"),
  ahr_12m  = analyse_ahr(type="AHR", max_time=m2d(12), alternative = "one.sided"),
  ahr_360m  = analyse_ahr(type="AHR", max_time=m2d(360), alternative = "one.sided"),
  milestone_6m = analyse_milestone_survival(times=m2d(6), alternative = "one.sided"),
  milestone_12m = analyse_milestone_survival(times=m2d(12), alternative = "one.sided"),
  cox = analyse_coxph(alternative = "one.sided"),
  descriptive = analyse_describe()
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  # estimation
  milestone_6m = summarise_estimator(est=milestone_surv_ratio, real= milestone_survival_trt_6m/ milestone_survival_ctrl_6m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
  milestone_12m = summarise_estimator(est=milestone_surv_ratio, real=milestone_survival_trt_12m/milestone_survival_ctrl_12m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
  ahr_6m  = summarise_estimator(est=AHR, real=AHRoc_6m, lower=AHR_lower, upper=AHR_upper, null=1),
  ahr_12m = summarise_estimator(est=AHR, real=AHRoc_12m, lower=AHR_lower, upper=AHR_upper, null=1),
  cox_gahr_360m = summarise_estimator(est=hr, real=gAHR_360m, lower=hr_lower, upper=hr_upper, null=1),
  cox_ahr_360m = summarise_estimator(est=hr, real=AHR_360m, lower=hr_lower, upper=hr_upper, null=1),
  cox_ahroc_360m = summarise_estimator(est=hr, real=AHRoc_360m, lower=hr_lower, upper=hr_upper, null=1),
  cox_gahr_fup = summarise_estimator(est=hr, real=gAHR_360m, lower=hr_lower, upper=hr_upper, null=1),
  cox_ahr_fup = summarise_estimator(est=hr, real=AHR_360m, lower=hr_lower, upper=hr_upper, null=1),
  cox_ahroc_fup = summarise_estimator(est=hr, real=AHRoc_360m, lower=hr_lower, upper=hr_upper, null=1),
  cox = summarise_estimator(est=hr, real=AHRoc_360m, lower=hr_lower, upper=hr_upper, null=1),
  # tests
  peto_peto = summarise_test(alpha),
  pw_exp_3 = summarise_test(2*alpha),
  pw_exp_12 = summarise_test(2*alpha),
  # tests group sequential
  peto_peto_gs = summarise_group_sequential(),
  descriptive = summarise_describe()
)

