# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  rmst_diff_6m = analyse_rmst_diff(max_time = m2d(6), alternative = "one.sided"),
  ahr_6m   = analyse_ahr(type="AHR", max_time=m2d(6), alternative = "one.sided"),
  gahr_6m  = analyse_ahr(type="gAHR", max_time=m2d(6), alternative = "one.sided"),
  milestone_6m = analyse_milestone_survival(times=m2d(6), alternative = "one.sided"),

  rmst_diff_12m = analyse_rmst_diff(max_time = m2d(12), alternative = "one.sided"),
  ahr_12m   = analyse_ahr(type="AHR", max_time=m2d(12), alternative = "one.sided"),
  gahr_12m  = analyse_ahr(type="gAHR", max_time=m2d(12), alternative = "one.sided"),
  milestone_12m = analyse_milestone_survival(times=m2d(12), alternative = "one.sided"),

  rmst_diff_24m = analyse_rmst_diff(max_time = m2d(24), alternative = "one.sided"),
  ahr_24m   = analyse_ahr(type="AHR", max_time=m2d(24), alternative = "one.sided"),
  gahr_24m  = analyse_ahr(type="gAHR", max_time=m2d(24), alternative = "one.sided"),
  milestone_24m = analyse_milestone_survival(times=m2d(24), alternative = "one.sided"),

  rmst_diff_36m = analyse_rmst_diff(max_time = m2d(36), alternative = "one.sided"),
  ahr_36m   = analyse_ahr(type="AHR", max_time=m2d(36), alternative = "one.sided"),
  gahr_36m  = analyse_ahr(type="gAHR", max_time=m2d(36), alternative = "one.sided"),
  milestone_36m = analyse_milestone_survival(times=m2d(36), alternative = "one.sided"),

  descriptive = analyse_describe()
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  # estimation
  rmst_diff_6m = summarise_estimator(est=rmst_diff, real=rmst_trt_6m - rmst_ctrl_6m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  ahr_6m  = summarise_estimator(est=AHR, real=AHR_6m, lower=AHR_lower, upper=AHR_upper, null=1),
  gahr_6m  = summarise_estimator(est=gAHR, real=gAHR_6m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  milestone_6m = summarise_estimator(est=milestone_surv_ratio, real= milestone_survival_ctrl_6m/ milestone_survival_trt_6m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),

  rmst_diff_12m = summarise_estimator(est=rmst_diff, real=rmst_trt_12m - rmst_ctrl_12m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  ahr_12m  = summarise_estimator(est=AHR, real=AHR_12m, lower=AHR_lower, upper=AHR_upper, null=1),
  gahr_12m  = summarise_estimator(est=gAHR, real=gAHR_12m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  milestone_12m = summarise_estimator(est=milestone_surv_ratio, real= milestone_survival_ctrl_12m/ milestone_survival_trt_12m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),

  rmst_diff_24m = summarise_estimator(est=rmst_diff, real=rmst_trt_24m - rmst_ctrl_24m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  ahr_24m  = summarise_estimator(est=AHR, real=AHR_24m, lower=AHR_lower, upper=AHR_upper, null=1),
  gahr_24m  = summarise_estimator(est=gAHR, real=gAHR_24m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  milestone_24m = summarise_estimator(est=milestone_surv_ratio, real= milestone_survival_ctrl_24m/ milestone_survival_trt_24m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),

  rmst_diff_36m = summarise_estimator(est=rmst_diff, real=rmst_trt_36m - rmst_ctrl_36m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  ahr_36m  = summarise_estimator(est=AHR, real=AHR_36m, lower=AHR_lower, upper=AHR_upper, null=1),
  gahr_36m  = summarise_estimator(est=gAHR, real=gAHR_36m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  milestone_36m = summarise_estimator(est=milestone_surv_ratio, real= milestone_survival_ctrl_36m/ milestone_survival_trt_36m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),

  descriptive = summarise_describe()
)

