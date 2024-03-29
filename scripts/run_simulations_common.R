# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  ahr_6m   = analyse_ahr(type="AHR", max_time=m2d(6), alternative = "one.sided"),
  ahr_12m  = analyse_ahr(type="AHR", max_time=m2d(12), alternative = "one.sided"),
  gahr_6m  = analyse_ahr(type="gAHR", max_time=m2d(6), alternative = "one.sided"),
  gahr_12m = analyse_ahr(type="gAHR", max_time=m2d(12), alternative = "one.sided"),
  median_surv = analyse_diff_median_survival(alternative = "one.sided"),
  milestone_6m = analyse_milestone_survival(times=m2d(6), alternative = "one.sided"),
  milestone_12m = analyse_milestone_survival(times=m2d(12), alternative = "one.sided"),
  rmst_diff_6m  = analyse_rmst_diff(max_time = m2d(6), alternative = "one.sided"),
  rmst_diff_12m = analyse_rmst_diff(max_time = m2d(12), alternative = "one.sided"),
  cox = analyse_coxph(alternative = "one.sided"),
  weighted_cox_6m  = analyse_coxph_weighted(type="G", max_time=m2d(6), alternative = "one.sided"),
  weighted_cox_12m = analyse_coxph_weighted(type="G", max_time=m2d(12), alternative = "one.sided"),
  aft_weibull = analyse_aft(dist="weibull", alternative = "one.sided"),
  aft_lognormal = analyse_aft(dist="lognormal", alternative = "one.sided"),
  diff_med_weibull = analyse_weibull(alternative = "one.sided"),
  # tests
  pw_exp_3  = analyse_piecewise_exponential(cuts=m2d(seq(0, 240, by= 3)), testing_only=TRUE),
  pw_exp_12 = analyse_piecewise_exponential(cuts=m2d(seq(0, 240, by=12)), testing_only=TRUE),
  peto_peto = analyse_gehan_wilcoxon(alternative = "one.sided"),
  fh_0_0 = analyse_logrank_fh_weights(rho=0, gamma=0, alternative = "one.sided"),
  fh_0_1 = analyse_logrank_fh_weights(rho=0, gamma=1, alternative = "one.sided"),
  fh_1_0 = analyse_logrank_fh_weights(rho=1, gamma=0, alternative = "one.sided"),
  fh_1_1 = analyse_logrank_fh_weights(rho=1, gamma=1, alternative = "one.sided"),
  logrank = analyse_logrank(alternative = "one.sided"),
  max_combo = analyse_maxcombo(alternative = "one.sided"),
  modest_6 = analyse_modelstly_weighted(t_star=m2d(6)),
  modest_8 = analyse_modelstly_weighted(t_star=m2d(8)),
  # tests group sequential
  peto_peto_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_gehan_wilcoxon(alternative = "one.sided")
  ),
  fh_gs_0_0 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=0, gamma=0, alternative = "one.sided")
  ),
  fh_gs_0_1 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=0, gamma=1, alternative = "one.sided")
  ),
  fh_gs_1_0 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=1, gamma=0, alternative = "one.sided")
  ),
  fh_gs_1_1 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=1, gamma=1, alternative = "one.sided")
  ),
  logrank_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank(alternative = "one.sided")
  ),
  max_combo_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_maxcombo(alternative = "one.sided")
  ),
  modest_gs_6 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_modelstly_weighted(t_star=m2d(6))
  ),
  modest_gs_8 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_modelstly_weighted(t_star=m2d(8))
  ),
  descriptive = analyse_describe()
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  # estimation
  ahr_6m  = summarise_estimator(est=AHR, real=AHR_6m, lower=AHR_lower, upper=AHR_upper, null=1),
  ahr_12m = summarise_estimator(est=AHR, real=AHR_12m, lower=AHR_lower, upper=AHR_upper, null=1),
  gahr_6m  = summarise_estimator(est=gAHR, real=gAHR_6m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  gahr_12m = summarise_estimator(est=gAHR, real=gAHR_12m, lower=gAHR_lower, upper=gAHR_upper, null=1),
  median_surv = summarise_estimator(est=diff_Q, real=median_survival_trt-median_survival_ctrl, lower=diff_Q_lower, upper=diff_Q_upper, null=0),
  milestone_6m = summarise_estimator(est=milestone_surv_ratio, real= milestone_survival_ctrl_6m/ milestone_survival_trt_6m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
  milestone_12m = summarise_estimator(est=milestone_surv_ratio, real=milestone_survival_ctrl_12m/milestone_survival_trt_12m, lower=milestone_surv_ratio_lower, upper=milestone_surv_ratio_upper, null=1),
  rmst_diff_6m  = summarise_estimator(est=rmst_diff, real=rmst_trt_6m -rmst_ctrl_6m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  rmst_diff_12m = summarise_estimator(est=rmst_diff, real=rmst_trt_12m-rmst_ctrl_12m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  cox = summarise_estimator(est=hr, real=gAHR_12m, lower=hr_lower, upper=hr_upper, null=1),
  weighted_cox_6m  = summarise_estimator(est=hr, real=gAHR_6m, lower=hr_lower, upper=hr_upper, null=1),
  weighted_cox_12m = summarise_estimator(est=hr, real=gAHR_12m, lower=hr_lower, upper=hr_upper, null=1),
  aft_weibull = summarise_estimator(est=coef, real=NA_real_, lower=lower, upper=upper, null=0),
  aft_lognormal = summarise_estimator(est=coef, real=NA_real_, lower=lower, upper=upper, null=0),
  diff_med_weibull = summarise_estimator(est=diff_med_est, real=median_survival_trt-median_survival_ctrl,lower=diff_med_lower, upper=diff_med_upper, null=0),
  # tests
  pw_exp_3 = summarise_test(2*alpha),
  pw_exp_12 = summarise_test(2*alpha),
  peto_peto = summarise_test(alpha),
  fh_0_0 = summarise_test(alpha),
  fh_0_1 = summarise_test(alpha),
  fh_1_0 = summarise_test(alpha),
  fh_1_1 = summarise_test(alpha),
  logrank = summarise_test(alpha),
  max_combo = summarise_test(alpha),
  modest_6 = summarise_test(alpha),
  modest_8 = summarise_test(alpha),
  # tests group sequential
  peto_peto_gs = summarise_group_sequential(),
  fh_gs_0_0 = summarise_group_sequential(),
  fh_gs_0_1 = summarise_group_sequential(),
  fh_gs_1_0 = summarise_group_sequential(),
  fh_gs_1_1 = summarise_group_sequential(),
  logrank_gs = summarise_group_sequential(),
  max_combo_gs = summarise_group_sequential(),
  modest_gs_6 = summarise_group_sequential(),
  modest_gs_8 = summarise_group_sequential(),
  descriptive = summarise_describe()
)

