# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  ahr_6m   = analyse_ahr(type="AHR", max_time=m2d(6)),
  ahr_12m  = analyse_ahr(type="AHR", max_time=m2d(12)),
  gahr_6m  = analyse_ahr(type="gAHR", max_time=m2d(6)),
  gahr_12m = analyse_ahr(type="gAHR", max_time=m2d(12)),
  median_surv = analyse_diff_median_survival(),
  milestone = analyse_milestone_survival(times=m2d(c(6, 12))),
  rmst_diff_6m  = analyse_rmst_diff(max_time = m2d(6)),
  rmst_diff_12m = analyse_rmst_diff(max_time = m2d(12)),
  cox = analyse_coxph(),
  weighted_cox_6m  = analyse_coxph_weighted(type="G", max_time=m2d(6)),
  weighted_cox_12m = analyse_coxph_weighted(type="G", max_time=m2d(12)),
  # tests
  pw_exp_3  = analyse_piecewise_exponential(cuts=m2d(seq(0, 240, by= 3)), testing_only=TRUE),
  pw_exp_12 = analyse_piecewise_exponential(cuts=m2d(seq(0, 240, by=12)), testing_only=TRUE),
  peto_peto = analyse_gehan_wilcoxon(),
  fh_0_0 = analyse_logrank_fh_weights(rho=0, gamma=0),
  fh_0_1 = analyse_logrank_fh_weights(rho=0, gamma=1),
  fh_1_0 = analyse_logrank_fh_weights(rho=1, gamma=0),
  fh_1_1 = analyse_logrank_fh_weights(rho=1, gamma=1),
  logrank = analyse_logrank(),
  max_combo = analyse_maxcombo(),
  modest_6 = analyse_modelstly_weighted(t_star=m2d(6)),
  modest_8 = analyse_modelstly_weighted(t_star=m2d(8)),
  # tests group sequential
  peto_peto_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_gehan_wilcoxon()
  ),
  fh_gs_0_0 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=0, gamma=0)
  ),
  fh_gs_0_1 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=0, gamma=1)
  ),
  fh_gs_1_0 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=1, gamma=0)
  ),
  fh_gs_1_1 = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank_fh_weights(rho=1, gamma=1)
  ),
  logrank_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_logrank()
  ),
  max_combo_gs = analyse_group_sequential(
    followup = c(condition$interim_events, condition$final_events),
    followup_type = c("event", "event"),
    alpha = nominal_alpha,
    analyse_functions = analyse_maxcombo()
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
  ahr_6m  = summarise_estimator(est=AHR, real=AHR_6m),
  ahr_12m = summarise_estimator(est=AHR, real=AHR_12m),
  gahr_6m  = summarise_estimator(est=gAHR, real=gAHR_6m),
  gahr_12m = summarise_estimator(est=gAHR, real=gAHR_12m),
  median_surv = summarise_estimator(est=diff_Q, real=median_survival_trt/median_survival_ctrl),
  milestone = summarise_estimator(est=milestone_surv_ratio[1], real= milestone_survival_ctrl_6m/ milestone_survival_trt_6m, name="milestone_6" ),
  milestone = summarise_estimator(est=milestone_surv_ratio[2], real=milestone_survival_ctrl_12m/milestone_survival_trt_12m, name="milestone_12"),
  rmst_diff_6m  = summarise_estimator(est=rmst_diff, real=rmst_trt_6m -rmst_ctrl_6m),
  rmst_diff_12m = summarise_estimator(est=rmst_diff, real=rmst_trt_12m-rmst_ctrl_12m),
  cox = summarise_estimator(est=hr, real=gAHR_12m, lower=hr_lower, upper=hr_upper),
  weighted_cox_6m  = summarise_estimator(est=hr, real=gAHR_6m, lower=hr_lower, upper=hr_upper),
  weighted_cox_12m = summarise_estimator(est=hr, real=gahr_6m, lower=hr_lower, upper=hr_upper),
  # tests
  pw_exp_3 = summarise_test(alpha),
  pw_exp_12 = summarise_test(alpha),
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
