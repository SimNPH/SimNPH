# define analyse functions ------------------------------------------------

my_analyse <- list(
  # estimation
  rmst_diff_24m = analyse_rmst_diff(max_time = m2d(24), alternative = "one.sided"),
  rmst_diff_36m = analyse_rmst_diff(max_time = m2d(36), alternative = "one.sided"),
  descriptive = analyse_describe()
)

my_analyse <- wrap_all_in_trycatch(my_analyse)

# define summaries --------------------------------------------------------

my_summarise <- create_summarise_function(
  # estimation
  rmst_diff_24m = summarise_estimator(est=rmst_diff, real=rmst_trt_24m - rmst_ctrl_24m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  rmst_diff_36m = summarise_estimator(est=rmst_diff, real=rmst_trt_36m - rmst_ctrl_36m, lower=rmst_diff_lower, upper=rmst_diff_upper, null=0),
  descriptive = summarise_describe()
)

