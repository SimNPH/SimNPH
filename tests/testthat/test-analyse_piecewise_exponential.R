test_that("analyse_piecewise_exponential outputs plausible data.frame for delayed effect", {
  capture_output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      head(1)
  )
  dat <- generate_delayed_effect(condition)
  my_analyse <- analyse_piecewise_exponential(cuts=c(30))
  res <- my_analyse(condition, dat)


  expect_named(res, c("hr_s", "p_s", "hr_s_lower", "hr_s_upper", "N_pat", "N_evt", "interval_table"), ignore.order = TRUE)
  expect(all( (res$p_s >= 0) %in% c(TRUE, NA)), "all p values >= 0 or NA")
  expect(all( (res$p_s <= 1) %in% c(TRUE, NA)), "all p values <= 1 or NA")
})
