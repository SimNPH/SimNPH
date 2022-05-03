test_that("analyse_rmst outputs plausible data.frame for delayed effect", {
  capture_output(
    condition <- desing_skeleton_delayed_effect() |>
      head(1)
  )
  dat <- generate_delayed_effect(condition)
  res <- analyse_rmst(condition, dat)

  expect_true(hasName(res, "p"), label="result has the column p")
  expect_lte(res$p, 1, label="p value le 1")
  expect_gte(res$p, 0, label="p value ge 0")
})
