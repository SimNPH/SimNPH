test_that("creating summarise function from many functions works", {

  Summarise <- create_summarise_function(
    method1 = `attr<-`(function(condition, results, fixed_objects){
      mean(results$value1)
    }, "name", "mean_val_1"),
    method2 = function(condition, results, fixed_objects){
      mean(results$value2)
    },
    method1 = function(condition, results, fixed_objects){
      median(results$value1)
    },
    method3 = function(condition, results, fixed_objects){
      data.frame(x=2, y=2)
    }
  )

  condition <- numeric(0)

  results <- list(
    list(method1=list(value1= 1), method2=list(value2= 1)),
    list(method1=list(value1= 1), method2=list(value2= 1)),
    list(method1=list(value1=10), method2=list(value2=10))
  )

  summary <- Summarise(condition, results)

  expect_type(Summarise, "closure")
  expect_equal(as.numeric(summary[1, ]), c(4,4,1,NA_real_))
  expect_type(summary, "list")
  expect_s3_class(summary, "data.frame")
  expect_named(summary, c("method1.mean_val_1", "method2", "method1", "method3"))
})

test_that("creating a summarise function for an estimator works", {

  capture.output(
    # generate the design matrix and append the true summary statistics
    condition <- desing_skeleton_delayed_effect() |>
      tail(4) |>
      head(1) |>
      true_summary_statistics_delayed_effect(cutoff_stats = 15)
  )

  # create some summarise functions
  summarise_all <- create_summarise_function(
    coxph=summarise_estimator(hr, gAHR, hr_lower, hr_upper),
    coxph=summarise_estimator(hr, hazard_trt/hazard_ctrl, hr_lower, hr_upper),
    coxph=summarise_estimator(exp(coef), gAHR),
    coxph=summarise_estimator(hr, NA_real_)
  )

  # runs simulations
  capture.output(
    type="output",
    capture.output(
      type="message",
      withr::with_seed(1, {
        sim_results <- runSimulation(
          design=condition,
          replications=10,
          generate=generate_delayed_effect,
          analyse=list(
            coxph=analyse_coxph
          ),
          summarise = summarise_all,
          save = FALSE
        )
      })
    )
  )

  expected_names <- expand.grid("coxph.", c("bias", "var", "mse", "mae", "coverage", "width"), c("", paste0(1:3, "."))) |>
    subset(select=c(1,3,2)) |>
    apply(1, paste, collapse="") |>
    unname()

  expected_names <- c(names(condition), expected_names, c("REPLICATIONS", "SIM_TIME", "COMPLETED", "SEED"))

  expect_named(sim_results, expected_names)

  expect(all(is.na(sim_results[, c("coxph.3.mse", "coxph.3.mae", "coxph.3.bias", "coxph.3.coverage")])), "summary results depending on the true value should be missing when the true value is not given")
  expect(all(is.na(sim_results[, c("coxph.2.coverage", "coxph.2.width")])), "summary results depending on the CI should be missing if no CI boundaries are given")

})
