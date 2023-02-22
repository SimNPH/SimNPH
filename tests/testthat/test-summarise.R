test_that("creating summarise function from many functions works", {

  Summarise <- create_summarise_function(
    method1 = `attr<-`(function(condition, results, fixed_objects){
      data.frame(val=mean(results$value1))
    }, "name", "mean_val_1"),
    method2 = function(condition, results, fixed_objects){
      data.frame(val=mean(results$value2))
    },
    method1 = function(condition, results, fixed_objects){
      data.frame(val=median(results$value1))
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
  expect_equal(as.numeric(summary[1, ]), c(4,4,1))
  expect_type(summary, "list")
  expect_s3_class(summary, "data.frame")
  expect_named(summary, c("method1.mean_val_1.val", "method2.val", "method1.val"))


  Summarise_err <- create_summarise_function(
    method1 = `attr<-`(function(condition, results, fixed_objects){
      data.frame(val=mean(results$value1))
    }, "name", "mean_val_1"),
    method2 = function(condition, results, fixed_objects){
      data.frame(val=mean(results$value2))
    },
    method1 = function(condition, results, fixed_objects){
      stop("test")
    }
  )

  summary_2 <- Summarise_err(condition, results)
  expect_named(summary_2, c("method1.mean_val_1.val", "method2.val", "method1.err"))
  expect_equal(summary_2$method1.err, "test")
})

test_that("creating a summarise function for an estimator works", {

  capture.output(
    # generate the design matrix and append the true summary statistics
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      tail(4) |>
      head(1) |>
      true_summary_statistics_delayed_effect(cutoff_stats = 15)
  )

  # create some summarise functions
  summarise_all <- create_summarise_function(
    coxph=summarise_estimator(hr, gAHR_15, hr_lower, hr_upper, name="gAHR"),
    coxph=summarise_estimator(hr, hazard_trt/hazard_ctrl, hr_lower, hr_upper, name="hr"),
    coxph=summarise_estimator(exp(coef), gAHR_15),
    coxph=summarise_estimator(hr, NA_real_)
  )

  # runs simulations
  capture_warnings(
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
              coxph=analyse_coxph()
            ),
            summarise = summarise_all,
            save = FALSE
          )
        })
      )
    )
  )

  expected_names <- expand.grid("coxph.", c("bias", "sd_bias", "var", "mse", "sd_mse", "N_missing", "N", "mae", "sd_mae", "coverage", "width", "sd_width", "N_missing_CI"), c("gAHR.", "hr.", "", "1.")) |>
    subset(select=c(1,3,2)) |>
    apply(1, paste, collapse="") |>
    unname()

  expected_names <- c(names(condition), expected_names, c("REPLICATIONS", "SIM_TIME", "COMPLETED", "SEED"))

  expect_named(sim_results, expected_names, ignore.order = TRUE)

  expect(all(is.na(sim_results[, c("coxph.1.mse", "coxph.1.mae", "coxph.1.bias", "coxph.1.coverage")])), "summary results depending on the true value should be missing when the true value is not given")
  expect(all(is.na(sim_results[, c("coxph.coverage", "coxph.width")])), "summary results depending on the CI should be missing if no CI boundaries are given")

})

test_that("generic summarise for tests works", {
  capture.output(
    condition <- merge(
      assumptions_delayed_effect(),
      design_fixed_followup(),
      by=NULL
    ) |>
      tail(4) |>
      head(1)
  )
  summarise_all <- create_summarise_function(
    logrank=summarise_test(alpha=c(0.95, 0.99)),
    logrank=summarise_test(alpha=c(0.9), name="innovative")
  )

  # runs simulations
  capture.output(
    suppressMessages(
      sim_results <- runSimulation(
        design=condition,
        replications=10,
        generate=generate_delayed_effect,
        analyse=list(
          logrank=analyse_logrank()
        ),
        summarise = summarise_all
      )
    )
  )

  expect(
    all(hasName(sim_results, c(
      "logrank.rejection_0.95", "logrank.rejection_0.99", "logrank.innovative.rejection_0.9",
      "logrank.N_missing_0.95", "logrank.N_missing_0.99", "logrank.innovative.N_missing_0.9",
      "logrank.N"
      ))),
    "expected names not present in sim_results"
  )

  expect_gte(sim_results$logrank.rejection_0.95,           0)
  expect_gte(sim_results$logrank.rejection_0.99,           0)
  expect_gte(sim_results$logrank.innovative.rejection_0.9, 0)
  expect_lte(sim_results$logrank.rejection_0.95,           1)
  expect_lte(sim_results$logrank.rejection_0.99,           1)
  expect_lte(sim_results$logrank.innovative.rejection_0.9, 1)
})

test_that("missings are treated correctly for summarise estimator", {
  my_summarise <- summarise_estimator(est, real, lower, upper)

  condition_and_results <- tibble::tribble(
    ~real,     ~est,    ~lower,  ~upper,
        0,      0.1,       -1,        1,
        0,        0,        2,        4,
        0,     -0.1,       -1, NA_real_,
        0, NA_real_, NA_real_,        1,
  )

  my_results <- my_summarise(condition_and_results, condition_and_results)

  tmp <- c(0.1, 0, -0.1)

  expect_equal(my_results$bias, 0)
  expect_equal(my_results$sd_bias, sd(tmp))
  expect_equal(my_results$var, var(tmp))
  expect_equal(my_results$mse, mean(tmp^2))
  expect_equal(my_results$sd_mse, sd(tmp^2))
  expect_equal(my_results$mae, mean(abs(tmp)))
  expect_equal(my_results$sd_mae, sd(abs(tmp)))
  expect_equal(my_results$N_missing, 1)
  expect_equal(my_results$N, 4)
  expect_equal(my_results$coverage, 0.5)
  expect_equal(my_results$width, 2)
  expect_equal(my_results$sd_width, 0)
  expect_equal(my_results$N_missing_CI, 2)
})

test_that("missings are treated correctly for summarise test", {
  my_summarise <- summarise_test(alpha = c(0.05, 0.01))

  condition_and_results <- tibble::tribble(
    ~real,       ~p,
        0,    0.001,
        0,    0.04 ,
        0,    0.1  ,
        0, NA_real_,
  )

  my_results <- my_summarise(condition_and_results, condition_and_results)

  expect_equal(my_results$rejection_0.05, 2/3)
  expect_equal(my_results$rejection_0.01, 1/3)
  expect_equal(my_results$N_missing_0.05, 1)
  expect_equal(my_results$N_missing_0.01, 1)
  expect_equal(my_results$N, 4)
})
