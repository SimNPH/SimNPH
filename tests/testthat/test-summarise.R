test_that("creating summarise function works", {

  Summarise <- create_summarise_function(list(
    method1 = function(condition, results, fixed_objects){
      mean(results$value1)
    },
    method2 = function(condition, results, fixed_objects){
      mean(results$value2)
    }
  ))


  condition <- numeric(0)

  results <- list(
    list(method1=list(value1=1), method2=list(value2=1)),
    list(method1=list(value1=3), method2=list(value2=3))
  )

  summary <- Summarise(condition, results)

  expect_type(Summarise, "closure")
  expect_equal(summary[1, ], c(2,2))
})
