#' Create a summarise function from a named list of functions
#'
#' @param summarise_functions named list of functions
#'
#' @return a function with arguments condition, results, fixed objects
#' @export
#'
#' @details the names of the list of functions correspond to the names in the
#'   list of analyse functions, each summarise function is applied to the
#'   results of the analyse function of the same name, names not present in both
#'   lists are ommitted in either list.
#'
#'   The functions in the list should have the arguments `condition`, `results`
#'   and `fixed_objects`. `results` is a list of lists. The outer list has one
#'   element for each replication, the inner list has one entry for each Analyse
#'   function. (Analyse functions have to return lists for this to work,
#'   otherwise the results are simplified to data.frames. Analyse functions from
#'   the SimNPH package all return lists.)
#'
#'   The individual summarise functions have to return data.frames, which are
#'   concatendated column-wise to give one row per condition. The names of the
#'   analyse methods are prepended to the respective coumn names, if the
#'   functions have a "name" attribute this is appended to the column names of
#'   the output. Column names not unique after that are appended numbers by
#'   `make.unique`.
#'
#' @examples
#' Summarise <- create_summarise_function(
#'   maxcombo = function(condition, results, fixed_objects=NULL){
#'     data.frame("rejection"=mean(results$p < alpha))
#'   },
#'   logrank  = function(condition, results, fixed_objects=NULL){
#'     data.frame("rejection"=mean(results$p < alpha))
#'   }
#' )
create_summarise_function <- function(...){
  summarise_functions <- list(...)

  function(condition, results, fixed_objects = NULL){
    # transpose lists, so each analyse function can be summarised in a map
    results_t <- results |>
      purrr::transpose() |>
      purrr::map(function(x){
        purrr::map(x, tibble::as_tibble) |>
          dplyr::bind_rows()
        })

    names_funs <- purrr::map_chr(summarise_functions, \(x){
      name <- attr(x, "name")
      if(is.null(name)){
        ""
      } else {
        paste0(".", name)
      }
    })

    aggregate_results <- purrr::imap(summarise_functions, function(f,i){
      if(i %in% names(results_t)){
        f(condition, results_t[[i]], fixed_objects)
      } else {
        NA_real_
      }
    })

    names(aggregate_results) <- paste0(names(aggregate_results), names_funs)
    names(aggregate_results) <- make.unique(names(aggregate_results))

    aggregate_results <- do.call(cbind, aggregate_results) |>
      as.data.frame()

    aggregate_results
  }
}


#' Generic Summarise function for esitmators
#'
#' @param est estimator, expression evaluated in results
#' @param real real summary statistic, expression evaluated in condition
#' @param lower lower CI, expression evaluated in results
#' @param upper upper CI, expression evaluated in results
#' @param name = NULL name for the summarise function, appended to the name of the analysis method in the final results
#'
#' @return
#'
#' A data frame with the columns
#'
#' * bias the bias with respect to the true value
#' * var the variance
#' * mse the mean squared error with respect to the true value
#' * mae the mean absolute error with respect to the true value
#' * coverage the coverage of the confidence interval
#' * width the mean width of the confidence interval
#'
#' @export
#'
#' @examples
#' # generate the design matrix and append the true summary statistics
#' condition <- desing_skeleton_delayed_effect() |>
#'   tail(4) |>
#'   head(1) |>
#'   true_summary_statistics_delayed_effect(cutoff_stats = 15)
#'
#' # create some summarise functions
#' summarise_all <- create_summarise_function(
#'   coxph=summarise_estimator(hr, gAHR, hr_lower, hr_upper, name="aAHR"),
#'   coxph=summarise_estimator(hr, hazard_trt/hazard_ctrl, hr_lower, hr_upper, name="HR"),
#'   coxph=summarise_estimator(exp(coef), gAHR),
#'   coxph=summarise_estimator(hr, NA_real_)
#' )
#'
#' # runs simulations
#' sim_results <- runSimulation(
#'   design=condition,
#'   replications=10,
#'   generate=generate_delayed_effect,
#'   analyse=list(
#'     coxph=analyse_coxph
#'   ),
#'   summarise = summarise_all
#' )
#'
#' # mse is missing for the summarise function in which the real value was NA
#' sim_results[, names(sim_results) |> grepl(pattern="mse")]
#' # but the variance can be estimated in all cases
#' sim_results[, names(sim_results) |> grepl(pattern="var")]
summarise_estimator <- function(est, real, lower=NULL, upper=NULL, name=NULL){

  est <- substitute(est)
  real <- substitute(real)
  lower <- substitute(lower)
  upper <- substitute(upper)

  res <- function(condition, results, fixed_objects){
    est    <- eval(est,   envir = results)
    real   <- eval(real,  envir = condition)
    lower  <- eval(lower, envir = results)
    upper  <- eval(upper, envir = results)

    results_tmp <- data.frame(
      bias     = mean(est - real),
      var      = var(est),
      mse      = mean((est-real)^2),
      mae      = mean(abs(est-real)),
      coverage = NA_real_,
      width    = NA_real_
    )

    if(!is.null(lower) && !is.null(upper)){
      results_tmp$coverage <- mean( (lower <= real) & (upper >= real) )
      results_tmp$width    <- mean( abs(upper - lower) )
    }

    results_tmp
  }

  attr(res, "name") <- name

  res
}
