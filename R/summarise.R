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

  # prepend name of the list element to the name attribute
  # then call make unique
  unique_tags <- purrr::imap_chr(summarise_functions, function(x, i){
    if(is.null(attr(x, "name"))){
      i
    } else {
      paste(i, attr(x, "name"), sep=".")
    }
  }) |>
    make.unique()

  # change the name attribute to the new unique names
  summarise_functions <- purrr::map2(summarise_functions, unique_tags, \(x,y){
    attr(x, "name") <- y
    x
  })

  function(condition, results, fixed_objects = NULL){
    # transpose list of results,
    # so each analyse function can be summarised in a map
    results_t <- results |>
      purrr::transpose() |>
      purrr::map(function(x){
        purrr::map(x, tibble::as_tibble) |>
          dplyr::bind_rows()
        })

    # call summarise functions on the results of the corresponding analysis methods
    aggregate_results <- purrr::imap_dfc(summarise_functions, function(f,i){
      if(i %in% names(results_t)){
        res <- f(condition, results_t[[i]], fixed_objects)
        names(res) <- paste(attr(f, "name"), names(res), sep=".")
        res
      } else {
        NULL
      }
    })

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
#' A function that can be used in Summarise that returns a data frame with the columns
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
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
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

#' Generic summarise function for tests
#'
#' @param alpha the significance level(s)
#' @param name = NULL name for the summarise function, appended to the name of the analysis method in the final results
#'
#' @return
#'
#' A function that can be used in Summarise that returns a data frame with the columns
#'
#' * rejection_X
#' * rejection_Y
#' * ...
#'
#' Where X, Y, ... are the alpha levels given in the argument
#'
#' @export
#'
#' @examples
#' condition <- merge(
#'   assumptions_delayed_effect(),
#'   design_fixed_followup(),
#'   by=NULL
#' ) |>
#'   tail(4) |>
#'   head(1)
#'
#' summarise_all <- create_summarise_function(
#'   logrank=summarise_test(alpha=c(0.5, 0.9, 0.95, 0.99))
#' )
#'
#' # runs simulations
#' sim_results <- runSimulation(
#'   design=condition,
#'   replications=100,
#'   generate=generate_delayed_effect,
#'   analyse=list(
#'     logrank=analyse_logrank
#'   ),
#'   summarise = summarise_all
#' )
#'
#' sim_results[, grepl("rejection", names(sim_results))]
summarise_test <- function(alpha, name=NULL){
  res <- function(condition, results, fixed_objects){
    results_tmp <- outer(results$p, 1-alpha, FUN=`<`) |>
      colMeans() |>
      as.list() |>
      as.data.frame() |>
      setNames(paste0("rejection_", alpha))

    results_tmp
  }

  attr(res, "name") <- name

  res
}
