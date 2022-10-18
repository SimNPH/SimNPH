#' Create Analyse Functions for Group Sequential Design
#'
#' @param followup followup events or time
#' @param followup_type "events" or "time"
#' @param alpha nominal alpha at each stage
#' @param analyse_functions analyse function or list of analyse functions
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
#' @details
#' `followup`, `followup_type` and `alpha` are evaluated for every simulated
#' dataset, i.e. the arguments to the Analyse function are available,
#' expressions like `followup=c(condition$interim, condition$max_followup)` are
#' valid arguments.
#'
#' analyse_functions should take arguments condition, dataset and fixed_objects
#' and return a list conatining p-value, number of patients and number of event
#' in the columsn `p`, `N_pat` and `N_evt`.
#'
#' @examples
#' # create a function to analyse after interim_events and maximum followup time
#' # given in the condition row of the design data.frame with given
#' # nominal alpha
#' analyse_maxcombo_sequential <- analyse_group_sequential(
#'   followup = c(condition$interim_events, condition$followup),
#'   followup_type = c("event", "time"),
#'   alpha = c(0.025, 0.05),
#'   analyse_functions = analyse_maxcombo
#' )
analyse_group_sequential <- function(followup, followup_type, alpha, analyse_functions){

  if(length(analyse_functions) == 1 && is(analyse_functions, "function")){
    analyse_functions <- list(analyse_functions)
  }

  followup <- substitute(followup)
  followup_type <- substitute(followup_type)
  alpha <- substitute(alpha)

  function(condition, dat, fixed_objects = NULL){

    followup <- eval(followup)
    followup_type <- eval(followup_type)
    alpha <- eval(alpha)

    datasets_stages <- mapply(
      function(followup, followup_type){
        if(followup_type == "event"){
          admin_censoring_events(dat, followup)
        } else if(followup_type == "time"){
          admin_censoring_time(dat, followup)
        }
      },
      followup, followup_type,
      SIMPLIFY = FALSE
    )

    results_stages <- mapply(
      function(dataset, analyse_function){
        analyse_function(condition, dataset, fixed_objects=fixed_objects)
      },
      datasets_stages, analyse_functions,
      SIMPLIFY = FALSE
    )

    p_values_stages <- sapply(results_stages, function(x){
      x[["p"]]
    })

    rejected_at_stage <- suppressWarnings(min(which(p_values_stages < alpha)))
    max_stage <- min(rejected_at_stage, length(followup))

    result_tmp <- list(
      rejected_at_stage = rejected_at_stage,
      N_pat = results_stages[[max_stage]][["N_pat"]],
      N_evt = results_stages[[max_stage]][["N_evt"]],
      followup = attr(datasets_stages[[max_stage]], "followup"),
      results_stages = list(list(results_stages))
    )

    result_tmp
  }
}


#' Summarise Output from Analyse Functions for Group Sequential Design
#'
#' @param condition condition row passed from runSimulation
#' @param results list of results from analyse functions
#' @param fixed_objects fixed_objects passed from runSimulation
#'
#' @describeIn analyse_group_sequential Summarise Output from Analyse Functions for Group Sequential Design
#'
#' @return
#' Returns a `data.frame` with the columns
#'  * `rejection` the empirical rejection rate
#'  * `n_pat` the mean number of patients recruited
#'  * `n_evt` the mean number of observed events
#'  * `followup` the mean followup time
#'
#' @export
#'
#' @examples
#' Summarise <- create_summarise_function(
#'   maxcombo_seq = summarise_group_sequential,
#'   logrank_seq = summarise_group_sequential
#' )
summarise_group_sequential <- function(condition, results, fixed_objects=NULL){
  data.frame(
    "rejection" = mean(is.finite(results$rejected_at_stage)),
    "n_pat" = mean(results$N_pat),
    "n_evt" = mean(results$N_evt),
    "followup" = mean(results$followup)
  )
}
