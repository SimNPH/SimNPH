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
#'   analyse methods are prepended to the respective coumn names.
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

    common_names <- intersect(names(results_t), names(summarise_functions))

    aggregate_results <- purrr::map(
      common_names,
      function(i){
        res <- summarise_functions[[i]](condition, results_t[[i]], fixed_objects)
        names(res) <- paste0(i, ".", names(res))
        res
      }
    )

    aggregate_results <- do.call(cbind, aggregate_results) |>
      as.data.frame()

    aggregate_results
  }
}
