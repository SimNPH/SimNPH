#' Wrap all functions in a list in tryCatch calls
#'
#' @param list_of_functions the list of functions to be wrapped
#' @param error the error function in the tryCatch call
#'
#' @return a list of functions
#' @export
#'
#' @details SimDesign redraws data if one analysis function fails. This is not
#'   only highly inefficient for large studies, but failure of a method is
#'   informative and might be of interest. Moreover redrawing of data might
#'   introduce bias if the failure of the method is not independent of the
#'   parameter value, which would be a strong assumption.
#'
#'   To avoid redrawing data, we can catch all errors the analysis methods could
#'   throw and return `NA` instead.
#'
#'   This is handled well by the summarise functions generated with
#'   `create_summarise_function` other summarise functions might throw errors
#'   when trying to `rbind` a data.frame to a scalar `NA` value. In this case
#'   add another `error` argument. For example `\(e){NULL}` could work in some
#'   cases, in other cases you'll have to give a function that returns a
#'   data.frame with the same columns as the analyse functions and only NA
#'   values.
#'
#' @examples
#' Analyse <- list(
#' A = \(condition, dat, fixed_objects=NULL){
#'   if( runif(1) >= 0.5){
#'     counter_11(1)
#'     list(x=1, y=1)
#'   } else {
#'     counter_12(1)
#'     stop("asdf")
#'   }
#' },
#' B = \(condition, dat, fixed_objects=NULL){
#'   counter_13(1)
#'   list(z=1)
#' }
#' )
#'
#' Analyse2 <- wrap_all_in_trycatch(Analyse)
#' Analyse2 <- wrap_all_in_trycatch(Analyse, error=\(e){NULL})
wrap_all_in_trycatch <- function(list_of_functions, error=\(e){warning(e$message); NA}){
  lapply(
    list_of_functions,
    \(f){
      function(...){
        tryCatch(
          f(...),
          error = error
        )
      }
    }
  )
}
