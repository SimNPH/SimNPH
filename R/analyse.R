#' Wrappers around Analyse Functions
#'
#' @param list_of_functions the list of functions to be wrapped
#' @param error the error function in the tryCatch call
#'
#' @return a list of functions
#'
#' @describeIn wrap_all_in_trycatch Wrap all functions in a list in tryCatch calls
#'
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
#' funs1 <- list(\(){stop("test")}, \(){1})
#' funs2 <- wrap_all_in_trycatch(funs1)
#' \dontrun{
#' lapply(funs1, \(f){f()})
#' }
#' lapply(funs2, \(f){f()})
#'
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

#' @describeIn wrap_all_in_trycatch wrap all functions in `withr::with_preserve_seed`
#'
#' @export
#'
#' @details Analysis functions might use random numbers. If simulations should
#'   be replicated this can interfere with the RNG state of other analysis
#'   functions. To avoid this you can wrap all analysis function in a
#'   `withr::with_preserve_seed` call, so that the RNG state is reset after each
#'   analysis function is called. This way adding, removing or changing one
#'   analysis function has no effect on the other analysis functions, even if
#'   the analysis functions use random numbers.
#'
#' @examples
#' funs1 <- list(\(){rnorm(1)})
#' funs2 <- list(\(){runif(1)}, \(){rnorm(1)})
#' funs3 <- funs2 |> wrap_all_in_preserve_seed()
#' set.seed(1)
#' lapply(funs1, \(f){f()})
#' set.seed(1)
#' lapply(funs2, \(f){f()})
#' set.seed(1)
#' lapply(funs3, \(f){f()})
wrap_all_in_preserve_seed <- function(list_of_functions){
  lapply(list_of_functions,
  \(f){
    function(...){
      withr::with_preserve_seed(f(...))
    }
  })
}
