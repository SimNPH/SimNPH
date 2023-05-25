#' @describeIn merge_additional_results Update or add Rows and Columns
#'
#' @param x left data.frame
#' @param y right data.frame
#' @param by columns to match by
#'
#' @return a data.frame
#' @export
#'
#' @details updates columns in x with values from matched rows in y and add
#'   joins columns from y not present in x. Calls `rows_upsert` and then
#'   `full_join`.
#'
#' @examples
#' a <- data.frame(x=5:2, y=5:2, a=5:2)
#' b <- data.frame(x=1:4, y=1:4+10, b=1:4*10)
#' upsert_merge(a, b, by="x")
upsert_merge <- function(x, y, by){
  y_upsert <- y |>
    dplyr::select(any_of(colnames(x)))

  y_merge <- y |>
    dplyr::select(any_of(by), any_of(setdiff(colnames(y), colnames(x))))

  res <- x |>
    dplyr::rows_upsert(y_upsert, by=by) |>
    dplyr::full_join(y_merge, by=by)

  res
}

#' Merge results from additional or updated simulations
#'
#' @param old old results
#' @param new new/additional results
#' @param design_names names of the paramterst
#' @param descriptive_regex regular expression for columns of descriptive statistics
#'
#' @return a data.frame of the merged simulation results
#' @export
#'
#' @details if `design_names` is omitted its value is taken from the
#'   `design_names` attribute of the simulation results.
#'
#'   If `descriptive_regex` is given, columns matching the regular expression in
#'   both datasets are compared, a warning is given, if the values of those
#'   columns do not match. This is intended to compare descriptive statistics or
#'   results of unchanged analysis methods to ensure, that both results stem
#'   from an exact replication of the simulation results.
#'
#' @examples
#' \dontrun{
#' old <- readRDS("old_results.Rds")
#' new <- readRDS("new_results.Rds")
#' combined <- merge_additional_results(old, new)
#' }
merge_additional_results <- function(old, new, design_names=NULL, descriptive_regex=NULL){
  if( is.null(design_names) ){
    # check if the datasets have the same design columns
    if( !setequal(attr(new, "design_names")$design, attr(old, "design_names")$design) ){
      warning(gettext("Different design columns in old and new simulations."))
    }
    design_names <- intersect(attr(old, "design_names")$design, attr(new, "design_names")$design)
  }

  if(!is.null(descriptive_regex)){
    compare_old <- old |>
      dplyr::semi_join(new, by=design_names) |>
      dplyr::select(all_of(design_names), matches(descriptive_regex)) |>
      dplyr::arrange(!!!design_names)

    compare_new <- new |>
      dplyr::semi_join(compare_old, by=design_names) |>
      dplyr::select(all_of(design_names), matches(descriptive_regex)) |>
      dplyr::arrange(!!!design_names)

    common_names <- intersect(names(compare_old), names(compare_new)) |>
      setdiff(design_names)

    for(i in common_names){
      if( any( dplyr::pull(compare_old, i) != dplyr::pull(compare_new, i)) ){
        warning(gettext('Different values in descriptive statistics, column "', i, '" no exact replication.'))
      }
    }
  }


  combined <- upsert_merge(old, new, by=design_names)

  # update design names attribute
  attr(combined, "design_names") <- mapply(
    union,
    attr(old, "design_names"),
    attr(new, "design_names"),
    SIMPLIFY = FALSE
    )
  attributes(combined)[c("ERROR_msg", "WARNING_msg", "extra_info")] <- NULL

  combined
}

#' Rename Columns in Simulation Results and Update Attributes
#'
#' @param results `SimDesign` object
#' @param rename named vector of new names
#'
#' @return `SimDesign` object with updated column names
#' @export
#'
#' @describeIn rename_results_column Rename Columns in Simulation Results
#'
#' @examples
#' \dontrun{
#' results2 <- rename_results_column(results, c("old_name"="new_name"))
#' }
rename_results_column <- function(results, rename){
  rename_helper <- function(x, rename){
    x[x %in% names(rename)] <- rename[x[x %in% names(rename)]]
    x
  }

  names(results) <- rename_helper(names(results), rename)
  attr(results, "design_names") <- lapply(
    attr(results, "design_names"),
    rename_helper,
    rename=rename
  )

  results
}

#' @param pattern regexp pattern as understood by `stringr::str_replace_all`
#' @param replacement replacement as understood by `stringr::str_replace_all`
#'
#' @export
#'
#' @describeIn rename_results_column Rename Columns in Simulation Results by Pattern
#'
#' @examples
#' \dontrun{
#' results2 <- rename_results_column_pattern(results, "old_name", "new_name")
#' }
rename_results_column_pattern <- function(results, pattern, replacement){
  names(results) <- stringr::str_replace_all(names(results), pattern, replacement)

  attr(results, "design_names") <- lapply(
    attr(results, "design_names"),
    stringr::str_replace_all,
    pattern=pattern,
    replacement=replacement
  )

  results
}
