#' Create a Function for Descriptive Statistics of a Dataset
#'
#' @return an analyse function that returns a list with the elements
#' * `followup` follow up time
#' * `events` table of events vs. treatment
#' * `ice` if column ice is present, table of intercurrent events, events, treatment
#' * `subgroup` if column subgroup is present, table of subgroup, events, treatment
#'
#' @export
#'
#' @examples
#' condition <- merge(
#'     assumptions_delayed_effect(),
#'     design_fixed_followup(),
#'     by=NULL
#'   ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_describe()(condition, dat)
analyse_describe <- function(){

  tabulate_helper <- function(dat, var){
    tmp <- list(
      sum(dat[, var]),
      sum(dat[dat$trt==0, var]),
      sum(dat[dat$trt==1, var])
    )

    names(tmp) <- c(var, paste0(var, "_ctrl"), paste0(var, "_trt"))
    tmp
  }

  function(condition, dat, fixed_objects = NULL){
    result <- list(
      n_pat = nrow(dat),
      max_followup = max(dat$t[dat$evt])
    )

    result <- c(result, tabulate_helper(dat, "evt"))

    if(!is.null(attr(dat, "followup"))){
      result$study_time <- attr(dat, "followup")
    } else {
      result$study_time <- NA_real_
    }

    if(hasName(dat, "ice")){
      result <- c(result, tabulate_helper(dat, "ice"))
    }

    if(hasName(dat, "subgroup")){
      result <- c(result, tabulate_helper(dat, "subgroup"))
    }

    result
  }
}


#' @describeIn analyse_describe Summarise Descriptive Statistics
#' @param name name for the summarise function, appended to the name of the analysis method in the final results
#'
#' @return
#'
#' A function that can be used in Summarise that returns a data frame with
#' columns with means and standard deviations for every variable in the
#' description.
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
#'   describe=summarise_describe()
#' )
#'
#' # runs simulations
#' sim_results <- runSimulation(
#'   design=condition,
#'   replications=100,
#'   generate=generate_delayed_effect,
#'   analyse=list(
#'     describe=analyse_describe()
#'   ),
#'   summarise = summarise_all
#' )
#'
#' # study time is missing, since there was no admin. censoring
#' sim_results[, 9:16]
summarise_describe <- function(name=NULL){
  res <- function(condition, results, fixed_objects=NULL){
    means <- results |>
      apply(2, mean, na.rm=TRUE) |>
      t() |>
      as.data.frame() |>
      setNames(names(results))

    sds <- results |>
      apply(2, sd, na.rm=TRUE) |>
      t() |>
      as.data.frame() |>
      setNames(paste0("sd_", names(results)))

    cbind(means, sds)
  }

  attr(res, "name") <- name

  res
}




