#' functions for plotting results
#'
#' @describeIn results_pivot_longer pivot simulation results into long format
#'
#' @param data for results_pivot_longer: simulation result as retured by SimDesign
#' @param exclude_from_methods=c("descriptive") "methods" that should not be pivoted into long format
#'
#' @return dataset in long format with one row per method and scenario and one
#'   column per metric
#'
#' @details With `exclude_from_methods` descriptive statistics or results of
#'   reference methods can be kept as own columns and used like the columns of
#'   the simulation parameters.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   plot_data <- simulation_results |>
#'   results_pivot_longer
#' }
results_pivot_longer <- function(data, exclude_from_methods=c("descriptive")){
  # delete potentially huge attributes that are not needed for plots
  attr(data, "ERROR_msg")    <- NULL
  attr(data, "WARNING_msg")  <- NULL
  attr(data, "extra_info")   <- NULL

  methods <- attr(data, "design_names") |>
    getElement("sim") |>
    str_extract(".*(?=\\.[^\\d])")

  summaries <- attr(data, "design_names") |>
    getElement("sim") |>
    str_remove(str_c(methods, "."))

  include <- !(methods %in% exclude_from_methods)

  pivot_spec <- tibble(
    .name=attr(data, "design_names")$sim[include],
    .value=summaries[include],
    method=methods[include]
    )

  result <- data |>
    rename(
      n_pat_design = n_pat
    ) |>
    pivot_longer_spec(pivot_spec)
}

order_combine_xvars <- function(data, xvars){
  data |>
    arrange(across(all_of(xvars))) |>
    unite(x, !!!xvars) |>
    mutate(
      x = fct_inorder(x)
    )
}

#' @param data for plot functions: simulation results in long format
#' @param xvars variables to be displayed on the x axis
#' @param yvar variables to be displayed on the y axis
#' @param yvar_sd=NULL variables that contain the standard deviations of the y-variables
#'
#' @describeIn results_pivot_longer nested loop plot of simulation results
#'
#' @return a ggplot object
#' @export
#'
#' @details A nested loop plot is a plot in which all variables on the x-axis
#'   are displayed together by sorting by the first, the second, ... variable
#'   and then creating a factor containing all variable levels.In this way a
#'   metric can be plotted for all simulation parameters in one plot.
#'
#'   `xvars` is typically a vector of all simulation parameters. `yvars` is
#'   typically one metric, for example bias, mse, ...
#'
#' @examples
#' \dontrun{
#' plotdata |>
#'  filter(
#'    method %in% c("ahr", "gahr", "cox", "weighted_cox")
#'  ) |>
#'  nested_loop_plot(
#'    xvars=c("recruitment", "n_pat_design", "interim_events", "final_events", "n_ctrl",
#'            "n_trt", "hazard_ctrl", "prog_prop_trt", "prog_prop_ctrl", "hr_before_after",
#'            "censoring_prop", "effect_size_ph"),
#'    yvar="bias",
#'    yvar_sd="sd_bias"
#'  ) +
#'  geom_hline(yintercept=0)
#' }
nested_loop_plot <- function(data, xvars, yvar, yvar_sd=NULL){
  yvar <- sym(yvar)

  data <- data |>
    order_combine_xvars(xvars)


  if(is.null(yvar_sd)){
    data <- data |>
      mutate(
        lower=NA_real_,
        upper=NA_real_
      )
  } else {
    yvar_sd <- sym(yvar_sd)
    data <- data |>
      mutate(
        sd = !!yvar_sd,
        lower = !!yvar - 2*sd/sqrt(REPLICATIONS),
        upper = !!yvar + 2*sd/sqrt(REPLICATIONS)
      )
  }

  ggplot(data, aes(x=x, y=!!yvar, group=method, colour=method, ymin=lower, ymax=upper, fill=method)) +
    geom_line() +
    geom_ribbon(alpha=0.3)
}


#' @describeIn results_pivot_longer nested loop plot of simulation results
#'
#' @param trellis_var variable for faccets
#' @param trellis_var_y=NULL optional additional variable for faccets
#'
#' @details A trellis plot is a series of plot in which one simulation parameter
#'   is plotted against one metric. The plot is split into the levels of one or
#'   two additional simulation parameters. This function also allows more then
#'   one x variable to give a hybrid trellis nested loop plot. (TODO: add
#'   reference)
#'
#'   If `trellis_var_y` is omitted, a plot using `facet_wrap` facceted by
#'   `trellis_var` is created. If `trellis_var_y` is given, then the facets use
#'   `trellis_var` for the columns and `trellis_var_y` for the rows of the plot
#'   matrix.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plotdata |>
#'  filter(
#'    method %in% c("ahr", "gahr", "cox", "weighted_cox")
#'  ) |>
#'  trellis_plot(
#'    xvars="effect_size_ph",
#'    yvar="bias",
#'    yvar_sd="sd_bias",
#'    trellis_var="prog_prop_trt"
#'  )
#'
#'plotdata |>
#'  filter(
#'    method %in% c("ahr", "gahr", "cox", "weighted_cox")
#'  ) |>
#'  trellis_plot(
#'    xvars="effect_size_ph",
#'    yvar="bias",
#'    yvar_sd="sd_bias",
#'    trellis_var="prog_prop_trt",
#'    trellis_var_y="prog_prop_ctrl"
#'  )
#' }
trellis_plot <- function(data, xvars, yvar, trellis_var, yvar_sd=NULL, trellis_var_y=NULL){
  yvar <- sym(yvar)

  data <- data |>
    order_combine_xvars(xvars)

  if(is.null(yvar_sd)){
    data$sd <- NA_real_
  } else {
    yvar_sd <- sym(yvar_sd)
    data <- data |>
      mutate(
        sd = !!yvar_sd
      )
  }

  data <- data |> mutate(
    lower = !!yvar - 2*sd/sqrt(REPLICATIONS),
    upper = !!yvar + 2*sd/sqrt(REPLICATIONS)
  )

  gg <- ggplot(data, aes(x=x, y=!!yvar, group=method, colour=method, ymin=lower, ymax=upper, fill=method)) +
    geom_line() +
    geom_ribbon(alpha=0.3)

  if(is.null(trellis_var_y)){
    gg <- gg +
      facet_wrap(
        trellis_var,
        labeller = label_both
      )
  } else {
    trellis_var <- sym(trellis_var)
    trellis_var_y <- sym(trellis_var_y)

    gg <- gg +
      facet_grid(
        rows=vars(!!trellis_var_y),
        cols=vars(!!trellis_var),
        labeller = label_both
      )
  }

  gg
}


