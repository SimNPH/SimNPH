#' Tabulate Parameters
#'
#' Return a table with selected parameters and corresponding unique values
#'
#' @param data tibble with simulation results (or parameter design) in long format
#' @param parameters which parameters to select for tabulation
#' @param .names column title for variables in table
#'
#' @return tibble containing selected parameters and corresponding unique values
#' @export
#'
#' @examples
#' \dontrun{
#' tabulate_parameters(sim_results_long, c("bias", "mse"))
#' }
tabulate_parameters <- function(data, parameters, .names="Parameter"){
  summarise(data,across(all_of(parameters),~paste(unique(.x),collapse=', '))) |>
    pivot_longer(everything(),names_to=.names,values_to = "Unique Values")
}


#' Tabulate Methods
#'
#' Return a table of all methods investigated in simulation results that
#' categorizes whether methods provide tests, estimators, and/or CIs
#'
#' @param data tibble with simulation results in long format
#'
#' @return tibble containing methods and what they provide
#' @export
#'
#' @examples
#' \dontrun{
#' tabulate_methods(delayed_long)
#' }
tabulate_methods <- function(data){
  mutate(data,
         test = !{is.na(rejection)&is.na(`rejection_0.025`)&is.na(null_lower)},
         estimator = !{is.na(bias)},
         ci = !is.na(coverage)) |>
    group_by(method) |>
    summarise(across(all_of(c('test','estimator','ci')),~.x[1])) |>
    mutate(across(2:4,~ifelse(.x,'X','')))
}

#' Nested Loop Plot of Simulation Results
#'
#' Generates a nested loop plot of simulation results
#'
#' Parameters selected for x-axis, rows, cols, and steps need to be such that together with filters all remaining scenarios uniquely identify a point on the y-axis.
#' rel_logrank will plot the results relative to the logrank test, this requires that the corresponding y_parameter is available for the logrank test and practically is only sensible for plots of power
#'
#' @param data tibble with simulation results in long format
#' @param methods methods for which operating characteristics should be plotted
#' @param parameter_x parameter to show on x-axis
#' @param parameter_y characteristic to show on y-axis
#' @param filters character vector defining filters to apply to dataset (see Details)
#' @param parameter_row parameter to define row facets
#' @param parameter_col parameter to defin column facets
#' @param parameters_steps parameters across which x-scale is looped
#' @param rel_logrank should results be ploted relative to the logrank test (see Details)
#' @param ... additional options to looplot
#'
#' @return
#' @export
#'
#' @examples
plot_sim_results <- function(data,
                             methods,
                             parameter_x,
                             parameter_y,
                             filters = c("abs(hazard_ctrl - nph::m2r(12))<1e-5",
                                         "recruitment == 18",
                                         "censoring_prop == 0"),
                             parameter_row = "effect_size_ph",
                             parameter_col = "n_pat_design",
                             parameters_steps = NULL,
                             rel_logrank = FALSE,
                             ...){
  ## Limit to interesting scenarios, take from shiny input
  filters = rlang::parse_exprs(filters)
  columns = c("method",parameter_x,parameter_y,parameter_row,parameter_col,parameters_steps)
  if(any((columns %in% names(data))==FALSE)){
    stop(paste("Column missing:",columns[which(!(columns %in% names(data)))]))
  }
  if(rel_logrank){
    if(!'logrank' %in% methods){
      stop("Error: cannot plot results relative to logrank in absence of logrank results")
    }
  }

  #  ggplot2::scale_colour_discrete()$palette(n=length(unique(metadata$category))) |>
  pal <- rep(ggsci::pal_jco()(10),2)
  my_colours <- pal[1:length(unique(metadata$category))] |>
    setNames(unique(metadata$category))
  my_colours <- my_colours[metadata$category] |> setNames(metadata$method)

  my_shapes <- metadata |> group_by(category) |>
    mutate(symbol = 1:n()) |> pull(symbol) |>
    setNames(metadata$method)


  plot_data <- filter(data,
                      method %in% methods) |>
    filter(!!!filters) |>
    select(all_of(columns))

    plot_rows <- plot_data[[parameter_row]] |> unique() |> length()

  if(rel_logrank){
    plot_data <- plot_data |>
      pivot_wider(names_from = method,values_from = all_of(parameter_y)) |>
      mutate(across(all_of(methods),~.x-logrank)) |>
      select(-logrank) |>
      pivot_longer(any_of(methods),values_to = parameter_y,names_to = "method")
  }

  gg <- plot_data |>
    combined_plot(yvar = parameter_y,
                  methods = methods,
                  facet_x_vars = parameter_col,
                  facet_y_vars = parameter_row,
                  xvars = c(parameters_steps,parameter_x),
                  use_colours = my_colours,
                  use_shapes = my_shapes,
                  grid_level = 2,
                  scale_stairs = .5,
                  heights_plots = c(4,1),
                  ...)
   gg[[1]] <- gg[[1]] + scale_shape_manual(values = c(as.character(c(1:9)),letters))
   gg[[1]] <- gg[[1]] |> labs_from_labels()
   gg
}


#' Plot Survival Curve
#'
#' Plot a survival curve for a given set of parameters.
#'
#' @param data A data frame containing the parameters needed to generate the survival curve. The required columns depend on the type of survival curve to generate.
#' @param type A character string indicating the type of survival curve to generate. Valid values are 'delayed', 'progression', 'subgroup', and 'crossing'.
#'
#' @return A survival curve plot.
#'
#' @importFrom nph pchaz subpop_pchaz pop_pchaz plot_shhr
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   hazard_trt = 0.01,
#'   hazard_ctrl = 0.008,
#'   delay = 30,
#'   hazard_after_prog = 0.005,
#'   prog_rate_trt = 0.2,
#'   prog_rate_ctrl = 0.15,
#'   prevalence = 0.5
#' )
#' plot_nph_curves(data = data, type = 'delay')
#' }
#' @export
plot_nph_curves <- function(data, type) {
  if (!(type %in% c("delay", "progression", "subgroup", "crossing"))) {
    stop("Error: 'type' argument must be one of 'delay', 'progression', 'subgroup', or 'crossing'")
  }

  # Determine t_max
  t_max <- if (hasName(data, "descriptive.max_followup")) {
    ## data$descriptive.max_followup # - this leads to different x-axis ranges
    log(100) / data$hazard_ctrl
  } else if (hasName(data, "hazard_ctrl")) {
    log(100) / data$hazard_ctrl
  } else {
    1095.75
  }

  # Create plot
  if (type == "delay") {
    if (data$delay == 0) {
      pch_trt <- nph::pchaz(Tint = c(0, t_max), lambda = c(data$hazard_trt))
    } else {
      pch_trt <- nph::pchaz(Tint = c(0, data$delay, t_max), lambda = c(data$hazard_ctrl, data$hazard_trt))
    }

    pch_ctrl <- nph::pchaz(Tint = c(0, t_max), lambda = c(data$hazard_ctrl))

  } else if (type == "progression") {
    pch_trt <- nph::subpop_pchaz(
      Tint=c(0, t_max),
      lambda1 = data$hazard_trt,
      lambda2 = data$hazard_after_prog,
      lambdaProg = data$prog_rate_trt,
      discrete_approximation = TRUE
    )

    pch_ctrl <- nph::subpop_pchaz(
      Tint=c(0, t_max),
      lambda1 = data$hazard_ctrl,
      lambda2 = data$hazard_after_prog,
      lambdaProg = data$prog_rate_ctrl,
      discrete_approximation = TRUE
    )

  } else if (type == "subgroup") {
    pch_trt <- nph::pop_pchaz(
      Tint = c(0, t_max),
      lambdaMat1 = matrix(c(data$hazard_trt, data$hazard_subgroup), 2, 1),
      lambdaMat2 = matrix(0, 2, 1),
      lambdaProgMat = matrix(0, 2, 1),
      p = c((1 - data$prevalence), data$prevalence)
    )

    pch_ctrl <- nph::pchaz(
      Tint = c(0, t_max),
      lambda = data$hazard_ctrl
    )

  } else if (type == "crossing") {
    if (data$crossing == 0) {
      pch_trt <- nph::pchaz(Tint=c(0, t_max), lambda = c(data$hazard_trt_after))
    } else {
      pch_trt <- nph::pchaz(Tint=c(0, data$crossing, t_max), lambda = c(data$hazard_trt_before, data$hazard_trt_after))
    }

    pch_ctrl <- nph::pchaz(Tint=c(0, t_max), lambda = c(data$hazard_ctrl))
  }

  shhr_gg(pch_ctrl,pch_trt)
#  legend(y=1,x=0, legend = c("control","treatment"), lty = 1, col = c("black", "green"), cex = 1.2, xjust=.3,yjust = -2,lwd=4)
}

#' Interactive Filters
#'
#' Utility function to define interactive (or not filters)
#'
#' If input is set to NULL a character vector with static default filters will be returned,
#' this is useful for debugging without the need to rebuild the shiny document. Also input is not
#' available to the function unless it is defined within the markdown document (even if only sourced)
#' so you need to pass it from within the report.
#'
#' @param prefix prefix to identify shiny input panel
#' @param input pass the shiny input object
#'
#' @return
#' @export
#'
#' @examples
interactive_filters <- function(prefix,input=NULL) {
  if(!is.null(input)) {
    c(str_c("abs(hazard_ctrl - nph::m2r(", input[[str_replace("prefix_median_survival_ctrl", "prefix_", prefix)]],"))<1e-6"),
      str_c("recruitment == ", input[[str_replace("prefix_recruitment", "prefix_", prefix)]]),
      str_c("censoring_prop == ", input[[str_replace("prefix_censoring", "prefix_", prefix)]]))
  } else {
    c("abs(hazard_ctrl - nph::m2r(12))<1e-6",
      "recruitment == 18",
      "censoring_prop == 0.1")
  }
}

#' Generate a sidebar panel with select inputs for clinical trial parameters.
#'
#' @param prefix A character string to use as a prefix for the input IDs.
#' @param default_mts The default value for the median survival time in the control arm.
#' @param default_recruitment The default value for the time to recruit target number of subjects.
#' @param default_censoring The default value for the proportion of subjects with random censoring.
#'
#' @return A sidebar panel with select inputs for clinical trial parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' # generate the sidebar panel with default values
#' sidebar_panel <- input_trial_inputs(
#'   "trial1_",
#'   default_mts = 12,
#'   default_recruitment = 18,
#'   default_censoring = 0
#' )
#' }
input_trial_inputs <- function(prefix, default_mts = 12, default_recruitment = 18, default_censoring = 0) {
  panel <- sidebarPanel(
    selectInput(paste0(prefix,"median_survival_ctrl"), label = "Median Survival in Control arm (mts):",
                choices = c(36, 12, 6), selected = default_mts),
    selectInput(paste0(prefix,"recruitment"), label = "Time to recruit target number of subjects (mts):",
                choices = c(18, 30), selected = default_recruitment),
    selectInput(paste0(prefix,"censoring"), label = "Proportion of subjects with random censoring:",
                choices = c(0, 0.1, 0.3), selected = default_censoring)
  )
  return(panel)
}

