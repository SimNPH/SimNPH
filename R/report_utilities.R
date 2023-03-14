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
tabulate_parameters <- function(data,parameters,.names="Parameter"){
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
tabulate_methods <- function(data){
  mutate(data,
         test = !{is.na(rejection)&is.na(`rejection_0.025`)},
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
#'
#' @param data tibble with simulation results in long format
#' @param methods methods for which operating characteristics should be plotted
#' @param parameter_x parameter to show on x-axis
#' @param parameter_y characteristic to show on y-axis
#' @param filters character vector defining filters to apply to dataset (see Details)
#' @param parameter_row parameter to define row facets
#' @param parameter_col parameter to defin column facets
#' @param parameters_steps parameters across which x-scale is looped
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
                             filters = c("median_survival_ctrl == 12",
                                         "recruitment == 18",
                                         "censoring_prop == 0"),
                             parameter_row = "effect_size_ph",
                             parameter_col = "n_pat_design",
                             parameters_steps = NULL,...){
  ## Limit to interesting scenarios, take from shiny input
  filters = rlang::parse_exprs(filters)
  columns = c("method",parameter_x,parameter_y,parameter_row,parameter_col,parameters_steps)
  if(any((columns %in% names(data))==FALSE)){
    stop(paste("Column missing:",columns[which(!(columns %in% names(data)))]))
  }
  filter(data,
         method %in% methods) |>
    filter(!!!filters) |>
    select(all_of(columns)) |>
    pivot_wider(names_from = "method",values_from = parameter_y) |> #View()
    looplot::nested_loop_plot(x=parameter_x,
                              grid_rows=parameter_row,
                              grid_cols=parameter_col,
                              steps=parameters_steps,
                              design_type = 'partial',
                              grid_scale="free_y",
                              spu_x_shift = .2, #spase between groups
                              # steps_y_base=0.5,
                              # steps_y_height=0.1,
                              # steps_y_shift = 0.5,
                              colors = scales::brewer_pal(palette = "Set1"),
                              steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                              post_processing = list(
                                add_custom_theme = list(
                                  axis.text.x = element_text(angle = -90,
                                                             vjust = 0.5,
                                                             size = 8)
                                )),...)
}

#' Plot parameter scenarios
#'
#' Nested loop plot of parameter scenarios
#'
#' @param data tibble with simulation results in long format
#' @param parameters which parameter values to show on y-axis
#' @param parameter_x parameter to show on x-axis
#' @param filters character vector defining filters to apply to dataset (see Details)
#' @param parameter_row parameter to define row facets
#' @param parameter_col parameter to defin column facets
#' @param parameters_steps parameters across which x-scale is looped
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_parameter_scenario <- function(data,
                                    parameters,
                                    parameter_x,
                                    filters = c("median_survival_ctrl == 12",
                                                "recruitment == 18",
                                                "censoring_prop == 0"),
                                    parameter_row = "effect_size_ph",
                                    parameter_col = "n_pat",
                                    parameters_steps = NULL,...){
  ## Limit to interesting scenarios, take from shiny input
  filters = rlang::parse_exprs(filters)
  columns = c(parameter_x,parameters,parameter_row,parameter_col,parameters_steps)
  if(any((columns %in% names(data))==FALSE)){
    stop(paste("Column missing:",columns[which(!(columns %in% names(data)))]))
  }
  filter(data,!!!filters) |>
    select(all_of(columns)) |>
    #    pivot_wider(names_from = "method",values_from = parameter_y) |> #View()
    looplot::nested_loop_plot(x=parameter_x,
                              grid_rows=parameter_row,
                              grid_cols=parameter_col,
                              steps=parameters_steps,
                              design_type = 'partial',
                              grid_scale="free_y",
                              spu_x_shift = .2, #space between groups
                              # steps_y_base=0.5,
                              # steps_y_height=0.1,
                              # steps_y_shift = 0.5,
                              colors = scales::brewer_pal(palette = "Set1"),
                              steps_values_annotate = TRUE, steps_annotation_size = 2.5,
                              post_processing = list(
                                add_custom_theme = list(
                                  axis.text.x = element_text(angle = -90,
                                                             vjust = 0.5,
                                                             size = 8)
                                )),...)
}


