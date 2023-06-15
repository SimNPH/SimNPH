#' Functions for Plotting and Reporting Results
#'
#' @describeIn results_pivot_longer pivot simulation results into long format
#'
#' @param data for results_pivot_longer: simulation result as retured by SimDesign
#' @param exclude_from_methods "methods" that should not be pivoted into long format
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
    stringr::str_extract(".*(?=\\.[^\\d])")

  summaries <- attr(data, "design_names") |>
    getElement("sim") |>
    stringr::str_remove(str_c(methods, "."))

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
    tidyr::pivot_longer_spec(pivot_spec)
}

order_combine_xvars <- function(data, xvars, facet_vars=c(), height_x_axis=0.8, grid_level=2){

  result <- data |>
    arrange(!!!xvars) |>
    unite(x, !!!xvars, remove=FALSE) |>
    mutate(
      x = fct_inorder(x)
    )

  x_axis <- result |>
    select(x, !!!xvars, !!!facet_vars) |>
    unique() |>
    pivot_longer(cols=c(-x, -any_of(facet_vars))) |>
    group_by(name) |>
    mutate(
      y=(as.integer(as.factor(value))-1) / (length(unique(value))-1)
    ) |>
    ungroup() |>
    mutate(
      level = match(name, xvars),
      y = level - (y*height_x_axis) - (0.5 * (1-height_x_axis))
    )

  x_axis_labels <- x_axis |>
    group_by(name) |>
    summarise(
      level=first(level),
      label=str_c(first(name), ": ", str_c(unique(value), collapse=", "))
    )

  x_axis_breaks <- result |>
    select(!!!xvars[1:grid_level], x) |>
    group_by(!!!xvars[1:grid_level]) |>
    filter(1:n() == 1) |>
    pull(x)

  attr(result, "x_axis") <- x_axis
  attr(result, "x_labels") <- x_axis_labels
  attr(result, "x_axis_breaks") <- x_axis_breaks
  result
}

#' @describeIn results_pivot_longer Nested Loop Plot with optional Facets
#'
#' @param data for combined_plto simulation results in long format, as returned by `results_pivot_longer`.
#' @param methods methods to include in the plot
#' @param xvars orderd vector of variable names to display on the x axis
#' @param yvar variable name of the variable to be displayed on the y axis (metric)
#' @param facet_x_vars vector of variable names to create columns of facets
#' @param facet_y_vars vector of variable names to create rows of facets
#' @param split_var index of xvars along groups of which the plot should be split
#' @param heights_plots relative heights of the main plot and the stairs on the bottom
#' @param scale_stairs height of the stairs for each variable between 0 and 1
#' @param grid_level depth of loops for which the grid-lines are drawn
#' @param scales passed on to facet_grid
#' @param hlines position of horizontal lines, passed as `yintercept` to
#'   `geom_hline`
#' @param use_colours optional named vector of colours used in `scale_colour_manual`
#' @param use_shapes optional named vector of shapes used in `scale_shape_manual`
#'
#' @return a ggplot/patchwork object conatining the plots
#' @export
#'
#' @details `use_colours` and `use_shapes` both use the `method` variable in their respective aesthetics.
#'
#' @examples
#' \dontrun{
#' # plot the rejection rate of two methods
#' combined_plot(
#'   results_long,
#'   c("logrank", "max_combo"),
#'   c("effect_size_ph", "delay", "hazard_ctrl", "n_pat_design", "recruitment", "censoring_prop"),
#'   "rejection_0.025",
#'   grid_level=2
#' )
#'
#' # use custom colour and shape scales
#' # this can be used to group methods by shape or colour
#' # this is also helpful if methods should have the same aesthetics across plots
#' my_colours <- c(
#'   logrank="black",
#'   max_combo="blue",
#'   fh_1_0="green",
#'   fh_0_1="red"
#' )
#'
#' my_shapes <- c(
#'   logrank=1,
#'   max_combo=2,
#'   fh_1_0=3,
#'   fh_0_1=3
#' )
#'
#' combined_plot(
#'   results_long,
#'   c("logrank", "max_combo", "fh_1_0", "fh_0_1"),
#'   c("effect_size_ph", "delay", "hazard_ctrl", "n_pat_design", "recruitment", "censoring_prop"),
#'   "rejection_0.025",
#'   grid_level=2,
#'   use_colours = my_colours,
#'   use_shapes = my_shapes
#' )
#'
#' combined_plot(
#'   results_long,
#'   c("logrank", "max_combo"),
#'   c("effect_size_ph", "delay", "hazard_ctrl", "n_pat_design", "recruitment", "censoring_prop"),
#'   "rejection_0.025",
#'   grid_level=2,
#'   use_colours = my_colours,
#'   use_shapes = my_shapes
#' )
#'
#' # if one has a dataset of metadata with categories of methods
#' # one could uses those two definitions
#' # colours for methods, same shapes for methods of same category
#' my_colours <- ggplot2::scale_colour_discrete()$palette(n=nrow(metadata)) |>
#'   sample() |>
#'   setNames(metadata$method)
#'
#' my_shapes <- metadata$category |>
#'   as.factor() |>
#'   as.integer() |>
#'   setNames(metadata$method)
#'
#'   combined_plot(
#' results_long,
#' c("logrank", "max_combo", "fh_1_0", "fh_0_1"),
#' c("effect_size_ph", "delay", "hazard_ctrl", "n_pat_design", "recruitment", "censoring_prop"),
#' "rejection_0.025",
#' grid_level=2,
#' use_colours = my_colours,
#' use_shapes = my_shapes
#' )
#'
#' combined_plot(
#'   results_long,
#'   c("logrank", "max_combo"),
#'   c("effect_size_ph", "delay", "hazard_ctrl", "n_pat_design", "recruitment", "censoring_prop"),
#'   "rejection_0.025",
#'   grid_level=2,
#'   use_colours = my_colours,
#'   use_shapes = my_shapes
#' )
#' }
combined_plot <- function(
    data,
    methods,
    xvars,
    yvar,
    facet_x_vars=c(),
    facet_y_vars=c(),
    split_var = 1,
    heights_plots = c(3,1),
    scale_stairs = 0.75,
    grid_level = 2,
    scales = "fixed",
    hlines = numeric(0),
    use_colours = NULL,
    use_shapes  = NULL
    ){

  if(!all(c("ggplot2", "patchwork") %in% .packages())){
    stop(gettext("Packages ggplot2 and patchwork need to be loaded for the functionality provided in the plot functions."))
  }

  facet_vars_y_sym <- syms(facet_y_vars)
  facet_vars_x_sym <- syms(facet_x_vars)
  xvars <- syms(xvars)
  yvar  <- sym(yvar)

  data <- data |>
    filter(method %in% methods)

  # remove facets in which all y values are empty
  # dont remove empty y-values in facets where there are some y-values
  # (so gaps in lines remain gaps in each facet and only completely facets are
  # dropped)
  data <- data |>
    ungroup() |>
    group_by(!!!facet_vars_x_sym) |>
    filter(!all(is.na(!!yvar))) |>
    ungroup() |>
    group_by(!!!facet_vars_y_sym) |>
    filter(!all(is.na(!!yvar))) |>
    ungroup()

  ## split lines

  len_x <- length(xvars)
  if(len_x > 1){
    lastvar <- xvars[[len_x]]
    splitvar <- xvars[[split_var]]
    data <- data |>
      group_by(method,!!!facet_vars_y_sym,!!!facet_vars_x_sym,!!splitvar) |>
      group_modify(~add_row(.x,.before = 1)) |>
      #    mutate(!!lastvar := ifelse(is.na(!!yvar),!!lastvar + .0,!!lastvar)) |>
      fill(!!!xvars[-split_var],.direction='up') |>
      ungroup()
  }


  data <- data |>
    order_combine_xvars(xvars, facet_vars=facet_x_vars, height_x_axis=scale_stairs, grid_level=grid_level)


  plot_2 <- lapply(xvars, \(xx){
      ggplot(data, aes(x=x, y=factor(format(!!xx, digits=3)), group=method)) +
        geom_step(linewidth=0.25) +
        theme_void(
          base_size = 9
        ) +
        theme(
          axis.text.y = element_text(),
          axis.title.y = element_text(angle=75),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.grid.major.y = element_line(
            linewidth = 0.125,
            colour="lightgray"
          )
        ) +
        ylab(as.character(xx))  +
        facet_grid(
          cols = vars(!!!facet_vars_x_sym)
        )
    })

  plot_2 <- Reduce(`/`, plot_2)

  plot_1 <- ggplot(data, aes(x=x, y=!!yvar, group=method, colour=method, shape=method)) +
    geom_line() +
    geom_point(size=4) +
    scale_x_discrete(breaks = attr(data, "x_axis_breaks")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )  +
    facet_grid(
      cols = vars(!!!facet_vars_x_sym),
      rows = vars(!!!facet_vars_y_sym),
      labeller = label_both,
      scales = scales
    ) +
    geom_hline(yintercept=hlines)

  if(!is.null(use_colours)){
    plot_1 <- plot_1 +
      scale_colour_manual(values=use_colours)
  }

  if(!is.null(use_shapes)){
    plot_1 <- plot_1 +
      scale_shape_manual(values=use_shapes)
  }
  (plot_1 / plot_2) + plot_layout(heights=heights_plots)
}


#' Add ggplot axis labels from labels attribute
#'
#' @param gg a ggplot object
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' test <- mtcars
#' # add a label attribute
#' attr(test$cyl, "label") <- "cylinders"
#'
#' # plot witht the variable names as axis titles
#' gg1 <- ggplot(test, aes(x=wt, y=cyl)) +
#'   geom_point()
#' gg1
#'
#' # add labels where defined in the attribute
#' gg2 <- ggplot(test, aes(x=wt, y=cyl)) +
#'   geom_point()
#'
#' gg2 <- labs_from_labels(gg2)
#' gg2
labs_from_labels <- function(gg){
  new_labels <- gg$mapping |>
    purrr::map(rlang::as_name) |>
    purrr::map(\(i){
      attr(gg$data[[i]], "label")
    }) |>
    unlist()

  gg + ggplot2::labs(!!!new_labels)
}
