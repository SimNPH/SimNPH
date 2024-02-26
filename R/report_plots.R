#' Functions for Plotting and Reporting Results
#'
#' @describeIn results_pivot_longer pivot simulation results into long format
#'
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
#' \donttest{
#' data("combination_tests_delayed")
#'
#' combination_tests_delayed |>
#'   results_pivot_longer() |>
#'   head()
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
    stringr::str_remove(stringr::str_c(methods, "."))

  include <- !(methods %in% exclude_from_methods)

  pivot_spec <- tibble::tibble(
    .name=attr(data, "design_names")$sim[include],
    .value=summaries[include],
    method=methods[include]
    )

  result <- data |>
    dplyr::rename_with(
      .fn = \(name){rep("n_pat_design", length(name))},
      .cols=dplyr::any_of("n_pat")
    ) |>
    tidyr::pivot_longer_spec(pivot_spec)
}

#' @describeIn results_pivot_longer Nested Loop Plot with optional Facets
#'
#' @param data for results_pivot_longer: simulation result as retured by SimDesign, for combined_plot: simulation results in long format, as returned by `results_pivot_longer`.
#' @param methods methods to include in the plot
#' @param xvars orderd vector of variable names to display on the x axis
#' @param yvar variable name of the variable to be displayed on the y axis (metric)
#' @param facet_x_vars vector of variable names to create columns of facets
#' @param facet_y_vars vector of variable names to create rows of facets
#' @param split_var index of xvars along groups of which the plot should be split
#' @param heights_plots relative heights of the main plot and the stairs on the bottom
#' @param scale_stairs this argument is deprecated and will be ignored
#' @param grid_level depth of loops for which the grid-lines are drawn
#' @param scales passed on to facet_grid
#' @param hlines position of horizontal lines, passed as `yintercept` to
#'   `geom_hline`
#' @param use_colours optional named vector of colours used in `scale_colour_manual`
#' @param use_shapes optional named vector of shapes used in `scale_shape_manual`
#'
#' @return a ggplot/patchwork object containing the plots
#' @export
#'
#' @details `use_colours` and `use_shapes` both use the `method` variable in their respective aesthetics.
#'
#' @examples
#' \donttest{
#' library("ggplot2")
#' library("patchwork")
#' data("combination_tests_delayed")
#'
#' results_long <- results_pivot_longer(combination_tests_delayed)
#'
#' # plot the rejection rate of two methods
#' combined_plot(
#'   results_long,
#'   c("logrank", "mwlrt", "maxcombo"),
#'   c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
#'   "rejection_0.025",
#'   grid_level=2
#' )
#'
#' # use custom colour and shape scales
#' # this can be used to group methods by shape or colour
#' # this is also helpful if methods should have the same aesthetics across plots
#' my_colours <- c(
#'   logrank="black",
#'   mwlrt="blue",
#'   maxcombo="green"
#' )
#'
#' my_shapes <- c(
#'   logrank=1,
#'   mwlrt=2,
#'   maxcombo=2
#' )
#'
#' combined_plot(
#'   results_long,
#'   c("logrank", "mwlrt", "maxcombo"),
#'   c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
#'   "rejection_0.025",
#'   grid_level=2,
#'   use_colours = my_colours,
#'   use_shapes = my_shapes
#' )
#'
#' # if one has a dataset of metadata with categories of methods
#' # one could uses those two definitions
#' # colours for methods, same shapes for methods of same category
#' metadata <- data.frame(
#'   method = c("logrank", "mwlrt", "maxcombo"),
#'   method_name = c("logrank test", "modestly weighed logrank test", "maxcombo test"),
#'   category = c("logrank test", "combination test", "combination test")
#' )
#'
#' my_colours <- ggplot2::scale_colour_discrete()$palette(n=nrow(metadata)) |>
#'   sample() |>
#'   setNames(metadata$method)
#'
#' my_shapes <- metadata$category |>
#'   as.factor() |>
#'   as.integer() |>
#'   setNames(metadata$method)
#'
#' combined_plot(
#'   results_long,
#'   c("logrank", "mwlrt", "maxcombo"),
#'   c("hr", "n_pat_design", "delay", "hazard_ctrl", "recruitment"),
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
    scale_stairs = NULL,
    grid_level = 2,
    scales = "fixed",
    hlines = numeric(0),
    use_colours = NULL,
    use_shapes  = NULL
){

  if(!is.null(scale_stairs)){
    warning("scale_stairs argument is deprecated and will be ignored. This argument will be removed in a future release.")
  }

  stopifnot(split_var <= length(xvars))
  stopifnot(split_var > 0)

  stopifnot(grid_level > 0)
  grid_level <- min(grid_level, length(xvars))

  facet_vars_y_sym <- rlang::syms(facet_y_vars)
  facet_vars_x_sym <- rlang::syms(facet_x_vars)
  yvar  <- rlang::sym(yvar)

  data <- data |>
    subset(method %in% methods)

  # create combined x variable
  data$x <- do.call(interaction, c(data[, xvars], lex.order=TRUE, drop=TRUE)) |>
    as.integer()

  # remove rows with missing yvar
  data <- data |>
    subset(!is.na(get(yvar)))

  # split lines
  group_vars <- xvars[1:split_var]
  data$x_split <- do.call(interaction, c(data[, group_vars], lex.order=TRUE, drop=TRUE)) |>
    as.integer()

  # grid breaks
  grid_vars <- xvars[1:grid_level]
  data$x_grid <-  do.call(interaction, c(data[, grid_vars], lex.order=TRUE))
  grid_breaks <- data$x[order(data$x)][!duplicated(data$x_grid[order(data$x)])]

  plot_1 <- ggplot2::ggplot(data, ggplot2::aes(
    x=x,
    y=!!yvar,
    group=interaction(method, x_split),
    colour=method,
    shape=method
  )) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size=4) +
    ggplot2::scale_x_continuous(
      breaks = grid_breaks,
      minor_breaks = NULL,
      expand = ggplot2::expansion(0,0)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )  +
    ggplot2::facet_grid(
      cols = dplyr::vars(!!!facet_vars_x_sym),
      rows = dplyr::vars(!!!facet_vars_y_sym),
      labeller = ggplot2::label_both,
      scales = scales
    ) +
    ggplot2::geom_hline(yintercept=hlines)

  data_plot2 <- data[!duplicated(do.call(interaction, data[,c("x", facet_x_vars)])), ]

  plot_2 <- lapply(xvars, \(xx){
    data_plot2 <- data_plot2 |>
      within(
        tmp_yvar <- factor(format(get(xx), digits=3))
      )

    ggplot2::ggplot(data_plot2, ggplot2::aes(x=x, y=tmp_yvar, group=method)) +
      ggplot2::geom_step(linewidth=0.25) +
      # ggplot2::geom_point(shape=4) +
      ggplot2::scale_x_continuous(
        breaks = grid_breaks,
        minor_breaks = NULL,
        expand = ggplot2::expansion(0,0)
      ) +
      ggplot2::theme_void(
        base_size = 9
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(),
        axis.title.y = ggplot2::element_text(angle=75),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(
          linewidth = 0.125,
          colour="lightgray"
        )
      ) +
      ggplot2::ylab(as.character(xx))  +
      ggplot2::facet_grid(
        cols = dplyr::vars(!!!facet_vars_x_sym)
      )
  })

  plot_2 <- patchwork::wrap_plots(plot_2, ncol=1)

  if(!is.null(use_colours)){
    plot_1 <- plot_1 +
      ggplot2::scale_colour_manual(values=use_colours)
  }

  if(!is.null(use_shapes)){
    plot_1 <- plot_1 +
      ggplot2::scale_shape_manual(values=use_shapes)
  }
  result <- (plot_1 / plot_2) + patchwork::plot_layout(heights=heights_plots)
  result
}


#' Add ggplot axis labels from labels attribute
#'
#' @param gg a ggplot object
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' library("ggplot2")
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
#' }
labs_from_labels <- function(gg){
  new_labels <- gg$mapping |>
    purrr::map(rlang::as_name) |>
    purrr::map(\(i){
      attr(gg$data[[i]], "label")
    }) |>
    unlist()

  gg + ggplot2::labs(!!!new_labels)
}
