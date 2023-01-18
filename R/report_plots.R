results_pivot_longer <- function(data){
  attr(data, "ERROR_msg")    <- NULL
  attr(data, "WARNING_msg")  <- NULL
  attr(data, "extra_info")   <- NULL

  methods <- attr(data, "design_names") |>
    getElement("sim") |>
    str_extract("^[^\\.]*")

  summaries <- attr(data, "design_names") |>
    getElement("sim") |>
    str_extract("(?<=\\.).*")

  pivot_spec <- tibble(.name=attr(data, "design_names")$sim, .value=summaries, method=methods)

  data |>
    rename(
      n_pat_design = n_pat
    ) |>
    pivot_longer_spec(pivot_spec)
}

order_combine_xvars <- function(data, xvars){
  data |>
    arrange(!!!xvars) |>
    unite(x, !!!xvars) |>
    mutate(
      x = fct_inorder(x)
    )
}

nested_loop_plot <- function(data, xvars, yvar, yvar_sd=NULL){
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

  ggplot(data, aes(x=x, y=!!yvar, group=method, colour=method, ymin=lower, ymax=upper, fill=method)) +
    geom_line() +
    geom_ribbon(alpha=0.3)
}


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


# usage examples, load simulation results first:
#
# library(tidyverse)
#
# plotdata <- results_pivot_longer(results_progression)
#
# plotdata |>
#   filter(
#     method %in% c("ahr", "gahr", "cox", "weighted_cox")
#   ) |>
#   nested_loop_plot(
#     xvars=c("recruitment", "n_pat_design", "interim_events", "final_events", "n_ctrl",
#             "n_trt", "hazard_ctrl", "prog_prop_trt", "prog_prop_ctrl", "hr_before_after",
#             "censoring_prop", "effect_size_ph"),
#     yvar="bias",
#     yvar_sd="sd_bias"
#   ) +
#   geom_hline(yintercept=0)
#
#
# plotdata |>
#   filter(
#     method %in% c("ahr", "gahr", "cox", "weighted_cox")
#   ) |>
#   trellis_plot(
#     xvars="effect_size_ph",
#     yvar="bias",
#     yvar_sd="sd_bias",
#     trellis_var="prog_prop_trt"
#   )
#
# plotdata |>
#   filter(
#     method %in% c("ahr", "gahr", "cox", "weighted_cox")
#   ) |>
#   trellis_plot(
#     xvars="effect_size_ph",
#     yvar="bias",
#     yvar_sd="sd_bias",
#     trellis_var="prog_prop_trt",
#     trellis_var_y="prog_prop_ctrl"
#   )
