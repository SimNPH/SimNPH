#' Plot of survival, hazard and hazard ratio of two groups as a function of time using ggplot and patchwork
#'
#' @param A mixpch object for group 1 (reference)
#' @param B mixpch object for group 2
#' @param main Title for the overall plot
#' @param sub Subtitle for the overall plot
#' @param group_names Group Names
#' @param lab_time Title for the time axis
#' @param lab_group Title group legend
#' @param trafo_time Function to transform time
#'
#' @return a `patchwork` object as defined in the patchwork package
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(patchwork)
#' B <- pchaz(c(0, 10, 100), c(0.1, 0.05))
#' A <- pchaz(c(0, 100), c(0.1))
#' shhr_gg(A, B)
#' shhr_gg(A, B, lab_time="Months", trafo_time=d2m)
#' }
shhr_gg <- function(A, B, main=NULL, sub=NULL, group_names=c("control", "treatment"), lab_time="Days", lab_group="Group", trafo_time=identity){

  if( !(requireNamespace("ggplot2", quietly = TRUE) & requireNamespace("patchwork", quietly = TRUE)) ){
    message("Packages ggplot2 and patchwork required for plotting functionality.")
    return(invisible(NULL))
  }

  plotdata <- data.frame(
    t     = sort(unique(union(A$t, B$t)))
  )

  plotdata$surv_a <- A$S[match(plotdata$t, A$t)]
  plotdata$surv_b <- B$S[match(plotdata$t, B$t)]
  plotdata$haz_a <- A$haz[match(plotdata$t, A$t)]
  plotdata$haz_b <- B$haz[match(plotdata$t, B$t)]
  plotdata$hr <- plotdata$haz_b / plotdata$haz_a

  plotdata$t <- trafo_time(plotdata$t)

  gg_surv <- ggplot2::ggplot(plotdata, ggplot2::aes(x=t)) +
    ggplot2::geom_line(ggplot2::aes(y=surv_a, colour=group_names[1], lty=group_names[1]), linewidth=1.3) +
    ggplot2::geom_line(ggplot2::aes(y=surv_b, colour=group_names[2], lty=group_names[2]), linewidth=1.3) +
    ggplot2::labs(
      x=lab_time,
      y="Survival",
      colour=lab_group,
      lty=lab_group
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0,1)
    )

  gg_haz <- ggplot2::ggplot(plotdata, ggplot2::aes(x=t)) +
    ggplot2::geom_line(ggplot2::aes(y=haz_a, colour=group_names[1], lty=group_names[1]), linewidth=1.3) +
    ggplot2::geom_line(ggplot2::aes(y=haz_b, colour=group_names[2], lty=group_names[2]), linewidth=1.3) +
    ggplot2::labs(
      x=lab_time,
      y="Hazard",
      colour=lab_group,
      lty=lab_group
    ) +
    ggplot2::expand_limits(y=0)


  gg_hr <- ggplot2::ggplot(plotdata, ggplot2::aes(x=t, y=hr)) +
    ggplot2::geom_line(linewidth=1.3) +
    ggplot2::labs(
      x=lab_time,
      y="Hazard ratio"
    ) +
    ggplot2::expand_limits(y=1)

  tmp_colours <- palette()[c(1,3)]
  names(tmp_colours) <- group_names

  tmp_lty <- c(1, 3)
  names(tmp_lty) <- group_names

  patchwork::wrap_plots(gg_surv, gg_haz, gg_hr) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(main, subtitle = sub) &
    ggplot2::theme_bw() &
    ggplot2::theme(legend.position="bottom") &
    ggplot2::scale_colour_manual(
      values = tmp_colours
    ) &
    ggplot2::scale_linetype_manual(
      values = tmp_lty
    )
}
