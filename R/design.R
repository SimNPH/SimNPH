#' Create a data.frame with an example fixed design
#'
#' @param print print code to generate parameter set?
#'
#' @return For design_fixed_followup: a design tibble with default values invisibly
#'
#' @details design_fixed_followup generates a default design `data.frame` for
#'   use with generate_delayed_effect or other generate_... functions. If print
#'   is `TRUE` code to produce the template is also printed for copying, pasting
#'   and editing by the user. (This is the default when run in an interactive
#'   session.)
#'
#' @export
#' @describeIn design_fixed_followup generate default fixed design
#'
#' @examples
#' Design <- design_fixed_followup()
#' Design
design_fixed_followup <- function(print=interactive()){
  skel <- "expand.grid(
  n_trt=150,         # 150 patients in the treatment arm
  n_ctrl=150,        # 150 patients in the control arm
  followup=m2d(24),  # followup 2 years
  recruitment=m2d(6) # recruitment time 6 months

)
"

  if(print){
    cat(skel)
  }

  invisible(
    skel |>
      str2expression() |>
      eval()
  )
}


#' Create a data.frame with an example group sequential design
#'
#' @param print print code to generate parameter set?
#'
#' @return For design_group_sequential: a design tibble with default values invisibly
#'
#' @details design_group_sequential generates a default design `data.frame` for
#'   use with generate_delayed_effect or other generate_... functions. If print
#'   is `TRUE` code to produce the template is also printed for copying, pasting
#'   and editing by the user. (This is the default when run in an interactive
#'   session.)
#'
#' @export
#' @describeIn design_group_sequential generate default group sequential design
#'
#' @examples
#' Design <- design_group_sequential()
#' Design
design_group_sequential <- function(print=interactive()){
  skel <- "expand.grid(
  n_trt=200,          # 200 patients in the treatment arm
  n_ctrl=200,         # 200 patients in the control arm
  followup=m2d(48),   # maximum followup time 4 years
  recruitment=m2d(6), # recruitment time 6 months
  interim_events=150, # interim analysis after 150 events
  final_events=300    # final analysis after 300 events
)
"

  if(print){
    cat(skel)
  }

  invisible(
    skel |>
      str2expression() |>
      eval()
  )
}
