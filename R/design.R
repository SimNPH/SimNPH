#' Create a data.frame with an example fixed design
#'
#' @return For design_fixed_followup: a design tibble with default values invisibly
#'
#' @details design_fixed_followup prints the code to generate a default
#'   design tibble for use with generate_delayed_effect or other generate_...
#'   functions and returns the evaluated code invisibly. This function is
#'   intended to be used to copy paste the code and edit the parameters.
#'
#' @export
#' @describeIn design_fixed_followup generate default fixed design
#'
#' @examples
#' Design <- design_fixed_followup()
#' Design
design_fixed_followup <- function(){
  skel <- "expand.grid(
  n_trt=150,         # 150 patients in the treatment arm
  n_ctrl=150,        # 150 patients in the control arm
  followup=m2d(24),  # followup 2 years
  recruitment=m2d(6) # recruitment time 6 months

)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}


#' Create a data.frame with an example group sequential design
#'
#' @return For design_group_sequential: a design tibble with default values invisibly
#'
#' @details design_group_sequential prints the code to generate a default
#'   design tibble for use with generate_delayed_effect or other generate_...
#'   functions and returns the evaluated code invisibly. This function is
#'   intended to be used to copy paste the code and edit the parameters.
#'
#' @export
#' @describeIn design_group_sequential generate default group sequential desing
#'
#' @examples
#' Design <- design_group_sequential()
#' Design
design_group_sequential <- function(){
  skel <- "expand.grid(
  n_trt=200,          # 200 patients in the treatment arm
  n_ctrl=200,         # 200 patients in the control arm
  followup=m2d(48),   # maximum followup time 4 years
  recruitment=m2d(6), # recruitment time 6 months
  interim_events=150, # interim analysis after 150 events
  final_events=300    # final analysis after 300 events
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}
