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
  n_trt=50,      # 100 patients in the treatment arm
  n_ctrl=50,     # 100 patients in the control arm
  followup=200,  # fixed followup time
  recruitment=50 # recruitment time

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
  n_trt=50,          # 100 patients in the treatment arm
  n_ctrl=50,         # 100 patients in the control arm
  followup=200,      # maximum followup time
  recruitment=50,    # recruitment time
  interim_events=25, # events interim analysis
  final_events=50    # events final analysis
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}
