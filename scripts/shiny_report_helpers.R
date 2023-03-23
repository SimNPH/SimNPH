#' Shiny interface for selection of secondary parameters
#'
#' Returns shiny inputPanel to select filter values for secondary parameters: Median Survival in Control Arm, Recruitment Speed, and Proportion of Subjects with random Censoring
#'
#' @param prefix character to uniquely identify panel
#'
#' @return shiny inputPanel
#' @export
#'
#' @examples
interactive_input_secondary_parameters <- function(prefix=""){
  inputPanel(
    selectInput(paste0(prefix,"median_survival_ctrl"),label="Median Survival in Control arm (mts):",
                choices = c(36,12,6),selected=12),
    selectInput(paste0(prefix,"recruitment"),label="Time to recruit target number of subjects (mts):",
                choices = c(18,30),selected=18),
    selectInput(paste0(prefix,"censoring"),label="Proportion of subjects with random censoring:",
                choices = c(0,0.1,0.3),selected=0)
  )
}

#' Interactive Filters
#'
#' Utility function to define interactive (or not filters)
#'
#' If interactive is set to FALSE a character vector with static default filters will be returned, this is useful for debugging without the need to rebuild the shiny document
#'
#' @param prefix prefix to identify shiny input panel
#'
#' @return
#' @export
#'
#' @examples
interactive_filters <- function(prefix,...) {
  if (exists("session") && exists("input")) {
    c(str_c("abs(hazard_ctrl - nph::m2r(", input[[str_replace("prefix_median_survival_ctrl", "prefix_", prefix)]],"))<1e-6"),
      str_c("recruitment == ", input[[str_replace("prefix_recruitment", "prefix_", prefix)]]),
      str_c("censoring_prop == ", input[[str_replace("prefix_censoring", "prefix_", prefix)]]))
  } else {
    c("abs(hazard_ctrl - nph::m2r(12))<1e-6",
      "recruitment == 18",
      "censoring_prop == 0.0")
  }
}
