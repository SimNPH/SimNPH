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
  fluidPage(
    inputPanel(
      selectInput(paste0(prefix,"median_survival_ctrl"),label="Median Survival in Control arm (mts):",
                  choices = c(36,12,6),selected=12),
      selectInput(paste0(prefix,"recruitment"),label="Time to recruit target number of subjects (mts):",
                  choices = c(18,30),selected=18),
      selectInput(paste0(prefix,"censoring"),label="Proportion of subjects with random censoring:",
                  choices = c(0,0.1,0.3),selected=0)
    ),
    plotOutput("plot")
  )
}
