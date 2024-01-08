#' SimNPH: Simulate Non Proportional Hazards
#'
#' This package provides several functions to simulate survival data with non
#' proportional hazards using the general purpose simulation package SimDesign.
#'
#' @docType package
#' @name SimNPH
#' @import SimDesign
#' @import survival
#' @importFrom grDevices palette
#' @importFrom methods is
#' @importFrom stats anova coefficients confint convolve glm integrate median na.omit pchisq pnorm poisson qnorm rbinom rexp rmultinom runif sd setNames uniroot
#' @importFrom utils hasName tail
NULL
#> NULL

# declaring varibles to avoid R CMD check notes.
# those are mostly column names that occur in with, within, subset functions,
# dplyr verbs and ggplot calls.
globalVariables(c(
  "interval", "trt", "method", "x", "name", "value", "level", "y", "n_pat",
  "surv_a", "surv_b", "haz_a", "haz_b", "hr", "x_split"
))


