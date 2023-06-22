#' SimNPH: Simulate Non Proportional Hazards
#'
#' This package provides several functions to simulate survival data with non
#' proportional hazards using the general purpose simulation package SimDesign.
#'
#' @docType package
#' @name SimNPH
#' @import SimDesign
#' @import survival
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices palette
#' @importFrom methods is
#' @importFrom stats anova coefficients confint convolve filter glm integrate median na.omit pchisq pnorm poisson qnorm rbinom rexp rmultinom runif sd setNames uniroot
#' @importFrom utils hasName tail
#' @useDynLib SimNPH, .registration=TRUE
NULL
#> NULL


