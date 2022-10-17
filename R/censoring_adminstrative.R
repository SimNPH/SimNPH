#' Add recruitment time to Dataset
#'
#' @param dat a simulated dataset
#' @param recruitment_until time of end of recruitment
#' @param recruitment_from time of start of recruitment (defaults to 0)
#'
#' @return
#' Returns the dataset with added recruitment times.
#' @export
#'
#' @describeIn admin_censoring_time add recruitment time
#'
#' @examples
#' dat <- data.frame(t=c(0, 1, 2), trt=c(FALSE, FALSE, TRUE))
#' recruitment_uniform(dat, 7, 0)
recruitment_uniform <- function(dat, recruitment_until, recruitment_from=0){
  dat$rec_time <- runif(nrow(dat), recruitment_from, recruitment_until)
  dat
}

#' Apply Administrative Censoring After Fixed Time
#'
#' @param dat a simulated dataset
#' @param followup followup time
#'
#' @return
#' Returns the dataset with administrative censoring after `follwup`, adds the
#' attribute `followup` with the followup time to the dataset.
#'
#' @export
#' @details
#'
#' The Dataset hast to include a column `rec_time` containing the recruitment
#' time as well as the columns with the event times `t` and a column with the
#' event indicator `evt`.
#'
#' Times and event indicaotrs for patients recruited after followup are set to
#' `NA`.
#'
#' @describeIn admin_censoring_time apply administrative censoring after fixed time
#'
#' @examples
#' dat <- data.frame(
#'   t = 1:10,
#'   rec_time = rep(1:5, each=2),
#'   trt = rep(c(TRUE, FALSE), times=5),
#'   evt = rep(TRUE, times=10)
#' )
#' dat
#'
#' admin_censoring_time(dat, 4)
#' admin_censoring_time(dat, 4, keep_non_recruited = TRUE)
#'
#' dat_censored <- admin_censoring_time(dat, 5)
#' attr(dat_censored, "followup")
admin_censoring_time <- function(dat, followup, keep_non_recruited = FALSE){
  calendar_time <- dat$t + dat$rec_time
  dat$t <- pmin(calendar_time, followup) - dat$rec_time
  dat$evt[calendar_time > followup] <- FALSE

  if(keep_non_recruited){
    dat$evt[dat$t < 0] <- NA
    dat$t[dat$t < 0] <- NA_real_
  } else {
    dat <- dat[dat$t >= 0, ]
  }

  attr(dat, "followup") <- followup
  dat
}

#' Apply Administrative Censoring After Fixed Number of Events
#'
#' @param dat a simulated dataset
#' @param events number of events after which the dataset is analyzed
#'
#' @return
#' Returns the dataset with administrative censoring after `events` events, adds
#' the attribute `followup` with the followup time to the dataset.
#'
#' @export
#' @details
#'
#' The Dataset hast to include a column `rec_time` containing the recruitment
#' time as well as the columns with the event times `t` and a column with the
#' event indicator `evt`.
#'
#' Times and event indicaotrs for patients recruited after followup are set to
#' `NA`.
#'
#' @describeIn admin_censoring_time apply administrative censoring after fixed number of events
#'
#' @examples
#' dat <- data.frame(
#'   t = 1:10,
#'   rec_time = rep(2*(1:5), each=2),
#'   trt = rep(c(TRUE, FALSE), times=5),
#'   evt = rep(TRUE, times=10)
#' )
#' dat
#'
#' admin_censoring_events(dat, 4)
#' admin_censoring_events(dat, 4, keep_non_recruited = TRUE)
#'
#' dat_censored <- admin_censoring_events(dat, 4)
#' attr(dat_censored, "followup")
admin_censoring_events <- function(dat, events, keep_non_recruited = FALSE){

  # find time at which number of events is reached
  tmp <- data.frame(
    calendar_time = dat$t + dat$rec_time,
    evt           = dat$evt
  )
  tmp <- tmp[order(tmp$calendar_time), ]
  followup <- tmp$calendar_time[min(which(cumsum(tmp$evt %in% TRUE) >= events))]

  admin_censoring_time(dat, followup, keep_non_recruited=keep_non_recruited)
}

