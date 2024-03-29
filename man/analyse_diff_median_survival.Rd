% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyse_diff_median_survival.R
\name{analyse_diff_median_survival}
\alias{analyse_diff_median_survival}
\title{Analyse the dataset using differnce in median survival}
\usage{
analyse_diff_median_survival(
  quant = 0.5,
  level = 0.95,
  alternative = "two.sided"
)
}
\arguments{
\item{quant}{quantile for which the difference should be calculated, defaults to the median}

\item{level}{confidence level for CI computation}

\item{alternative}{alternative hypothesis for the tests "two.sided" or "one.sieded"}
}
\value{
Returns an analysis function, that can be used in runSimulations
}
\description{
Analyse the dataset using differnce in median survival
}
\details{
The implementation from the nph package is used, see the documentation there
for details.

The data.frame returned by the created function includes the follwing
columns:
\itemize{
\item \code{p} p value of the test, see Details
\item \code{alternative} the alternative used
\item \code{diff_Q} estimated differnce in quantile of the suvivla functions
\item \code{diff_Q_lower} unadjusted lower bound of the confidence interval for the differnce in quantile of the suvivla functions
\item \code{diff_Q_upper} unadjusted upper bound of the confidence interval for the differnce in quantile of the suvivla functions
\item \code{CI_level} the CI level used
\item \code{quantile} quantile used for extimation
\item \code{N_pat} number of patients
\item \code{N_evt} number of events
}
}
\examples{
condition <- merge(
  assumptions_delayed_effect(),
  design_fixed_followup(),
  by = NULL
) |>
  head(1)
dat <- generate_delayed_effect(condition)
analyse_diff_median_survival()(condition, dat)
}
\seealso{
\link[nph:nphparams]{nph::nphparams}
}
