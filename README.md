# SimNPH

**Sim**ulate **N**on **P**roportional **H**azards

This package provides several functions to simulate survival data with non
proportional hazards using the general purpose simulation package
[SimDesign](https://cran.r-project.org/package=SimDesign).

This package follows the structure of SimDesing and provides functions that
can readily be used as `generate`, `analyse` and `summarise` arguments in 
`SimDesign`'s `runSimulation` function. 

## Usage

**Please be aware, that this package is in an early stage of development.
Neither is it feature complete nor is the interface completely stable, expect
bugs, changes and missing functionality!**

### Installation

The current unstable version can be installed with:

```
remotes::install_git("https://github.com/SimNPH/SimNPH.git")
```

### Getting Started

Some examples of data generation, testing and estimation can be found in this
[vignette](https://simnph.github.io/SimNPH/doc/simple_example.html).

Documentation for all functions can be found in the respective help topics.

## Contributing

### Project Structure

The folder structure follows the structure of a typical R package. Documentation
and the NAMESPACE file are generated with roxygen, so don't edit those files
manually. Unit tests using
[testthat](https://cran.r-project.org/package=testthat) are in `tests/testthat`
with filenames corresponding to the filenames of the functions being tested. 

This package uses the devtools, with, testthat and roxygen2 packages to
facilitate package development, all of those packages are conveniently
integrated within rstudio.

### Coding Style

Please roughly keep to the [tidyverse styleguide](https://style.tidyverse.org/)

When using functions from other packages (`nph`, `cmprsk`, ...) refer to them 
by namespace (for example `nph::rSurv_fun`). In the `DESCRIPTION` file: don't 
add the package under `Depends`, do add the package under `Suggests` or 
`Imports`.

### Adding functions

#### Data Generating Models

Functions to generate datasets should be named beginning with `generate_` and a 
descriptive name of the scenario they should simulate. 

In order to work with the SimDesign package, the function has to take the
arguments `condition` and `fixed_objects`. `condition` is a one row `data.frame`
that contains all parameters for the scenario that is simulated. `fixed_objects`
is a list in which variables, that are the same across all secenarios are
passed.

Generate functions should come with an `assumptions_...` function that outputs
code to generate a default dataset with all the parameters relevant for the date
generation model, which can be used to construct a `Design` dataset. The
`assumptions_...` function should also return the dataset created by the printed
code invisibly. Please implement this function in the same R file as the
`generate_...` function and document it in the documentation of the
`generate_...` function using the roxygen tag `@describein`.

Here's boilerplate for a generate and assumptions functions, you can just copy
it and fill in your own code:
``` r
#' Generate Dataset ...
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condidtion has to contain the following columns:
#'
#'   * X ...
#'   * Y ...
#'   * ...
#'
#' @return
#' For generate_x: A dataset with the columns t (time) and trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_x simulates a dataset with ...
#'
#' @examples
generate_x <- function(condition, fixed_objects=NULL){

  data.frame(t=numeric(0), trt=integer(0), evt=logical(0))
}

#' Create an empty assumtions data.frame for generate_x
#'
#' @return For assumptions_x: a design tibble with default values invisibly
#'
#' @details assumptions_x prints the code to generate a default
#'   design tibble for use with generate_delayed_effect and returns the
#'   evaluated code invisibly. This function is intended to be used to copy
#'   paste the code and edit the parameters.
#'
#' @export
#' @describeIn generate_x generate default design tibble
#'
#' @examples
#' Design <- assumptions_x()
#' Design
assumptions_delayed_effect <- function(){
  skel <- "expand.grid(
  X = 1:3,
  Y = c("a", "b")
)
"

cat(skel)
invisible(
  skel |>
    str2expression() |>
    eval()
)
}

#' Calculate true summary statistics for scenarios with delayed treatment effect
#'
#' @param Design Design data.frame for x
#' @param cutoff_stats Cutoff time for rmst and average hazard ratios
#' @param fixed_objects fixed objects not used for now
#'
#' @return For true_summary_statistics_x: the design data.frame
#'   passed as argument with the additional columns:
#' * `rmst_trt` rmst in the treatment group
#' * `median_surv_trt` median survival in the treatment group
#' * `rmst_ctrl` rmst in the control group
#' * `median_surv_ctrl` median survial in the control group
#' * `gAHR` geometric average hazard ratio
#' * `AHR` average hazard ratio
#'
#' @export
#'
#' @describeIn generate_x  calculate true summary statistics for ...
#'
#' @examples
true_summary_statistics_x <- function(Design, cutoff_stats=10, fixed_objects=NULL){

  true_summary_statistics_x_rowwise <- function(condition, cutoff_stats){
    res <- data.frame(
      rmst_trt = NA_real_,
      medial_surv_trt = NA_real_,
      rmst_ctrl = NA_real_,
      median_surv_ctrl = NA_real_,
      gAHR = NA_real_,
      AHR = NA_real_
    )
    
    res
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    mapply(FUN=true_summary_statistics_x_rowwise, cutoff_stats = cutoff_stats, SIMPLIFY = FALSE)

  Design <- do.call(rbind, Design)

  Design
}
```

#### Analysis Functions

Analysis functions should take parameters regarding the analysis method (for
example weights of a weighted logrank test) and return a function that can be
used in SimDesigns `runSimulation`. To use a function as analysis function in
`runSimulation` it has to take the arguments `condition`, `dat`,
`fixed_objects`. Condition and fixed_objects are the same as in the
`generate_...` functions. `dat` is the Dataset generated by the `generate_...`
function. `dat` contains at least the columns `t`, `trt`, `evt` for time,
treatment and event.

In order to work with the summarise functions in this package, analysis
functions have to return a list. (Unfortunately vectors and data.frames are
concatenated by SimDesign in a not very useful way, lists are passed on as is.)

If your function requires summaries that are more complex than
the general functions for estimators and tests please provide a summarise
function in the same file as the analysis method.

``` r
#' Create Analyse Functions for ...
#'
#' @param X
#'
#' @return an analyse function that can be used in runSimulation
#' @export
#'
analyse_X <- function(X){
  function(condition, dat, fixed_objects = NULL){
    result_tmp <- list(
      A = NA_real_,
      N_pat=nrow(dat),
      N_evt=sum(dat$evt)
    )

    result_tmp
  }
}


#' Summarise Output from Analyse Functions for ...
#'
#' @param X
#'
#' @describeIn analyse_X Summarise Output from Analyse X
#'
#' @return
#' Returns a function with the arguments:
#'  * condition
#'  * results
#'  * fixed objects
#'
#' that can be passed to create_summarise_function or to
#' SimDesign::runSimulation and that returns a `data.frame` with the columns
#'  * `Y` ...
#'  * ...
#'
#' @export
#'
#' @examples
summarise_x <- function(name=NULL){
  res <- function(condition, results, fixed_objects=NULL){
    data.frame(
      "Y"=NA_real_
    )
  }

  attr(res, "name") <- name

  res
}
```

#### Functions for Censoring and Experiment Designs

Different study setups are also implemented via functions that manipulate the
simulated dataset. For example to apply administrative censoring after a fixed
followup time, first a recruitment time is sampled from a uniform distribution
using the function `recruitment_uniform` then the observations are censored
according to their survival time, the recruitment time and the follwup time
using the function `admin_censoring_time`. This can also be combined with random
censoring with a fixed rate, ...

These functions should take the simulated dataset as the first argument and
further arguments at later positions so that they can be used in pipes.

``` r
#' Modify a simulated dataset by ...
#'
#' @param dat a simulated dataset
#' @param X ...
#' @param ...
#'
#' @return
#' Returns the dataset with ...
#' @export
#'
#' @examples
censoring_function_x <- function(dat, X){
  dat$X <- X
  dat
}
```

#### Documentation and Tests

When adding functions provide documentation for all functions exported in the
package namespace and make sure `devtools::document()` runs without errors or
warnings.

Please also add tests that at least test if the functions can be called without 
error and produce datasets with the right columns (name and datatype). To add a 
tests file with the correct name in the correct place use `usethis::use_test()`.

