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
remotes::install_github("SimNPH/SimNPH")
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

### Project Setup

1. clone this repository
2. install needed packages by running `install_dependencies.R`
3. do everything else (run tests, build documentation, build the package,
   build vignettes, ...) using `devtools` or RStudio

### Coding Style

Please roughly keep to the [tidyverse styleguide](https://style.tidyverse.org/)

When using functions from other packages (`nph`, `cmprsk`, ...) refer to them 
by namespace (for example `nph::rSurv_fun`). In the `DESCRIPTION` file: don't 
add the package under `Depends`, do add the package under `Suggests`.

### Adding functions

Functions to generate datasets should be named beginning with `generate_` and a 
descriptive name of the scenario they should simulate. Functions to analyse 
should be named like `analyse_` followed by the name of the test or model they
implement.

Generate functions should come with a `design_skeleton_...` function that
outputs code and returns a `tibble` for an example design dataset to be used in
`SimDesign` so users can copy paste the needed variable names and just fill in
their values. Please add this function to the same file as the `generate_`
function and document it in the documentation of the `generate_` function using
the roxygen tag `@describein`. As a reference, look at
`generate_delayed_effect.R`.

Analyse functions have to return a list in order to work with the summarise
functions. Lists can also include more complex objects, for example, the
`analyse_group_sequential` function output the whole results of all stages in
its return value. If your function requires summaries that are more complex than
the general functions for estimators and tests please provide a summarise
function in the same file as the analysis method, see again
`analyse_group_squential.R` for an example.

When adding functions provide documentation for all functions exported in the
package namespace and make sure `devtools::document()` runs without errors or
warnings.

Please also add tests that at least test if the functions can be called without 
error and produce datasets with the right columns (name and datatype). To add a 
tests file with the correct name in the correct place use `usethis::use_test()`.
