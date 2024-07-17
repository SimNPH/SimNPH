# SimNPH

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/SimNPH/SimNPH/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SimNPH/SimNPH?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/SimNPH)](https://CRAN.R-project.org/package=SimNPH)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/SimNPH/SimNPH/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SimNPH/SimNPH/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Sim**ulate **N**on **P**roportional **H**azards

This package provides several functions to simulate survival data with non
proportional hazards using the general purpose simulation package
[SimDesign](https://cran.r-project.org/package=SimDesign).

This package follows the structure of SimDesing and provides functions that
can readily be used as `generate`, `analyse` and `summarise` arguments in 
`SimDesign`'s `runSimulation` function. 

## Usage

### Installation

The current development version can be installed with:

```
remotes::install_git("https://github.com/SimNPH/SimNPH.git")
```

### Getting Started

Documentation for all functions can be found in the respective help topics in
the package after installation or
[here](https://simnph.github.io/SimNPH/reference/index.html)

Some examples of data generation, testing and estimation can be found in this
[vignette](https://simnph.github.io/SimNPH/articles/vignettes_prebuild/simple_example.html).

## Results of the Simulation Study

The results of the simulation study done by the CONFIRMS consortium in the
course of the EMA tender will be published in an upcoming paper 
([preprint on arXiv](https://arxiv.org/abs/2310.05622)) and are presented in a 
[shinylive App](https://simnph.github.io/SimResultsShinylive/about.html). 

## Reproducing the Results of the Simulation study

The parameters for the simulations were set using the following scripts:

* `scripts/set_parameters_delayed_effect.R`
* `scripts/set_parameters_crossing_hazards.R`
* `scripts/set_parameters_subgroup.R`
* `scripts/set_parameters_progression.R`

The simulations were then run using the following scripts:

* `scripts/run_simulations_delayed_effect.R`
* `scripts/run_simulations_crossing_hazards.R`
* `scripts/run_simulations_subgroup.R`
* `scripts/run_simulations_progression.R`

Which in turn use `scripts/run_simulations_common.R`.

All scripts are executed with the base directory of this repository as working
directory.

Changes since the conduct of the simulation study:

* The version number of SimNPH was increased, to reproduce the simulation study, 
  remove the check from `run_simulations_...`
* Simulation parameters and results were previously stored in the folder `data`.
  This folder was renamed to `data_sim_study` because `data` is reserved for the
  use in R packages. Change the paths accordingly.
* The `run_simulations_...` scripts use the parameters from a specified data. If
  you want to use parameter-values you generated yourself, update the path 
  accordingly.
  
## Re-Using the Scenarios

If you want to re-use the scenarios from the simulation study to investigate
additional methods, use the scripts and documentation provided in the
[SimTemplate](https://github.com/SimNPH/SimTemplate) github repository. Also
refer to the documentation on
[contributing](https://simnph.github.io/SimNPH/CONTRIBUTING.html). If you do so
please cite the simulation study, you can get a bibtex entry with
`citation(package="SimNPH")`.

