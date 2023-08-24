# SimNPH

**Sim**ulate **N**on **P**roportional **H**azards

This package provides several functions to simulate survival data with non
proportional hazards using the general purpose simulation package
[SimDesign](https://cran.r-project.org/package=SimDesign).

This package follows the structure of SimDesing and provides functions that
can readily be used as `generate`, `analyse` and `summarise` arguments in 
`SimDesign`'s `runSimulation` function. 

## Usage

### Installation

The current unstable version can be installed with:

```
remotes::install_git("https://github.com/SimNPH/SimNPH.git")
```

Some functions require functions from the `nph` package, that have an error in
the current CRAN version. A patched version can be installed from this
repository with the code below. The errors should be fixed in `nph` in an
upcoming CRAN release.

```
install.packages(
  "https://github.com/SimNPH/SimNPH/raw/sims_week12/nph_2.1-1.tar.gz",
  repos=NULL
)
```

### Getting Started

Documentation for all functions can be found in the respective help topics in
the package after installation or
[here](https://simnph.github.io/SimNPH/reference/index.html)

Some examples of data generation, testing and estimation can be found in this
[vignette](https://simnph.github.io/SimNPH/articles/vignettes_prebuild/simple_example.html).

## Results of the Simulation Study

The results of the simulation study done by the CONFIRMS consortium in the
course of the EMA tender will be published in an upcoming paper and are
presented in a [Shiny App](https://sny.cemsiis.meduniwien.ac.at/~mp314/rsnph/)
the source code of which is also on
[github](https://github.com/SimNPH/ShinySimNPH)
