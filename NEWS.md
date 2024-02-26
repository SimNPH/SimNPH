# SimNPH 0.5.5

* Sped up calculation of treatment arm hazard from effect size and calculation
  of censoring rate from censoring proportion for progression secenarios.

# SimNPH 0.5.4

* Reworked the `combined_plot` function. Fixed a bug where the `split_var`
  argument was not working. Deprecated the defunct argument `scale_stairs`,
  providing this argument now gives a warning. Improved speed and removed
  dependencies.

# SimNPH 0.5.3

* Changed the calculation of the real summary measures for the progression
  scenario: corrected one mistake in the calculation of the gAHR and the
  functions now use the new, faster implementation from the miniPCH package.

# SimNPH 0.5.2

* Added option to suppress output in design and assumption templates.
* Changed `\dontrun` to `\donttest` in examples.
* Improved examples in documentation.
* Cleaned up package metadata.

# SimNPH 0.5.1

* Added references to the preprint of the paper to description and citation. 
  Added references to the preprint and the shiny app to the readme.
* First CRAN submission.

# SimNPH 0.5.0

* Moved compiled functions to new package "miniPCH" to clarify structure and to
  be compliant with Rcpp licence. 

# SimNPH 0.4.3

* Added a workaround for a bug in the `nph` package, so the package no longer
  requires a custom patched version of `nph` to work.

# SimNPH 0.4.2

* Built first version of pkdown Site.
