# SimNPH

Simulate *N*on *P*roportional *H*azards

This package provides several functions to simulate survival data with non
proportional hazards using the general purpose simulation package
[SimDesign](https://cran.r-project.org/package=SimDesign).

This package follows the structure of SimDesing in that it provides seperate
functions to generate, analyse and summarise.

## Contributing

### Project Structure

The folder structure follows the structure of a typical R package. Documentation
and the NAMESPACE file are generated with roxygen, so don't edit those files
manually. Unit tests using
[testthat](https://cran.r-project.org/package=testthat) are in `tests/testthat`
with filenames corresponding to the filenames of the functions being tested.

### Coding Style

Please roughly keep to the [tidyverse styleguide](https://style.tidyverse.org/)
