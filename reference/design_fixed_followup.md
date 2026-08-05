# Create a data.frame with an example fixed design

Create a data.frame with an example fixed design

## Usage

``` r
design_fixed_followup(print = interactive())
```

## Arguments

- print:

  print code to generate parameter set?

## Value

For design_fixed_followup: a design tibble with default values invisibly

## Details

design_fixed_followup generates a default design `data.frame` for use
with generate_delayed_effect or other generate\_... functions. If print
is `TRUE` code to produce the template is also printed for copying,
pasting and editing by the user. (This is the default when run in an
interactive session.)

## Functions

- `design_fixed_followup()`: generate default fixed design

## Examples

``` r
Design <- design_fixed_followup()
Design
#>   n_trt n_ctrl followup recruitment
#> 1   150    150    730.5     182.625
```
