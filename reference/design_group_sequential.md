# Create a data.frame with an example group sequential design

Create a data.frame with an example group sequential design

## Usage

``` r
design_group_sequential(print = interactive())
```

## Arguments

- print:

  print code to generate parameter set?

## Value

For design_group_sequential: a design tibble with default values
invisibly

## Details

design_group_sequential generates a default design `data.frame` for use
with generate_delayed_effect or other generate\_... functions. If print
is `TRUE` code to produce the template is also printed for copying,
pasting and editing by the user. (This is the default when run in an
interactive session.)

## Functions

- `design_group_sequential()`: generate default group sequential design

## Examples

``` r
Design <- design_group_sequential()
Design
#>   n_trt n_ctrl followup recruitment interim_events final_events
#> 1   200    200     1461     182.625            150          300
```
