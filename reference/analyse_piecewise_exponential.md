# Create Analyse function for piecewise exponential model

Create Analyse function for piecewise exponential model

## Usage

``` r
analyse_piecewise_exponential(cuts, testing_only = FALSE)
```

## Arguments

- cuts:

  interval boundaries for the piecewise exponential model

- testing_only:

  if set to `TRUE` omits all statistics in the intervals and just
  returns the p value of the global test.

## Value

an analyse function that can be used in runSimulation

## Details

If there's any time interval no patients ever enter, NA is returned for
all time intervals. This behavior will likely change in future package
versions.

## Examples

``` r
condition <- merge(
   assumptions_delayed_effect(),
   design_fixed_followup(),
   by=NULL
 ) |>
 head(1)
dat <- generate_delayed_effect(condition)
analyse_piecewise_exponential(cuts=c(90, 360))(condition, dat)
#> $p
#> [1] 0.0005337612
#> 
#> $hr_s
#> # A tibble: 1 × 3
#>   `trt:intervalI1` `trt:intervalI2` `trt:intervalI3`
#>              <dbl>            <dbl>            <dbl>
#> 1            0.515            0.575            0.636
#> 
#> $p_s
#> # A tibble: 1 × 3
#>   `trt:intervalI1` `trt:intervalI2` `trt:intervalI3`
#>              <dbl>            <dbl>            <dbl>
#> 1            0.129           0.0432         0.000815
#> 
#> $hr_s_lower
#> # A tibble: 1 × 3
#>   `trt:intervalI1` `trt:intervalI2` `trt:intervalI3`
#>              <dbl>            <dbl>            <dbl>
#> 1            -1.57            -1.10           -0.716
#> 
#> $hr_s_upper
#> # A tibble: 1 × 3
#>   `trt:intervalI1` `trt:intervalI2` `trt:intervalI3`
#>              <dbl>            <dbl>            <dbl>
#> 1            0.169          -0.0238           -0.186
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
#> $interval_table
#>     
#>        0   1
#>   I1 277  23
#>   I2 221  56
#>   I3   0 221
#> 
```
