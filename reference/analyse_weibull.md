# Analyse Dataset with Weibull Regression

Analyse Dataset with Weibull Regression

## Usage

``` r
analyse_weibull(level = 0.95, alternative = "two.sided")
```

## Arguments

- level:

  confidence level for CI computation

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

an analysis function that returns a data.frame

## Details

the columns in the return are the two-sided p-value for the test of
equal medians. The estimated medians in the treatment and control group
and the estimated difference in median survival with confidence
intervals.

The estimates and tests are comstructed by fitting seperate Weibull
regression models in the treatment and control groups and then
estimating the medians and respective variances with the delta-method.

## Examples

``` r
condition <- merge(
    assumptions_delayed_effect(),
    design_fixed_followup(),
    by=NULL
  ) |>
  head(3) |>
  tail(1)
dat <- generate_delayed_effect(condition)
analyse_weibull()(condition, dat)
#> $p
#> [1] 0.07112782
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $med_trt_est
#> [1] 778.092
#> 
#> $med_trt_lower
#> [1] 631.5979
#> 
#> $med_trt_upper
#> [1] 924.5862
#> 
#> $med_ctrl_est
#> [1] 1011.487
#> 
#> $med_ctrl_lower
#> [1] 804.6256
#> 
#> $med_ctrl_upper
#> [1] 1218.348
#> 
#> $diff_med_est
#> [1] 233.3947
#> 
#> $diff_med_lower
#> [1] -20.08522
#> 
#> $diff_med_upper
#> [1] 486.8746
#> 
#> $CI_level
#> [1] 0.95
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
```
