# Analyse Dataset with accelarated failure time models

Analyse Dataset with accelarated failure time models

## Usage

``` r
analyse_aft(level = 0.95, dist = "weibull", alternative = "two.sided")
```

## Arguments

- level:

  confidence level for CI computation

- dist:

  passed to survival::survreg

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

an analyse function that returns a list with the elements

- `p` p value of the score test (two.sided) or the Wald test (one.sided)

- `alternative` the alternative used

- `coef` coefficient for `trt`

- `lower` lower 95% confidence intervall boundary for the coefficient

- `upper`lower 95% confidence intervall boundary for the coefficient

- `CI_level` the CI level used

- `N_pat` number of patients

- `N_evt` number of events

## Details

`alternative` can be "two.sided" for a two sided test of equality of the
summary statistic or "one.sided" for a one sided test testing H0:
treatment has equal or shorter survival than control vs. H1 treatment
has longer survival than control.

## Examples

``` r
condition <- merge(
  assumptions_delayed_effect(),
  design_fixed_followup(),
  by = NULL
) |>
  head(1)
dat <- generate_delayed_effect(condition)
analyse_aft()(condition, dat)
#> $p
#> [1] 0.01392478
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $coef
#>       trt 
#> 0.2875961 
#> 
#> $lower
#> [1] 0.05838397
#> 
#> $upper
#> [1] 0.5168082
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
analyse_aft(dist="lognormal")(condition, dat)
#> $p
#> [1] 0.03025141
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $coef
#>       trt 
#> 0.3155178 
#> 
#> $lower
#> [1] 0.0301161
#> 
#> $upper
#> [1] 0.6009195
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
