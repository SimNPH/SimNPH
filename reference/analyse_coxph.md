# Analyse Dataset with the Cox Protportional Hazards Model

Analyse Dataset with the Cox Protportional Hazards Model

## Usage

``` r
analyse_coxph(level = 0.95, alternative = "two.sided")
```

## Arguments

- level:

  confidence level for CI computation

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

an analyse function that returns a list with the elements

- `p` p value of the score test (two.sided) or the Wald test (one.sided)

- `alternative` the alternative used

- `coef` coefficient for `trt`

- `hr` hazard ratio for `trt`

- `hr_lower` lower 95% confidence intervall boundary for the hazard
  ratio for `trt`

- `hr_upper`lower 95% confidence intervall boundary for the hazard ratio
  for `trt`

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
analyse_coxph()(condition, dat)
#> $p
#> [1] 1.470751e-07
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $coef
#>       trt 
#> -0.627749 
#> 
#> $hr
#>      trt 
#> 0.533792 
#> 
#> $hr_lower
#> [1] 0.4209495
#> 
#> $hr_upper
#> [1] 0.6768838
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
