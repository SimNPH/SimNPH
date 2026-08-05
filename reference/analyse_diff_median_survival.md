# Analyse the dataset using differnce in median survival

Analyse the dataset using differnce in median survival

## Usage

``` r
analyse_diff_median_survival(
  quant = 0.5,
  level = 0.95,
  alternative = "two.sided"
)
```

## Arguments

- quant:

  quantile for which the difference should be calculated, defaults to
  the median

- level:

  confidence level for CI computation

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

Returns an analysis function, that can be used in runSimulations

## Details

The implementation from the nph package is used, see the documentation
there for details.

The data.frame returned by the created function includes the follwing
columns:

- `p` p value of the test, see Details

- `alternative` the alternative used

- `diff_Q` estimated differnce in quantile of the suvivla functions

- `diff_Q_lower` unadjusted lower bound of the confidence interval for
  the differnce in quantile of the suvivla functions

- `diff_Q_upper` unadjusted upper bound of the confidence interval for
  the differnce in quantile of the suvivla functions

- `CI_level` the CI level used

- `quantile` quantile used for extimation

- `N_pat` number of patients

- `N_evt` number of events

## See also

[nph::nphparams](https://rdrr.io/pkg/nph/man/nphparams.html)

## Examples

``` r
condition <- merge(
  assumptions_delayed_effect(),
  design_fixed_followup(),
  by = NULL
) |>
  head(1)
dat <- generate_delayed_effect(condition)
analyse_diff_median_survival()(condition, dat)
#> $p
#> [1] 0.0005225646
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $diff_Q
#> [1] 449
#> 
#> $diff_Q_lower
#> [1] 195.3115
#> 
#> $diff_Q_upper
#> [1] 702.6885
#> 
#> $CI_level
#> [1] 0.95
#> 
#> $quantile
#> [1] 0.5
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
```
