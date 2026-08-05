# Analyse the Dataset using the difference in RMST

Analyse the Dataset using the difference in RMST

## Usage

``` r
analyse_rmst_diff(max_time = NA, level = 0.95, alternative = "two.sided")
```

## Arguments

- max_time:

  time for which the RMST is calculated

- level:

  confidence level for CI computation

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

Returns an analysis function, that can be used in runSimulations

## Details

The implementation from the nph package is used, see the documentation
there for details.

`alternative` can be "two.sided" for a two sided test of equality of the
summary statistic or "one.sided" for a one sided test testing H0:
treatment has equal or shorter survival than control vs. H1 treatment
has longer survival than control.

The data.frame returned by the created function includes the follwing
columns:

- `p` p value of the test, see Details

- `alternative` the alternative used

- `rmst_diff` estimated differnce in RMST

- `rmst_diff_lower` unadjusted lower bound of the confidence interval
  for differnce in RMST

- `rmst_diff_upper` unadjusted upper bound of the confidence interval
  for differnce in RMST

- `CI_level` the CI level used

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
analyse_rmst_diff()(condition, dat)
#> $p
#> [1] 0.02254084
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $rmst_diff
#> [1] 297.9007
#> 
#> $rmst_diff_lower
#> [1] 41.94202
#> 
#> $rmst_diff_upper
#> [1] 553.8594
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
