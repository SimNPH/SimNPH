# Analyse the Dataset using difference or quotient of milestone survival

Analyse the Dataset using difference or quotient of milestone survival

## Usage

``` r
analyse_milestone_survival(
  times,
  what = "quot",
  level = 0.95,
  alternative = "two.sided"
)
```

## Arguments

- times:

  followup times at which the the survival should be compared

- what:

  "quot" for quotient and "diff" for differnce of surival probabilities

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

- `milestone_surv_ratio` / `milestone_surv_diff` ratio or differnce of
  survival probabilities

- `times` followup times at which the the survival are compared

- `N_pat` number of patients

- `N_evt` number of events

- `p` p value for the H0 that the ratios are 1 or the differnce is 0
  respectively

- `alternative` the alternative used

- `milestone_surv_ratio_lower` / `milestone_surv_diff_lower` upper/lower
  CI for the estimate

- `milestone_surv_ratio_upper` / `milestone_surv_diff_upper` upper/lower
  CI for the estimate

- `CI_level` the CI level used

## See also

[nph::nphparams](https://rdrr.io/pkg/nph/man/nphparams.html)

## Examples

``` r
condition <- merge(
    assumptions_delayed_effect(),
    design_fixed_followup(),
    by=NULL
  ) |>
  head(1)
dat <- generate_delayed_effect(condition)
analyse_milestone_survival(3:5)(condition, dat)
#> Warning: replacement has length zero
#> $p
#> [1] NA
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $milestone_surv_ratio
#> [1] NA
#> 
#> $milestone_surv_ratio_lower
#> [1] NA
#> 
#> $milestone_surv_ratio_upper
#> [1] NA
#> 
#> $CI_level
#> [1] 0.95
#> 
#> $times
#> [1] 3 4 5
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
analyse_milestone_survival(3:5, what="diff")(condition, dat)
#> Warning: replacement has length zero
#> $p
#> [1] NA
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $milestone_surv_diff
#> [1] NA
#> 
#> $milestone_surv_diff_lower
#> [1] NA
#> 
#> $milestone_surv_diff_upper
#> [1] NA
#> 
#> $CI_level
#> [1] 0.95
#> 
#> $times
#> [1] 3 4 5
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
```
