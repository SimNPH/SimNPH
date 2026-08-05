# Analyse Dataset with the Logrank Test

Analyse Dataset with the Logrank Test

## Usage

``` r
analyse_logrank(alternative = "two.sided")
```

## Arguments

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

an analysis function that returns a data.frame with the columns

- `p` p-value of the logrank test

- `alternative` the alternative used

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
analyse_logrank()(condition, dat)
#> $p
#> [1] 4.028494e-07
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
```
