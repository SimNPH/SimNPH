# Analyse Dataset with the Maxcombo Test

Analyse Dataset with the Maxcombo Test

## Usage

``` r
analyse_maxcombo(alternative = "two.sided")
```

## Arguments

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

an analyse function that returns a data.frame with the combined p-value
of the max combo test in the column p

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
analyse_maxcombo()(condition, dat)
#> $p
#> [1] 9.762037e-06
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
```
