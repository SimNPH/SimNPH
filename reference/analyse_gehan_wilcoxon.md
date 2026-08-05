# Create Analyse function for Gehan Wilcoxon test

Create Analyse function for Gehan Wilcoxon test

## Usage

``` r
analyse_gehan_wilcoxon(alternative = "two.sided")
```

## Arguments

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

an analyse function that can be used in runSimulation

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
analyse_gehan_wilcoxon()(condition, dat)
#> $p
#> [1] 0.009155326
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
