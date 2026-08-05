# Analyse Dataset with the Fleming Harrington weighted Logrank Test

Analyse Dataset with the Fleming Harrington weighted Logrank Test

## Usage

``` r
analyse_logrank_fh_weights(rho, gamma, alternative = "two.sided")
```

## Arguments

- rho:

  rho for the rho-gamma family of weights

- gamma:

  gamma for the rho-gamma family of weights

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

a function with the arguments condition, dat and fixed_objects that
returns a dataframe with the p-value of the weighted logrank test in the
column p. See ?SimDesign::Analyse for details on the arguments
condition, dat, fixed_arguments.

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
# create two functions with different weights
analyse_01 <- analyse_logrank_fh_weights(rho = 0, gamma = 1)
analyse_10 <- analyse_logrank_fh_weights(rho = 1, gamma = 0)
# run the tests created before
analyse_01(condition, dat)
#> $p
#> [1] 5.437706e-06
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
analyse_10(condition, dat)
#> $p
#> [1] 5.61106e-05
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
