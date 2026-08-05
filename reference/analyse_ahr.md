# Analyse the dataset using extimators for the the average hazard ratio

Analyse the dataset using extimators for the the average hazard ratio

## Usage

``` r
analyse_ahr(
  max_time = NA,
  type = "AHR",
  level = 0.95,
  alternative = "two.sided"
)
```

## Arguments

- max_time:

  time for which the AHR is calculated

- type:

  "AHR" for average hazard ratio "gAHR" for geometric average hazard
  ratio

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

- `AHR`/`gAHR` estimated (geometric) average hazard ratio

- `AHR_lower`/`gAHR_lower` unadjusted lower bound of the confidence
  interval for the (geometric) average hazard ratio

- `AHR_upper`/`gAHR_upper` unadjusted upper bound of the confidence
  interval for the (geometric) average hazard ratio

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
analyse_ahr()(condition, dat)
#> $p
#> [1] 0.0002724291
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $AHR
#> [1] 0.6076783
#> 
#> $AHR_lower
#> [1] 0.4647279
#> 
#> $AHR_upper
#> [1] 0.7946001
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
analyse_ahr(type = "gAHR")(condition, dat)
#> $p
#> [1] 1.171566e-05
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $gAHR
#> [1] 0.5910285
#> 
#> $gAHR_lower
#> [1] 0.4671682
#> 
#> $gAHR_upper
#> [1] 0.7477279
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
analyse_ahr(max_time = 50, type = "AHR")(condition, dat)
#> $p
#> [1] 0.8071216
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $AHR
#> [1] 0.8881195
#> 
#> $AHR_lower
#> [1] 0.34261
#> 
#> $AHR_upper
#> [1] 2.302199
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
analyse_ahr(max_time = 50, type = "gAHR")(condition, dat)
#> $p
#> [1] 0.8049128
#> 
#> $alternative
#> [1] "two.sided"
#> 
#> $gAHR
#> [1] 0.8870164
#> 
#> $gAHR_lower
#> [1] 0.3425746
#> 
#> $gAHR_upper
#> [1] 2.29672
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
