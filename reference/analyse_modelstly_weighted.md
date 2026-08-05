# Create Analyse function for the modestly weighted logrank test

Create Analyse function for the modestly weighted logrank test

## Usage

``` r
analyse_modelstly_weighted(t_star)
```

## Arguments

- t_star:

  parameter t\* of the modestly weighted logrank test

## Value

an analyse function that can be used in runSimulation

## Examples

``` r
condition <- merge(
   assumptions_delayed_effect(),
   design_fixed_followup(),
   by=NULL
 ) |>
 head(1)
dat <- generate_delayed_effect(condition)
analyse_modelstly_weighted(20)(condition, dat)
#> $p
#> [1] 1.205258e-05
#> 
#> $N_pat
#> [1] 300
#> 
#> $N_evt
#> [1] 300
#> 
```
