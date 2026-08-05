# Generic Summarise function for esitmators

Generic Summarise function for esitmators

## Usage

``` r
summarise_estimator(
  est,
  real,
  lower = NULL,
  upper = NULL,
  null = NULL,
  est_sd = NULL,
  name = NULL
)
```

## Arguments

- est:

  estimator, expression evaluated in results

- real:

  real summary statistic, expression evaluated in condition

- lower:

  lower CI, expression evaluated in results

- upper:

  upper CI, expression evaluated in results

- null:

  parameter value under the null hypothesis

- est_sd:

  standard deviation estimated by the method, evaluated in results

- name:

  name for the summarise function, appended to the name of the analysis
  method in the final results

## Value

A function that can be used in Summarise that returns a data frame with
summary statistics of the performance measures in the columns.

## Details

The different parameters are evaluated in different envionments, `est`,
`lower`, `upper`, `est_sd` refer to output of the method and are
evaluated in the results dataset. `real` refers to a real value of a
summary statistic in this scenario and is therefore evaluated in the
condition dataset. `null` and `name` are constants and directly
evaluated when the function is defined. The argument `null`, the
parameter value under the null hypothesis is used to output the
rejection rate based on the confidence intervall. Which is output in the
column `null_cover`

## Examples

``` r
# \donttest{
# generate the design matrix and append the true summary statistics
condition <- merge(
  assumptions_delayed_effect(),
  design_fixed_followup(),
  by=NULL
) |>
  tail(4) |>
  head(1) |>
  true_summary_statistics_delayed_effect(cutoff_stats = 15)

# create some summarise functions
summarise_all <- create_summarise_function(
  coxph=summarise_estimator(hr, gAHR_15, hr_lower, hr_upper, name="gAHR"),
  coxph=summarise_estimator(hr, hazard_trt/hazard_ctrl, hr_lower, hr_upper, name="HR"),
  coxph=summarise_estimator(hr, NA_real_, name="NA")
)

# runs simulations
sim_results <- runSimulation(
  design=condition,
  replications=10,
  generate=generate_delayed_effect,
  analyse=list(
    coxph=analyse_coxph()
  ),
  summarise = summarise_all
)
#> 
#> Replications: 10;   RAM Used: 213 Mb;   
#>  Conditions: delay=121., hzrd_c=0.0009, hzrd_t=0.0006, rndm_w=0.0001, n_trt=150, n_ctrl=150, follwp=730., rcrtmn=182., mdn_srvvl_t=1034, mdn_srvvl_c=730., rmst_t_15=14.8, rmst_c_15=14.8, gAHR_1=1, AHR_15=1, AHRc_15=1, AHR__1=1
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%

# mse is missing for the summarise function in which the real value was NA
sim_results[, names(sim_results) |> grepl(pattern="\\.mse$")]
#> # A tibble: 1 × 3
#>   coxph.gAHR.mse coxph.HR.mse coxph.NA.mse
#>            <dbl>        <dbl>        <dbl>
#> 1        0.10149     0.025430          NaN
# but the standard deviation can be estimated in all cases
sim_results[, names(sim_results) |> grepl(pattern="\\.sd_est$")]
#> # A tibble: 1 × 3
#>   coxph.gAHR.sd_est coxph.HR.sd_est coxph.NA.sd_est
#>               <dbl>           <dbl>           <dbl>
#> 1           0.15869         0.15869         0.15869
# }
```
