# Generate Dataset with crossing hazards

Generate Dataset with crossing hazards

Create an empty assumtions data.frame for generate_crossing_hazards

Calculate hr after crossing the hazard functions

Calculate true summary statistics for scenarios with crossing hazards

## Usage

``` r
generate_crossing_hazards(condition, fixed_objects = NULL)

assumptions_crossing_hazards(print = interactive())

hr_after_crossing_from_PH_effect_size(
  design,
  target_power_ph = NA_real_,
  final_events = NA_real_,
  target_alpha = 0.025
)

cen_rate_from_cen_prop_crossing_hazards(design)

true_summary_statistics_crossing_hazards(
  Design,
  cutoff_stats = NULL,
  milestones = NULL,
  fixed_objects = NULL
)
```

## Arguments

- condition:

  condition row of Design dataset

- fixed_objects:

  additional settings, see details

- print:

  print code to generate parameter set?

- design:

  design data.frame

- target_power_ph:

  target power under proportional hazards

- final_events:

  target events for inversion of Schönfeld Formula, defaults to
  `condition$final_events`

- target_alpha:

  target one-sided alpha level for the power calculation

- Design:

  Design data.frame for crossing hazards

- cutoff_stats:

  (optionally named) cutoff time, see details

- milestones:

  (optionally named) vector of times at which milestone survival should
  be calculated

## Value

For generate_crossing_hazards: A dataset with the columns t (time) and
trt (1=treatment, 0=control), evt (event, currently TRUE for all
observations)

For assumptions_crossing_hazards: a design tibble with default values
invisibly

For hr_after_crossing_from_PH_effect_size: the design data.frame passed
as argument with the additional column hazard_trt.

for cen_rate_from_cen_prop_crossing_hazards: design data.frame with the
additional column random_withdrawal

For true_summary_statistics_crossing_hazards: the design data.frame
passed as argument with additional columns,

## Details

Condidtion has to contain the following columns:

- n_trt number of paitents in treatment arm

- n_ctrl number of patients in control arm

- crossing time of crossing of the hazards

- hazard_ctrl hazard in the control arm = hazard before onset of
  treatment effect

- hazard_trt_before hazard in the treatment arm before onset of
  treatment effect

- hazard_trt_after hazard in the treatment arm afert onset of treatment
  effect

If fixed_objects is given and contains an element `t_max`, then this is
used as the cutoff for the simulation used internally. If t_max is not
given in this way the 1-(1/10000) quantile of the survival distribution
in the control or treatment arm is used (which ever is larger).

assumptions_crossing_hazards generates a default design `data.frame` for
use with generate_crossing_hazards If print is `TRUE` code to produce
the template is also printed for copying, pasting and editing by the
user. (This is the default when run in an interactive session.)

`hr_after_crossing_from_PH_effect_size` calculates the hazard ratio
after crossing of hazards as follows: First, the hazard ratio needed to
archive the desired power under proportional hazards is calculated by
inverting Schönfeld's sample size formula. Second the median survival
times for both arm under this hazard ratio and proportional hazards are
calculated. Finally the hazard rate of the treatment arm after crossing
of hazards is set such that the median survival time is the same as the
one calculated under proportional hazards.

This is a heuristic and to some extent arbitrary approach to calculate
hazard ratios that correspond to reasonable and realistic scenarios.

cen_rate_from_cen_prop_crossing_hazards takes the proportion of censored
patients from the column `censoring_prop`. This column describes the
proportion of patients who are censored randomly before experiencing an
event, without regard to administrative censoring.

`cutoff_stats` are the times used to calculate the statistics like
average hazard ratios and RMST, that are only calculated up to a certain
point.

## Functions

- `generate_crossing_hazards()`: simulates a dataset with crossing
  hazards

- `assumptions_crossing_hazards()`: generate default assumptions
  `data.frame`

- `hr_after_crossing_from_PH_effect_size()`: Calculate hr after crossing
  of the hazards from PH effect size

- `cen_rate_from_cen_prop_crossing_hazards()`: calculate censoring rate
  from censoring proportion

- `true_summary_statistics_crossing_hazards()`: calculate true summary
  statistics for crossing hazards

## Examples

``` r
one_simulation <- merge(
    assumptions_crossing_hazards(),
    design_fixed_followup(),
    by=NULL
  ) |>
  head(1) |>
  generate_crossing_hazards()
head(one_simulation)
#>      t trt  evt
#> 1  862   1 TRUE
#> 2  298   1 TRUE
#> 3   41   1 TRUE
#> 4 2273   1 TRUE
#> 5 1927   1 TRUE
#> 6 1050   1 TRUE
tail(one_simulation)
#>        t trt  evt
#> 295  259   0 TRUE
#> 296  245   0 TRUE
#> 297  788   0 TRUE
#> 298 1348   0 TRUE
#> 299 1905   0 TRUE
#> 300  878   0 TRUE
Design <- assumptions_crossing_hazards()
Design
#>   crossing  hazard_ctrl hazard_trt_before hazard_trt_after random_withdrawal
#> 1    0.000 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 2   60.875 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 3  121.750 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 4  182.625 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 5  243.500 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 6  304.375 0.0009488668       0.001265156     0.0006325779      0.0001897734
my_design <- merge(
    assumptions_crossing_hazards(),
    design_fixed_followup(),
    by=NULL
  )

my_design$final_events <- ceiling((my_design$n_trt + my_design$n_ctrl)*0.75)
my_design$hazard_trt <- NA
my_design <- hr_after_crossing_from_PH_effect_size(my_design, target_power_ph=0.9)
my_design
#>   crossing  hazard_ctrl hazard_trt_before hazard_trt_after random_withdrawal
#> 1    0.000 0.0009488668       0.001265156     0.0006158887      0.0001897734
#> 2   60.875 0.0009488668       0.001265156     0.0005787618      0.0001897734
#> 3  121.750 0.0009488668       0.001265156     0.0005371313      0.0001897734
#> 4  182.625 0.0009488668       0.001265156     0.0004901248      0.0001897734
#> 5  243.500 0.0009488668       0.001265156     0.0004366293      0.0001897734
#> 6  304.375 0.0009488668       0.001265156     0.0003752012      0.0001897734
#>   n_trt n_ctrl followup recruitment final_events hazard_trt target_median_trt
#> 1   150    150    730.5     182.625          225         NA          1125.442
#> 2   150    150    730.5     182.625          225         NA          1125.442
#> 3   150    150    730.5     182.625          225         NA          1125.442
#> 4   150    150    730.5     182.625          225         NA          1125.442
#> 5   150    150    730.5     182.625          225         NA          1125.442
#> 6   150    150    730.5     182.625          225         NA          1125.442
#>   target_hr
#> 1 0.6490782
#> 2 0.6490782
#> 3 0.6490782
#> 4 0.6490782
#> 5 0.6490782
#> 6 0.6490782
design <- data.frame(
  crossing = c(2, 4, 6),
  hazard_ctrl = c(0.05, 0.05, 0.05),
  hazard_trt_before = c(0.025, 0.025, 0.025),
  hazard_trt_after = c(0.1, 0.1, 0.1),
  censoring_prop = c(0.1, 0.3, 0.2),
  n_trt = c(50, 50, 50),
  n_ctrl = c(50, 50, 50),
  followup = c(200, 200, 200),
  recruitment = c(50, 50, 50)
)
cen_rate_from_cen_prop_crossing_hazards(design)
#>   crossing hazard_ctrl hazard_trt_before hazard_trt_after censoring_prop n_trt
#> 1        2        0.05             0.025              0.1            0.1    50
#> 2        4        0.05             0.025              0.1            0.3    50
#> 3        6        0.05             0.025              0.1            0.2    50
#>   n_ctrl followup recruitment random_withdrawal
#> 1     50      200          50       0.007469191
#> 2     50      200          50       0.029373101
#> 3     50      200          50       0.016885532
my_design <- merge(
    assumptions_crossing_hazards(),
    design_fixed_followup(),
    by=NULL
  )
my_design$follwup <- 15
my_design <- true_summary_statistics_crossing_hazards(my_design)
my_design
#>   crossing  hazard_ctrl hazard_trt_before hazard_trt_after random_withdrawal
#> 1    0.000 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 2   60.875 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 3  121.750 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 4  182.625 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 5  243.500 0.0009488668       0.001265156     0.0006325779      0.0001897734
#> 6  304.375 0.0009488668       0.001265156     0.0006325779      0.0001897734
#>   n_trt n_ctrl followup recruitment follwup median_survival_trt
#> 1   150    150    730.5     182.625      15            1095.750
#> 2   150    150    730.5     182.625      15            1034.875
#> 3   150    150    730.5     182.625      15             974.000
#> 4   150    150    730.5     182.625      15             913.125
#> 5   150    150    730.5     182.625      15             852.250
#> 6   150    150    730.5     182.625      15             791.375
#>   median_survival_ctrl
#> 1                730.5
#> 2                730.5
#> 3                730.5
#> 4                730.5
#> 5                730.5
#> 6                730.5
```
