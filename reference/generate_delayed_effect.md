# Generate Dataset with delayed effect

Generate Dataset with delayed effect

Create an empty assumtions data.frame for generate_delayed_effect

Calculate hr after onset of treatment effect

Calculate true summary statistics for scenarios with delayed treatment
effect

## Usage

``` r
generate_delayed_effect(condition, fixed_objects = NULL)

assumptions_delayed_effect(print = interactive())

hr_after_onset_from_PH_effect_size(
  design,
  target_power_ph = NA_real_,
  final_events = NA_real_,
  target_alpha = 0.025
)

cen_rate_from_cen_prop_delayed_effect(design)

true_summary_statistics_delayed_effect(
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

  target events for inversion of Schönfeld Formula defaults to
  `condition$final_events`

- target_alpha:

  target one-sided alpha level for the power calculation

- Design:

  Design data.frame for delayed effect

- cutoff_stats:

  (optionally named) cutoff times, see details

- milestones:

  (optionally named) vector of times at which milestone survival should
  be calculated

## Value

For generate_delayed_effect: A dataset with the columns t (time) and trt
(1=treatment, 0=control), evt (event, currently TRUE for all
observations)

For assumptions_delayed_effect: a design tibble with default values
invisibly

For hr_after_onset_from_PH_effect_size: the design data.frame passed as
argument with the additional column hazard_trt.

for cen_rate_from_cen_prop_delayed_effect: design data.frame with the
additional column random_withdrawal

For true_summary_statistics_delayed_effect: the design data.frame passed
as argument with additional columns

## Details

Condidtion has to contain the following columns:

- n_trt number of paitents in treatment arm

- n_ctrl number of patients in control arm

- delay time until onset of effect

- hazard_ctrl hazard in the control arm = hazard before onset of
  treatment effect

- hazard_trt hazard in the treatment arm afert onset of treatment effect

If fixed_objects is given and contains an element `t_max`, then this is
used as the cutoff for the simulation used internally. If t_max is not
given in this way the 1-(1/10000) quantile of the survival distribution
in the control or treatment arm is used (which ever is larger).

assumptions_delayed_effect generates a default design `data.frame` for
use with generate_delayed_effect. If print is `TRUE` code to produce the
template is also printed for copying, pasting and editing by the user.
(This is the default when run in an interactive session.)

`hr_after_onset_from_PH_effect_size` calculates the hazard ratio after
onset of treatment effect as follows: First, the hazard ratio needed to
archive the desired power under proportional hazards is calculated by
inverting Schönfeld's sample size formula. Second the median survival
times for both arm under this hazard ratio and proportional hazards are
calculated. Finally the hazard rate of the treatment arm after onset of
treatment effect is set such that the median survival time is the same
as the one calculated under proportional hazards.

This is a heuristic and to some extent arbitrary approach to calculate
hazard ratios that correspond to reasonable and realistic scenarios.

cen_rate_from_cen_prop_delayed_effect takes the proportion of censored
patients from the column `censoring_prop`. This column describes the
proportion of patients who are censored randomly before experiencing an
event, without regard to administrative censoring.

`cutoff_stats` are the times used to calculate the statistics like
average hazard ratios and RMST, that are only calculated up to a certain
point.

## Functions

- `generate_delayed_effect()`: simulates a dataset with delayed
  treatment effect

- `assumptions_delayed_effect()`: generate default assumptions
  `data.frame`

- `hr_after_onset_from_PH_effect_size()`: Calculate hr after onset of
  treatment effect of the hazards from PH effect size

- `cen_rate_from_cen_prop_delayed_effect()`: calculate censoring rate
  from censoring proportion

- `true_summary_statistics_delayed_effect()`: calculate true summary
  statistics for delayed effect

## Examples

``` r
one_simulation <- merge(
    assumptions_delayed_effect(),
    design_fixed_followup(),
    by=NULL
  ) |>
  head(1) |>
  generate_delayed_effect()
head(one_simulation)
#>      t trt  evt
#> 1 5593   1 TRUE
#> 2 3411   1 TRUE
#> 3  955   1 TRUE
#> 4 1498   1 TRUE
#> 5 6771   1 TRUE
#> 6  421   1 TRUE
tail(one_simulation)
#>        t trt  evt
#> 295  158   0 TRUE
#> 296 1428   0 TRUE
#> 297  578   0 TRUE
#> 298  964   0 TRUE
#> 299 2495   0 TRUE
#> 300  572   0 TRUE
Design <- assumptions_delayed_effect()
Design
#>     delay  hazard_ctrl   hazard_trt random_withdrawal
#> 1   0.000 0.0009488668 0.0006325779      0.0001897734
#> 2  60.875 0.0009488668 0.0006325779      0.0001897734
#> 3 121.750 0.0009488668 0.0006325779      0.0001897734
#> 4 182.625 0.0009488668 0.0006325779      0.0001897734
#> 5 243.500 0.0009488668 0.0006325779      0.0001897734
#> 6 304.375 0.0009488668 0.0006325779      0.0001897734
my_design <- merge(
  assumptions_delayed_effect(),
  design_fixed_followup(),
  by=NULL
)

my_design$hazard_ctrl <- 0.05
my_design$final_events <- ceiling((my_design$n_trt + my_design$n_ctrl)*0.75)
my_design$hazard_trt <- NA
my_design <- hr_after_onset_from_PH_effect_size(my_design, target_power_ph=0.9)
#> Warning: Median survival is shorter than delay of treatment effect, calculation not possible
#> Warning: Median survival is shorter than delay of treatment effect, calculation not possible
#> Warning: Median survival is shorter than delay of treatment effect, calculation not possible
#> Warning: Median survival is shorter than delay of treatment effect, calculation not possible
#> Warning: Median survival is shorter than delay of treatment effect, calculation not possible
my_design
#>     delay hazard_ctrl hazard_trt random_withdrawal n_trt n_ctrl followup
#> 1   0.000        0.05 0.03245391      0.0001897734   150    150    730.5
#> 2  60.875        0.05         NA      0.0001897734   150    150    730.5
#> 3 121.750        0.05         NA      0.0001897734   150    150    730.5
#> 4 182.625        0.05         NA      0.0001897734   150    150    730.5
#> 5 243.500        0.05         NA      0.0001897734   150    150    730.5
#> 6 304.375        0.05         NA      0.0001897734   150    150    730.5
#>   recruitment final_events target_median_trt target_hr
#> 1     182.625          225          21.35789 0.6490782
#> 2     182.625          225          21.35789 0.6490782
#> 3     182.625          225          21.35789 0.6490782
#> 4     182.625          225          21.35789 0.6490782
#> 5     182.625          225          21.35789 0.6490782
#> 6     182.625          225          21.35789 0.6490782
design <- expand.grid(
  delay=seq(0, 10, by=5),            # delay of 0, 1, ..., 10 days
  hazard_ctrl=0.2,                   # hazard under control and before treatment effect
  hazard_trt=0.02,                   # hazard after onset of treatment effect
  censoring_prop=c(0.1, 0.25, 0.01), # 10%, 25%, 1% random censoring
  followup=100,                      # followup of 100 days
  n_trt=50,                          # 50 patients treatment
  n_ctrl=50                          # 50 patients control
)
cen_rate_from_cen_prop_delayed_effect(design)
#>   delay hazard_ctrl hazard_trt censoring_prop followup n_trt n_ctrl
#> 1     0         0.2       0.02           0.10      100    50     50
#> 2     5         0.2       0.02           0.10      100    50     50
#> 3    10         0.2       0.02           0.10      100    50     50
#> 4     0         0.2       0.02           0.25      100    50     50
#> 5     5         0.2       0.02           0.25      100    50     50
#> 6    10         0.2       0.02           0.25      100    50     50
#> 7     0         0.2       0.02           0.01      100    50     50
#> 8     5         0.2       0.02           0.01      100    50     50
#> 9    10         0.2       0.02           0.01      100    50     50
#>   random_withdrawal
#> 1      0.0043517713
#> 2      0.0047198579
#> 3      0.0050796441
#> 4      0.0150805823
#> 5      0.0162247001
#> 6      0.0173312396
#> 7      0.0003698016
#> 8      0.0004022565
#> 9      0.0004341304
my_design <- merge(
    assumptions_delayed_effect(),
    design_fixed_followup(),
    by=NULL
  )
my_design <- true_summary_statistics_delayed_effect(my_design)
my_design
#>     delay  hazard_ctrl   hazard_trt random_withdrawal n_trt n_ctrl followup
#> 1   0.000 0.0009488668 0.0006325779      0.0001897734   150    150    730.5
#> 2  60.875 0.0009488668 0.0006325779      0.0001897734   150    150    730.5
#> 3 121.750 0.0009488668 0.0006325779      0.0001897734   150    150    730.5
#> 4 182.625 0.0009488668 0.0006325779      0.0001897734   150    150    730.5
#> 5 243.500 0.0009488668 0.0006325779      0.0001897734   150    150    730.5
#> 6 304.375 0.0009488668 0.0006325779      0.0001897734   150    150    730.5
#>   recruitment median_survival_trt median_survival_ctrl
#> 1     182.625           1095.7500                730.5
#> 2     182.625           1065.3125                730.5
#> 3     182.625           1034.8750                730.5
#> 4     182.625           1004.4375                730.5
#> 5     182.625            974.0000                730.5
#> 6     182.625            943.5625                730.5
```
