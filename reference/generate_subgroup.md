# Generate Dataset with different treatment effect in subgroup

Generate Dataset with different treatment effect in subgroup

Create an empty assumtions data.frame for generate_subgroup

Calculate true summary statistics for scenarios with differential
treatment effect in subgroup

Calculate hazards in treatment arm in subgroup and compliment

## Usage

``` r
generate_subgroup(condition, fixed_objects = NULL)

assumptions_subgroup(print = interactive())

true_summary_statistics_subgroup(
  Design,
  cutoff_stats = NULL,
  milestones = NULL,
  fixed_objects = NULL
)

hazard_subgroup_from_PH_effect_size(
  design,
  target_power_ph = NA_real_,
  final_events = NA_real_,
  target_alpha = 0.025
)

cen_rate_from_cen_prop_subgroup(design)
```

## Arguments

- condition:

  condition row of Design dataset

- fixed_objects:

  additional settings, see details

- print:

  print code to generate parameter set?

- Design:

  Design data.frame for subgroup

- cutoff_stats:

  (optionally named) cutoff times, see details

- milestones:

  (optionally named) vector of times at which milestone survival should
  be calculated

- design:

  design data.frame

- target_power_ph:

  target power under proportional hazards

- final_events:

  target events for inversion of Schönfeld Formula, defaults to
  `condition$final_events`

- target_alpha:

  target one-sided alpha level for the power calculation

## Value

For generate_subgroup: A dataset with the columns t (time) and trt
(1=treatment, 0=control), evt (event, currently TRUE for all
observations)

For assumptions_subgroup: a design tibble with default values invisibly

For true_summary_statistics_subgroup: the design data.frame passed as
argument with the additional columns

For hazard_subgroup_from_PH_effect_size: the design data.frame passed as
argument with the additional columns hazard_trt and hazard_subgroup.

for cen_rate_from_cen_prop_subgroup: design data.frame with the
additional column random_withdrawal

## Details

Condidtion has to contain the following columns:

- n_trt number of paitents in treatment arm

- n_ctrl number of patients in control arm

- hazard_ctrl hazard in the control arm

- hazard_trt hazard in the treatment arm for not cured patients

- hazard_subgroup hazard in the subgroup in the treatment arm

- prevalence proportion of cured patients

assumptions_subgroup generates a default design `data.frame` for use
with generate_subgroup If print is `TRUE` code to produce the template
is also printed for copying, pasting and editing by the user. (This is
the default when run in an interactive session.)

`cutoff_stats` are the times used to calculate the statistics like
average hazard ratios and RMST, that are only calculated up to a certain
point.

`hazard_subgroup_from_PH_effect_size` calculates the hazard rate in the
subgroup and the compliment of the subgroup in the treatment arm as
follows: First, the hazard ratio needed to archive the desired power
under proportional hazards is calculated by inverting Schönfeld's sample
size formula. Second the median survival times for both arms under this
hazard ratio and proportional hazards are calculated. Finally the hazard
rate of the treatment arm in the subgroup and its complement are set
such that the median survival time is the same as the one calculated
under proportional hazards.

This is a heuristic and to some extent arbitrary approach to calculate
hazard ratios that correspond to reasonable and realistic scenarios.

cen_rate_from_cen_prop_subgroup takes the proportion of censored
patients from the column `censoring_prop`. This column describes the
proportion of patients who are censored randomly before experiencing an
event, without regard to administrative censoring.

## Functions

- `generate_subgroup()`: simulates a dataset with a mixture of cured
  patients

- `assumptions_subgroup()`: generate default assumptions `data.frame`

- `true_summary_statistics_subgroup()`: calculate true summary
  statistics for subgroup

- `hazard_subgroup_from_PH_effect_size()`: Calculate hazards in
  treatement arm

- `cen_rate_from_cen_prop_subgroup()`: calculate censoring rate from
  censoring proportion

## Examples

``` r
one_simulation <- merge(
    assumptions_subgroup(),
    design_fixed_followup(),
    by=NULL
  ) |>
  head(1) |>
  generate_subgroup()
head(one_simulation)
#>      t trt  evt subgroup
#> 1  428   1 TRUE        0
#> 2  132   1 TRUE        0
#> 3 3087   1 TRUE        0
#> 4  734   1 TRUE        0
#> 5  578   1 TRUE        0
#> 6  192   1 TRUE        0
tail(one_simulation)
#>        t trt  evt subgroup
#> 295 6391   0 TRUE        0
#> 296 2128   0 TRUE        0
#> 297  762   0 TRUE        0
#> 298 2131   0 TRUE        1
#> 299 1323   0 TRUE        0
#> 300  840   0 TRUE        0
Design <- assumptions_subgroup()
Design
#>    hazard_ctrl   hazard_trt hazard_subgroup prevalence random_withdrawal
#> 1 0.0009488668 0.0006325779    9.488668e-05        0.2      0.0001897734
#> 2 0.0009488668 0.0006325779    9.488668e-05        0.4      0.0001897734
#> 3 0.0009488668 0.0006325779    9.488668e-05        0.6      0.0001897734
#> 4 0.0009488668 0.0006325779    9.488668e-05        0.8      0.0001897734
my_design <- merge(
    assumptions_subgroup(),
    design_fixed_followup(),
    by=NULL
  )
my_design <- true_summary_statistics_subgroup(my_design)
my_design
#>    hazard_ctrl   hazard_trt hazard_subgroup prevalence random_withdrawal n_trt
#> 1 0.0009488668 0.0006325779    9.488668e-05        0.2      0.0001897734   150
#> 2 0.0009488668 0.0006325779    9.488668e-05        0.4      0.0001897734   150
#> 3 0.0009488668 0.0006325779    9.488668e-05        0.6      0.0001897734   150
#> 4 0.0009488668 0.0006325779    9.488668e-05        0.8      0.0001897734   150
#>   n_ctrl followup recruitment median_survival_trt median_survival_ctrl
#> 1    150    730.5     182.625            1422.745                730.5
#> 2    150    730.5     182.625            2001.270                730.5
#> 3    150    730.5     182.625            3143.728                730.5
#> 4    150    730.5     182.625            5119.929                730.5

my_design <- merge(
  assumptions_subgroup(),
  design_fixed_followup(),
  by=NULL
)

my_design$hazard_trt <- NA
my_design$hazard_subgroup <- NA
my_design$hr_subgroup_relative <- 0.9
my_design$final_events <- ceiling((my_design$n_ctrl + my_design$n_trt) * 0.75)
my_design <- hazard_subgroup_from_PH_effect_size(my_design, target_power_ph=0.9)
my_design
#>    hazard_ctrl   hazard_trt hazard_subgroup prevalence random_withdrawal n_trt
#> 1 0.0009488668 0.0006288263    0.0005659437        0.2      0.0001897734   150
#> 2 0.0009488668 0.0006421335    0.0005779202        0.4      0.0001897734   150
#> 3 0.0009488668 0.0006558156    0.0005902340        0.6      0.0001897734   150
#> 4 0.0009488668 0.0006698768    0.0006028891        0.8      0.0001897734   150
#>   n_ctrl followup recruitment hr_subgroup_relative final_events
#> 1    150    730.5     182.625                  0.9          225
#> 2    150    730.5     182.625                  0.9          225
#> 3    150    730.5     182.625                  0.9          225
#> 4    150    730.5     182.625                  0.9          225
#>   target_median_trt
#> 1          1125.442
#> 2          1125.442
#> 3          1125.442
#> 4          1125.442
design <- expand.grid(
  hazard_ctrl=0.2,                   # hazard under control and before treatment effect
  hazard_trt=0.02,                   # hazard after onset of treatment effect
  hazard_subgroup=0.01,              # hazard in the subgroup in treatment
  prevalence = c(0.2, 0.5),           # subgroup prevalence
  censoring_prop=c(0.1, 0.25, 0.01), # 10%, 25%, 1% random censoring
  followup=100,                      # followup of 100 days
  n_trt=50,                          # 50 patients treatment
  n_ctrl=50                          # 50 patients control
)
cen_rate_from_cen_prop_subgroup(design)
#>   hazard_ctrl hazard_trt hazard_subgroup prevalence censoring_prop followup
#> 1         0.2       0.02            0.01        0.2           0.10      100
#> 2         0.2       0.02            0.01        0.5           0.10      100
#> 3         0.2       0.02            0.01        0.2           0.25      100
#> 4         0.2       0.02            0.01        0.5           0.25      100
#> 5         0.2       0.02            0.01        0.2           0.01      100
#> 6         0.2       0.02            0.01        0.5           0.01      100
#>   n_trt n_ctrl random_withdrawal
#> 1    50     50      0.0030440737
#> 2    50     50      0.0026408014
#> 3    50     50      0.0108976367
#> 4    50     50      0.0095627716
#> 5    50     50      0.0002558513
#> 6    50     50      0.0002211651
```
