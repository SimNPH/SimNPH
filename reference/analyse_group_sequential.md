# Create Analyse Functions for Group Sequential Design

Create Analyse Functions for Group Sequential Design

Summarise Output from Analyse Functions for Group Sequential Design

## Usage

``` r
analyse_group_sequential(followup, followup_type, alpha, analyse_functions)

summarise_group_sequential(name = NULL)
```

## Arguments

- followup:

  followup events or time

- followup_type:

  "events" or "time"

- alpha:

  nominal alpha at each stage

- analyse_functions:

  analyse function or list of analyse functions

- name:

  name attribute of the returned closure

## Value

an analyse function that can be used in runSimulation

Returns a function with the arguments:

- condition

- results

- fixed objects

that can be passed to create_summarise_function or to
SimDesign::runSimulation and that returns a `data.frame`.

## Details

`followup`, `followup_type` and `alpha` are evaluated for every
simulated dataset, i.e. the arguments to the Analyse function are
available, expressions like
`followup=c(condition$interim, condition$max_followup)` are valid
arguments.

analyse_functions should take arguments condition, dataset and
fixed_objects and return a list conatining p-value, number of patients
and number of event in the columsn `p`, `N_pat` and `N_evt`.

## Functions

- `summarise_group_sequential()`: Summarise Output from Analyse
  Functions for Group Sequential Design

## Examples

``` r
# create a function to analyse after interim_events and maximum followup time
# given in the condition row of the design data.frame with given
# nominal alpha
analyse_maxcombo_sequential <- analyse_group_sequential(
  followup = c(condition$interim_events, condition$followup),
  followup_type = c("event", "time"),
  alpha = c(0.025, 0.05),
  analyse_functions = analyse_maxcombo()
)
Summarise <- create_summarise_function(
  maxcombo_seq = summarise_group_sequential(),
  logrank_seq = summarise_group_sequential(name="logrank")
)
```
