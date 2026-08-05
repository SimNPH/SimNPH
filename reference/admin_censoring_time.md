# Add recruitment time to Dataset

Add recruitment time to Dataset

Apply Administrative Censoring After Fixed Time

Apply Administrative Censoring After Fixed Number of Events

## Usage

``` r
recruitment_uniform(
  dat,
  recruitment_until,
  recruitment_from = 0,
  discrete = TRUE
)

admin_censoring_time(dat, followup, keep_non_recruited = FALSE)

admin_censoring_events(
  dat,
  events,
  keep_non_recruited = FALSE,
  on_incomplete = "ignore"
)
```

## Arguments

- dat:

  a simulated dataset

- recruitment_until:

  time of end of recruitment

- recruitment_from:

  time of start of recruitment (defaults to 0)

- discrete:

  should the recruitment time be rounded to full days?

- followup:

  followup time

- keep_non_recruited:

  should patients recruited after end of study be kept

- events:

  number of events after which the dataset is analyzed

- on_incomplete:

  what to do if there are fewer events than planned
  "ignore","warn","stop"

## Value

Returns the dataset with added recruitment times.

Returns the dataset with administrative censoring after `follwup`, adds
the attribute `followup` with the followup time to the dataset.

Returns the dataset with administrative censoring after `events` events,
adds the attribute `followup` with the followup time to the dataset.

## Details

The Dataset hast to include a column `rec_time` containing the
recruitment time as well as the columns with the event times `t` and a
column with the event indicator `evt`.

Times and event indicaotrs for patients recruited after followup are set
to `NA`.

The Dataset hast to include a column `rec_time` containing the
recruitment time as well as the columns with the event times `t` and a
column with the event indicator `evt`.

Times and event indicaotrs for patients recruited after followup are set
to `NA`.

If there are less events than planned for study end `on_incomplete`
defines what should be done. "ignore" simply returns the dataset with
the maximum of the observed times as followup. "warn" does the same but
gives a warning. "stop" stopps with an error.

## Functions

- `recruitment_uniform()`: add recruitment time

- `admin_censoring_time()`: apply administrative censoring after fixed
  time

- `admin_censoring_events()`: apply administrative censoring after fixed
  number of events

## Examples

``` r
dat <- data.frame(t=c(0, 1, 2), trt=c(FALSE, FALSE, TRUE))
recruitment_uniform(dat, 7, 0)
#>   t   trt rec_time
#> 1 0 FALSE        1
#> 2 1 FALSE        6
#> 3 2  TRUE        5
dat <- data.frame(
  t = 1:10,
  rec_time = rep(1:5, each=2),
  trt = rep(c(TRUE, FALSE), times=5),
  evt = rep(TRUE, times=10)
)
dat
#>     t rec_time   trt  evt
#> 1   1        1  TRUE TRUE
#> 2   2        1 FALSE TRUE
#> 3   3        2  TRUE TRUE
#> 4   4        2 FALSE TRUE
#> 5   5        3  TRUE TRUE
#> 6   6        3 FALSE TRUE
#> 7   7        4  TRUE TRUE
#> 8   8        4 FALSE TRUE
#> 9   9        5  TRUE TRUE
#> 10 10        5 FALSE TRUE

admin_censoring_time(dat, 4)
#>   t rec_time   trt   evt
#> 1 1        1  TRUE  TRUE
#> 2 2        1 FALSE  TRUE
#> 3 2        2  TRUE FALSE
#> 4 2        2 FALSE FALSE
#> 5 1        3  TRUE FALSE
#> 6 1        3 FALSE FALSE
#> 7 0        4  TRUE FALSE
#> 8 0        4 FALSE FALSE
admin_censoring_time(dat, 4, keep_non_recruited = TRUE)
#>     t rec_time   trt   evt
#> 1   1        1  TRUE  TRUE
#> 2   2        1 FALSE  TRUE
#> 3   2        2  TRUE FALSE
#> 4   2        2 FALSE FALSE
#> 5   1        3  TRUE FALSE
#> 6   1        3 FALSE FALSE
#> 7   0        4  TRUE FALSE
#> 8   0        4 FALSE FALSE
#> 9  NA        5  TRUE    NA
#> 10 NA        5 FALSE    NA

dat_censored <- admin_censoring_time(dat, 5)
attr(dat_censored, "followup")
#> [1] 5
dat <- data.frame(
  t = 1:10,
  rec_time = rep(2*(1:5), each=2),
  trt = rep(c(TRUE, FALSE), times=5),
  evt = rep(TRUE, times=10)
)
dat
#>     t rec_time   trt  evt
#> 1   1        2  TRUE TRUE
#> 2   2        2 FALSE TRUE
#> 3   3        4  TRUE TRUE
#> 4   4        4 FALSE TRUE
#> 5   5        6  TRUE TRUE
#> 6   6        6 FALSE TRUE
#> 7   7        8  TRUE TRUE
#> 8   8        8 FALSE TRUE
#> 9   9       10  TRUE TRUE
#> 10 10       10 FALSE TRUE

admin_censoring_events(dat, 4)
#>   t rec_time   trt   evt
#> 1 1        2  TRUE  TRUE
#> 2 2        2 FALSE  TRUE
#> 3 3        4  TRUE  TRUE
#> 4 4        4 FALSE  TRUE
#> 5 2        6  TRUE FALSE
#> 6 2        6 FALSE FALSE
#> 7 0        8  TRUE FALSE
#> 8 0        8 FALSE FALSE
admin_censoring_events(dat, 4, keep_non_recruited = TRUE)
#>     t rec_time   trt   evt
#> 1   1        2  TRUE  TRUE
#> 2   2        2 FALSE  TRUE
#> 3   3        4  TRUE  TRUE
#> 4   4        4 FALSE  TRUE
#> 5   2        6  TRUE FALSE
#> 6   2        6 FALSE FALSE
#> 7   0        8  TRUE FALSE
#> 8   0        8 FALSE FALSE
#> 9  NA       10  TRUE    NA
#> 10 NA       10 FALSE    NA

dat_censored <- admin_censoring_events(dat, 4)
attr(dat_censored, "followup")
#> [1] 8
```
