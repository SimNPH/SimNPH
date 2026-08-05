# Functions to Convert Between Days and Months and Medians and Rates

Some functions to convert between days and months and rates and medians.

## Usage

``` r
r2m(lambda)

m2r(med)

m2d(mon)

d2m(day)
```

## Arguments

- lambda:

  hazard rate

- med:

  median in months

- mon:

  time in months

- day:

  time in days

## Value

median survival time in months (`r2m`)

hazard rate per day (`m2r`)

time in days (`m2d`)

time in months (`d2m`)

## Functions

- `r2m()`: daily rate to median in months

- `m2r()`: median to months to daily rate

- `m2d()`: months to days

- `d2m()`: days to months

## Examples

``` r
r2m(0.002)
#> [1] 11.3864
m2r(12)
#> [1] 0.001897734
m2d(1)
#> [1] 30.4375
d2m(31)
#> [1] 1.01848
```
