# Create a summarise function from a named list of functions

Create a summarise function from a named list of functions

## Usage

``` r
create_summarise_function(...)
```

## Arguments

- ...:

  summarise function

## Value

a function with arguments condition, results, fixed objects

## Details

the names of the list of functions correspond to the names in the list
of analyse functions, each summarise function is applied to the results
of the analyse function of the same name, names not present in both
lists are ommitted in either list.

The functions in the list should have the arguments `condition`,
`results` and `fixed_objects`. `results` is a list of lists. The outer
list has one element for each replication, the inner list has one entry
for each Analyse function. (Analyse functions have to return lists for
this to work, otherwise the results are simplified to data.frames.
Analyse functions from the SimNPH package all return lists.)

The individual summarise functions have to return data.frames, which are
concatendated column-wise to give one row per condition. The names of
the analyse methods are prepended to the respective coumn names, if the
functions have a "name" attribute this is appended to the column names
of the output. Column names not unique after that are appended numbers
by `make.unique`.

## Examples

``` r
Summarise <- create_summarise_function(
  maxcombo = function(condition, results, fixed_objects=NULL){
    data.frame("rejection"=mean(results$p < alpha))
  },
  logrank  = function(condition, results, fixed_objects=NULL){
    data.frame("rejection"=mean(results$p < alpha))
  }
)
```
