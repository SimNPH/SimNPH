# Add ggplot axis labels from labels attribute

Add ggplot axis labels from labels attribute

## Usage

``` r
labs_from_labels(gg)
```

## Arguments

- gg:

  a ggplot object

## Value

a ggplot object

## Examples

``` r
# \donttest{
library("ggplot2")
test <- mtcars
# add a label attribute
attr(test$cyl, "label") <- "cylinders"

# plot witht the variable names as axis titles
gg1 <- ggplot(test, aes(x=wt, y=cyl)) +
  geom_point()
gg1


# add labels where defined in the attribute
gg2 <- ggplot(test, aes(x=wt, y=cyl)) +
  geom_point()

gg2 <- labs_from_labels(gg2)
gg2

# }
```
