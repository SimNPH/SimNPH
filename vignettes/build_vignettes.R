# if vignettes have images: move images to vignettes folder
knitr::knit(
  "vignettes/simple_example.Rmd.orig",
  output = "vignettes/simple_example.Rmd"
  )
