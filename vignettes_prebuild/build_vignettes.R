

rmarkdown::render(
  "vignettes_prebuild/simple_example.Rmd",
  output_file = "simple_example.html",
  output_dir = "vignettes",
  knit_root_dir = here::here()
  )

knitr::purl(
  "vignettes_prebuild/simple_example.Rmd",
  output = "vignettes/simple_example.R"
)

file.copy(
  "vignettes_prebuild/simple_example.html.asis",
  "vignettes/",
  overwrite = TRUE
  )

