
folder_name <- "sims_2"
git_tag <- "sims_2"

zip_name <- paste0(folder_name, ".zip")
package_url <- paste0("https://github.com/SimNPH/SimNPH/archive/refs/tags/", git_tag, ".zip")

download.file(package_url, destfile = zip_name)
dir.create(folder_name)
unzip(zip_name, exdir = folder_name)
setwd(paste0(folder_name, "/", git_tag))
devtools::install(".", build=TRUE, quick=TRUE, dependencies=TRUE)
