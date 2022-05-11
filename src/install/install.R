#!/usr/bin/env Rscript

options(
  Ncpus = parallel::detectCores() / 2,
  Ncpu  = parallel::detectCores() / 2
)

library(tidyverse)
library(devtools)

get_requirements <- function(path) {
	path %>%
		readr::read_lines() %>%
		purrr::map_chr(~ .x %>% stringr:::str_remove("#.*$|^$") %>% stringr::str_trim()) %>%
  		purrr::discard(~ .x %>% stringr::str_detect("^$"))
}

options(repos = c(CRAN = "https://packagemanager.rstudio.com/all/__linux__/focal/2021-02-25+Y3JhbiwyOjc4NDM0OTs4NzIzN0YzQg"), download.file.method = 'libcurl')
options(BioC_mirror = 'https://packagemanager.rstudio.com/bioconductor')


# resolve special dependencies
# these packages need to be installed first, otherwise installing other packages crashes

# install from CRAN
install.packages("Rfast", type="source")

# install from BiocConductor & CRAN
print(get_requirements("src/install/packages.txt"))
pkgs <- get_requirements("src/install/packages.txt")
BiocManager::install(pkgs, update = FALSE, Ncpus = getOption("Ncpus", 1L))


