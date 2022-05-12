#!/usr/bin/env R

# Load libraries
library(tidyverse)

c(
  "src/requirements/cran.txt",
  "src/requirements/github.txt",
  "src/requirements/bioconductor.txt"
) %>%
  purrr::map(readr::read_lines) %>%
  purrr::simplify() %>%
  purrr::discard(~ .x %>% stringr:::str_detect("^#|^$")) %>%
  # remove user name and slash from github packages
  purrr::map_chr(~ .x %>% str_extract("[A-z]+$")) %>%
  purrr::walk(require, character.only = TRUE)

# load tidyverse again to use as default namespace
library(tidyverse)

# log session
dir.create("log", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "log/R_sessionInfo.txt")

# set working directory to the current directory
setwd(here::here())

# Source all R files inside the src directory
# list.files(path = "src", pattern = "\\.R$", full.names = TRUE, recursive = TRUE) %>%
#   purrr::discard(~ .x %>% str_detect("^src/Rprofile")) %>%
#   setdiff(c("src/defaults.R", "src/install.R", "src/main.R")) %>%
#   purrr::walk(source)

# Global options
options(
  mc.cores = parallel::detectCores()
)

# Default ggplot2 theme
theme_my <-
  ggplot2::theme_minimal(base_size = 20) +
  ggplot2::theme(
    axis.line.x = ggplot2::element_line(size = 0.8),
    axis.line.y = ggplot2::element_line(size = 0.8),
    axis.ticks = ggplot2::element_line(colour = "black", size = 0.8),
    axis.text = ggplot2::element_text(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )

ggplot2::theme_set(theme_my)
