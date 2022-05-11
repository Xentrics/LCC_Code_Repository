#!/usr/bin/env R

#
# .Rprofile for both run and env
#
# Do not load libraries here. Only basic configurations
#

# Global options
options(
  mc.cores = parallel::detectCores()
)

# force R to use analysis directory as home
Sys.setenv(HOME = "/project/analysis")

# log session
dir.create("log", showWarnings = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), "log/R_sessionInfo.txt")

# set working directory to the current directory
setwd(here::here())

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
