#!/usr/bin/env R

#
# .Rprofile for interaktive analysis with make env
#

message("Calling .Rprofile for interactive analysis")
source("/project/analysis/src/Rprofile/all.R")

# needed for future to source external site libraries
.libPaths(c("/analysis/R-site-library/", .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths()[1], Sys.getenv("R_LIBS"), sep = .Platform$path.sep))

# needed for RStudio server: otherwise overwritten
options(repos = Sys.getenv("CRAN"))

# Auto start project
setHook("rstudio.sessionInit", function(newSession) {
  if (newSession && is.null(rstudioapi::getActiveProject()))
    rstudioapi::openProject(".")
}, action = "append")
