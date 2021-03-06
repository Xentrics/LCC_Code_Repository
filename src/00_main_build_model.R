#!/usr/bin/env Rscript

# In Rstudio, the can find 'Run' on the top-right border of this text editor!
# Alternatively, you can type 'ctrl+enter' to execute a line

# Load R packages
library("tidyverse")
library("magrittr")
library("drake")

#== SIAMCAT - Augmentation ====

###== Minimal workflow ====

# this code section will train 2 models for each, Candida genus and Candida glabrata
# tested on a 12 core machine with 64GB of RAM, but only a few GB are required
# this should finish within 5-30 minutes

# visualize workflow
g <- drake::r_drake_ggraph(source = "src/01_augment_train_siamcat_MINIMAL.HL.HU2.make.R")
g$data %<>%
  mutate(label = if_else(type == "object", label, ""))
g + ggraph::geom_node_text(aes(label = label), angle = 45, size = 2)

# compute workflow
drake::r_make(source = "src/01_augment_train_siamcat_MINIMAL.HL.HU2.make.R")


#
###== Validate models in additional test samples ====

# visualize workflow
g <- drake::r_drake_ggraph(source = "src/02_augment_test_siamcat.HL.make.R")
g$data %<>%
  mutate(label = if_else(type == "object", label, ""))
g + ggraph::geom_node_text(aes(label = label), angle = 45, size = 2)

# this will test prediction in independent test samples
# the final results will be a table 'performance_test_summary.xls'
drake::r_make(source = "src/02_augment_test_siamcat.HL.make.R")


#
###== TRAIN - FULL Workflow ====

# train models across different settings
# THIS CAN TAKE SEVERAL HOURS OR DAYS DEPENDING ON YOUR HARDWARE

# visualize workflow
g <- drake::r_drake_ggraph(source = "src/01_augment_train_siamcat.HL.HU2.make.R")
g$data %<>%
  mutate(label = if_else(type == "object", label, ""))
g + ggraph::geom_node_text(aes(label = label), angle = 45, size = 1.5)

# compute workflow
drake::r_make(source = "src/01_augment_train_siamcat.HL.HU2.make.R")






