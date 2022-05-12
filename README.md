# Supplementary Code for Publication

## Primary information about the software environment
- This repository contains code to train [SIAMCAT](https://siamcat.embl.de) models with augmented samples
- The code is written in [R](https://www.r-project.org) 4.0.4 and executed in a [docker](https://www.docker.com) container (based on [rocker/tidyverse](https://hub.docker.com/r/rocker/verse))
- Caching and workflow is control using [Drake](https://books.ropensci.org/drake/) and [tidyverse](https://www.tidyverse.org)
- [gnu make](https://www.gnu.org/software/make/) is used to start the container
- Smaller models from the publication are shipped with this repository, but full models can be build from the supplied code

## Install
- this program only requires [docker](https://www.docker.com) and [gnu make](https://www.gnu.org/software/make/) on the system
  - in ubuntu or similar, this can be done using **sudo apt install make dockerd**
  - this program was only tested on linux machines, but should run on any x86 or x64 system, **but NOT M1 processors**

## Run
- start an r-studio server within this repository: **make env**
  - this command will download the docker image & start the container
  - it will also inform you about the **url** at which the container has been started
  - **type the corresponding url into a web-browser**
- you can send execute the scripts within **src/** in order
  - **src/00_main_build.R** is the main entry point, including further comments on the individual scripts
  - it also contains a minimal example that will build models based on **bacterial abundances** as used in the manuscript
  - the resulting SIAMCAT models will be saved in **rds** format, which is a native data type for **R** data objects
 
