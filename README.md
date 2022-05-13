# Supplementary Code for Publication

## Overview over Software Environment
- This repository contains code to train [SIAMCAT](https://siamcat.embl.de) models with augmented samples
- The code is written in [R](https://www.r-project.org) and executed in a [docker](https://www.docker.com) container (based on [rocker/tidyverse](https://hub.docker.com/r/rocker/verse))
- Caching and workflow is control using [Drake](https://books.ropensci.org/drake/) and [tidyverse](https://www.tidyverse.org)
- [gnu make](https://www.gnu.org/software/make/) is used to start the container
- Smaller models from the publication are shipped with this repository, but full models can be build from the supplied code

## Installation
- this program only requires [docker](https://www.docker.com) and [gnu make](https://www.gnu.org/software/make/) on the system
  - in Ubuntu or similar, this can be done using **sudo apt install make dockerd**
  - this program was only tested on linux machines, but should run on any **x86 or x64 system**, but **NOT M1** processors
- this process tested on SLES15 and Ubuntu 20.04
- Installation time depends on your internet connection
  - The size of the **docker image** is **2.15GB compressed** (download size) and **6.15GB decompressed** (final size on hard-drive)

## Run
- start an r-studio server within this repository: **make env**
  - this command will download the docker image & start the container
  - it will also inform you about the **url** at which the container has been started ([localhost:4444](http://localhost:4444) or [lcc_rserver:4444](http://lcc_rserver:4444))
  	- if the port 4444 is already in use, you need to change it to some other number within the **Makefile**
  - **type the corresponding url into a web-browser**
- The **credentials** are:
  - **user** : rstudio
  - **password**: tesseract
- you can send execute the scripts within **src/** in order
  - **src/00_main_build.R** is the main entry point, including further comments on the individual scripts
  - it also contains a minimal example that will build models based on **bacterial abundances** as used in the manuscript
  - the resulting **SIAMCAT models** will be saved in **rds** format, which is a native data type for **R** data objects
- the **example_output** directory pre-computed results as a result of testing this repository

## Stop
- once you are down, you can simply stop and/or remove the container
  - identify the **ID** container using **docker ps** (e.g. 09d17d227d)
  - docker rm -f ID

