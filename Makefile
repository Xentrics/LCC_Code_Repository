#!/usr/bin/make -f

# This is the main entry point for the analysis
# All results must be created by making the rule 'analysis'

.PHONY: build run env status analyse

USER_ID := $(shell id -u)
GROUP_ID := $(shell id -g)
PROJECT_ROOT := $(shell pwd)

CUSTOM_NAME := lcc_siamcat
# ongoing numbering of env instances
DOCKER_ENV_NR := $(shell docker container ls -f name=^/$(CUSTOM_NAME)_env$ | wc -l)
DOCKER_ENV_PORT := "4444"
DOCKER_ENV_NAME := $(CUSTOM_NAME)_env_$(DOCKER_ENV_NR)
DOCKER_IMAGE := "xentrics/lcc_code_repository:latest"


define DOCKER_RUN_ARGS
	--volume ${PWD}/../:/project \
	--volume ${PWD}:/project/analysis \
	--hostname $(CUSTOM_NAME)
endef

# Run docker container. This will also run the analysis by default
run: build
	docker run \
		--user $(USER_ID):$(GROUP_ID) \
		-v $(PROJECT_ROOT)/src/Rprofile/run.R:/usr/local/lib/R/etc/Rprofile.site \
		$(DOCKER_RUN_ARGS) $(DOCKER_NAME) make env

# Start docker container for interactive analysis
# - detached RStudio Server
# - user rstudio has the same UID than the host user
# - user rstudio as the analysis root directory as home
env:
	docker run -d --rm \
		--user root:root \
		-e PASSWORD=tesseract \
		-p $(DOCKER_ENV_PORT):8787 \
		--name $(DOCKER_ENV_NAME) \
		$(DOCKER_RUN_ARGS) $(DOCKER_IMAGE) \
		bash -c " \
			groupadd -g $(GROUP_ID) GROUP ; \
			usermod -u $(USER_ID) -g $(GROUP_ID) -d /project/analysis rstudio && \
			env >> /usr/local/lib/R/etc/Renviron && \
			/init" && \
		echo SUCCESS: Container $(DOCKER_ENV_NAME) started at http://$(shell hostname):$(DOCKER_ENV_PORT)

# Status of DVC, git and make
status:
	@echo "----------- DVC status -----------"
	dvc status
	@echo "\n\n----------- git status -----------"
	git status
	@echo "\n\n---------- make status -----------"
	make -n

pull:
	docker pull $(DOCKER_IMAGE)

