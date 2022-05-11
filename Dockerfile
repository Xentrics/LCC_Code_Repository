# Please consider a version stable image e.g. do not use the tag latest
FROM sbi-registry.hki-jena.de/seelbind/dockerfiles:rserver_analyse_2021_03

USER root
WORKDIR /build
RUN apt-get update && apt-get install -y make && apt-get clean

# update default mirrors
ENV CRAN="https://packagemanager.rstudio.com/all/__linux__/focal/2021-02-25+Y3JhbiwyOjc4NDM0OTs4NzIzN0YzQg"
ENV BioC_mirror="https://packagemanager.rstudio.com/bioconductor"
ENV R_HOME="${R_HOME:-/usr/local/lib/R}"

RUN echo "options(repos = c(CRAN = '${CRAN}'), download.file.method = 'libcurl')" >> "${R_HOME}/etc/Rprofile.site" && \
    echo "options(BioC_mirror = 'https://packagemanager.rstudio.com/bioconductor')" >> "${R_HOME}/etc/Rprofile.site"

COPY src/install/install.R     ./
COPY src/install/packages.txt  ./src/install/packages.txt
COPY src/install/requirements/ ./src/install/requirements/
RUN Rscript install.R

USER rstudio
EXPOSE 8787
WORKDIR /analysis
CMD make analyse

