FROM rocker/r-ver:4.3.2

LABEL maintainer="Christopher Nelson <christopher.s.nelson.01@gmail.com>"
LABEL description="R environment for SLE methylation pipeline (ComBat-meth, PCA, Houseman deconvolution)"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN install2.r --error --skipinstalled --ncpus -1 \
    quadprog \
    ggplot2 \
    data.table \
    dplyr \
    optparse \
    matrixStats \
    tidyr \
    RColorBrewer \
    remotes \
    betareg \
    aod \
    RcppArmadillo \
    BiocManager

# Install Bioconductor sva
RUN R -e "BiocManager::install('sva', ask = FALSE, update = FALSE)"

# Install ComBatMet from GitHub (pre-compiled, no runtime compilation needed)
RUN R -e "remotes::install_github('JmWangBio/ComBatMet', upgrade = 'never')"

# Verify all key packages load
RUN R -e "\
    library(quadprog); \
    library(ggplot2); \
    library(data.table); \
    library(dplyr); \
    library(optparse); \
    library(matrixStats); \
    library(tidyr); \
    library(RColorBrewer); \
    library(betareg); \
    library(aod); \
    library(RcppArmadillo); \
    library(sva); \
    library(ComBatMet); \
    message('All packages loaded successfully')"

ENV PATH="/usr/local/bin:${PATH}"
