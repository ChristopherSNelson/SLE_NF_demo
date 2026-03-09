FROM mambaorg/micromamba:1.5-jammy

LABEL maintainer="Christopher Nelson <christopher.s.nelson.01@gmail.com>"
LABEL description="Python ML environment for SLE methylation pipeline (region detection, NMF stratification)"

COPY envs/python_ml.yml /tmp/env.yml

RUN micromamba install -y -n base -f /tmp/env.yml \
    && micromamba clean --all --yes

# Nextflow needs aws CLI for S3 file staging
RUN micromamba install -y -n base -c conda-forge awscli \
    && micromamba clean --all --yes

# Verify key packages
RUN micromamba run -n base python3 -c "\
    import sklearn; import numpy; import pandas; \
    import scipy; import matplotlib; import umap; \
    print('All packages loaded successfully')"

USER root
