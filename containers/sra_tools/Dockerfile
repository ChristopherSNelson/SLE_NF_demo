FROM mambaorg/micromamba:1.5-jammy

LABEL maintainer="Christopher Nelson <christopher.s.nelson.01@gmail.com>"
LABEL description="SRA tools for SLE methylation pipeline (fasterq-dump for SRA data download)"

COPY envs/sra_tools.yml /tmp/env.yml

RUN micromamba install -y -n base -f /tmp/env.yml \
    && micromamba clean --all --yes

# Nextflow needs aws CLI for S3 file staging
RUN micromamba install -y -n base -c conda-forge awscli \
    && micromamba clean --all --yes

# Verify
RUN micromamba run -n base fasterq-dump --version 2>&1 | head -1 \
    && echo "sra-tools verified"

USER root
