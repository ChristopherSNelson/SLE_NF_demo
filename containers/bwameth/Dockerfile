FROM mambaorg/micromamba:1.5-jammy

LABEL maintainer="Christopher Nelson <christopher.s.nelson.01@gmail.com>"
LABEL description="Combined bioinformatics tools for SLE methylation pipeline (fastp, bwa-meth, picard, samtools, MethylDackel)"

# Install all bioinformatics tools into base env
COPY envs/bwameth.yml envs/picard.yml envs/methyldackel.yml envs/fastp.yml /tmp/

RUN micromamba install -y -n base \
    -c bioconda -c conda-forge -c defaults \
    bwameth=0.2.7 bwa samtools=1.20 \
    picard=3.1.1 \
    methyldackel=0.6.1 \
    fastp=0.23.4 \
    && micromamba clean --all --yes

# Nextflow needs aws CLI for S3 file staging
RUN micromamba install -y -n base -c conda-forge awscli \
    && micromamba clean --all --yes

# Verify key tools are available
RUN micromamba run -n base bash -c "\
    bwameth.py --help > /dev/null 2>&1 && \
    samtools --version > /dev/null && \
    picard MarkDuplicates --version 2>&1 | head -1 && \
    MethylDackel --version 2>&1 | head -1 && \
    fastp --version 2>&1 && \
    echo 'All tools verified'"

USER root
