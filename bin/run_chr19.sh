#!/usr/bin/env bash
# run_chr19.sh — Chr19 subsample + full pipeline run
# Run inside tmux: tmux new -s pipeline
# Survives lid close? No — caffeinate -i pauses on lid close. Keep lid open.

set -e
cd /Users/chrisnelson/Desktop/SLE_NF_demo

# Step 1: Extract chr19 (skip if already done)
if [ ! -f chr19.fa ]; then
    echo ">>> Extracting chr19..."
    /opt/homebrew/bin/samtools faidx GRCh38.primary_assembly.genome.fa chr19 > chr19.fa
    /opt/homebrew/bin/samtools faidx chr19.fa
    echo ">>> chr19.fa ready"
else
    echo ">>> chr19.fa already exists, skipping"
fi

# Step 2: Subsample 10M reads per file in parallel (skip if already done)
if [ ! -d fastq_chr19 ] || [ "$(ls fastq_chr19/*.fastq.gz 2>/dev/null | wc -l)" -lt 12 ]; then
    echo ">>> Subsampling 10M reads per file..."
    mkdir -p fastq_chr19
    for f in fastq_downsampled/*.fastq.gz; do
        seqtk sample -s42 "$f" 10000000 | gzip > "fastq_chr19/$(basename "$f")" &
    done
    wait
    echo ">>> Subsampling done"
    ls -lh fastq_chr19/
else
    echo ">>> fastq_chr19/ already populated, skipping"
fi

# Step 3: Run Nextflow pipeline
echo ">>> Starting Nextflow..."
caffeinate -i nextflow run main.nf -profile local,conda \
    --sample_sheet samples_chr19.csv \
    --genome chr19.fa \
    --min_depth 1 \
    --outdir results_chr19 \
    -resume 2>&1 | tee chr19_run.log
