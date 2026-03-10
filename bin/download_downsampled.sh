#!/bin/bash
#
# Download downsampled WGBS reads from 6 SRP410780 samples.
# Uses fastq-dump -X to stream only the first N reads from NCBI (no prefetch).
# ~50M reads/sample = ~2x coverage. Enough for methylation calling at min_depth=2.
#
set -euo pipefail

NUM_READS=50000000  # 50M reads (~2x WGS coverage, ~5.5% of full depth)
OUTDIR="downsampled_fastq"
CONDA_ENV="work/conda/env-83b9441f58bbf195-81269c29b376af0d44fd374274f6f1f3"

# 6 samples: 3 SLE + 3 Control, balanced across 2 batches
# (meets ComBatMet ≥3/batch and dmrseq ≥2/condition requirements)
SAMPLES=(
    "SRR22476697"   # SLE, batch1
    "SRR22476698"   # SLE, batch1
    "SRR22476700"   # SLE, batch2
    "SRR22476701"   # Control, batch2
    "SRR22476704"   # Control, batch1
    "SRR22476705"   # Control, batch1
)

mkdir -p "$OUTDIR"

echo "========================================="
echo " Downloading ${NUM_READS} reads per sample"
echo " Samples: ${#SAMPLES[@]}"
echo " Output: $OUTDIR/"
echo "========================================="

for SRR in "${SAMPLES[@]}"; do
    if [ -f "$OUTDIR/${SRR}_1.fastq.gz" ] && [ -f "$OUTDIR/${SRR}_2.fastq.gz" ]; then
        echo "[$SRR] Already downloaded, skipping"
        continue
    fi

    echo ""
    echo "[$SRR] Downloading ${NUM_READS} reads..."
    START=$(date +%s)

    # Stream directly from NCBI — no prefetch, no temp files
    # Use the pipeline's cached conda env for sra-tools
    if [ -d "$CONDA_ENV" ]; then
        conda run -p "$CONDA_ENV" \
            fastq-dump -X "$NUM_READS" --split-files --gzip \
            --outdir "$OUTDIR" "$SRR"
    else
        # Fallback: try system sra-tools
        fastq-dump -X "$NUM_READS" --split-files --gzip \
            --outdir "$OUTDIR" "$SRR"
    fi

    END=$(date +%s)
    ELAPSED=$(( END - START ))
    SIZE=$(du -sh "$OUTDIR/${SRR}_1.fastq.gz" 2>/dev/null | cut -f1)
    echo "[$SRR] Done in ${ELAPSED}s (R1: $SIZE)"
done

echo ""
echo "========================================="
echo " Download complete"
echo "========================================="
du -sh "$OUTDIR"
ls -lh "$OUTDIR"
