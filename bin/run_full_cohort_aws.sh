#!/usr/bin/env bash
# Two-phase AWS run for full 9-sample cohort.
#
# Phase 1: Align the 3 Sjogren's samples (downloads from SRA, runs through METHYLDACKEL)
# Phase 2: Consolidate all 9 bedGraphs to a single S3 prefix, then run downstream
#
# Prerequisites:
#   - AWS credentials configured (aws configure)
#   - Nextflow + micromamba on PATH
#   - TOWER_ACCESS_TOKEN exported (optional, for Seqera monitoring)

set -euo pipefail

BUCKET=sle-methylation-pipeline
REGION=us-east-2
GENOME=s3://${BUCKET}/input/GRCh38.primary_assembly.genome.fa
GENOME_FAI=s3://${BUCKET}/input/GRCh38.primary_assembly.genome.fa.fai
BEDGRAPHS_S3=s3://${BUCKET}/input/bedgraphs_9sample

# -----------------------------------------------------------------------
# PHASE 1 — Align Sjogren's samples (alignment_only, downloads from SRA)
# -----------------------------------------------------------------------
echo "=== Phase 1: Sjogren's alignment ==="
nextflow run main.nf \
  -profile aws,conda \
  --sample_sheet sampleSheets/samples_sjogrens_align.csv \
  --genome        ${GENOME} \
  --genome_fai    ${GENOME_FAI} \
  --alignment_only true \
  --outdir        s3://${BUCKET}/results_sjogrens \
  -w              s3://${BUCKET}/work \
  -resume

echo ""
echo "=== Phase 1 complete. Consolidating bedGraphs to ${BEDGRAPHS_S3} ==="

# -----------------------------------------------------------------------
# CONSOLIDATE — copy bedGraphs from both runs into a single S3 prefix
# -----------------------------------------------------------------------

# 6 SLE/Control bedGraphs from the previous 6-sample run
aws s3 cp s3://${BUCKET}/results_6sample/methyldackel/ ${BEDGRAPHS_S3}/ \
    --exclude "*" --include "*_CpG.bedGraph" --recursive --region ${REGION}

# 3 Sjogren's bedGraphs from Phase 1
aws s3 cp s3://${BUCKET}/results_sjogrens/methyldackel/ ${BEDGRAPHS_S3}/ \
    --exclude "*" --include "*_CpG.bedGraph" --recursive --region ${REGION}

echo "Consolidated bedGraphs:"
aws s3 ls ${BEDGRAPHS_S3}/ --region ${REGION} | grep bedGraph

# -----------------------------------------------------------------------
# PHASE 2 — Full downstream analysis (skip alignment, all 9 samples)
# -----------------------------------------------------------------------
echo ""
echo "=== Phase 2: Full 9-sample downstream analysis ==="
nextflow run main.nf \
  -profile aws,conda \
  --sample_sheet    sampleSheets/samples_9sample_downstream.csv \
  --genome          ${GENOME} \
  --genome_fai      ${GENOME_FAI} \
  --skip_alignment  true \
  --bedgraph_dir    ${BEDGRAPHS_S3} \
  --outdir          s3://${BUCKET}/results_9sample \
  -w                s3://${BUCKET}/work \
  -resume

echo ""
echo "=== Pipeline complete. Results at s3://${BUCKET}/results_9sample ==="
