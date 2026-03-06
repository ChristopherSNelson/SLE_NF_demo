# SLE Methylation Pipeline

End-to-end Nextflow DSL2 pipeline for Systemic Lupus Erythematosus (SLE) bisulfite sequencing analysis.

**FASTQ → alignment → methylation calling → batch correction → cell deconvolution → DMR detection → patient stratification**

## Pipeline Stages

| Stage | Tool | Description |
|-------|------|-------------|
| QC | FastQC | Pre- and post-trim quality reports |
| Trimming | Trim Galore | Adapter and quality trimming |
| Alignment | bwa-meth | Bisulfite-aware alignment to hg38 |
| Deduplication | Picard | Mark and remove PCR duplicates |
| Methylation calling | MethylDackel | Per-CpG methylation extraction |
| Batch correction | ComBat-meth | Methylation-specific batch correction on M-values |
| PCA | R/ggplot2 | Raw and corrected PCA plots by batch and condition |
| Cell deconvolution | Houseman (quadprog) | Blood cell type fraction estimation |
| Region detection | Python/sklearn | Sub-population-aware DMR candidate identification |
| Patient stratification | NMF | Unsupervised subtype discovery with rank selection |

## Prerequisites

### Nextflow (requires Java 11+)

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/   # or /usr/local/bin/
```

### Micromamba

```bash
# macOS
brew install micromamba

# Linux
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Nextflow will use micromamba to create per-process conda environments automatically. You do **not** need to create any environments manually.

### Reference genome

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38.fa
samtools faidx hg38.fa
```

## Quick Start

### Run with test data (simulated chr22, ~5 min)

```bash
nextflow run main.nf -profile test,conda
```

### Run with real data

1. Create a sample sheet CSV:

```csv
sample,fastq_1,fastq_2,condition,batch,age
SRR22476697,/path/to/SRR22476697_1.fq.gz,/path/to/SRR22476697_2.fq.gz,SLE,batch1,51
SRR22476701,/path/to/SRR22476701_1.fq.gz,/path/to/SRR22476701_2.fq.gz,control,batch1,21
```

2. Run:

```bash
nextflow run main.nf \
    -profile conda \
    --sample_sheet samplesheet.csv \
    --genome /path/to/hg38.fa \
    --outdir /path/to/results
```

## Output Structure

```
results/
├── fastqc/                 # QC HTML reports
├── trim_galore/            # Trimmed FASTQs
├── alignments/             # Sorted, deduplicated BAMs
├── methyldackel/           # CpG bedGraphs + M-bias plots
├── combat_meth/            # Batch-corrected M-value and beta matrices
├── pca/                    # PCA plots (raw/corrected x batch/condition)
├── houseman/               # Cell type fraction estimates
├── region_detect/          # Candidate DMR regions + manhattan plot
├── nmf/                    # Patient clusters, W/H matrices, rank selection plot, UMAP
└── pipeline_info/          # Nextflow timeline, report, trace, DAG
```

## Dataset

Current validation cohort: [SRP410780](https://www.ncbi.nlm.nih.gov/sra/?term=SRP410780) — 11 WGBS samples (4 SLE, 3 Sjogren's, 4 healthy controls).

This dataset exercises all pipeline code paths but is underpowered for biological discovery (n=11). The future goal is to expand to 50+ samples across diverse ethnicities, reflecting the disproportionate impact of SLE on Black, Hispanic, and Asian populations.

## Key Design Choices

- **ComBat-meth** (not plain ComBat) for methylation-appropriate batch correction
- **M-values** (logit of beta) for all statistical and ML analyses
- **Houseman via quadprog** — avoids Bioconductor/minfi compilation issues
- **NMF rank selection** uses cophenetic correlation + dispersion + silhouette together
- **Per-module conda environments** managed automatically by Nextflow

## License

TBD
