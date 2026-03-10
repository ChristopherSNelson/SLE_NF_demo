# SLE Methylation Pipeline

End-to-end Nextflow DSL2 pipeline for Systemic Lupus Erythematosus (SLE) bisulfite sequencing analysis.

**FASTQ → alignment → methylation calling → batch correction → cell deconvolution → DMR detection → patient stratification**

## Pipeline Stages

| Stage | Tool | Description |
|-------|------|-------------|
| QC + trimming | fastp | Adapter trimming, quality filtering, QC reports |
| Alignment | bwa-meth | Bisulfite-aware alignment to hg38 |
| Deduplication | Picard | Mark and remove PCR duplicates |
| CRAM conversion | samtools | BAM → CRAM (40–60% smaller) |
| Methylation calling | MethylDackel | Per-CpG methylation extraction |
| Batch correction | ComBat-meth | Methylation-specific batch correction on M-values |
| PCA | R/ggplot2 | Raw and corrected PCA plots by batch and condition |
| Cell deconvolution | Houseman (quadprog) | Blood cell type fraction estimation |
| DMR detection | dmrseq | Count-level DMR detection with permutation p-values |
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

Validation cohort: [SRP410780](https://www.ncbi.nlm.nih.gov/sra/?term=SRP410780) (GEO: [GSE213478](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213478)) — 11 whole-genome bisulfite sequencing samples (4 SLE, 3 Sjogren's, 4 healthy controls).

### Downsampled runs

The current workflow uses pre-subsampled local FASTQs — the pipeline no longer auto-fetches from SRA (`FETCH_SRA` is disabled in `main.nf`). Reads are subsampled to ~50 M read pairs per sample (~2× CpG coverage, sufficient at `--min_depth 2`) using `fastq-dump -X 50000000 --split-files --gzip` (sra-tools; Leinonen et al., *Nucleic Acids Research* 2011) and stored in `fastq_downsampled/`. To re-enable SRA streaming, uncomment `FETCH_SRA` in `main.nf`.

The 6-sample downsampled subset (`samples_downsampled.csv`) uses 3 SLE + 3 Control samples balanced across 2 batches (SRR22476697–698, SRR22476700–701, SRR22476704–705).

This dataset exercises all pipeline code paths but is underpowered for biological discovery (n=11). The future goal is to expand to 50+ samples across diverse ethnicities, reflecting the disproportionate impact of SLE on Black, Hispanic, and Asian populations.

## Key Design Choices

- **ComBat-meth** (not plain ComBat) for methylation-appropriate batch correction
- **M-values** (logit of beta) for all statistical and ML analyses
- **Houseman via quadprog** — avoids Bioconductor/minfi compilation issues
- **NMF rank selection** uses cophenetic correlation + dispersion + silhouette together
- **VDJ risk annotation** on DMRs — immune receptor loci (IGH, IGK, IGL, TRA/TRD, TRB, TRG) undergo somatic recombination that mimics differential methylation. DMRs overlapping these regions are flagged with `vdj_risk=TRUE`
- **Per-module conda environments** managed automatically by Nextflow

## Companion: ImmuneMethylTools

[ImmuneMethylTools](https://github.com/ChristopherSNelson/ImmuneMethylTools) is a complementary pre-alignment QC framework for immune-cell methylation studies. It provides upstream quality gates that run before this pipeline's alignment stage:

- **Bisulfite conversion QC** — flags samples with non-CpG methylation rate >1% (incomplete conversion)
- **Sample contamination detection** — identifies muddy beta distributions via Sarle's Bimodality Coefficient
- **Batch × disease confound check** — Cramér's V test to detect batch-confounded designs before batch correction absorbs biological signal
- **VDJ clonal artifact masking** — detects and masks CpGs in immune receptor loci where somatic recombination creates artifactual methylation signals, preserving samples for downstream analysis rather than dropping them

Together, the two tools form a full QC → analysis workflow for autoimmune WGBS cohorts.

## License

TBD
