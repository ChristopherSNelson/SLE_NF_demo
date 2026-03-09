# SLE Methylation Nextflow DSL2 Pipeline

## Overview
End-to-end bisulfite sequencing pipeline for Systemic Lupus Erythematosus (SLE) cohort analysis.
FASTQ → alignment → methylation calling → batch correction → cell deconvolution → DMR detection → patient stratification.

## Prerequisites & Setup

### 1. Install Nextflow (requires Java 11+)
```bash
curl -s https://get.nextflow.io | bash
# Move to a directory on your PATH
mv nextflow ~/bin/   # or /usr/local/bin/
```

### 2. Install Micromamba
```bash
# macOS
brew install micromamba
# Linux
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```
Nextflow will use micromamba to create per-process conda environments automatically.
You do NOT need to create any environments manually.

### 3. Download reference genome
```bash
# hg38 primary assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38.fa
samtools faidx hg38.fa
```

### 4. Download SRP410780 data (optional — test profile uses simulated data)
```bash
# Requires sra-tools
for SRR in SRR22476695 SRR22476696 SRR22476697 SRR22476698 SRR22476699 \
           SRR22476700 SRR22476701 SRR22476702 SRR22476703 SRR22476704 SRR22476705; do
    fasterq-dump --split-files --threads 4 ${SRR}
    gzip ${SRR}_1.fastq ${SRR}_2.fastq
done
```

### How Nextflow manages dependencies
- Each process declares a `conda` directive pointing to a YAML in `envs/`
- When you run with `-profile conda`, Nextflow creates and caches envs in `work/conda/`
- Works with either **conda** or **micromamba** — set `conda.useMicromamba` in nextflow.config
  - `conda.useMicromamba = true` → uses micromamba (faster solving, smaller footprint)
  - `conda.useMicromamba = false` → uses standard conda (must be on PATH)
- No manual env creation needed; just have Nextflow + one of conda/micromamba installed

## Current Dataset: SRP410780 (Demo/Validation Cohort)

11 samples, 3 conditions, whole-genome bisulfite sequencing:

| Condition | N | Samples |
|-----------|---|---------|
| SLE | 4 | SRR22476697–700 |
| Sjogren's syndrome | 3 | SRR22476695–696, 703 |
| Healthy control | 4 | SRR22476701–702, 704–705 |

**Note:** SLE and Sjogren's are clinically overlapping autoimmune diseases — some patients have both.
This makes the NMF stratification biologically interesting but statistically underpowered at n=11.

### Limitations at n=11
- NMF feasible for k=2..4 only; consensus matrix is coarse with 11 entries
- Silhouette scores noisy (~3-4 samples per cluster)
- Sufficient to validate all code paths and produce expected outputs
- NOT sufficient for biological discovery or subtype claims

### Future Goal
- Expand to 50+ samples across diverse ethnicities
- SLE disproportionately affects Black, Hispanic, and Asian populations — cohort must reflect this
- Additional SRA WGBS/RRBS cohorts can be incorporated as they become available
- At 50+ samples, NMF rank selection and patient stratification become statistically meaningful

## Run Commands

### Local (default)
```bash
nextflow run main.nf -profile local,conda --sample_sheet samples.csv --genome hg38.fa --outdir results
```

### Test data
```bash
python3 bin/generate_test_data.py --outdir test_data
nextflow run main.nf -profile test,conda
```

### Single-sample alignment test (skip cohort steps)
```bash
nextflow run main.nf -profile local,conda \
  --sample_sheet samples_single.csv \
  --genome GRCh38.primary_assembly.genome.fa \
  --alignment_only true \
  --outdir results_single
```

### HPC (SLURM)
```bash
nextflow run main.nf -profile slurm,conda --sample_sheet samples.csv --genome hg38.fa
```

### AWS Batch
```bash
nextflow run main.nf -profile aws,docker --sample_sheet s3://bucket/samples.csv --genome s3://bucket/hg38.fa
```

### Google Cloud Life Sciences
```bash
nextflow run main.nf -profile gcloud,docker --sample_sheet gs://bucket/samples.csv --genome gs://bucket/hg38.fa
```

**Note:** Cloud profiles require editing `nextflow.config` to set your project/bucket/queue names.

## Executor, Resume & Crash Recovery

### Resume is on by default
- `resume = true` in nextflow.config — completed tasks are skipped on re-run
- Each task is fingerprinted by its inputs + script hash; if nothing changed, the cached result is reused
- To force a full re-run: `nextflow run main.nf -resume false`
- Resume works across executor changes (e.g., start local, move to SLURM)

### Crash recovery
- **Automatic retry:** Processes retry up to 3 times on exit codes 137 (OOM kill), 143 (SIGTERM), 247 (out of memory)
- **Dynamic resource scaling:** Each retry doubles memory and adds CPUs: `memory = { 8.GB * task.attempt }`
- **`errorStrategy = 'finish'`** for non-transient errors — lets running tasks complete before stopping, so you don't waste already-done compute
- **Preemption-safe:** Cloud spot/preemptible instances trigger exit 143, which is automatically retried

### Executor profiles
| Profile | Executor | Use case |
|---------|----------|----------|
| `local` | local | Laptop/workstation, auto-detects available CPUs |
| `slurm` | SLURM | HPC cluster, submits each process as a job |
| `aws` | AWS Batch | Cloud, uses spot instances, requires S3 workDir |
| `gcloud` | Google Life Sciences | Cloud, uses preemptible VMs, requires GCS workDir |

### Combining profiles
Profiles are composable with commas. Always pair an executor with an environment:
```bash
-profile slurm,conda        # HPC + conda envs
-profile aws,docker          # AWS + Docker containers
-profile local,singularity   # Local + Singularity containers
```

### Work directory
- Local: `./work/` (default)
- Cloud: must set `workDir` to `s3://` or `gs://` bucket in the profile
- Clean up after successful run: `nextflow clean -f -before <run_name>`

### Trace files & logging (enabled by default)
All enabled in `nextflow.config` — no flags needed at runtime:
- **Timeline** → `${outdir}/pipeline_info/timeline.html` — Gantt chart of process execution
- **Report** → `${outdir}/pipeline_info/report.html` — resource usage summary per process
- **Trace** → `${outdir}/pipeline_info/trace.txt` — TSV with per-task metrics (CPU%, peak RSS/VMem, I/O, duration)
- **DAG** → `${outdir}/pipeline_info/dag.svg` — process dependency graph (requires Graphviz)
- **Nextflow log** → `.nextflow.log` in the launch directory — full debug log, rotated automatically
- Per-task logs in work dir: `.command.sh`, `.command.err`, `.command.out`, `.command.log`, `.exitcode`

## Pipeline Architecture (12 processes)

| # | Process | Tool | Inputs | Key Outputs |
|---|---------|------|--------|-------------|
| 0 | FETCH_SRA | fasterq-dump | SRR accession | Paired-end FASTQ.gz (skipped for local paths) |
| 1 | FASTP | fastp | Raw FASTQs | Trimmed FASTQs + HTML/JSON QC reports |
| 2 | BWAMETH_INDEX | bwameth index | Genome FASTA | C2T converted index (6 files) |
| 3 | BWAMETH_ALIGN | bwameth | Trimmed FASTQs + index | Sorted BAMs + BAI + flagstat (not published) |
| 5 | MARK_DUPLICATES | picard MarkDuplicates | Sorted BAMs | Dedup BAMs + BAI + metrics (only metrics published) |
| 5b | BAM_TO_CRAM | samtools view -C | Dedup BAMs + genome + .fai | CRAMs + CRAI (40-60% smaller than BAM) |
| 6 | METHYLDACKEL | MethylDackel extract | CRAMs + genome + .fai | CpG bedGraphs + M-bias reports |
| 7 | COMBAT_METH | R (ComBat-meth) | bedGraphs + sample sheet | mvalues_corrected.tsv, beta_matrix.rds, cpg_manifest.tsv |
| 8 | PCA_PLOT | R (ggplot2) | M-value matrices + sample sheet | 4 PNGs (raw/corrected x batch/condition) + pca_variance.tsv |
| 9 | HOUSEMAN_DECONV | R (quadprog) | Beta matrix + reference panel | cell_fractions.tsv + cell_fractions.png |
| 10 | REGION_DETECT | Python (sklearn) | Corrected M-values + sample sheet | candidate_dmrs.bed, window_results.tsv, dmr_manhattan.png |
| 11 | NMF_STRATIFY | Python (sklearn, NMF) | Corrected M-values + sample sheet + cell fractions + DMRs | nmf_clusters.tsv, W/H matrices, rank_selection.png, nmf_umap.png, stability_loo.tsv |

### Sample sheet format
The sample sheet CSV supports two modes:

**Local FASTQs** — provide paths to existing files:
```csv
sample_id,fastq_1,fastq_2,condition,batch
SLE_01,/data/SLE_01_1.fastq.gz,/data/SLE_01_2.fastq.gz,SLE,batch1
```

**SRA accessions** — set fastq_1 to the SRR ID, leave fastq_2 empty:
```csv
sample_id,fastq_1,fastq_2,condition,batch
SRR22476697,SRR22476697,,SLE,batch1
```

Both modes can be mixed in the same sample sheet. SRA rows are auto-detected by the `SRR` prefix and routed through `FETCH_SRA`; local rows go straight to QC.

### Data flow
```
SRR accessions ──► FETCH_SRA ──┐
                                ├──► FASTQs ──► FASTP
Local FASTQs ──────────────────┘                  │
                                                  ▼
                                            BWAMETH_ALIGN ──► MARK_DUPLICATES ──► BAM_TO_CRAM ──► METHYLDACKEL
                                                 ▲                                                      │
                                          BWAMETH_INDEX                                                 ▼
                                                                                                  COMBAT_METH
                                                                                                   │    │   │
                                                                                                   ▼    ▼   ▼
                                                                                                 PCA  HOUSEMAN  REGION_DETECT
                                                                                                         │         │
                                                                                                         ▼         ▼
                                                                                                        NMF_STRATIFY
```

## File Structure
```
SLE_NF_demo/
├── main.nf                     # Entry workflow, channel wiring
├── nextflow.config             # Params, profiles (test, conda, aws), resource defaults
├── samples_single.csv          # Single-sample CSV for alignment testing
├── samples_full.csv            # Full 11-sample SRP410780 cohort
├── modules/
│   ├── fetch_sra.nf
│   ├── fastp.nf                # QC + trimming (replaces fastqc + trim_galore)
│   ├── fastqc.nf               # (legacy, kept for reference)
│   ├── trim_galore.nf          # (legacy, kept for reference)
│   ├── bwameth_index.nf
│   ├── bwameth_align.nf
│   ├── mark_duplicates.nf
│   ├── bam_to_cram.nf
│   ├── methyldackel.nf
│   ├── combat_meth.nf
│   ├── pca_plot.nf
│   ├── houseman_deconv.nf
│   ├── region_detect.nf
│   └── nmf_stratify.nf
├── bin/
│   ├── combat_meth.R
│   ├── pca_plot.R
│   ├── houseman_deconv.R
│   ├── region_detect.py
│   ├── nmf_stratify.py
│   ├── generate_test_data.py
│   └── cleanup_work.sh        # Purge work dir, keep conda cache
├── envs/
│   ├── sra_tools.yml           # sra-tools (fasterq-dump)
│   ├── fastp.yml               # fastp (QC + trimming)
│   ├── bwameth.yml             # bwa-meth, bwa, samtools
│   ├── picard.yml              # picard, samtools
│   ├── methyldackel.yml        # methyldackel
│   ├── r_methylation.yml       # R shared env: quadprog, ggplot2, data.table, dplyr, optparse, matrixStats
│   └── python_ml.yml           # Python shared env: scikit-learn, numpy, pandas, umap-learn, matplotlib
├── containers/
│   └── r_methylation/
│       └── Dockerfile          # Pre-built ComBatMet + all R deps
├── docs/
│   └── aws_setup.md            # Step-by-step AWS Batch deployment guide
└── conf/
    └── test.config             # Test profile: simulated chr22 data, 6 samples
```

## Dependency Strategy

### Per-module conda YAMLs with micromamba
- Each bioinformatics tool gets its own conda env YAML in `envs/`
- R processes (ComBat, PCA, Houseman) share `r_methylation.yml` to avoid redundant R installs
- Python processes (region detect, NMF) share `python_ml.yml`
- Set `conda.useMicromamba = true` in nextflow.config (avoids conda PATH issues)
- All channels: bioconda, conda-forge, defaults

### R environment (`r_methylation.yml`) must include:
- r-base, r-quadprog, r-ggplot2, r-data.table, r-dplyr, r-optparse, r-matrixstats, r-tidyr

### Python environment (`python_ml.yml`) must include:
- python, scikit-learn, numpy, pandas, scipy, matplotlib, seaborn, umap-learn

## Design Decisions

### CRITICAL: Use ComBatMet, NEVER plain ComBat for methylation data
- **DO NOT call `sva::ComBat()` directly.** It is designed for gene expression, not methylation.
- `bin/combat_meth.R` uses the **ComBatMet** R package (Wang, 2025, NAR Genomics & Bioinformatics):
  - GitHub: https://github.com/JmWangBio/ComBatMet
  - Uses beta regression to estimate batch-free distributions
  - Operates directly on beta values (bounded [0,1]), not M-values
  - Realigns quantiles to corrected counterparts — proper for methylation data
- The package is installed from GitHub via `remotes::install_github()` at first run
- Main function: `ComBat_met(beta_matrix, batch = batch_vec, group = group_vec, full_mod = TRUE)`
- After correction, beta values are logit-transformed to M-values for downstream stats/ML
- `r-remotes` is included in `envs/r_methylation.yml` for GitHub installation

### M-values (not beta) for all statistics and ML
- M-value = logit(beta) = log2(beta / (1 - beta))
- Beta values are bounded [0,1] and heteroscedastic; M-values are approximately normal
- Beta matrix is still saved for Houseman deconvolution (reference panels use beta scale)

### Houseman deconvolution via quadprog (NOT minfi)
- minfi requires Bioconductor compilation which fails frequently in conda
- Implement constrained least-squares projection using quadprog directly
- Reference panel: use pre-built methylation profiles for major blood cell types

### Region detection: sub-population-aware approach
- Sliding window over CpGs, not naive t-test
- Use mixed-effects or stratified test to handle patient sub-population correlation
- Low-signal case-control associations are expected in SLE — method must be sensitive to subtle effects
- Output: candidate DMR BED file with window-level statistics

### NMF patient stratification
- Sweep k=2..8
- 30-50 NMF runs per rank for consensus clustering
- Three complementary metrics for rank selection:
  - **Cophenetic correlation**: factorization stability across runs (pick elbow, not max)
  - **Dispersion**: consensus matrix crispness (companion to cophenetic)
  - **Silhouette**: geometric cluster quality on H matrix
- All three plotted together in rank_selection.png
- UMAP visualization of H matrix colored by assigned cluster

#### NMF enhancements (robustness)
- **Cell-type regression**: Houseman cell fractions are regressed out of M-values before NMF.
  Prevents clusters driven by cell-proportion differences rather than epigenetic subtypes.
  NMF now depends on HOUSEMAN_DECONV output (changed pipeline DAG).
- **DMR-based feature selection**: CpGs within candidate DMR regions are used instead of
  top-variance CpGs. Falls back to variance-based selection if < 500 DMR CpGs found.
  NMF now depends on REGION_DETECT output (changed pipeline DAG).
- **LOO stability**: Leave-one-out analysis runs NMF n_samples times, removing one sample
  each iteration. Reports cluster assignment consistency across iterations. Critical for
  small cohorts (n=11) where a single sample can drive cluster structure.
- **Clinical metadata correlation** (optional): If `--clinical_metadata` TSV is provided,
  correlates H-matrix factor weights with numeric clinical variables (e.g., SLEDAI, age).
  Outputs `clinical_correlations.tsv`. Code is ready for when metadata becomes available.

### PCA plots
- 4 plots: {raw, corrected} x {colored by batch, colored by condition}
- PNG only, no SVG output
- Use `aes(.data[[col]])` not deprecated `aes_string()`

## Known Bugs to Avoid (from Bio-omni validation)

These were all encountered and fixed during initial testing. They are baked into the pipeline:

1. **Use `workflow.manifest.version`** not `manifest.version` in nextflow.config
2. **Set `conda.useMicromamba = true`** — conda alone may not be on PATH
3. **Stage .fai as explicit input** in METHYLDACKEL — declare `path fasta_fai` in input block
4. **No `--ignore-flags` for MethylDackel** — flag is invalid; duplicates excluded by default after dedup
5. **Full R package list in conda spec** — include r-optparse, r-data.table, r-dplyr, r-matrixstats
6. **Skip bedGraph track header** — use `skip=1` in R `read.table()` calls for bedGraph input
7. **No `aes_string()` in ggplot2** — deprecated; use `aes(.data[[variable]])` instead
8. **No SVG outputs anywhere** — remove all SVG emit blocks from module output declarations
9. **Houseman via quadprog** — do NOT use minfi; Bioconductor compilation fails in conda
10. **Include scikit-learn in conda spec** — needed by region_detect.py
11. **Cap NMF k_max** — `k_max = min(k_max, n_samples - 1)` to prevent silhouette crash
12. **ComBatMet needs ≥3 samples per batch** — with only 2 per batch, beta regression can't estimate within-batch variance (subscript out of bounds)
13. **gcc version mismatch in conda** — R built with gcc 13.3.0 but conda installs 13.4.0; `combat_meth.R` auto-creates a symlink to fix this
14. **cxx-compiler + gfortran + llvm-openmp required** in `r_methylation.yml` — without these, `updog` (ComBatMet dep) fails to compile from source
15. **`Math.min()` doesn't work on Nextflow MemoryUnit** — use `[4.GB * task.attempt, 16.GB].min()` (Groovy collection min) instead of `Math.min(4.GB * task.attempt, 16.GB)`
16. **Picard MarkDuplicates needs JVM memory cap** — without `-Xmx`, Picard grabs all available heap and OOMs on 16 GB machines. Use `picard -Xmx${avail_mem}g` + `MAX_RECORDS_IN_RAM=500000`
17. **BH correction: don't depend on statsmodels** — implement Benjamini-Hochberg directly with numpy to avoid adding another conda dependency

## Nextflow DSL2 Gotchas

Common pitfalls when writing or modifying this pipeline. These apply to Nextflow DSL2 (>=22.10).

### Channel consumption
- **Channels are consumed on first use.** A queue channel (e.g. `Channel.fromPath(...)`) can only be read by one process. If a second process needs the same data, use `.first()` to convert it to a value channel, or assign it to a variable and reuse.
- **Value channels vs queue channels:** `channel.of(x)` and `.first()` produce value channels (reusable, never consumed). `Channel.fromPath(...)` and `.splitCsv()` produce queue channels (consumed once). Most process inputs that are "reference data" (genome, sample sheet) should be value channels via `.first()`.
- **`.collect()` produces a value channel** containing a list. This is reusable but blocks until all upstream items arrive.

### Process input/output
- **`emit:` labels are mandatory for multi-output processes** if you want to access them by name (`PROC.out.name`). Without labels, outputs are positional only (`PROC.out[0]`).
- **`tuple val(x), path(y)` order matters.** The order in the `input:` block must match the channel structure exactly. A `tuple(sample_id, file)` channel cannot feed a `path(file)` input without `.map { it[1] }` first.
- **`path` inputs are staged (symlinked) into the work dir.** If a tool expects a co-located index file (e.g., `.bam` + `.bai`, `.fa` + `.fai`), both must be declared as separate `path` inputs so they are staged together.
- **Glob outputs like `path("*.html")` fail if zero files match.** Use `optional: true` or ensure the command always produces at least one match.

### Channel operators
- **`.join()` matches on the first element (key) of tuples.** `ch_a.join(ch_b)` merges `[key, a1, a2]` with `[key, b1]` into `[key, a1, a2, b1]`. Keys must match exactly (type and value).
- **`.map {}` vs `.flatMap {}`:** Use `.map` for 1-to-1 transforms; `.flatMap` when each input produces multiple output items.
- **`.combine()` creates a cross product** — use only when you genuinely need every combination.
- **`.mix()` interleaves channels** and does NOT guarantee order.

### Workflow wiring
- **`include` imports are file-level, not directory-level.** `include { FASTQC } from './modules/fastqc'` — the `.nf` extension is optional.
- **Process calls in `workflow {}` are not sequential.** Nextflow infers execution order from data dependencies. If process B needs output from A, wire `A.out.x` into B's input — don't rely on code order.
- **`params` are frozen at startup.** You cannot modify `params.x` inside a process or workflow block.
- **`workflow.manifest.version`** not `manifest.version` — the latter silently resolves to null.

### Script blocks
- **Dollar signs in bash need escaping.** Inside `script: """..."""`, `$variable` is Nextflow interpolation. Use `\$variable` for bash variables, or use `script: '...'` (single-quoted) for pure bash with `!{nf_var}` syntax for Nextflow variables.
- **`task.cpus`, `task.memory`** are available inside script blocks. Use `${task.cpus}` not `${cpus}`.
- **`bin/` scripts are auto-added to PATH.** Any executable in `bin/` can be called by name in script blocks without `./bin/` prefix.

### Conda/environment
- **`conda.useMicromamba = true`** must be set in `nextflow.config` — plain `conda` may not be on PATH (especially on macOS).
- **Conda env creation is per-unique-YAML.** If two processes point to the same YAML, they share one env. Changing the YAML triggers a rebuild.
- **`-profile conda` is required at runtime** to activate conda. Without it, processes try to run in the host environment.

### Common runtime errors
- **"Missing output file"** — the glob pattern in `output:` didn't match any files. Check the exact filename produced by the tool.
- **"No such variable"** — usually a Nextflow vs bash variable collision. Escape bash `$` as `\$`.
- **"Process terminated with an error exit status"** — check `.command.err` and `.command.log` in the work directory for the actual error.
- **Work directory structure:** Each task runs in `work/xx/xxxxxxxxx.../`. Files there: `.command.sh` (the script), `.command.run` (wrapper), `.command.err`, `.command.out`, `.command.log`, `.exitcode`.

## Coding Conventions

- Nextflow DSL2 syntax throughout
- Process names in UPPER_SNAKE_CASE
- All scripts in `bin/` with proper shebangs (`#!/usr/bin/env Rscript`, `#!/usr/bin/env python3`)
- R scripts use optparse for CLI argument parsing
- Python scripts use argparse for CLI argument parsing
- PublishDir pattern: `"${params.outdir}/<stage_name>"`
- Resource labels: `label 'process_low'`, `label 'process_medium'`, `label 'process_high'`

## Build Progress & Next Steps

### COMPLETED
- [x] `CLAUDE.md` — full project spec
- [x] `README.md` — user-facing docs with install/setup
- [x] `nextflow.config` — params, profiles, resource labels, micromamba
- [x] `envs/fastqc.yml`
- [x] `envs/trim_galore.yml`
- [x] `envs/bwameth.yml`
- [x] `envs/picard.yml`
- [x] `envs/methyldackel.yml`
- [x] `envs/r_methylation.yml`
- [x] `envs/python_ml.yml`

### COMPLETED
- [x] All 12 Nextflow modules (modules/*.nf) including FETCH_SRA
- [x] All 5 bin scripts (combat_meth.R, pca_plot.R, houseman_deconv.R, region_detect.py, nmf_stratify.py)
- [x] main.nf — entry workflow with SRA/local branching
- [x] conf/test.config — test profile
- [x] bin/generate_test_data.py — simulated chr22 data generator
- [x] Executor profiles (local/SLURM/AWS/GCloud) with crash recovery
- [x] Fixed conda env specs for osx-arm64 compatibility
- [x] Implemented proper ComBat-meth (no sva dependency)

### COMPLETED
- [x] Pre-made bedGraphs bypass alignment for downstream testing (`--skip_alignment`)
- [x] 6 test samples (3 SLE + 3 Control, 2 batches) for ComBatMet variance estimation
- [x] Fixed ComBatMet installation: `cxx-compiler`, `gfortran`, `llvm-openmp` in conda
- [x] Auto-symlink for gcc version mismatch in `combat_meth.R` (13.3.0 → 13.4.0)
- [x] All 5 downstream processes validated: ComBat-meth, PCA, Houseman, Region Detect, NMF

### Current status: Full downstream pipeline passes with test profile

Test command: `nextflow run main.nf -profile test,conda`
- Skips alignment (uses pre-made bedGraphs)
- Runs ComBat-meth → PCA → Houseman → Region Detect → NMF
- NMF cleanly separates SLE from Control (silhouette=0.953 at k=2)

### COMPLETED
- [x] R methylation Dockerfile with pre-built ComBatMet (`containers/r_methylation/Dockerfile`)
- [x] Replaced FASTQC + Trim Galore with single FASTP process
- [x] `--alignment_only` mode for single-sample alignment testing
- [x] Picard JVM memory cap (`-Xmx` + `MAX_RECORDS_IN_RAM`)
- [x] Fixed `Math.min()` on Nextflow MemoryUnit types
- [x] BH multiple testing correction (removed statsmodels dependency)
- [x] BAM→CRAM conversion (40-60% storage savings)
- [x] AWS Batch profile with Spot instances and raised resource caps
- [x] AWS setup guide (`docs/aws_setup.md`)
- [x] Full 11-sample CSV (`samples_full.csv`)
- [x] Work directory cleanup script (`bin/cleanup_work.sh`)

### COMPLETED
- [x] NMF cell-type regression (regress out Houseman fractions before clustering)
- [x] NMF DMR-based feature selection (use DMR CpGs, fallback to top-variance if < 500)
- [x] NMF leave-one-out stability analysis (outputs stability_loo.tsv)
- [x] Optional `--clinical_metadata` for H-matrix correlation (code ready, waiting for data)
- [x] Pipeline DAG updated: NMF depends on HOUSEMAN_DECONV + REGION_DETECT

### Current status: Overnight single-sample real-data run in progress (SRR22476697)

### TODO — Next steps (in priority order)

#### 1. Validate overnight run results
- Check CRAM, bedGraph, metrics from SRR22476697
- Review peak memory from trace.txt
- Confirm MethylDackel produces non-empty bedGraphs from real data

#### 2. Run full 11-sample cohort on AWS
- Build and push 4 Docker images to ECR (see `docs/aws_setup.md`)
- Set up AWS Batch infrastructure (IAM, S3, compute env, job queue)
- Run `samples_full.csv` with `-profile aws`
- Estimated cost: $44-58 on Spot instances

#### 3. Build remaining Docker images
- `sle-bwameth` (combined bioinformatics tools)
- `sle-python-ml` (scikit-learn, numpy, pandas, umap-learn)
- `sle-sra-tools` (fasterq-dump)
- All need `aws` CLI for S3 file staging on AWS Batch

#### 4. Polish
- Add MultiQC summary step
- Add input validation (check sample sheet columns, genome file exists)
- CI/CD: GitHub Actions workflow for test profile

## Git Commit Conventions

All commits must include both co-authors:
```
Co-Authored-By: Chris Nelson <christopher.s.nelson.01@gmail.com>
Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
```
