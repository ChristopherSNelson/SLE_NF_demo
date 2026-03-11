# HANDOFF — 2026-03-10 Evening Session

## Active Pipeline Run
- tmux session: `pipeline` — attach with `tmux attach -t pipeline`
- seqtk subsampling still in progress (12 parallel jobs, ~20 processes) as of session close
- Nextflow auto-starts after `wait` completes — no action needed
- Nextflow command running:
  ```
  nextflow run main.nf -profile local,conda \
    --sample_sheet samples_chr19.csv \
    --genome chr19.fa \
    --min_depth 1 \
    --outdir results_chr19 \
    -resume 2>&1 | tee chr19_run.log
  ```
- Monitor: `tail -f chr19_run.log`
- Expected completion: midnight–2am (bottleneck is BWAMETH_ALIGN, 6 samples sequential at 6 CPUs)
- Laptop must stay open — lid close forces sleep even with `caffeinate -i`

## What Was Done This Session
- Killed stalled full-genome single-sample alignment (PIDs 16286/16287 and children)
- Pivoted strategy: chr19-only run to get real figures before Thursday 3pm interview
- Installed `seqtk` 1.5 and `tmux` via brew
- Extracted `chr19.fa` (57 MB) from GRCh38, indexed as `chr19.fa.fai`
- Created `samples_chr19.csv` pointing to `fastq_chr19/` subsampled FASTQs
- Launched parallel seqtk subsampling (10M reads/file, 12 files, bash `&` + `wait` in tmux)
- Queued Nextflow step 3 in same tmux session — auto-starts after subsampling

## What's Left (priority order)
1. Wait for chr19 pipeline to finish — verify `results_chr19/` outputs
2. Harvest figures: PCA, NMF UMAP, DMR manhattan, cell fractions, fastp QC, M-bias
3. SLIDES — interview Thursday Mar 12 3pm, not started yet
4. (Optional) delete `fastq_downsampled/` (32 GB) after figures confirmed good

## Key Decisions
- Chr19 over chr22: biologically richer (KIR cluster, gene-dense, relevant for autoimmunity)
- 10M reads per sample: ~300K map to chr19 (~0.8x), `--min_depth 1` to compensate
- Parallel subsampling via bash `&` + `wait` in tmux: all 12 files simultaneously
- `caffeinate -i` wraps full pipeline to prevent idle sleep on battery

## Potential Blockers
- Coverage on chr19 may be thin (0.8x) — if ComBatMet or dmrseq fails, re-run seqtk with 20M reads
  - ComBatMet needs ≥3 samples per batch
  - Original `fastq_downsampled/` (32 GB) preserved for re-subsampling if needed
- Conda `base` has ancient samtools 0.1.19 shadowing brew's 1.23 — use `/opt/homebrew/bin/samtools` in shell

## New Files This Session
- `chr19.fa` + `chr19.fa.fai` — chr19-only reference
- `samples_chr19.csv` — 6-sample sheet pointing to `fastq_chr19/`
- `fastq_chr19/` — subsampled FASTQs (in progress at session close)
- `chr19_run.log` — Nextflow stdout/stderr (created when Nextflow starts)
