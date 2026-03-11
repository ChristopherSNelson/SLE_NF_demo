# HANDOFF — 2026-03-10 Late Evening Session

## Active Pipeline Run
- tmux session: `pipeline` — attach with `tmux attach -t pipeline`
- Subsampling in progress: 4/12 files done as of session close, 8 remaining
- Nextflow will auto-start after subsampling completes
- Nextflow command:
  ```
  nextflow run main.nf -profile local,conda \
    --sample_sheet samples_chr19.csv \
    --genome chr19.fa \
    --min_depth 1 \
    --outdir results_chr19 \
    -resume 2>&1 | tee chr19_run.log
  ```
- Monitor subsampling: `ls -lh fastq_chr19/`
- Monitor Nextflow: `tail -f chr19_run.log`
- Expected subsampling done: ~10:30pm
- Expected Nextflow done: ~2–3am (bottleneck is BWAMETH_ALIGN x6 samples)
- Laptop must stay open — lid close forces sleep even with `caffeinate -i`

## What Was Done This Session
- Killed stalled all-parallel seqtk jobs (I/O contention, 15+ hrs estimated)
- Switched subsampling to `gzip -dc | head -n 40000000 | gzip`, 2 jobs at a time
- Discovered macOS `zcat` silently produces empty output on `.gz` files — first run produced empty FASTQs, pipeline "completed" in minutes with zero aligned reads
- Fixed script to use `gzip -dc`; wiped bad results (`rm -rf results_chr19/`, `nextflow clean -f`)
- Re-ran subsampling — confirmed files now have real content (~300–940 MB each)
- Updated `bin/run_chr19.sh` — idempotent, skips non-empty files

## What's Left (priority order)
1. Wait for chr19 pipeline to finish — verify `results_chr19/` outputs
2. Harvest figures: PCA, NMF UMAP, DMR manhattan, cell fractions, fastp QC, M-bias
3. SLIDES — interview Thursday Mar 12 3pm, not started yet
4. (Optional) delete `fastq_downsampled/` (32 GB) after figures confirmed good

## Key Decisions
- chr19 over chr22: biologically richer (KIR cluster, gene-dense, relevant for autoimmunity)
- 10M reads per sample: ~200–300K map to chr19, `--min_depth 1` to compensate
- `gzip -dc | head -n 40000000`, 2 at a time: avoids both reservoir-sampling overhead and I/O contention
- `(gzip -dc "$f" | head -n 40000000; true) | gzip` — the `; true` suppresses broken-pipe SIGPIPE from head exiting early

## Potential Blockers
- Coverage on chr19 may be thin — if ComBatMet or dmrseq fails, re-run with 20M reads (`head -n 80000000`) and wipe cache
- ComBatMet needs ≥3 samples per batch
- Original `fastq_downsampled/` (32 GB) preserved for re-subsampling if needed

## New/Modified Files This Session
- `bin/run_chr19.sh` — fixed: `gzip -dc` instead of `zcat`, 2-at-a-time throttle, skip non-empty
- `fastq_chr19/` — 4/12 files complete at session close
