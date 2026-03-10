# HANDOFF — 2026-03-10

## Active pipeline run

Single-sample alignment is running in a separate terminal:
- PIDs: 16286 (JVM), 16287 (caffeinate wrapper)
- Command: `nextflow run main.nf -profile local,conda --sample_sheet samples_single_downsampled.csv --genome GRCh38.primary_assembly.genome.fa --alignment_only true --min_depth 2 --outdir results_single_downsampled -resume`
- Expected duration: 6-10h (genome index cached, 6 CPUs)
- Monitor: `tail -f nf_run.log` or check `.nextflow.log`
- Kill if needed: `kill -9 16286 16287`

## What was done this session

- Fixed `process_high` cpus: 6 min / 8 max (was always 8; machine has 8 physical cores)
- Added VDJ risk annotation to `bin/region_detect.R` — flags DMRs overlapping IGH/IGK/IGL/TRA_TRD/TRB/TRG with `vdj_risk` + `vdj_locus` columns
- Fixed IGH hg38 start coordinate: 105586437 → 105586937 (verified against IMGT)
- Dry-run validated VDJ annotation logic with synthetic DMRs (6/6 loci flagged, clean control correct)
- Updated README: added pipeline DAG image, ImmuneMethylTools companion section, removed inline bold
- Generated `docs/dag.png` via `nextflow run main.nf -preview -with-dag`
- Added to CLAUDE.md: guardrail alerts table, verification loop, handoff discipline, mistakes log, model selection nuance, kill -9 tip, no-bold rule
- Ported relevant sections from bcherny-claude CLAUDE.md
- Created `/shutdown` slash command at `.claude/commands/shutdown.md` and `~/.claude/commands/shutdown.md`
- Pushed private GitHub repo: https://github.com/ChristopherSNelson/SLE_NF_demo
- All changes committed and pushed to main

## What's left (priority order)

1. Wait for single-sample alignment to finish — verify CRAM, bedGraph, M-bias plot, flagstat, trace.txt peak memory
2. Run full 6-sample pipeline: `nextflow run main.nf -profile local,conda --sample_sheet samples_downsampled.csv --genome GRCh38.primary_assembly.genome.fa --min_depth 2 --outdir results_downsampled`
3. Harvest figures for interview slides: PCA, NMF UMAP, DMR manhattan, cell fractions, fastp QC
4. AWS Batch: blocked on billing (no payment method); region us-east-2
5. Future polish: MultiQC, GitHub Actions CI, nf-core-style subway map diagram

## Key decisions

- FETCH_SRA stays disabled; pipeline uses `fastq_downsampled/` with `samples_downsampled.csv`
- VDJ annotation is post-hoc (no new Nextflow process needed) — just a flag column in DMR output
- ImmuneMethylTools positioned as complementary pre-alignment QC, not merged into this pipeline

## Blockers

- Single-sample run must complete before full 6-sample run (confirm alignment works end-to-end first)
- AWS Batch needs payment method added to account
- Interview Thursday Mar 12 — need real data figures before then
