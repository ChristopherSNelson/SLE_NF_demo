# Handoff — 2026-03-11 ~1:30am

## Active Run
- **tmux:** `fullrun` — `tmux attach -t fullrun`
- **PIDs:** 91471/91472 (kill with `kill -9 91471 91472` if needed)
- **Command:** `caffeinate -i nextflow run main.nf -profile local,conda --sample_sheet samples_fullgenome.csv --genome GRCh38.primary_assembly.genome.fa --outdir results_full -resume`
- **Status at handoff:** FASTP done (6/6 cached), BWAMETH_ALIGN just started (0/6)
- **Monitor:** `tail -f .nextflow.log`
- **Duration:** Unknown — no benchmark exists for bwameth on this machine at this data size. Could be 2h, could be 12h+.

## What Was Done This Session
- Chr19 pipeline run completed successfully — figures in `results_chr19/`
- Fixed ComBatMet crash on sparse data: `tryCatch` fallback to raw beta in `bin/combat_meth.R`
- Fixed MethylDackel mbias: `rsvg-convert` SVG→PNG, `librsvg` added to conda env, proper publishDir
- Converted existing chr19 mbias SVGs to PNG manually: `results_chr19/methyldackel/`
- Assessed chr19 figures for interview (see table below)
- Discovered `fastq_downsampled/` files are systemically corrupted — R1/R2 read count mismatches
- Created `samples_fullgenome.csv` using clean `fastq_chr19/` files (10M pairs) vs full genome
- Cleaned `work/` task dirs via `bin/cleanup_work.sh` — freed 95GB, conda envs at 4.5GB preserved
- Fixed batch assignments: 3/3 split, both conditions per batch — ComBatMet will run
- 3 atomic commits pushed

## Chr19 Figures — Usable for Slides
| Figure | Path | Use? |
|--------|------|------|
| NMF UMAP | `results_chr19/nmf/nmf_umap.png` | Yes |
| Rank selection | `results_chr19/nmf/rank_selection.png` | Yes |
| DMR manhattan | `results_chr19/region_detect/dmr_manhattan.png` | Yes — 502 regions |
| PCA raw batch | `results_chr19/pca/pca_raw_batch.png` | Yes — batch effect visible |
| M-bias PNGs | `results_chr19/methyldackel/*_mbias_OT/OB.png` | Yes — flat/clean |
| fastp HTML | `results_chr19/fastp/SRR22476697_fastp.html` | Yes |
| Cell fractions | `results_chr19/houseman/cell_fractions.png` | No — flat chr19 artifact |
| PCA corrected condition | `results_chr19/pca/pca_corrected_condition.png` | No — no separation |

## Next Steps (priority order)
1. Morning: check `tmux attach -t fullrun` — did alignment finish?
2. If `results_full/` complete: harvest PCA, cell fractions, NMF, DMR figures
3. BUILD SLIDES — interview Thursday Mar 12 2:30pm, nothing started yet
4. Slide structure suggestion:
   - Motivation: SLE disproportionately affects underrepresented populations
   - Pipeline architecture DAG
   - QC: fastp + M-bias (chr19 real data)
   - Batch correction: PCA before/after
   - Cell deconvolution: cell fractions
   - DMR detection: manhattan
   - NMF stratification: rank selection + UMAP
   - Future: 50+ samples, diverse cohort

## Key Decisions
- `fastq_chr19/` used for full-genome run — only clean balanced files available
- `fastq_downsampled/` abandoned — systemically corrupted
- 3/3 batch split in `samples_fullgenome.csv` — ComBatMet should run properly
- No SVG outputs — mbias converted to PNG via rsvg-convert

## Blockers
- Alignment time unknown — may not finish before Thursday afternoon
- Fallback: chr19 figures + test profile figures (labeled simulated) cover all slide slots
- Claude usage ~77% weekly limit — may need Gemini for slide session
