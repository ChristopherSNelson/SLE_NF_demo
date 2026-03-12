#!/usr/bin/env Rscript

# DMR detection using dmrseq (Korthauer et al., 2019)
# Works on count-level data from MethylDackel bedGraphs with batch as covariate.
# Produces variable-width DMRs with permutation-based p-values.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

# Auto-install dmrseq + bsseq from Bioconductor if not available
# (conda solver can't satisfy dmrseq's version-pinned dep chain on any platform)
if (!requireNamespace("dmrseq", quietly = TRUE) ||
    !requireNamespace("bsseq", quietly = TRUE)) {
  cat("Installing dmrseq + bsseq from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(c("dmrseq", "bsseq"), ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dmrseq)
  library(bsseq)
})

option_list <- list(
  make_option("--bedgraphs", type = "character", help = "Comma-separated bedGraph file paths"),
  make_option("--sample_sheet", type = "character", help = "Sample sheet CSV"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Load sample sheet ---
samples <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE)
bg_files <- strsplit(opt$bedgraphs, ",")[[1]]
bg_files <- trimws(bg_files)

cat("Loading", length(bg_files), "bedGraph files...\n")

# --- Parse bedGraphs into BSseq object ---
# MethylDackel format: chr start end beta% count_m count_u (with track header)
all_data <- list()

for (bg_file in bg_files) {
  # Extract sample_id from filename (e.g., SLE_01_CpG.bedGraph -> SLE_01)
  sid <- sub("_CpG\\.bedGraph$", "", basename(bg_file))

  dt <- fread(bg_file, skip = 1, header = FALSE,
              col.names = c("chr", "start", "end", "beta_pct", "count_m", "count_u"))

  dt$total <- dt$count_m + dt$count_u
  dt$sample_id <- sid
  all_data[[sid]] <- dt
}

# Find common CpG sites across all samples
cpg_keys <- lapply(all_data, function(dt) paste(dt$chr, dt$start, sep = ":"))
common_cpgs <- Reduce(intersect, cpg_keys)
cat("Found", length(common_cpgs), "common CpG sites\n")

if (length(common_cpgs) < 10) {
  cat("Too few common CpGs. Writing empty outputs.\n")
  write.table(data.frame(), file.path(opt$outdir, "candidate_dmrs.bed"),
              sep = "\t", row.names = FALSE, col.names = FALSE)
  write.table(data.frame(), file.path(opt$outdir, "window_results.tsv"),
              sep = "\t", row.names = FALSE)
  png(file.path(opt$outdir, "dmr_manhattan.png"), width = 800, height = 400)
  plot.new(); title("No DMRs detected (too few common CpGs)")
  dev.off()
  quit(save = "no")
}

# Build count matrices (CpGs x samples)
sample_ids <- names(all_data)
ref_dt <- all_data[[sample_ids[1]]]
ref_dt$key <- paste(ref_dt$chr, ref_dt$start, sep = ":")
ref_dt <- ref_dt[ref_dt$key %in% common_cpgs, ]
ref_dt <- ref_dt[order(ref_dt$chr, ref_dt$start), ]

M_mat <- matrix(0L, nrow = nrow(ref_dt), ncol = length(sample_ids))
Cov_mat <- matrix(0L, nrow = nrow(ref_dt), ncol = length(sample_ids))
colnames(M_mat) <- colnames(Cov_mat) <- sample_ids

for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  dt <- all_data[[sid]]
  dt$key <- paste(dt$chr, dt$start, sep = ":")
  dt <- dt[match(ref_dt$key, dt$key), ]
  M_mat[, i] <- as.integer(dt$count_m)
  Cov_mat[, i] <- as.integer(dt$total)
}

# Create BSseq object
bs <- BSseq(
  chr = ref_dt$chr,
  pos = ref_dt$start,
  M = M_mat,
  Cov = Cov_mat,
  sampleNames = sample_ids
)

cat("BSseq object:", nrow(bs), "CpGs x", ncol(bs), "samples\n")

# --- Match samples to metadata ---
# Align sample sheet to BSseq sample order
meta <- samples[match(sample_ids, samples$sample_id), ]

# Condition as factor (test SLE vs reference)
# Use the most common condition as reference, or "Control" if present
conditions <- meta$condition
if ("Control" %in% conditions) {
  conditions <- factor(conditions, levels = c("Control", setdiff(unique(conditions), "Control")))
} else {
  conditions <- factor(conditions)
}

cat("Conditions:", paste(levels(conditions), collapse = " vs "), "\n")
cat("Samples per condition:", paste(table(conditions), collapse = ", "), "\n")

# --- Run dmrseq ---
cat("\nRunning dmrseq...\n")

# Build design with batch as covariate if multiple batches exist
batches <- meta$batch
has_batches <- length(unique(batches)) > 1

tryCatch({
  if (has_batches) {
    cat("Including batch as covariate\n")
    pData(bs)$condition <- conditions
    pData(bs)$batch <- factor(batches)

    dmrs <- dmrseq(
      bs = bs,
      testCovariate = "condition",
      adjustCovariate = "batch",
      BPPARAM = BiocParallel::MulticoreParam(1)
    )
  } else {
    pData(bs)$condition <- conditions

    dmrs <- dmrseq(
      bs = bs,
      testCovariate = "condition",
      BPPARAM = BiocParallel::MulticoreParam(1)
    )
  }

  cat("Found", length(dmrs), "DMR candidates\n")

  # --- Write outputs ---
  if (length(dmrs) > 0) {
    # Extract fields from GRanges into a plain base data.frame
    # (as.data.frame on GRanges can return S4 DFrame which breaks write.table)
    mc <- GenomicRanges::mcols(dmrs)
    cat("GRanges mcols:", paste(colnames(mc), collapse = ", "), "\n")

    dmr_df <- data.frame(
      seqnames = as.character(GenomicRanges::seqnames(dmrs)),
      start    = as.integer(GenomicRanges::start(dmrs)),
      end      = as.integer(GenomicRanges::end(dmrs)),
      width    = as.integer(GenomicRanges::width(dmrs)),
      stringsAsFactors = FALSE
    )
    # Add metadata columns that exist
    mc_names <- colnames(mc)
    for (col_name in mc_names) {
      vals <- mc[[col_name]]
      if (is.numeric(vals) || is.integer(vals)) {
        dmr_df[[col_name]] <- as.numeric(vals)
      } else {
        dmr_df[[col_name]] <- as.character(vals)
      }
    }

    cat("DMR columns:", paste(colnames(dmr_df), collapse = ", "), "\n")

    # --- VDJ risk annotation ---
    # Human immune receptor loci (hg38) where somatic recombination creates
    # methylation artifacts that mimic DMRs in blood/PBMC samples.
    vdj_loci <- data.frame(
      chr   = c("chr14", "chr2",  "chr22", "chr14",    "chr7",  "chr7"),
      start = c(105586937, 88857361, 22026076, 21621904, 38240024, 142299011),
      end   = c(106879844, 90235368, 22922913, 22552132, 38368055, 142813287),
      locus = c("IGH",   "IGK",   "IGL",   "TRA_TRD",  "TRB",   "TRG"),
      stringsAsFactors = FALSE
    )

    dmr_df$vdj_risk <- FALSE
    dmr_df$vdj_locus <- NA_character_
    for (i in seq_len(nrow(vdj_loci))) {
      overlap <- dmr_df$seqnames == vdj_loci$chr[i] &
                 dmr_df$end >= vdj_loci$start[i] &
                 dmr_df$start <= vdj_loci$end[i]
      dmr_df$vdj_risk[overlap] <- TRUE
      dmr_df$vdj_locus[overlap] <- vdj_loci$locus[i]
    }
    n_vdj <- sum(dmr_df$vdj_risk)
    if (n_vdj > 0) cat("Flagged", n_vdj, "DMRs overlapping VDJ loci\n")

    # Full results
    write.table(dmr_df, file.path(opt$outdir, "window_results.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # BED file for significant DMRs (or top 20 by area stat if no p-values)
    # dmrseq uses "qval"/"pval"; handle both naming conventions
    n_top <- min(20, nrow(dmr_df))
    qcol <- intersect(c("qval", "qvalue"), colnames(dmr_df))[1]
    pcol <- intersect(c("pval", "pvalue"), colnames(dmr_df))[1]

    if (!is.na(qcol) && any(dmr_df[[qcol]] < 0.05, na.rm = TRUE)) {
      sig <- dmr_df[dmr_df[[qcol]] < 0.05, , drop = FALSE]
    } else if (!is.na(pcol)) {
      ord <- order(dmr_df[[pcol]])
      sig <- dmr_df[ord[seq_len(n_top)], , drop = FALSE]
    } else if ("areaStat" %in% colnames(dmr_df)) {
      ord <- order(abs(dmr_df$areaStat), decreasing = TRUE)
      sig <- dmr_df[ord[seq_len(n_top)], , drop = FALSE]
    } else {
      sig <- dmr_df[seq_len(n_top), , drop = FALSE]
    }

    bed_out <- data.frame(
      chr = sig$seqnames,
      start = sig$start,
      end = sig$end
    )
    write.table(bed_out, file.path(opt$outdir, "candidate_dmrs.bed"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    cat("Wrote", nrow(sig), "candidate DMRs to BED\n")

    # --- Manhattan plot ---
    if (!is.na(pcol)) {
      dmr_df$neg_log_p <- -log10(pmax(dmr_df[[pcol]], 1e-300))
    } else if ("areaStat" %in% colnames(dmr_df)) {
      dmr_df$neg_log_p <- abs(dmr_df$areaStat)
    } else {
      dmr_df$neg_log_p <- 1
    }
    dmr_df$midpoint <- (dmr_df$start + dmr_df$end) / 2
    dmr_df$chr_color <- factor(as.numeric(factor(dmr_df$seqnames)) %% 2)

    p <- ggplot(dmr_df, aes(x = midpoint, y = neg_log_p, color = chr_color)) +
      geom_point(size = 2, alpha = 0.6) +
      scale_color_manual(values = c("#1f77b4", "#aec7e8"), guide = "none") +
      facet_wrap(~seqnames, scales = "free_x", nrow = 1) +
      geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", alpha = 0.5) +
      labs(x = "Genomic Position", y = "-log10(p-value)",
           title = paste("DMR Manhattan Plot (dmrseq,", length(dmrs), "regions)")) +
      theme_minimal() +
      theme(axis.text.x = element_blank(), panel.spacing = unit(0.1, "lines"))

    ggsave(file.path(opt$outdir, "dmr_manhattan.png"), p,
           width = 14, height = 5, dpi = 150)

  } else {
    # No DMRs found — write empty outputs
    write.table(data.frame(), file.path(opt$outdir, "window_results.tsv"),
                sep = "\t", row.names = FALSE)
    write.table(data.frame(), file.path(opt$outdir, "candidate_dmrs.bed"),
                sep = "\t", row.names = FALSE, col.names = FALSE)
    png(file.path(opt$outdir, "dmr_manhattan.png"), width = 1400, height = 500, res = 150)
    plot.new(); title("DMR Manhattan Plot (no significant regions)")
    dev.off()
  }
}, error = function(e) {
  cat("dmrseq error:", conditionMessage(e), "\n")
  cat("Falling back to empty outputs\n")
  write.table(data.frame(), file.path(opt$outdir, "window_results.tsv"),
              sep = "\t", row.names = FALSE)
  write.table(data.frame(), file.path(opt$outdir, "candidate_dmrs.bed"),
              sep = "\t", row.names = FALSE, col.names = FALSE)
  png(file.path(opt$outdir, "dmr_manhattan.png"), width = 1400, height = 500, res = 150)
  plot.new(); title(paste("DMR detection failed:", conditionMessage(e)))
  dev.off()
})

cat("Region detection complete\n")
