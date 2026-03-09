#!/usr/bin/env Rscript

# Batch correction for DNA methylation data using ComBatMet
# (Wang, 2025, NAR Genomics and Bioinformatics)
# https://github.com/JmWangBio/ComBatMet
#
# ComBatMet uses beta regression to estimate batch-free distributions and
# realigns quantiles to their corrected counterparts. It operates directly
# on beta values (not M-values), which is the correct approach for bounded
# methylation data.
#
# DO NOT replace with sva::ComBat() — that assumes unbounded expression data.

# Install ComBatMet from GitHub if not already available
if (!requireNamespace("ComBatMet", quietly = TRUE)) {
  cat("Installing ComBatMet from GitHub...\n")

  # Fix gcc version mismatch in conda: R may expect a different minor version
  # of gfortran libs than what's installed. Create a symlink if needed.
  r_home <- R.home()
  gcc_dir <- file.path(dirname(dirname(r_home)), "lib", "gcc", "arm64-apple-darwin20.0.0")
  if (dir.exists(gcc_dir)) {
    existing <- list.dirs(gcc_dir, recursive = FALSE, full.names = FALSE)
    makeconf <- readLines(file.path(r_home, "etc", "Makeconf"), warn = FALSE)
    expected <- regmatches(makeconf, regexpr("gcc/arm64-apple-darwin20\\.0\\.0/[0-9.]+", makeconf))
    if (length(expected) > 0) {
      expected_ver <- sub(".*/", "", expected[1])
      if (!expected_ver %in% existing && length(existing) > 0) {
        target <- file.path(gcc_dir, expected_ver)
        source <- file.path(gcc_dir, existing[1])
        cat(sprintf("Symlinking gcc %s -> %s for Fortran libs\n", expected_ver, existing[1]))
        file.symlink(source, target)
      }
    }
  }

  remotes::install_github("JmWangBio/ComBatMet", upgrade = "never", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(matrixStats)
  library(ComBatMet)
})

option_list <- list(
  make_option("--bedgraph_dir", type = "character", help = "Directory containing bedGraph files"),
  make_option("--sample_sheet", type = "character", help = "Sample sheet CSV"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Read sample sheet ---
samples <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE)
cat("Loaded sample sheet with", nrow(samples), "samples\n")

# --- Read bedGraph files and build beta matrix ---
bedgraph_files <- list.files(opt$bedgraph_dir, pattern = "_CpG\\.bedGraph$", full.names = TRUE)
cat("Found", length(bedgraph_files), "bedGraph files\n")

read_bedgraph <- function(f) {
  # MethylDackel bedGraphs have a track header line; skip it
  dt <- fread(f, skip = 1, header = FALSE,
              col.names = c("chr", "start", "end", "beta", "count_m", "count_u"))
  sid <- sub("_CpG\\.bedGraph$", "", basename(f))
  dt$sample_id <- sid
  dt$beta <- dt$count_m / (dt$count_m + dt$count_u)
  dt[, .(chr, start, end, beta, sample_id)]
}

all_bg <- rbindlist(lapply(bedgraph_files, read_bedgraph))

# Pivot to matrix: rows = CpG sites, columns = samples
all_bg[, cpg_id := paste(chr, start, sep = ":")]
beta_wide <- dcast(all_bg, cpg_id ~ sample_id, value.var = "beta")

# Extract CpG manifest
cpg_ids <- beta_wide$cpg_id
cpg_info <- data.frame(
  cpg_id = cpg_ids,
  chr = sub(":.*", "", cpg_ids),
  start = as.integer(sub(".*:", "", cpg_ids)),
  stringsAsFactors = FALSE
)

# Build numeric matrix
beta_mat <- as.matrix(beta_wide[, -1, with = FALSE])
rownames(beta_mat) <- cpg_ids

# Filter: remove sites with any NA or invariant sites
complete_rows <- complete.cases(beta_mat)
beta_mat <- beta_mat[complete_rows, ]
cpg_info <- cpg_info[complete_rows, ]

row_vars <- rowVars(beta_mat)
variable_rows <- row_vars > 1e-6
beta_mat <- beta_mat[variable_rows, ]
cpg_info <- cpg_info[variable_rows, ]

# Clamp beta to strict (0,1) — ComBatMet requires beta in open interval
beta_mat <- pmin(pmax(beta_mat, 0.001), 0.999)

cat("After filtering:", nrow(beta_mat), "CpG sites across", ncol(beta_mat), "samples\n")

# --- Raw M-values (for PCA comparison) ---
mvalues_raw <- log2(beta_mat / (1 - beta_mat))

# --- ComBatMet batch correction ---
# Match sample order between matrix columns and sample sheet rows
sample_order <- match(colnames(beta_mat), samples$sample_id)
if (any(is.na(sample_order))) {
  sample_order <- match(colnames(beta_mat), basename(samples$sample_id))
}
batch_info <- samples$batch[sample_order]

# Encode condition as numeric group vector for ComBat_met
condition_info <- samples$condition[sample_order]
group_numeric <- as.integer(factor(condition_info)) - 1L  # 0/1 encoding

if (length(unique(batch_info)) > 1) {
  cat("Applying ComBatMet correction across", length(unique(batch_info)), "batches\n")
  cat("  Batches:", paste(unique(batch_info), collapse = ", "), "\n")
  cat("  Groups:", paste(unique(condition_info), collapse = ", "), "\n")

  beta_corrected <- ComBat_met(
    beta_mat,
    batch = batch_info,
    group = group_numeric,
    full_mod = TRUE
  )

  # Enforce bounds after correction
  beta_corrected <- pmin(pmax(beta_corrected, 0.001), 0.999)
} else {
  cat("Only one batch detected — skipping ComBatMet correction\n")
  beta_corrected <- beta_mat
}

# --- Convert corrected beta to M-values ---
mvalues_corrected <- log2(beta_corrected / (1 - beta_corrected))

# --- Write outputs ---
write.table(mvalues_raw, file = file.path(opt$outdir, "mvalues_raw.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(mvalues_corrected, file = file.path(opt$outdir, "mvalues_corrected.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
saveRDS(beta_corrected, file = file.path(opt$outdir, "beta_matrix.rds"))
write.table(cpg_info, file = file.path(opt$outdir, "cpg_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("ComBatMet complete. Outputs written to", opt$outdir, "\n")
