#!/usr/bin/env Rscript

# ComBat-meth: batch correction for methylation data
# Operates on M-values (logit-transformed beta) to respect the bounded nature of methylation data.
# Unlike plain ComBat (designed for gene expression), ComBat-meth applies an empirical Bayes
# framework on M-values and back-transforms to produce corrected beta values.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(sva)
  library(matrixStats)
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

# Parse each bedGraph into a data.table with chr, start, end, beta
read_bedgraph <- function(f) {
  dt <- fread(f, skip = 1, header = FALSE,
              col.names = c("chr", "start", "end", "beta", "count_m", "count_u"))
  # Extract sample ID from filename
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

# Remove invariant sites
row_vars <- rowVars(beta_mat)
variable_rows <- row_vars > 1e-6
beta_mat <- beta_mat[variable_rows, ]
cpg_info <- cpg_info[variable_rows, ]

cat("After filtering:", nrow(beta_mat), "CpG sites across", ncol(beta_mat), "samples\n")

# --- Convert beta to M-values ---
# Clamp beta to avoid log(0) or log(Inf)
beta_clamped <- pmin(pmax(beta_mat, 0.001), 0.999)
mvalues_raw <- log2(beta_clamped / (1 - beta_clamped))

# --- ComBat-meth: batch correction on M-values ---
# Match sample order between matrix and sample sheet
sample_order <- match(colnames(beta_mat), samples$sample_id)
if (any(is.na(sample_order))) {
  # Try matching without path prefix
  sample_order <- match(colnames(beta_mat), basename(samples$sample_id))
}
batch_info <- samples$batch[sample_order]

# Design matrix with condition as covariate to preserve biological signal
mod <- model.matrix(~condition, data = samples[sample_order, ])

if (length(unique(batch_info)) > 1) {
  cat("Applying ComBat-meth correction across", length(unique(batch_info)), "batches\n")
  mvalues_corrected <- ComBat(
    dat = mvalues_raw,
    batch = batch_info,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
  )
} else {
  cat("Only one batch detected — skipping ComBat-meth correction\n")
  mvalues_corrected <- mvalues_raw
}

# --- Back-transform corrected M-values to beta ---
beta_corrected <- 2^mvalues_corrected / (1 + 2^mvalues_corrected)

# --- Write outputs ---
write.table(mvalues_raw, file = file.path(opt$outdir, "mvalues_raw.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(mvalues_corrected, file = file.path(opt$outdir, "mvalues_corrected.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
saveRDS(beta_corrected, file = file.path(opt$outdir, "beta_matrix.rds"))
write.table(cpg_info, file = file.path(opt$outdir, "cpg_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("ComBat-meth complete. Outputs written to", opt$outdir, "\n")
