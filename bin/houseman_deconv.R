#!/usr/bin/env Rscript

# Houseman cell-type deconvolution using quadprog (NOT minfi)
# Constrained least-squares projection of sample methylation profiles
# onto a reference panel of blood cell-type methylation signatures.
# Reference: Salas et al. 2018 IDOL library (WGBS-optimized, 450 CpGs, 6 blood cell types)

suppressPackageStartupMessages({
  library(optparse)
  library(quadprog)
  library(ggplot2)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--beta_matrix", type = "character", help = "Beta matrix RDS file"),
  make_option("--sample_sheet", type = "character", help = "Sample sheet CSV"),
  make_option("--ref_panel",   type = "character", help = "Salas 2018 IDOL reference panel CSV"),
  make_option("--outdir",      type = "character", default = ".", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Load data ---
beta_mat <- readRDS(opt$beta_matrix)
samples   <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE)

cat("Beta matrix:", nrow(beta_mat), "CpGs x", ncol(beta_mat), "samples\n")

# --- Load Salas 2018 IDOL reference panel ---
ref_raw <- fread(opt$ref_panel, skip = 1)  # skip title row
# Columns: V1=probe_id, chr, start, end, ..., CD8T, CD4T, NK, Bcell, Mono, Neu
cell_types <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
ref_raw[, cpg_id := paste(chr, start, sep = ":")]
ref_panel <- as.matrix(ref_raw[, ..cell_types])
rownames(ref_panel) <- ref_raw$cpg_id

cat("Reference panel:", nrow(ref_panel), "CpGs x", ncol(ref_panel), "cell types\n")
cat("Overlapping CpGs with beta matrix:", sum(rownames(ref_panel) %in% rownames(beta_mat)), "\n")

# --- Houseman projection via quadprog ---
# For each sample, solve: min ||y - X*w||^2 s.t. w >= 0, sum(w) = 1
houseman_deconv <- function(y, ref) {
  common <- intersect(names(y), rownames(ref))
  if (length(common) < 50) {
    warning("Fewer than 50 overlapping CpGs — fractions unreliable")
    return(rep(NA, ncol(ref)))
  }
  y_sub   <- y[common]
  ref_sub <- ref[common, ]

  valid   <- complete.cases(cbind(y_sub, ref_sub))
  y_sub   <- y_sub[valid]
  ref_sub <- ref_sub[valid, ]

  n_ct <- ncol(ref_sub)
  Dmat <- 2 * crossprod(ref_sub) + diag(1e-6, n_ct)
  dvec <- 2 * crossprod(ref_sub, y_sub)
  Amat <- cbind(rep(1, n_ct), diag(n_ct))
  bvec <- c(1, rep(0, n_ct))

  result <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
    error = function(e) { warning("QP failed: ", e$message); NULL }
  )

  if (is.null(result)) return(rep(NA, n_ct))
  w <- pmax(result$solution, 0)
  names(w) <- colnames(ref_sub)
  w
}

# --- Run deconvolution for each sample ---
fractions <- matrix(NA, nrow = ncol(beta_mat), ncol = length(cell_types))
rownames(fractions) <- colnames(beta_mat)
colnames(fractions) <- cell_types

for (j in seq_len(ncol(beta_mat))) {
  y <- beta_mat[, j]
  names(y) <- rownames(beta_mat)
  fractions[j, ] <- houseman_deconv(y, ref_panel)
}

# --- Write output ---
frac_df <- as.data.frame(fractions)
frac_df$sample_id <- rownames(fractions)
frac_df <- frac_df[, c("sample_id", cell_types)]

write.table(frac_df, file = file.path(opt$outdir, "cell_fractions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Stacked bar plot ---
frac_long <- tidyr::pivot_longer(frac_df, cols = all_of(cell_types),
                                  names_to = "cell_type", values_to = "fraction")
frac_long <- merge(frac_long, samples[, c("sample_id", "condition")], by = "sample_id")

p <- ggplot(frac_long, aes(x = sample_id, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~condition, scales = "free_x", space = "free_x") +
  labs(title = "Estimated Cell-Type Fractions (Salas 2018 IDOL reference)",
       x = "Sample", y = "Fraction", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(opt$outdir, "cell_fractions.png"), p, width = 10, height = 6, dpi = 150)

cat("Houseman deconvolution complete\n")
