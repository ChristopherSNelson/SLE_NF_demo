#!/usr/bin/env Rscript

# PCA visualization of methylation M-values
# Produces 4 plots: {raw, corrected} x {colored by batch, colored by condition}

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--mvalues_raw", type = "character", help = "Raw M-values TSV"),
  make_option("--mvalues_corrected", type = "character", help = "Corrected M-values TSV"),
  make_option("--sample_sheet", type = "character", help = "Sample sheet CSV"),
  make_option("--outdir", type = "character", default = ".", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Load data ---
samples <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE)
mval_raw <- as.matrix(read.table(opt$mvalues_raw, sep = "\t", header = TRUE,
                                  row.names = 1, check.names = FALSE))
mval_corr <- as.matrix(read.table(opt$mvalues_corrected, sep = "\t", header = TRUE,
                                   row.names = 1, check.names = FALSE))

# --- PCA function ---
run_pca <- function(mat) {
  # Use top 5000 most variable sites for speed
  n_sites <- min(5000, nrow(mat))
  vars <- apply(mat, 1, var, na.rm = TRUE)
  top_idx <- order(vars, decreasing = TRUE)[1:n_sites]
  pca <- prcomp(t(mat[top_idx, ]), center = TRUE, scale. = TRUE)
  pca
}

# --- Plot function ---
make_pca_plot <- function(pca, samples, color_col, title, filename) {
  var_explained <- summary(pca)$importance[2, 1:2] * 100
  df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  # Match sample metadata to PCA sample order
  sample_order <- match(rownames(pca$x), samples$sample_id)
  df[[color_col]] <- samples[[color_col]][sample_order]

  p <- ggplot(df, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
    geom_point(size = 3) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme_minimal() +
    theme(legend.position = "right")

  ggsave(filename, p, width = 7, height = 5, dpi = 150)
  var_explained
}

# --- Run PCA ---
pca_raw <- run_pca(mval_raw)
pca_corr <- run_pca(mval_corr)

# --- Generate 4 plots ---
v1 <- make_pca_plot(pca_raw, samples, "batch", "PCA: Raw M-values (batch)",
                    file.path(opt$outdir, "pca_raw_batch.png"))
v2 <- make_pca_plot(pca_raw, samples, "condition", "PCA: Raw M-values (condition)",
                    file.path(opt$outdir, "pca_raw_condition.png"))
v3 <- make_pca_plot(pca_corr, samples, "batch", "PCA: Corrected M-values (batch)",
                    file.path(opt$outdir, "pca_corrected_batch.png"))
v4 <- make_pca_plot(pca_corr, samples, "condition", "PCA: Corrected M-values (condition)",
                    file.path(opt$outdir, "pca_corrected_condition.png"))

# --- Variance table ---
var_df <- data.frame(
  dataset = c("raw", "raw", "corrected", "corrected"),
  PC = rep(c("PC1", "PC2"), 2),
  variance_pct = c(v1[1], v1[2], v3[1], v3[2])
)
write.table(var_df, file = file.path(opt$outdir, "pca_variance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("PCA plots complete\n")
