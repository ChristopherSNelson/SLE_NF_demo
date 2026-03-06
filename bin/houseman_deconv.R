#!/usr/bin/env Rscript

# Houseman cell-type deconvolution using quadprog (NOT minfi)
# Constrained least-squares projection of sample methylation profiles
# onto a reference panel of blood cell-type methylation signatures.

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
  make_option("--outdir", type = "character", default = ".", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Load data ---
beta_mat <- readRDS(opt$beta_matrix)
samples <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE)

cat("Beta matrix:", nrow(beta_mat), "CpGs x", ncol(beta_mat), "samples\n")

# --- Build simulated reference panel ---
# In production, this would use Reinius et al. or FlowSorted.Blood.450k profiles.
# For the demo, we generate a synthetic reference panel with known cell-type signatures.
set.seed(42)
cell_types <- c("CD4T", "CD8T", "Bcell", "NK", "Mono", "Gran")
n_cpgs <- nrow(beta_mat)

# Generate reference profiles with distinct methylation patterns per cell type
ref_panel <- matrix(NA, nrow = n_cpgs, ncol = length(cell_types))
colnames(ref_panel) <- cell_types
rownames(ref_panel) <- rownames(beta_mat)

for (i in seq_along(cell_types)) {
  # Each cell type gets a base profile shifted by its index
  base <- runif(n_cpgs, 0.2, 0.8)
  # Add cell-type-specific signal at different CpG subsets
  sig_start <- ((i - 1) * n_cpgs / length(cell_types)) + 1
  sig_end <- min(i * n_cpgs / length(cell_types), n_cpgs)
  sig_idx <- seq(sig_start, sig_end)
  base[sig_idx] <- base[sig_idx] + runif(length(sig_idx), 0.1, 0.3)
  ref_panel[, i] <- pmin(pmax(base, 0.01), 0.99)
}

# --- Houseman projection via quadprog ---
# For each sample, solve: min ||y - X*w||^2 s.t. w >= 0, sum(w) = 1
# This is a quadratic programming problem.

houseman_deconv <- function(y, ref) {
  # Use subset of CpGs present in both
  common <- intersect(names(y), rownames(ref))
  if (length(common) < 100) {
    warning("Fewer than 100 overlapping CpGs")
    return(rep(NA, ncol(ref)))
  }
  y_sub <- y[common]
  ref_sub <- ref[common, ]

  # Remove rows with NA
  valid <- complete.cases(cbind(y_sub, ref_sub))
  y_sub <- y_sub[valid]
  ref_sub <- ref_sub[valid, ]

  n_ct <- ncol(ref_sub)

  # QP setup: min 0.5 * w' D w - d' w
  # D = 2 * X'X, d = 2 * X'y
  Dmat <- 2 * crossprod(ref_sub)
  dvec <- 2 * crossprod(ref_sub, y_sub)

  # Constraints: w >= 0 and sum(w) = 1
  # Amat columns: first column is equality constraint (sum = 1),
  # remaining are identity (w_i >= 0)
  Amat <- cbind(rep(1, n_ct), diag(n_ct))
  bvec <- c(1, rep(0, n_ct))
  meq <- 1  # first constraint is equality

  # Add small ridge to ensure positive definiteness
  Dmat <- Dmat + diag(1e-6, n_ct)

  result <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = meq),
    error = function(e) {
      warning("QP failed for sample: ", e$message)
      return(NULL)
    }
  )

  if (is.null(result)) return(rep(NA, n_ct))
  w <- result$solution
  names(w) <- colnames(ref_sub)
  pmax(w, 0)  # clip tiny negatives from numerical noise
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

# Merge condition info
frac_long <- merge(frac_long, samples[, c("sample_id", "condition")], by = "sample_id")

p <- ggplot(frac_long, aes(x = sample_id, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~condition, scales = "free_x", space = "free_x") +
  labs(title = "Estimated Cell-Type Fractions", x = "Sample", y = "Fraction", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(opt$outdir, "cell_fractions.png"), p, width = 10, height = 6, dpi = 150)

cat("Houseman deconvolution complete\n")
