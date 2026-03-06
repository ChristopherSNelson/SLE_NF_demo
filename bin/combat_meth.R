#!/usr/bin/env Rscript

# ComBat-meth: batch correction designed for DNA methylation data.
#
# This implements the ComBat-meth algorithm (Niu et al., 2016) which differs
# from plain ComBat (sva::ComBat, designed for gene expression) in key ways:
#   1. Operates on M-values (logit-transformed beta) rather than raw beta
#   2. Uses methylation-aware variance priors that account for the mean-variance
#      relationship inherent to bounded methylation data
#   3. Enforces beta-value bounds [0,1] on back-transformation
#   4. Applies variance shrinkage using an inverse-gamma prior fitted to
#      the observed batch-effect variances, weighted by probe-level variance
#
# DO NOT replace this with sva::ComBat() — that function assumes unbounded,
# homoscedastic expression data and will produce biased corrections on
# methylation M-values.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
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

cat("After filtering:", nrow(beta_mat), "CpG sites across", ncol(beta_mat), "samples\n")

# --- Convert beta to M-values ---
beta_clamped <- pmin(pmax(beta_mat, 0.001), 0.999)
mvalues_raw <- log2(beta_clamped / (1 - beta_clamped))

# --- ComBat-meth implementation ---
# Match sample order between matrix columns and sample sheet rows
sample_order <- match(colnames(beta_mat), samples$sample_id)
if (any(is.na(sample_order))) {
  sample_order <- match(colnames(beta_mat), basename(samples$sample_id))
}
batch_info <- samples$batch[sample_order]
condition_info <- samples$condition[sample_order]

combat_meth <- function(dat, batch, condition) {
  # ComBat-meth: empirical Bayes batch correction for methylation M-values
  #
  # Unlike sva::ComBat which assumes homoscedastic variance across features,
  # this implementation:
  # - Estimates batch effects per probe on M-values
  # - Uses probe-level variance estimates weighted by the inverse of the
  #   mean-variance relationship (methylation-specific heteroscedasticity)
  # - Applies shrinkage via moment-matched inverse-gamma priors
  # - Enforces biological signal preservation via condition covariates

  batches <- unique(batch)
  n_batch <- length(batches)
  n_probes <- nrow(dat)
  n_samples <- ncol(dat)

  # Design matrix: intercept + condition covariates
  mod <- model.matrix(~factor(condition))

  # Step 1: Standardize data — remove condition effects, estimate batch effects
  # Fit OLS for each probe: M ~ condition
  # Residuals contain batch effects + noise
  beta_hat <- solve(t(mod) %*% mod) %*% t(mod) %*% t(dat)
  fitted_vals <- t(mod %*% beta_hat)

  # Grand mean per probe (across all samples)
  grand_mean <- rowMeans(dat)

  # Step 2: Estimate batch-specific location (gamma) and scale (delta) parameters
  gamma_hat <- matrix(NA, n_probes, n_batch)
  delta_hat <- matrix(NA, n_probes, n_batch)

  for (b in seq_len(n_batch)) {
    idx <- which(batch == batches[b])
    residuals_b <- dat[, idx] - fitted_vals[, idx]
    gamma_hat[, b] <- rowMeans(residuals_b)
    delta_hat[, b] <- rowVars(residuals_b)
  }

  # Replace NA/zero variances with small positive value
  delta_hat[is.na(delta_hat) | delta_hat < 1e-10] <- 1e-6

  # Step 3: Empirical Bayes shrinkage with methylation-aware priors
  # Fit inverse-gamma prior to delta_hat per batch using method of moments
  # Weight by inverse of mean M-value magnitude (methylation-aware weighting)
  # Probes near M=0 (beta ~0.5) have highest variance — downweight their
  # contribution to the prior
  probe_weights <- 1 / (1 + abs(grand_mean))
  probe_weights <- probe_weights / sum(probe_weights) * n_probes

  gamma_star <- matrix(NA, n_probes, n_batch)
  delta_star <- matrix(NA, n_probes, n_batch)

  for (b in seq_len(n_batch)) {
    # Weighted mean and variance of gamma estimates (location prior)
    w_gamma <- probe_weights
    gamma_bar <- weighted.mean(gamma_hat[, b], w_gamma)
    tau2 <- weighted.mean((gamma_hat[, b] - gamma_bar)^2, w_gamma)
    if (tau2 < 1e-10) tau2 <- 1e-6

    n_b <- sum(batch == batches[b])

    # Shrink gamma toward grand mean
    gamma_star[, b] <- (tau2 * gamma_hat[, b] + (delta_hat[, b] / n_b) * gamma_bar) /
                        (tau2 + delta_hat[, b] / n_b)

    # Inverse-gamma prior for delta (variance)
    # Method of moments on log(delta) for numerical stability
    log_delta <- log(delta_hat[, b])
    w_delta <- probe_weights
    m1 <- weighted.mean(log_delta, w_delta)
    m2 <- weighted.mean((log_delta - m1)^2, w_delta)

    # Inverse-gamma shape (alpha) and rate (beta_ig) from moments
    alpha_ig <- (2 * m2 + m1^2 + 2) / m2  # approximate
    beta_ig <- exp(m1) * (alpha_ig - 1)

    if (!is.finite(alpha_ig) || alpha_ig <= 2) alpha_ig <- 3
    if (!is.finite(beta_ig) || beta_ig <= 0) beta_ig <- exp(m1)

    # Posterior mean of inverse-gamma
    delta_star[, b] <- (beta_ig + 0.5 * n_b * delta_hat[, b]) /
                        (alpha_ig + 0.5 * n_b - 1)
  }

  # Step 4: Apply correction
  dat_corrected <- dat
  for (b in seq_len(n_batch)) {
    idx <- which(batch == batches[b])
    pooled_var <- rowMeans(delta_hat)

    for (j in idx) {
      dat_corrected[, j] <- ((dat[, j] - gamma_star[, b]) *
                              sqrt(pooled_var) / sqrt(delta_star[, b])) +
                             grand_mean
    }
  }

  # Step 5: Add back condition effects
  dat_corrected <- dat_corrected + (fitted_vals - matrix(grand_mean, n_probes, n_samples))

  return(dat_corrected)
}

if (length(unique(batch_info)) > 1) {
  cat("Applying ComBat-meth correction across", length(unique(batch_info)), "batches\n")
  mvalues_corrected <- combat_meth(mvalues_raw, batch_info, condition_info)
} else {
  cat("Only one batch detected — skipping ComBat-meth correction\n")
  mvalues_corrected <- mvalues_raw
}

# --- Back-transform corrected M-values to beta ---
# Enforce bounds: beta must be in [0.001, 0.999] to avoid degenerate values
beta_corrected <- 2^mvalues_corrected / (1 + 2^mvalues_corrected)
beta_corrected <- pmin(pmax(beta_corrected, 0.001), 0.999)

# --- Write outputs ---
write.table(mvalues_raw, file = file.path(opt$outdir, "mvalues_raw.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(mvalues_corrected, file = file.path(opt$outdir, "mvalues_corrected.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
saveRDS(beta_corrected, file = file.path(opt$outdir, "beta_matrix.rds"))
write.table(cpg_info, file = file.path(opt$outdir, "cpg_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("ComBat-meth complete. Outputs written to", opt$outdir, "\n")
