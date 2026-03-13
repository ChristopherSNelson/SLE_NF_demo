#!/usr/bin/env Rscript
# Plot all dmrseq regions regardless of significance, colored by direction
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
input_file  <- args[1]
output_file <- args[2]
title_text  <- if (length(args) > 2) args[3] else "SLE Differential Methylation (dmrseq)"

df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Compute position and significance
df <- df %>%
  mutate(
    midpoint_mb = ((start + end) / 2) / 1e6,
    neg_log_p   = -log10(pmax(pval, 1e-300)),
    direction   = ifelse(beta > 0, "Hyper (SLE)", "Hypo (SLE)"),
    sig_tier    = case_when(
      qval < 0.05 ~ "FDR < 0.05",
      pval < 0.05 ~ "p < 0.05",
      TRUE        ~ "ns"
    )
  )

# Top 5 peaks for labeling
top5 <- df %>% arrange(pval) %>% head(5) %>%
  mutate(label = paste0(round(midpoint_mb, 2), " Mb"))

p <- ggplot(df, aes(x = midpoint_mb, y = neg_log_p)) +
  # All points, colored by methylation direction
  geom_point(aes(color = direction, alpha = sig_tier), size = 2.2) +
  scale_color_manual(
    name   = "Direction",
    values = c("Hyper (SLE)" = "#e41a1c", "Hypo (SLE)" = "#377eb8")
  ) +
  scale_alpha_manual(
    name   = "Significance",
    values = c("FDR < 0.05" = 1.0, "p < 0.05" = 0.75, "ns" = 0.35),
    guide  = guide_legend(override.aes = list(size = 3))
  ) +
  # Nominal p = 0.05 threshold line
  geom_hline(yintercept = -log10(0.05), color = "gray40",
             linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = max(df$midpoint_mb) * 0.98,
           y = -log10(0.05) + 0.05, label = "p = 0.05",
           hjust = 1, size = 3, color = "gray40") +
  # Label top 5 peaks
  geom_text(data = top5, aes(label = label),
            vjust = -0.6, size = 3, fontface = "bold", color = "black") +
  labs(
    title    = title_text,
    subtitle = paste0(nrow(df), " tested regions  |  n = 6 samples (3 SLE, 3 Control)"),
    x        = "Genomic Position on chr19 (Mb)",
    y        = expression(-log[10](p-value))
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position   = "right",
    plot.subtitle     = element_text(size = 10, color = "gray40")
  )

ggsave(output_file, p, width = 12, height = 5, dpi = 150)
cat("Manhattan plot saved to:", output_file, "\n")
