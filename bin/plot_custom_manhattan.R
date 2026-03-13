#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Use basic commandArgs
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]
plot_title <- if(length(args) > 2) args[3] else "SLE Differential Methylation"

if (is.na(input_file) || is.na(output_dir)) {
    cat("Usage: Rscript bin/plot_custom_manhattan.R <input_file> <output_dir> [title]\n")
    quit(status=1)
}

# Load data
df <- read.table(input_file, header=TRUE, sep="\t")

# 1. Transform coordinates to Megabases (Mb)
df <- df %>% 
  mutate(midpoint_mb = ((start + end) / 2) / 1e6) %>%
  mutate(neg_log_p = -log10(pmax(pval, 1e-10)))

# 2. Identify top peaks for labeling
top_hits <- df %>% arrange(pval) %>% head(5)
top_hits$label <- paste0(round(top_hits$midpoint_mb, 2), " Mb")

# 3. Create the Plot (REGULAR STYLE)
# Using standard blue points, no complex scales or transparency
p <- ggplot(df, aes(x = midpoint_mb, y = neg_log_p)) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_text(data = top_hits, aes(label = label), vjust = -1, size = 3.5) +
  labs(title = plot_title,
       x = "Genomic Position (Mb)",
       y = "-log10 (p-value)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

# Save the plot
output_path <- file.path(output_dir, "dmr_manhattan_presentation.png")
ggsave(output_path, p, width = 10, height = 5, dpi = 150)
cat("Regular style Manhattan plot saved to:", output_path, "\n")
