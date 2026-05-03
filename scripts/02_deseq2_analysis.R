#!/usr/bin/env Rscript

# 02_deseq2_analysis.R
# Purpose: Perform differential expression analysis using DESeq2 with configurable parameters.

suppressPackageStartupMessages({
  library(DESeq2)
  library(airway)
  library(yaml)
})

# Load config
config <- yaml::read_yaml("config/config.yaml")

padj_cutoff <- config$padj_cutoff
log2fc_cutoff <- config$log2fc_cutoff

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("Loading airway dataset...")

data("airway")

dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

dds$dex <- relevel(dds$dex, ref = "untrt")

message("Running DESeq2 model...")

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("dex", "trt", "untrt"),
  alpha = padj_cutoff
)

# Convert to dataframe
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Order by significance
res_df <- res_df[order(res_df$padj), ]

# Filter significant genes
sig_df <- subset(
  res_df,
  !is.na(padj) & padj < padj_cutoff
)

# Upregulated
up_df <- subset(
  sig_df,
  log2FoldChange > log2fc_cutoff
)

# Downregulated
down_df <- subset(
  sig_df,
  log2FoldChange < -log2fc_cutoff
)

# Save objects
saveRDS(dds, "data/processed/deseq2_dds.rds")
saveRDS(res, "data/processed/deseq2_results.rds")

# Save tables
write.csv(res_df, "results/tables/deseq2_results_all.csv", row.names = FALSE)
write.csv(sig_df, "results/tables/deseq2_results_significant.csv", row.names = FALSE)
write.csv(up_df, "results/tables/deseq2_results_upregulated.csv", row.names = FALSE)
write.csv(down_df, "results/tables/deseq2_results_downregulated.csv", row.names = FALSE)

# Summary
message("DESeq2 analysis completed.")
message("Total genes tested: ", nrow(res_df))
message("Significant genes (padj < ", padj_cutoff, "): ", nrow(sig_df))
message("Upregulated genes: ", nrow(up_df))
message("Downregulated genes: ", nrow(down_df))
