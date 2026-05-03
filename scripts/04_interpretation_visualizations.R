#!/usr/bin/env Rscript

# 04_interpretation_visualizations.R
# Purpose: Generate QC, diagnostic, and interpretation figures using config-driven settings.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(yaml)
})

# Load config
config <- yaml::read_yaml("config/config.yaml")

top_n_heatmap <- config$top_n_heatmap
plot_width <- config$plot_width
plot_height <- config$plot_height
dpi <- config$dpi

# Output directories
dir.create("results/figures/png", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/pdf", recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width = plot_width, height = plot_height) {
  ggsave(
    filename = file.path("results/figures/png", paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  ggsave(
    filename = file.path("results/figures/pdf", paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
}

message("Loading DESeq2 objects...")

dds <- readRDS("data/processed/deseq2_dds.rds")
res <- readRDS("data/processed/deseq2_results.rds")

vsd <- vst(dds, blind = FALSE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "dex", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = dex)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of airway RNA-seq samples",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  )

save_plot(p, "pca_samples", width = 7, height = 5)

# MA plot
png("results/figures/png/ma_plot.png", width = 1000, height = 800, res = dpi)
plotMA(res, ylim = c(-5, 5), main = "DESeq2 MA plot")
dev.off()

pdf("results/figures/pdf/ma_plot.pdf", width = 7, height = 5)
plotMA(res, ylim = c(-5, 5), main = "DESeq2 MA plot")
dev.off()

# Dispersion plot
png("results/figures/png/dispersion_plot.png", width = 1000, height = 800, res = dpi)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
dev.off()

pdf("results/figures/pdf/dispersion_plot.pdf", width = 7, height = 5)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
dev.off()

# Adjusted p-value histogram
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]

p <- ggplot(res_df, aes(x = padj)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(
    title = "Distribution of adjusted p-values",
    x = "Adjusted p-value",
    y = "Gene count"
  )

save_plot(p, "padj_histogram", width = 7, height = 5)

# Size factors
sf_df <- data.frame(
  sample = colnames(dds),
  size_factor = sizeFactors(dds),
  dex = as.factor(dds$dex)
)

p <- ggplot(sf_df, aes(x = sample, y = size_factor, fill = dex)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "DESeq2 size factors",
    x = "Sample",
    y = "Size factor"
  )

save_plot(p, "size_factors", width = 7, height = 5)

# Top up/downregulated genes
top_up <- res_df[order(-res_df$log2FoldChange), ][1:15, ]
top_down <- res_df[order(res_df$log2FoldChange), ][1:15, ]

top_genes <- rbind(
  data.frame(top_up, direction = "Upregulated"),
  data.frame(top_down, direction = "Downregulated")
)

p <- ggplot(
  top_genes,
  aes(
    x = reorder(gene_id, log2FoldChange),
    y = log2FoldChange,
    fill = direction
  )
) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top up- and downregulated genes",
    x = "Gene",
    y = "log2 fold change"
  )

save_plot(p, "top_up_down_genes", width = 8, height = 7)

# Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

png("results/figures/png/sample_distance_heatmap.png", width = 1000, height = 900, res = dpi)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample-to-sample distance heatmap"
)
dev.off()

pdf("results/figures/pdf/sample_distance_heatmap.pdf", width = 8, height = 7)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  main = "Sample-to-sample distance heatmap"
)
dev.off()

# Top significant genes heatmap
top_genes_heatmap <- head(res_df[order(res_df$padj), "gene_id"], top_n_heatmap)

heatmap_matrix <- assay(vsd)[top_genes_heatmap, ]
heatmap_matrix <- heatmap_matrix - rowMeans(heatmap_matrix)

annotation_col <- as.data.frame(colData(dds)[, "dex", drop = FALSE])

png("results/figures/png/heatmap_top50.png", width = 1000, height = 1200, res = dpi)
pheatmap(
  heatmap_matrix,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = paste("Top", top_n_heatmap, "differentially expressed genes")
)
dev.off()

pdf("results/figures/pdf/heatmap_top50.pdf", width = 8, height = 10)
pheatmap(
  heatmap_matrix,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  main = paste("Top", top_n_heatmap, "differentially expressed genes")
)
dev.off()

# GO enrichment dot plot
if (file.exists("results/tables/go_enrichment.csv")) {
  go <- read.csv("results/tables/go_enrichment.csv")

  if (nrow(go) > 0) {
    go_top <- go[order(go$Adjusted.P.value), ][1:min(15, nrow(go)), ]
    go_top$minus_log10_padj <- -log10(go_top$Adjusted.P.value)
    go_top$overlap_count <- as.numeric(sub("/.*", "", go_top$Overlap))

    p <- ggplot(
      go_top,
      aes(
        x = minus_log10_padj,
        y = reorder(Term, minus_log10_padj),
        size = overlap_count
      )
    ) +
      geom_point() +
      theme_minimal() +
      labs(
        title = "Top enriched GO Biological Processes",
        x = "-log10 adjusted p-value",
        y = "GO term",
        size = "Gene overlap"
      )

    save_plot(p, "go_dotplot_overlap", width = 9, height = 6)
  }
}

message("Interpretation visualizations completed.")
