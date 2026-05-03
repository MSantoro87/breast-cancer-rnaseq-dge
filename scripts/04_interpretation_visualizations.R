# 04_interpretation_visualizations.R
#
# Purpose:
#   Generate interpretation-focused visualizations for the RNA-seq analysis.
#   These plots support quality control, sample-level interpretation,
#   differential expression diagnostics, and functional enrichment reporting.
#
# Input:
#   - Airway RNA-seq dataset from Bioconductor (`airway`)
#   - DESeq2 model using design:
#       ~ cell + dex
#   - Differential expression results generated internally from DESeq2
#   - GO enrichment results:
#       results/tables/go_enrichment.csv
#
# Output:
#   - Sample distance heatmap:
#       results/figures/sample_distance_heatmap.png
#   - Size factor plot:
#       results/figures/size_factors.png
#   - Top up/downregulated genes plot:
#       results/figures/top_up_down_genes.png
#   - Adjusted p-value histogram:
#       results/figures/padj_histogram.png
#   - GO enrichment dot plot with gene overlap:
#       results/figures/go_dotplot_overlap.png
#
# Method:
#   - Reconstruct DESeq2 object from the airway dataset
#   - Run DESeq2 differential expression analysis
#   - Apply variance-stabilizing transformation (VST)
#   - Generate sample-level and gene-level diagnostic plots
#   - Visualize enriched GO terms using adjusted p-values and gene overlap
#
# Notes:
#   - This script is designed to be run independently with Rscript
#   - It does not rely on saved R workspace objects
#   - Figures are saved as PNG files for easy GitHub rendering
#
# Example usage:
#   Rscript scripts/04_interpretation_visualizations.R



# 04_interpretation_visualizations.R

.libPaths(c("/home/mariano/R/library", .libPaths()))

library(DESeq2)
library(airway)
library(ggplot2)
library(pheatmap)

# Load data
data("airway")

# Recreate DESeq2 object
dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

dds$dex <- relevel(dds$dex, ref = "untrt")

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# -----------------------------
# Sample distance heatmap
# -----------------------------

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

png("results/figures/sample_distance_heatmap.png", width = 1000, height = 900)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  main = "Sample-to-sample distance heatmap"
)

pdf("results/figures/heatmap_top50.pdf")

# plot code
dev.off()

# -----------------------------
# Size factors plot
# -----------------------------

sf <- sizeFactors(dds)

sf_df <- data.frame(
  sample = colnames(dds),
  size_factor = sf,
  dex = as.factor(dds$dex)
)

p <- ggplot(sf_df, aes(x = sample, y = size_factor, fill = dex)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "DESeq2 size factors",
    y = "Size factor",
    x = "Sample"
  )

ggsave("results/figures/size_factors.png", p, width = 7, height = 5)

# -----------------------------
# Top up/downregulated genes
# -----------------------------

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]

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

ggsave("results/figures/top_up_down_genes.png", p, width = 8, height = 7)

# -----------------------------
# Adjusted p-value histogram
# -----------------------------

p <- ggplot(res_df, aes(x = padj)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(
    title = "Distribution of adjusted p-values",
    x = "Adjusted p-value",
    y = "Gene count"
  )

ggsave("results/figures/padj_histogram.png", p, width = 7, height = 5)

# -----------------------------
# GO enrichment dot plot
# -----------------------------

go <- read.csv("results/tables/go_enrichment.csv")

go_top <- go[order(go$Adjusted.P.value), ][1:15, ]

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

ggsave("results/figures/go_dotplot_overlap.png", p, width = 9, height = 6)

cat("Interpretation visualizations completed.\n")
