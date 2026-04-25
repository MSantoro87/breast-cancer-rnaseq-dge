# 02_deseq2_analysis.R
#
# Purpose:
#   Perform differential gene expression analysis on RNA-seq count data
#   using DESeq2. The script identifies genes that are significantly
#   regulated in response to dexamethasone treatment in airway cells
#   and generates a volcano plot for visualization.
#
# Input:
#   - Airway RNA-seq dataset from Bioconductor (`airway` package)
#   - Experimental design:
#       ~ cell + dex
#     where:
#       cell = donor / biological replicate
#       dex  = treatment condition (trt vs untrt)
#
# Output:
#   - Differential expression results:
#       results/tables/deseq2_results.csv
#   - Volcano plot:
#       results/figures/volcano_plot.pdf
#
# Method:
#   - DESeq2 pipeline:
#       1. Create DESeqDataSet
#       2. Normalize counts
#       3. Estimate dispersion
#       4. Fit negative binomial model
#       5. Perform statistical testing (Wald test)
#
# Notes:
#   - The reference level for treatment is set to "untrt"
#   - Results include log2 fold changes and adjusted p-values (FDR)
#   - Volcano plot highlights significantly regulated genes
#
# Example usage:
#   Rscript scripts/02_deseq2_analysis.R


# 02_deseq2_analysis.R
# Differential expression analysis using DESeq2 on airway dataset

# Set library path
.libPaths(c("/home/mariano/R/library", .libPaths()))

# Load libraries
library(DESeq2)
library(airway)

# Load dataset
data("airway")

# Create DESeq2 dataset
dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

# Set reference level
dds$dex <- relevel(dds$dex, ref = "untrt")

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Save results
write.csv(
  res_df,
  "results/tables/deseq2_results.csv",
  row.names = FALSE
)

# Volcano plot
library(EnhancedVolcano)

pdf("results/figures/volcano_plot.pdf", width = 8, height = 7)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Differential Expression: Dexamethasone Treatment",
  subtitle = "Airway RNA-seq dataset"
)

dev.off()

png("results/figures/volcano_plot.png", width = 1000, height = 800)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Differential Expression: Dexamethasone Treatment",
  subtitle = "Airway RNA-seq dataset"
)

dev.off()



# Heatmap (top 50 genes)
library(pheatmap)

resOrdered <- res[order(res$pvalue), ]
top_genes <- rownames(resOrdered)[1:50]

normalized_counts <- counts(dds, normalized = TRUE)
mat <- normalized_counts[top_genes, ]

pdf("results/figures/heatmap_top50.pdf", width = 8, height = 10)

pheatmap(
  mat,
  scale = "row",
  show_rownames = FALSE,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

cat("DESeq2 analysis completed. Results saved.\n")
