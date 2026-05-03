# 05_up_down_enrichment.R
# Purpose: Separate enrichment analysis for upregulated and downregulated genes
# Input: DESeq2 results from airway RNA-seq analysis
# Output: up/down gene tables, enrichment tables, summary statistics

.libPaths(c("/home/mariano/R/library", .libPaths()))

library(DESeq2)
library(airway)
library(enrichR)

data("airway")

dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

dds$dex <- relevel(dds$dex, ref = "untrt")
dds <- DESeq(dds)

res <- results(dds, contrast = c("dex", "trt", "untrt"))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]

gene_info <- as.data.frame(rowData(dds))[, c("gene_id", "symbol")]

res_annotated <- merge(
  res_df,
  gene_info,
  by = "gene_id"
)

res_annotated <- res_annotated[
  !is.na(res_annotated$symbol) & res_annotated$symbol != "",
]

sig_genes <- res_annotated[res_annotated$padj < 0.05, ]

up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]

write.csv(up_genes, "results/tables/upregulated_genes.csv", row.names = FALSE)
write.csv(down_genes, "results/tables/downregulated_genes.csv", row.names = FALSE)

up_symbols <- unique(up_genes$symbol)
down_symbols <- unique(down_genes$symbol)

selected_dbs <- c(
  "GO_Biological_Process_2021",
  "KEGG_2021_Human"
)

up_enrich <- enrichr(up_symbols, selected_dbs)
down_enrich <- enrichr(down_symbols, selected_dbs)

write.csv(
  up_enrich[["GO_Biological_Process_2021"]],
  "results/tables/go_enrichment_upregulated.csv",
  row.names = FALSE
)

write.csv(
  down_enrich[["GO_Biological_Process_2021"]],
  "results/tables/go_enrichment_downregulated.csv",
  row.names = FALSE
)

write.csv(
  up_enrich[["KEGG_2021_Human"]],
  "results/tables/kegg_enrichment_upregulated.csv",
  row.names = FALSE
)

write.csv(
  down_enrich[["KEGG_2021_Human"]],
  "results/tables/kegg_enrichment_downregulated.csv",
  row.names = FALSE
)

summary_stats <- data.frame(
  metric = c(
    "genes_tested",
    "significant_genes_padj_0.05",
    "upregulated_genes",
    "downregulated_genes"
  ),
  value = c(
    nrow(res_df),
    nrow(sig_genes),
    nrow(up_genes),
    nrow(down_genes)
  )
)

write.csv(
  summary_stats,
  "results/tables/summary_statistics.csv",
  row.names = FALSE
)

cat("Up/down enrichment completed.\n")
cat("Upregulated genes:", nrow(up_genes), "\n")
cat("Downregulated genes:", nrow(down_genes), "\n")

write.csv(up_enrichment, "results/tables/upregulated_enrichment.csv")
write.csv(down_enrichment, "results/tables/downregulated_enrichment.csv")
