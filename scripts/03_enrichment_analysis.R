#!/usr/bin/env Rscript

# 03_enrichment_analysis.R
# Purpose: Run GO and KEGG enrichment analysis on significant DESeq2 genes.

suppressPackageStartupMessages({
  library(enrichR)
  library(ggplot2)
  library(DESeq2)
  library(airway)
})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

message("Loading significant DESeq2 results...")

sig_genes <- read.csv("results/tables/deseq2_results_significant.csv")

if (nrow(sig_genes) == 0) {
  stop("No significant genes found in deseq2_results_significant.csv")
}

message("Loading airway annotation...")

data("airway")

dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

gene_info <- as.data.frame(rowData(dds))[, c("gene_id", "symbol")]

sig_genes <- merge(
  sig_genes,
  gene_info,
  by = "gene_id"
)

gene_symbols <- unique(sig_genes$symbol)
gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_symbols <- gene_symbols[gene_symbols != ""]

if (length(gene_symbols) == 0) {
  stop("No valid gene symbols found for enrichment analysis.")
}

message("Running enrichR analysis...")

selected_dbs <- c(
  "GO_Biological_Process_2021",
  "KEGG_2021_Human"
)

enrich_results <- enrichr(gene_symbols, selected_dbs)

go_results <- enrich_results[["GO_Biological_Process_2021"]]
kegg_results <- enrich_results[["KEGG_2021_Human"]]

write.csv(
  go_results,
  "results/tables/go_enrichment.csv",
  row.names = FALSE
)

write.csv(
  kegg_results,
  "results/tables/kegg_enrichment.csv",
  row.names = FALSE
)

if (!is.null(go_results) && nrow(go_results) > 0) {
  go_top <- go_results[order(go_results$Adjusted.P.value), ]
  go_top <- head(go_top, 10)

  p <- ggplot(
    go_top,
    aes(
      x = reorder(Term, Adjusted.P.value),
      y = -log10(Adjusted.P.value)
    )
  ) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Top Enriched GO Biological Processes",
      x = "GO term",
      y = "-log10 adjusted p-value"
    )

  ggsave(
    "results/figures/go_barplot.pdf",
    p,
    width = 9,
    height = 6
  )
}

message("Enrichment analysis completed.")
message("Significant genes used: ", nrow(sig_genes))
message("Unique gene symbols used: ", length(gene_symbols))
message("GO terms returned: ", nrow(go_results))
message("KEGG terms returned: ", nrow(kegg_results))
