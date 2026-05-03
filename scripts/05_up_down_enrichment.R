#!/usr/bin/env Rscript

# 05_up_down_enrichment.R
# Purpose: Run enrichment separately for upregulated and downregulated genes.

suppressPackageStartupMessages({
  library(enrichR)
  library(DESeq2)
  library(airway)
  library(yaml)
})

# Load config
config <- yaml::read_yaml("config/config.yaml")
enrichment_databases <- config$enrichment_databases

# Output directory
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

message("Loading DESeq2 results...")

up_genes <- read.csv("results/tables/deseq2_results_upregulated.csv")
down_genes <- read.csv("results/tables/deseq2_results_downregulated.csv")

if (nrow(up_genes) == 0 || nrow(down_genes) == 0) {
  stop("Upregulated or downregulated gene list is empty.")
}

message("Loading gene annotation...")

data("airway")

dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

gene_info <- as.data.frame(rowData(dds))[, c("gene_id", "symbol")]

# Merge gene symbols
up_genes <- merge(up_genes, gene_info, by = "gene_id")
down_genes <- merge(down_genes, gene_info, by = "gene_id")

# Clean gene symbols
clean_symbols <- function(x) {
  x <- unique(x)
  x <- x[!is.na(x)]
  x <- x[x != ""]
  return(x)
}

up_symbols <- clean_symbols(up_genes$symbol)
down_symbols <- clean_symbols(down_genes$symbol)

message("Running enrichment for UPREGULATED genes...")
up_results <- enrichr(up_symbols, enrichment_databases)

message("Running enrichment for DOWNREGULATED genes...")
down_results <- enrichr(down_symbols, enrichment_databases)

# Extract results
go_up <- up_results[["GO_Biological_Process_2021"]]
kegg_up <- up_results[["KEGG_2021_Human"]]

go_down <- down_results[["GO_Biological_Process_2021"]]
kegg_down <- down_results[["KEGG_2021_Human"]]

# Save results
write.csv(go_up, "results/tables/go_enrichment_upregulated.csv", row.names = FALSE)
write.csv(kegg_up, "results/tables/kegg_enrichment_upregulated.csv", row.names = FALSE)

write.csv(go_down, "results/tables/go_enrichment_downregulated.csv", row.names = FALSE)
write.csv(kegg_down, "results/tables/kegg_enrichment_downregulated.csv", row.names = FALSE)

message("Up/down enrichment analysis completed.")
message("Upregulated genes used: ", length(up_symbols))
message("Downregulated genes used: ", length(down_symbols))
