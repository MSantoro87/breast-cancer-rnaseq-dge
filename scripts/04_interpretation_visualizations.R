.libPaths(c("/home/mariano/R/library", .libPaths()))

## Sample distance heatmap

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

png("results/figures/sample_distance_heatmap.png", width = 1000, height = 900)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  main = "Sample-to-sample distance heatmap"
)
dev.off()

## Size factors plot

sf <- sizeFactors(dds)
sf_df <- data.frame(sample = names(sf), size_factor = sf, dex = dds$dex)

p <- ggplot(sf_df, aes(x = sample, y = size_factor, fill = dex)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "DESeq2 size factors", y = "Size factor", x = "Sample")

ggsave("results/figures/size_factors.png", p, width = 7, height = 5)

## Top up/downregulated genes

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- na.omit(res_df)

top_up <- res_df[order(-res_df$log2FoldChange), ][1:15, ]
top_down <- res_df[order(res_df$log2FoldChange), ][1:15, ]

top_genes <- rbind(
  data.frame(top_up, direction = "Upregulated"),
  data.frame(top_down, direction = "Downregulated")
)

p <- ggplot(top_genes, aes(x = reorder(gene_id, log2FoldChange), y = log2FoldChange, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top up- and downregulated genes", x = "Gene", y = "log2 fold change")

ggsave("results/figures/top_up_down_genes.png", p, width = 8, height = 7)

## Adjusted p-value histogram

p <- ggplot(res_df, aes(x = padj)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(title = "Distribution of adjusted p-values", x = "Adjusted p-value", y = "Gene count")

ggsave("results/figures/padj_histogram.png", p, width = 7, height = 5)

## Enrichment dot plot with gene overlap

go <- read.csv("results/tables/go_enrichment.csv")
go_top <- go[order(go$Adjusted.P.value), ][1:15, ]

go_top$minus_log10_padj <- -log10(go_top$Adjusted.P.value)
go_top$overlap_count <- as.numeric(sub("/.*", "", go_top$Overlap))

p <- ggplot(go_top, aes(
  x = minus_log10_padj,
  y = reorder(Term, minus_log10_padj),
  size = overlap_count
)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Top enriched GO Biological Processes",
    x = "-log10 adjusted p-value",
    y = "GO term",
    size = "Gene overlap"
  )

ggsave("results/figures/go_dotplot_overlap.png", p, width = 9, height = 6)
