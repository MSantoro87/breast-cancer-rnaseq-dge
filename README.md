# RNA-seq Differential Expression (Airway)

## Date
2026-04-24

## Objective
Analyze gene expression changes in airway cells after dexamethasone treatment.

## Dataset
- Bioconductor: `airway`
- 8 samples (trt vs untrt)

## Environment
- Ubuntu 20.04
- R 4.5.2
- Key packages: DESeq2, ggplot2, pheatmap, EnhancedVolcano, enrichR

## Steps Completed
1. Loaded dataset and inspected metadata
2. Built DESeq2 model: ~ cell + dex
3. Ran differential expression (DESeq2)
4. Generated plots:
   - Volcano plot → `results/figures/volcano_plot.pdf`
   - Heatmap (top 50 genes) → `results/figures/heatmap_top50.pdf`
5. Extracted significant genes (padj < 0.05)
6. Ran enrichment (GO + KEGG)
7. Created GO barplot → `results/figures/go_barplot.pdf`

## Partial results (draft)
- ~4000 significant genes (padj < 0.05)
- Clear transcriptional response to dexamethasone

## Results

### DE Analysis

Differential gene expression analysis using DESeq2 identified approximately 4000 genes with significant expression changes (adjusted p-value < 0.05) between dexamethasone-treated and untreated airway smooth muscle cells. These results indicate a broad transcriptional response to glucocorticoid treatment.

### Visualization

A volcano plot highlighted the distribution of significantly upregulated and downregulated genes, showing a balanced but substantial shift in gene expression upon treatment. A heatmap of the top 50 differentially expressed genes revealed clear clustering of samples according to treatment condition, confirming that dexamethasone induces a distinct transcriptional signature.

### Functional Enrichment Analysis

Gene Ontology (GO) enrichment analysis revealed that the differentially expressed genes are significantly associated with biological processes related to immune response, inflammatory signaling, and cellular stress pathways. These findings are consistent with the known anti-inflammatory and immunomodulatory effects of dexamethasone.

KEGG pathway analysis further supported these results, identifying pathways involved in signaling and regulatory mechanisms that are typically modulated by glucocorticoids.

### Interpretation

Overall, the analysis demonstrates that dexamethasone treatment leads to coordinated changes in gene expression affecting immune and inflammatory processes. This is biologically plausible given the role of glucocorticoids in suppressing inflammation and regulating immune responses.

These results validate the computational workflow and illustrate how RNA-seq analysis can be used to uncover biologically meaningful patterns in gene expression data.

## Notes / Issues
- Fixed PDF device issue (`dev.off()`)
- Switched to minimal package set (no tidyverse/clusterProfiler)

## Next Steps
- Interpret GO/KEGG results
- Write biological interpretation
- Optional: add Python ML step
