# RNA-seq Differential Expression Analysis: Dexamethasone Response in Airway Cells

## 📌 Overview
This project implements a **reproducible RNA-seq analysis workflow** to characterize the transcriptional response of human airway smooth muscle cells to **dexamethasone treatment**.

Beyond standard differential expression analysis, the study aims to define a **treatment-responsive gene program** and interpret its biological significance using functional enrichment and visualization.

---

## 🎯 Biological Question
How does dexamethasone treatment reshape the transcriptional landscape of airway smooth muscle cells?

Specifically:
- Which genes are significantly regulated?
- What biological processes are affected?
- Can a coherent **transcriptional program** be identified?

---

## 🧬 Dataset
- Source: Bioconductor `airway` dataset  
- Samples: 8 (treated vs untreated)  
- Data type: RNA-seq count data  
- Design: paired biological samples across different donors  

---

## ⚙️ Computational Environment
- OS: Ubuntu 20.04  
- R version: 4.5.2  

### Key packages
- DESeq2  
- ggplot2  
- pheatmap  
- EnhancedVolcano  
- enrichR  

---

## 🧪 Analysis Workflow

## 1. Experimental Design
A DESeq2 model was constructed to account for both:
- donor-specific variability (`cell`)
- treatment effect (`dex`)

```r
design = ~ cell + dex
```

## 2. Differential Expression Analysis
- Statistical testing performed using DESeq2
- Genes filtered based on adjusted p-value (padj < 0.05)

This defines a dexamethasone-responsive gene set.

## 3. Quality Control & Diagnostics
- PCA analysis of samples
- MA plot for expression changes
- Dispersion plot for model quality
- Sample-to-sample distance heatmap
- Size factor normalization assessment

## 4. Gene-Level & Program Analysis
- Volcano plot of differential expression
- Heatmap of top differentially expressed genes
- Counts plot for representative genes
- Identification of upregulated and downregulated gene sets

## 5. Functional Enrichment 
- GO and KEGG enrichment analysis
- Separate enrichment for upregulated and downregulated genes
- GO barplot and dotplot visualizations

## 6. Results Transcriptional Response

Approximately ~4000 genes were significantly differentially expressed between treated and untreated samples, indicating a strong transcriptional response to dexamethasone. 
This gene set defines a treatment-responsive transcriptional program.

## 7. Sample-Level Structure

PCA analysis shows clear separation between treated and untreated samples, confirming a strong treatment effect. Clustering also reflects donor-specific variability, validating the experimental design.

## 8. Model Diagnostics
-MA plot confirms expected distribution of fold changes
-Dispersion estimates indicate appropriate modeling of variability
-Sample distance heatmap shows clustering by biological condition
-Gene-Level Insights
-Volcano plot highlights both upregulated and downregulated genes
-Heatmap reveals distinct expression patterns driven by treatment
-Individual gene plots confirm strong differential expression signals
-Functional Interpretation

## 9. Partials (in brief)

Enrichment analysis reveals that the dexamethasone-responsive gene program is associated with:
-Immune response pathways
-Inflammatory signaling
-Cellular stress responses

-KEGG pathways further highlight:

regulatory signaling cascades
glucocorticoid-mediated biological processes

## 10. Biological Interpretation (staged partials)

- Dexamethasone induces a coordinated transcriptional program consistent with its role as a glucocorticoid regulator of inflammation.

- The observed downregulation of immune-related pathways and modulation of signaling processes suggest suppression of inflammatory activity, aligning with known pharmacological effects.

- This analysis demonstrates how RNA-seq data can be leveraged to uncover mechanistic insights into treatment response.

- PCA Plot

- Volcano Plot

- Heatmap

- MA Plot

- Dispersion Plot

- GO Enrichment (Dotplot)

### 11. Functional Interpretation

The dexamethasone-responsive transcriptional program can be decomposed into:

- **Upregulated genes** associated with regulatory and signaling processes
- **Downregulated genes** enriched in immune and inflammatory pathways
- **This separation** highlights distinct biological programs that are activated or suppressed in response to treatment.

## 12. Summary Statistics

- Total genes tested: ~33,000  
- Significant genes (padj < 0.05): ~4,000  
- Upregulated genes: X  
- Downregulated genes: Y  

📁 Project Structure
.
├── data/
├── environment/
├── metadata/
├── results/
│   ├── figures/
│   └── tables/
├── scripts/
└── README.md
🔁 Reproducibility

Run the full pipeline:

- ``Rscript scripts/02_deseq2_analysis.R``
- ``Rscript scripts/03_enrichment_analysis.R``
- ``Rscript scripts/04_interpretation_visualizations.R``

## 13. Future Work
- Integrate machine learning (Python) to classify treatment conditions
- Apply workflow to cancer RNA-seq datasets
- Perform gene module and network analysis
- Extend enrichment analysis with additional databases
