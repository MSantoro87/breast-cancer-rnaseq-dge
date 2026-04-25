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

### 1. Experimental Design
A DESeq2 model was constructed to account for both:
- donor-specific variability (`cell`)
- treatment effect (`dex`)

```r
design = ~ cell + dex


### 2. Differential Expression Analysis
- Statistical testing performed using DESeq2
- Genes filtered based on adjusted p-value (padj < 0.05)

This defines a dexamethasone-responsive gene set.

### 3. Quality Control & Diagnostics

To validate the robustness of the analysis:

- PCA: global structure of samples
- MA plot: fold-change vs expression relationship
- Dispersion plot: gene-wise variance modeling
- Sample distance heatmap: clustering of biological replicates

### 4. Gene-Level Visualization
- Volcano plot: global differential expression landscape
- Heatmap: top 50 differentially expressed genes
- Counts plot: expression pattern of representative genes

### 5. Functional Enrichment Analysis
Gene Ontology (GO) enrichment
KEGG pathway analysis

These analyses contextualize the transcriptional response in terms of biological processes and pathways.

### Results Transcriptional Response

Approximately ~4000 genes were significantly differentially expressed between treated and untreated samples, indicating a strong transcriptional response to dexamethasone.

This gene set defines a treatment-responsive transcriptional program.

### Sample-Level Structure

PCA analysis shows clear separation between treated and untreated samples, confirming a strong treatment effect. Clustering also reflects donor-specific variability, validating the experimental design.

### Model Diagnostics
-MA plot confirms expected distribution of fold changes
-Dispersion estimates indicate appropriate modeling of variability
-Sample distance heatmap shows clustering by biological condition
-Gene-Level Insights
-Volcano plot highlights both upregulated and downregulated genes
-Heatmap reveals distinct expression patterns driven by treatment
-Individual gene plots confirm strong differential expression signals
-Functional Interpretation

-Enrichment analysis reveals that the dexamethasone-responsive gene program is associated with:

-Immune response pathways
-Inflammatory signaling
-Cellular stress responses

-KEGG pathways further highlight:

regulatory signaling cascades
glucocorticoid-mediated biological processes

🧠 Biological Interpretation

Dexamethasone induces a coordinated transcriptional program consistent with its role as a glucocorticoid regulator of inflammation.

The observed downregulation of immune-related pathways and modulation of signaling processes suggest suppression of inflammatory activity, aligning with known pharmacological effects.

This analysis demonstrates how RNA-seq data can be leveraged to uncover mechanistic insights into treatment response.

- PCA Plot

- Volcano Plot

- Heatmap

- MA Plot

- Dispersion Plot

- GO Enrichment (Dotplot)

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

- Rscript scripts/02_deseq2_analysis.R
- Rscript scripts/03_enrichment_analysis.R
- Rscript scripts/04_interpretation_visualizations.R

🚀 Future Work
- Integrate machine learning (Python) to classify treatment conditions
- Apply workflow to cancer RNA-seq datasets
- Perform gene module and network analysis
- Extend enrichment analysis with additional databases
