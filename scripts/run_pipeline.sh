#!/usr/bin/env bash
set -euo pipefail

echo "Starting RNA-seq DGE pipeline..."

Rscript scripts/00_check_environment.R
Rscript scripts/01_download_dataset.R
Rscript scripts/02_deseq2_analysis.R
Rscript scripts/03_enrichment_analysis.R
Rscript scripts/04_interpretation_visualizations.R
Rscript scripts/05_up_down_enrichment.R

echo "Pipeline completed successfully!"
