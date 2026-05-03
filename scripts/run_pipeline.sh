#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/processed logs results/tables results/figures

echo "Running RNA-seq DGE pipeline..."

Rscript scripts/00_check_environment.R 2>&1 | tee logs/00_check_environment.log
Rscript scripts/01_download_dataset.R 2>&1 | tee logs/01_download_dataset.log
Rscript scripts/02_deseq2_analysis.R 2>&1 | tee logs/02_deseq2_analysis.log
Rscript scripts/03_enrichment_analysis.R 2>&1 | tee logs/03_enrichment_analysis.log
Rscript scripts/04_interpretation_visualizations.R 2>&1 | tee logs/04_interpretation_visualizations.log
Rscript scripts/05_up_down_enrichment.R 2>&1 | tee logs/05_up_down_enrichment.log

echo "Pipeline completed successfully."
