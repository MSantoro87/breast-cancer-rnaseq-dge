# 00_check_environment.R
#
# Purpose:
#   Verify that the required R packages for the RNA-seq analysis pipeline
#   are installed and can be loaded successfully. The script also records
#   session information to support reproducibility.
#
# Input:
#   - User-specific R library path:
#       /home/mariano/R/library
#   - Installed R packages:
#       DESeq2, ggplot2, pheatmap, EnhancedVolcano, enrichR
#
# Output:
#   - Printed R session information via sessionInfo()
#   - Intended output file when run from Bash:
#       environment/sessionInfo.txt
#
# Example usage:
#   Rscript scripts/00_check_environment.R > environment/sessionInfo.txt

.libPaths(c("/home/mariano/R/library", .libPaths()))

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(enrichR)

sessionInfo()
