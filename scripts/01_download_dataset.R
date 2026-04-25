# 01_download_dataset.R
#
# Purpose:
#   Download and store metadata for a publicly available RNA-seq dataset
#   from the Gene Expression Omnibus (GEO). This script retrieves the
#   dataset GSE183947 and saves it locally for downstream analysis.
#
# Input:
#   - GEO dataset identifier:
#       GSE183947
#   - Internet connection to access GEO repository
#
# Output:
#   - GEO metadata object saved as:
#       data/raw/GSE183947_geo_metadata.rds
#   - Printed summary of the dataset structure (ExpressionSet)
#
# Dependencies:
#   - GEOquery (Bioconductor)
#   - BiocManager (for installation if missing)
#
# Notes:
#   - This script retrieves metadata only, not processed count matrices.
#   - The saved RDS object can be reused without re-downloading the dataset.
#
# Example usage:
#   Rscript scripts/01_download_dataset.R


# 01_download_dataset.R
# Goal: Download and inspect GEO dataset GSE183947

.libPaths(c("/home/mariano/R/library", .libPaths()))

# Install GEOquery if missing
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", lib = "/home/mariano/R/library")
}

library(GEOquery)

# Download GEO metadata
gse <- getGEO("GSE183947", GSEMatrix = TRUE)

# Inspect object
print(gse)

# Save metadata object
saveRDS(gse, file = "data/raw/GSE183947_geo_metadata.rds")
