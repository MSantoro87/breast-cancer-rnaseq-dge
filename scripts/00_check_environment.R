#!/usr/bin/env Rscript

# 00_check_environment.R
# Purpose: Validate the R environment, required packages, configuration file,
# project directories, and optional Quarto installation for the RNA-seq pipeline.

required_packages <- c(
  "airway",
  "SummarizedExperiment",
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "enrichR",
  "yaml",
  "readr",
  "dplyr",
  "knitr",
  "rmarkdown"
)

required_dirs <- c(
  "config",
  "scripts",
  "data/processed",
  "results/tables",
  "results/figures/png",
  "results/figures/pdf",
  "reports",
  "environment",
  "logs"
)

message("Checking project directories...")

for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    message("Created missing directory: ", dir)
  }
}

message("Checking configuration file...")

if (!file.exists("config/config.yaml")) {
  stop("Missing config/config.yaml")
}

message("Checking required R packages...")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall missing CRAN packages with install.packages() and Bioconductor packages with BiocManager::install()."
  )
}

message("All required R packages are installed.")

message("Checking Quarto installation...")

quarto_path <- Sys.which("quarto")

if (quarto_path == "") {
  warning("Quarto was not found. Report rendering will not work until Quarto is installed.")
} else {
  message("Quarto found: ", quarto_path)
  system("quarto --version")
}

message("Writing session information...")

sink("environment/R_environment.txt")
cat("RNA-seq DGE pipeline environment report\n")
cat("Generated on: ", as.character(Sys.time()), "\n\n")
cat("R version:\n")
print(R.version.string)
cat("\nLibrary paths:\n")
print(.libPaths())
cat("\nSession info:\n")
print(sessionInfo())
sink()

message("Environment check completed successfully.")
message("Session information saved to environment/R_environment.txt")
