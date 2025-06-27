# This file ensures proper repository configuration for Bioconductor packages
# It should be sourced at the beginning of your app.R or server.R file

# Check if BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Set repositories to include Bioconductor
options(repos = BiocManager::repositories())

# Ensure S4Arrays can be found
if (!requireNamespace("S4Arrays", quietly = TRUE)) {
  BiocManager::install("S4Arrays", version = "3.20")
}