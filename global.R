# Set repositories to include both CRAN and Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  install.packages("BiocManager")
}

options(repos = BiocManager::repositories())

# IMPORTANT: renv::restore() disabled for Docker deployments
# Packages are installed during Docker image build, not at runtime
# If running locally (non-Docker), uncomment the line below:
# renv::restore()

# Docker runtime diagnostics (only runs in containerized environment)
if (file.exists("/.dockerenv")) {
  cat("\n=== miRQuest Docker Runtime ===\n")
  cat("Working directory:", getwd(), "\n")
  cat("Library paths:", paste(.libPaths(), collapse = "\n  "), "\n")
  cat("R version:", R.version.string, "\n")
  cat("Container startup time:", Sys.time(), "\n")
  cat("================================\n\n")
}

library(renv)
library(shiny)
library(vroom)
library(dplyr)
library(tidyverse)
library(miRNAtap) 
library(miRNAtap.db)
library(topGO) 
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(cowplot)
library(circlize)
library(biomaRt)
library(paletteer)
library(cowplot)
library(edgeR)
library(pheatmap)
library(limma)
library(clusterProfiler)
library(ggraph)
library(igraph)
library(reactome.db)
library(ReactomePA)
library(bslib)
library(here)
library(DT)
library(plotly)
library(visNetwork)

if (!requireNamespace("shinyjs", quietly = TRUE)) {
    install.packages("shinyjs")
}

library(shinyjs)

# Increase max upload size
options(shiny.maxRequestSize = 500 * 1024^2)

# Source all R scripts in the scripts directory
scripts_to_source <- 
    list.files(
        path = "Utility_Functions",
        pattern = "\\.R$",
        full.names = TRUE
    )

for (script in scripts_to_source) {
    cat("Sourcing script:", script, "\n")
    source(script)
}
