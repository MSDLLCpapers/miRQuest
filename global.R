# Set repositories to include both CRAN and Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  install.packages("BiocManager")
}

options(repos = BiocManager::repositories())

# renv::restore()

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

if (!requireNamespace("shinyjs", quietly = TRUE)) {
    install.packages("shinyjs")
}

library(shinyjs)

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