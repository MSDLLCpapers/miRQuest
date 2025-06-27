# Custom installation script for RStudio Connect

# Remove any system library paths to avoid conflicts
.libPaths(.libPaths()[1])

# Set environment to avoid using system libraries during compilation
Sys.setenv(R_LIBS_USER = .libPaths()[1])
Sys.setenv(R_LIBS = .libPaths()[1])

# Configure options
renv::settings$install.packages.compile.from.source("never")
options(install.packages.compile.from.source = "never")
options(pkgType = "binary")

# Restore packages
renv::restore()

# If SparseArray fails, try manual installation
if (!requireNamespace("SparseArray", quietly = TRUE)) {
  message("SparseArray not found after restore, attempting manual installation...")
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Force removal of any cached/partial installations
  try(remove.packages("SparseArray"), silent = TRUE)
  
  # Install with binary preference
  BiocManager::install("SparseArray", 
                       version = "3.20",
                       ask = FALSE,
                       force = TRUE,
                       type = "binary")
}

# Update snapshot
renv::snapshot(type = "all")