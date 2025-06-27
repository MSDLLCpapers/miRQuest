# Pre-installation script for Posit Connect
# This ensures all Bioconductor packages are properly installed

# Set up BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Configure repositories
options(repos = BiocManager::repositories())

# Get the correct library path for packrat
lib_path <- .libPaths()[1]
message("Installing to library: ", lib_path)

# Force installation to the packrat library only
.libPaths(lib_path)

# Install dependencies directly from Bioconductor URLs
# This bypasses any system library issues

# S4Vectors 0.44.0
tryCatch({
  install.packages("https://bioconductor.org/packages/3.20/bioc/src/contrib/S4Vectors_0.44.0.tar.gz",
                   repos = NULL,
                   type = "source",
                   lib = lib_path)
}, error = function(e) {
  message("S4Vectors installation failed: ", e$message)
})

# IRanges
tryCatch({
  install.packages("https://bioconductor.org/packages/3.20/bioc/src/contrib/IRanges_2.40.1.tar.gz",
                   repos = NULL,
                   type = "source",
                   lib = lib_path)
}, error = function(e) {
  message("IRanges installation failed: ", e$message)
})

# XVector
tryCatch({
  install.packages("https://bioconductor.org/packages/3.20/bioc/src/contrib/XVector_0.46.0.tar.gz",
                   repos = NULL,
                   type = "source",
                   lib = lib_path)
}, error = function(e) {
  message("XVector installation failed: ", e$message)
})

# S4Arrays
tryCatch({
  install.packages("https://bioconductor.org/packages/3.20/bioc/src/contrib/S4Arrays_1.6.0.tar.gz",
                   repos = NULL,
                   type = "source",
                   lib = lib_path)
}, error = function(e) {
  message("S4Arrays installation failed: ", e$message)
})

# Set environment variable to use correct include paths
Sys.setenv(PKG_CPPFLAGS = paste0("-I'", file.path(lib_path, "S4Vectors/include"), "' ",
                                  "-I'", file.path(lib_path, "IRanges/include"), "' ",
                                  "-I'", file.path(lib_path, "XVector/include"), "'"))

# Now install SparseArray
tryCatch({
  install.packages("https://bioconductor.org/packages/3.20/bioc/src/contrib/SparseArray_1.6.2.tar.gz",
                   repos = NULL,
                   type = "source",
                   lib = lib_path)
}, error = function(e) {
  message("SparseArray installation failed: ", e$message)
  # If it fails, try with BiocManager
  BiocManager::install("SparseArray", lib = lib_path, force = TRUE)
})