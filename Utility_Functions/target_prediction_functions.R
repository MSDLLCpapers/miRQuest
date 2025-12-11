### Helper functions for target prediction ---

# helper to run an expression and capture warnings without printing them
captureWarnings <- function(expr) {
  warnings <- character(0)
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings)
}


# Helpers
is_empty_df <- function(x) is.null(x) || (is.data.frame(x) && nrow(x) == 0)

normalize_columns <- function(df, org_col) {
  desired_cols <- c(
    "entrezgene_id",
    "rank_final",
    "ensembl_gene_id",
    "description",
    org_col,
    "miRNA",
    "miRTarBase.ID",
    "Experiments",
    "Support.Type",
    "References..PMID."
  )
  
  if (nrow(df) == 0) {
    # Build an empty schema safely
    return(make_empty_schema(org_col))
  }
  
  # For non-empty data frames, add missing columns with correct types
  for (nm in desired_cols) {
    if (!(nm %in% names(df))) {
      if (nm == "rank_final") {
        df[[nm]] <- as.numeric(NA)
      } else {
        df[[nm]] <- as.character(NA)
      }
    }
  }
  df <- df[, desired_cols, drop = FALSE]
  df
}

make_empty_schema <- function(org_col) {
  cols_chr <- c(
    "entrezgene_id",
    "ensembl_gene_id",
    "description",
    org_col,
    "miRNA",
    "miRTarBase.ID",
    "Experiments",
    "Support.Type",
    "References..PMID."
  )
  cols_num <- c("rank_final")
  
  # Create zero-row data frames with correct types
  df_chr <- as.data.frame(
    setNames(
      replicate(length(cols_chr), character(0), simplify = FALSE),
      cols_chr
    ),
    stringsAsFactors = FALSE
  )
  df_num <- as.data.frame(
    setNames(
      replicate(length(cols_num), numeric(0), simplify = FALSE),
      cols_num
    )
  )
  # Bind and order
  empty_df <- cbind(df_chr, df_num)
  empty_df <- empty_df[, c(
    "entrezgene_id",
    "rank_final",
    "ensembl_gene_id",
    "description",
    org_col,
    "miRNA",
    "miRTarBase.ID",
    "Experiments",
    "Support.Type",
    "References..PMID."
  ), drop = FALSE]
  empty_df
}