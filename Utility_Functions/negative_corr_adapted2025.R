#' This function has been minimally modified from https://rdrr.io/bioc/anamiR/src/R/negative_cor.R (GPL-2 license) and made to be compatible with miRQuest.
#' Please cite: Wang TT, Lee CY, Lai LC, Tsai MH, Lu TP, Chuang EY. anamiR: integrated analysis of MicroRNA and gene expression profiling. BMC Bioinformatics. 2019 May 14;20(1):239. doi: 10.1186/s12859-019-2870-x. PMID: 31088348; PMCID: PMC6518761.
#' @import stats
#' @export
negative_cor <- function(
  mrna_data,
  mirna_data,
  method = c("pearson",
             "kendall",
             "spearman"),
  cut.off = -0.5
) {

  method <- match.arg(method)

  common_column <- intersect(colnames(mrna_data),
                             colnames(mirna_data))

  mrna_data <- as.data.frame(mrna_data) %>% dplyr::select(all_of(common_column))

  mirna_data <- as.data.frame(mirna_data) %>% dplyr::select(all_of(common_column))

  cal_cor <- function(data_1, data_2, cor_cut) {
    n <- 1
    corr <- list()
    
    withProgress(message = 'Calculating correlations...', value = 0, {
      for (i in seq_len(nrow(data_1))) {
        for (j in seq_len(nrow(data_2))) {
          mrna <- as.numeric(data_1[i, 1:(ncol(mrna_data) - 5)])
          mirna <- as.numeric(data_2[j, 1:(ncol(mirna_data) - 5)])
          tmp <- stats::cor(mrna, mirna, method = method)
          
            
          print(paste0("The row indices are ", i, " for mrna, and ",j, " for mirna"))
          print(tmp)
          if (tmp < cor_cut) {
            corr[[n]] <- row.names(data_2)[j]
            corr[[n]][2] <- row.names(data_1)[i]
            corr[[n]][3] <- tmp
            corr[[n]][4] <- data_2[["log_ratio"]][j]
            corr[[n]][5] <- data_2[["P-adjust"]][j]
            corr[[n]][6] <- data_2[["mean_case"]][j]
            corr[[n]][7] <- data_2[["mean_control"]][j]
            corr[[n]][8] <- data_1[["log_ratio"]][i]
            corr[[n]][9] <- data_1[["P-adjust"]][i]
            corr[[n]][10] <- data_1[["mean_case"]][i]
            corr[[n]][11] <- data_1[["mean_control"]][i]
            n <- n + 1
          }
        
        }
        # Update progress based on the mrna_data row
        incProgress(1/ nrow(data_1), detail = paste("Correlating mRNA", i, "of", nrow(data_1), "with", nrow(data_2), "miRNAs"))
        
      }
    })
    return(corr)
  }

  corr <- cal_cor(mrna_data, mirna_data, cut.off)

  corr <- do.call(rbind, corr)

  if (is.null(corr)) {
    cut.off <- cut.off + 0.2
    corr <- cal_cor(mrna_data, mirna_data, cut.off)
    corr <- do.call(rbind, corr)
  }
  last_column <- paste("Correlation")
  colnames(corr) <- c("miRNA", "Gene", last_column,
                      "logratio_miRNA", "P-adjust(miRNA)",
                      "mean_case(miRNA)", "mean_control(miRNA)",
                      "logratio_gene", "P-adjust(gene)",
                      "mean_case(gene)", "mean_control(gene)")
  return(corr)
}


