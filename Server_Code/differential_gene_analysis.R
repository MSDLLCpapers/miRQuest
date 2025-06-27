## Lines 53-81 have been modified from https://rdrr.io/bioc/anamiR/src/R/differExp_discrete.R (AnamiR, GPL-2 license)
## Therefore, if you use this software, please also cite: Wang TT, Lee CY, Lai LC, Tsai MH, Lu TP, Chuang EY. anamiR: integrated analysis of MicroRNA and gene expression profiling. BMC Bioinformatics. 2019 May 14;20(1):239. doi: 10.1186/s12859-019-2870-x. PMID: 31088348; PMCID: PMC6518761.
# Module 6: Do Differential Gene Expression Analysis with DESEQ2

observeEvent(input$mrna_runDESeq, {
  req(input$designColumn)
  
  notification_id <- showNotification("Computing, please wait...", type = "message", duration = NULL)
  
  # Use pre-loaded data if available, otherwise use uploaded file
  if (!is.null(mrna_wideData())) {
    mRNAcountData <- mrna_wideData()
    mRNAcountData <- mRNAcountData %>% column_to_rownames("X")
  } else if (!is.null(input$mrna_countTable)) {
    mRNAcountData <- read.csv(input$mrna_countTable$datapath, header = TRUE)
    mRNAcountData <- mRNAcountData %>% column_to_rownames("X")
  } else {
    showNotification("Please upload mRNA count data or use pre-loaded demo data", type = "error")
    removeNotification(notification_id)
    return()
  }
  
  metadataData <- metadata()
  
  metadataData <- metadataData %>% filter(SampleID %in% colnames(mRNAcountData))
  mRNAcountData <- mRNAcountData[,metadataData$SampleID]
  metadataData$SampleID <- factor(metadataData$SampleID, levels =colnames(mRNAcountData))
  metadataData <- metadataData[order(metadataData$SampleID), ]
  
  all(metadataData$SampleID==colnames(mRNAcountData))
  
  designFormula <- as.formula(paste("~", input$designColumn))
  dds_mRNA <- DESeqDataSetFromMatrix(countData = mRNAcountData, 
                                     colData = metadataData, 
                                     design = designFormula)
  
  dds_mRNA <- DESeq(dds_mRNA)
  res_mRNA <- results(dds_mRNA)
  
  res_mRNA <- as.data.frame(res_mRNA)
  res_significant_mRNA <- res_mRNA %>% 
    filter(padj < input$padj_filter) %>%
    rownames_to_column("mRNA") %>%
    arrange(padj)
  
  unique_groups <- unique(metadataData[[input$designColumn]])
  
  if (length(unique_groups) < 2) {
    showNotification("Not enough unique groups in the selected design column to compare.", type = "error")
    return()
  }
  
  gp1 <- which(metadataData[[input$designColumn]] == unique_groups[1])
  gp2 <- which(metadataData[[input$designColumn]] == unique_groups[2])
  
  p_value <- res_mRNA[["pvalue"]]
  p_adjust <- res_mRNA[["padj"]]
  FC <- res_mRNA[["log2FoldChange"]]
  
  idx <- which(p_adjust < input$padj_filter)
  DE_data <- mRNAcountData[idx, ]
  
  mean_gp1 <- if (length(gp1) == 1) {
    mRNAcountData[, gp1]
  } else {
    apply(mRNAcountData[, gp1], 1, mean)
  }
  
  mean_gp2 <- if (length(gp2) == 1) {
    mRNAcountData[, gp2]
  } else {
    apply(mRNAcountData[, gp2], 1, mean)
  }
  
  DE_data <- cbind(DE_data, FC[idx], p_value[idx], p_adjust[idx], mean_gp1[idx], mean_gp2[idx])
  
  len_col <- ncol(DE_data)
  colnames(DE_data)[(len_col - 4):len_col] <- c("log_ratio", "P-Value", "P-adjust", "mean_case", "mean_control")
  
  FC_rows <- abs(DE_data[, len_col - 4])
  DE_data <- DE_data[FC_rows > 0.5, ]
  
  gene_names <- row.names(DE_data)
  mrna_data <- as.matrix(sapply(DE_data, as.numeric))
  row.names(mrna_data) <- gene_names
  print(length(res_significant_mRNA$mRNA))
  print(length(intersect(row.names(DE_data), res_significant_mRNA$mRNA)))
  
  mrna_data_reactive(mrna_data)
  
  output$mrna_DE_table <- renderTable({
    req(mrna_data_reactive())
    as.data.frame(mrna_data_reactive())
  })
  
  mrna_res_significant(res_significant_mRNA)
  
  output$mrna_resSignificantTable <- renderTable({
    req(mrna_res_significant())  
    if (!is.null(mrna_res_significant())) {
      filtered_data <- mrna_res_significant() %>%
        filter(padj <= input$padj_filter)  
      
      output$deg_count <- renderText({
        paste("Number of DEGs meeting the criteria:", nrow(filtered_data))
      })
      
      if (nrow(filtered_data) == 0) {
        output$NoDEGmessage <- renderText("No DEGs meet the filtering criteria.")
        return()
      } else {
        return(head(as.data.frame(filtered_data)))
      }
    } else {
      output$NoDEGmessage <- renderText("No DEGs were identified")
      return()
    }
  })
  
  removeNotification(notification_id)
  
  print(dim(mrna_res_significant))
})

# Download handler for DEG table (outside observeEvent)
output$download_DEG_Table <- downloadHandler(
  filename = function() {
    padj_threshold <- input$padj_filter
    paste("Differentially_expressed_genes", 
          "padj_lt",
          padj_threshold, Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(mrna_res_significant())
    data <- mrna_res_significant()
    
    # Handle NULL or empty results
    if (is.null(data) || (is.data.frame(data) && nrow(data) == 0)) {
      # Create empty data frame with appropriate columns
      data <- data.frame(
        mRNA = character(),
        baseMean = numeric(),
        log2FoldChange = numeric(),
        lfcSE = numeric(),
        stat = numeric(),
        pvalue = numeric(),
        padj = numeric()
      )
    }
    
    write.csv(as.data.frame(data), file, row.names = FALSE, quote = FALSE)
  }
)
