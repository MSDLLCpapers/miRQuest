## Lines 102-131 have been adapted from https://rdrr.io/bioc/anamiR/src/R/differExp_discrete.R (AnamiR, GPL-2 license)
## Therefore, if you use this software, please also cite: Wang TT, Lee CY, Lai LC, Tsai MH, Lu TP, Chuang EY. anamiR: integrated analysis of MicroRNA and gene expression profiling. BMC Bioinformatics. 2019 May 14;20(1):239. doi: 10.1186/s12859-019-2870-x. PMID: 31088348; PMCID: PMC6518761.
# Module 3: Do Differential miRNA Expression Analysis with DESEQ2

# Reactive values for storing heatmap and volcano plot functions
significant_heatmap_reactive <- reactiveVal()
volcano_plot_reactive <- reactiveVal()

observeEvent(input$runDESeq, {
  req(wideData(), 
      metadata(), 
      input$designColumn)
  
  notification_id <- showNotification("Computing, please wait...", type = "message", duration = NULL)
  
  countData <- wideData()
  metadataData <- metadata()
  
  metadataData <- metadataData %>% filter(SampleID %in% colnames(countData))
  
  countData <- countData[,metadataData$SampleID]
  metadataData$SampleID <- factor(metadataData$SampleID, levels =colnames(countData))
  metadataData <- metadataData[order(metadataData$SampleID), ]
  all(metadataData$SampleID==names(countData))
  
  if (!all(metadataData$SampleID == names(countData))) {
    showNotification("Sample IDs do not match between metadata and count data", type = "error")
    return()
  }
  
  designFormula <- as.formula(paste("~", input$designColumn))
  
  dds <- DESeqDataSetFromMatrix(countData = countData, 
                                colData = metadataData, 
                                design = designFormula)
  
  dds <- DESeq(dds)
  DESEQ_obj(dds)
  
  norm_cts <- vst(dds, blind = TRUE, sum( rowMeans( counts(dds, normalized=TRUE)) > 5 )) %>% assay()
  norm_cts_reactive(norm_cts)
  heatmap_annotation <- as.data.frame(colData(dds)) %>%
    dplyr::select(c(input$designColumn))
  heatmap_annotation_reactive(heatmap_annotation)
  print(rownames(norm_cts))
  
  res <- results(dds)
  
  res <- as.data.frame(res)
  DE_miRNA_results_reactive(res)
  
  res_significant_data <- res %>% 
    filter(padj< input$DEM_padj_filter) %>%
    rownames_to_column("miRNA") %>% 
    arrange(padj)
  
  res_significant(res_significant_data)
  
  output$resSignificantTable <- renderTable({
    req(res_significant())  
    if (!is.null(res_significant())) {
      filtered_data <- res_significant() %>%
        filter(padj < input$DEM_padj_filter)  
      
      output$dem_count <- renderText({
        paste("Number of DE MiRNAs meeting the threshold criteria:", nrow(filtered_data))
      })
      
      if (nrow(filtered_data) == 0) {
        output$NoDEGmessage <- renderText("No DE MiRNAs meet the filtering criteria.")
        return()
      } else {
        return(head(as.data.frame(filtered_data)))
      }
    } else {
      output$NoDEGmessage <- renderText("No DE MiRNAs were identified")
      return()
    }
  })
  
  removeNotification(notification_id)
  
  # Download handler for differential miRNA
  output$download_significant_miRNA <- downloadHandler(
    filename = function() {
      paste("Differentially_Expressed_miRNA", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(as.data.frame(res_significant()), file, row.names = FALSE, quote = FALSE)
    }
  )
  
  ### Prepare the counts for negative correlation module ---
  
  unique_groups <- unique(metadataData[[input$designColumn]])
  
  if (length(unique_groups) < 2) {
    showNotification("Not enough unique groups in the selected design column to compare.", type = "error")
    return()
  }
  
  gp1 <- which(metadataData[[input$designColumn]] == unique_groups[1])
  gp2 <- which(metadataData[[input$designColumn]] == unique_groups[2])
  
  p_value <- res[["pvalue"]]
  p_adjust <- res[["padj"]]
  FC <- res[["log2FoldChange"]]
  
  idx <- which(p_adjust < input$DEM_padj_filter)
  
  DE_data <- countData[idx, ]
  
  mean_gp1 <- if (length(gp1) == 1) {
    countData[, gp1]
  } else {
    apply(countData[, gp1], 1, mean)
  }
  
  mean_gp2 <- if (length(gp2) == 1) {
    countData[, gp2]
  } else {
    apply(countData[, gp2], 1, mean)
  }
  
  DE_data <- cbind(DE_data, FC[idx], p_value[idx], p_adjust[idx], mean_gp1[idx], mean_gp2[idx])
  
  len_col <- ncol(DE_data)
  colnames(DE_data)[(len_col - 4):len_col] <- c("log_ratio", "P-Value", "P-adjust", "mean_case", "mean_control")
  
  FC_rows <- abs(DE_data[, len_col - 4])
  DE_data <- DE_data[FC_rows > 0.5, ]
  
  gene_names <- row.names(DE_data)
  mirna_data <- as.data.frame(DE_data)
  row.names(mirna_data) <- gene_names
  
  mirna_data_reactive(mirna_data)
  
  # future_expansion_code: DE miRNA data table for debugging/correlation analysis
  # output$mirna_DE_table <- renderTable({
  #   req(mirna_data_reactive())
  #   as.data.frame(mirna_data_reactive())
  # })
  # 
  # # Download handler for miRNA DE table (for correlation analysis)
  # output$download_mirna_DE_table <- downloadHandler(
  #   filename = function() {
  #     paste("miRNA_DE_data_for_correlation_", Sys.Date(), ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     write.csv(as.data.frame(mirna_data_reactive()), file, row.names = TRUE)
  #   }
  # )
  
  ### Show a Select DE miRNA as a Boxplot ---
  
  updateSelectInput(session, "selectedMiRNA", choices = res_significant()$miRNA)
  updateSelectInput(session, "DEMiRNA", choices = res_significant()$miRNA)
  
  grouping_variable <- input$designColumn
  levels_group <- levels(as.factor(metadataData[[grouping_variable]]))
  
  if (length(levels_group) < 2) {
    showNotification("Not enough levels in the selected grouping variable for comparison.", type = "error")
    return()
  }
  
  output$boxplot <- renderPlot({
    req(input$groupBy_Boxplot, 
        input$selectedMiRNA, 
        DESEQ_obj())
    grouping_variable <- input$groupBy_Boxplot
    
    dds <- DESEQ_obj()
    selected_miRNA <- input$selectedMiRNA
    
    create_boxplot(dds, selected_miRNA, grouping_variable)
  })
  
  output$downloadBoxplot <- downloadHandler(
    filename = function() {
      paste("boxplot_", input$selectedMiRNA, ".png", sep = "")
    },
    content = function(file) {
      png(file)
      dds <- DESEQ_obj()
      selected_miRNA <- input$selectedMiRNA
      grouping_variable <- input$groupBy_Boxplot
      
      print(create_boxplot(dds, selected_miRNA, grouping_variable))
      dev.off()
    }
  )
})

observeEvent(input$significant_heatmap_button, {
  req(res_significant(), norm_cts_reactive())
  
  res_significant_data <- res_significant()
  norm_cts <- norm_cts_reactive()
  heatmap_annotation <- heatmap_annotation_reactive()
  
  # Handle both DESeq2 (log2FoldChange) and limma (logFC) column names
  if ("log2FoldChange" %in% colnames(res_significant_data)) {
    select_significant <- res_significant_data %>% 
      mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
      filter(abs_log2FoldChange >= input$DEM_log2FoldChange_filter) %>% 
      pull(miRNA)
  } else if ("logFC" %in% colnames(res_significant_data)) {
    select_significant <- res_significant_data %>% 
      mutate(abs_log2FoldChange = abs(logFC)) %>%
      filter(abs_log2FoldChange >= input$DEM_log2FoldChange_filter) %>% 
      pull(miRNA)
  } else {
    # If neither column exists, show all significant miRNAs
    select_significant <- res_significant_data$miRNA
  }
  
  print(intersect(rownames(norm_cts),select_significant))
  
  subset_norm_cts <- norm_cts[rownames(norm_cts) %in% select_significant, ]
  
  # Store heatmap function for download
  heatmap_significant <- pheatmap(subset_norm_cts, 
                                  annotation_col = heatmap_annotation,
                                  scale= "row",
                                  cluster_rows=TRUE, 
                                  show_rownames=input$significant_heatmap_showRowLabels,
                                  show_colnames = input$significant_heatmap_showColLabels,
                                  cluster_cols=TRUE)
  
  output$significant_miRNA_heatmap <- renderPlot({
    heatmap_significant
  })
  
  # Store heatmap in reactive for download handler
  significant_heatmap_reactive(heatmap_significant)
})

observeEvent(input$Volcano_Plot_Button, {  
  res <- DE_miRNA_results_reactive()
  
  # Create volcano plot
  output$volcanoPlot <- renderPlot({
    EnhancedVolcano(res,
                    lab = rownames(res),
                    labSize = input$volcano_plot_label_size,
                    x = 'log2FoldChange',
                    y = 'padj',
                    ylab = bquote(~-Log[10]~ '(p-adjusted)'),
                    pCutoff = input$DEM_padj_filter,
                    title = NULL,
                    subtitle = NULL,
                    FCcutoff = input$DEM_log2FoldChange_filter,
                    max.overlaps = Inf) +
      theme_cowplot(12) +
      theme(legend.position = "top",
            plot.subtitle = element_text(hjust = 0.5),  
            legend.justification = "center") +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  # Store volcano plot in reactive for download handler
  volcano_plot_reactive(function() {
    EnhancedVolcano(res,
                    lab = rownames(res),
                    labSize = input$volcano_plot_label_size,
                    x = 'log2FoldChange',
                    y = 'padj',
                    ylab = bquote(~-Log[10]~ '(p-adjusted)'),
                    pCutoff = input$DEM_padj_filter,
                    title = NULL,
                    subtitle = NULL,
                    FCcutoff = input$DEM_log2FoldChange_filter,
                    max.overlaps = Inf) +
      theme_cowplot(12) +
      theme(legend.position = "top",
            plot.subtitle = element_text(hjust = 0.5),  
            legend.justification = "center") +
      theme(plot.title = element_text(hjust = 0.5))
  })
})

### Do Differential miRNA Expression Analysis with limma ---

# future_expansion_code: limma::voom analysis
# observeEvent(input$runVoom, {
#   req(wideData(), 
#       metadata(), 
#       input$designColumn)
#   
#   countData <- wideData()
#   metadataData <- metadata()
#   
#   metadataData <- metadataData %>% filter(SampleID %in% colnames(countData))
#   countData <- countData[, metadataData$SampleID]
#   metadataData$SampleID <- factor(metadataData$SampleID, levels = colnames(countData))
#   metadataData <- metadataData[order(metadataData$SampleID), ]
#   
#   countData <- DGEList(counts=countData,genes=rownames(countData))
#   
#   design <- model.matrix(~ metadataData[[input$designColumn]])
#   
#   v <- voom(countData, design, plot = TRUE)
#   
#   fit <- lmFit(v, design)
#   
#   fit <- eBayes(fit)
#   
#   res <- topTable(fit, adjust = "BH", sort.by = "P", number = Inf)
#   
#   res_significant_data <- res %>% 
#     filter(adj.P.Val < 0.05) %>%
#     rownames_to_column("miRNA")
#   
#   res_significant(res_significant_data)
#   
#   output$resSignificantTable <- renderTable({
#     req(res_significant())
#     as.data.frame(res_significant())
#   })
#   
#   updateSelectInput(session, "selectedMiRNA", choices = res_significant()$miRNA)
#   
#   output$volcanoPlot <- renderPlot({
#     EnhancedVolcano(res,
#                     lab = rownames(res),
#                     labSize = 6,
#                     x = 'logFC',
#                     y = 'adj.P.Val',
#                     ylab = bquote(~-Log[10]~ '(adjusted p-value)'),
#                     pCutoff = 0.05,
#                     title = "Differential Expression Analysis",
#                     subtitle = "limma::voom",
#                     FCcutoff = 1) +
#       theme_cowplot(12) +
#       theme(legend.position = "top",
#             plot.subtitle = element_text(hjust = 0.5),
#             legend.justification = "center") +
#       theme(plot.title = element_text(hjust = 0.5))
#   })
#   
#   # Download handler for limma volcano plot
#   output$downloadVolcanoLimma <- downloadHandler(
#     filename = function() {
#       paste("Volcano_Plot_limma_voom_", Sys.Date(), ".png", sep = "")
#     },
#     content = function(file) {
#       png(file, width = 800, height = 600)
#       print(
#         EnhancedVolcano(res,
#                         lab = rownames(res),
#                         labSize = 6,
#                         x = 'logFC',
#                         y = 'adj.P.Val',
#                         ylab = bquote(~-Log[10]~ '(adjusted p-value)'),
#                         pCutoff = 0.05,
#                         title = "Differential Expression Analysis",
#                         subtitle = "limma::voom",
#                         FCcutoff = 1) +
#           theme_cowplot(12) +
#           theme(legend.position = "top",
#                 plot.subtitle = element_text(hjust = 0.5),
#                 legend.justification = "center") +
#           theme(plot.title = element_text(hjust = 0.5))
#       )
#       dev.off()
#     }
#   )
# })
