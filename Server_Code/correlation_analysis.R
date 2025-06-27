# Module 7: miRNA-mRNA Correlation Analysis

mrna_no_zero_reactive <- reactive({
  req(mrna_data_reactive())
  mrna_data <- as.data.frame(mrna_data_reactive())
  last_five_columns <- (ncol(mrna_data) - 4):ncol(mrna_data)
  mrna_data[rowSums(mrna_data[, -last_five_columns, drop = FALSE]) != 0, ] %>% 
    filter(abs(log_ratio) >= input$filter_DEG_correlation)
})

mirna_no_zero_reactive <- reactive({
  req(mirna_data_reactive())
  mirna_data <- as.data.frame(mirna_data_reactive())
  last_five_columns <- (ncol(mirna_data) - 4):ncol(mirna_data)
  mirna_data[rowSums(mirna_data[, -last_five_columns, drop = FALSE]) != 0, ] %>% 
    filter(abs(log_ratio) >= input$filter_DEM_correlation)
})

output$DEG_DEM_filter_message <- renderText({
  mrna_no_zero <- mrna_no_zero_reactive()
  mirna_no_zero <- mirna_no_zero_reactive()
  
  num_mrna_remaining <- nrow(mrna_no_zero)
  num_mirna_remaining <- nrow(mirna_no_zero)
  
  paste("This will correlate", num_mrna_remaining, "DEGs with",
        num_mirna_remaining, "DE miRNAs")
})

observeEvent(input$runCorrelationAnalysis, {
  if (is.null(mirna_data_reactive()) || is.null(mrna_data_reactive())) {
    showNotification(
      HTML("Please complete the following steps first:<br>
            1. Run Module 3 (Differential miRNA Expression Analysis)<br>
            2. Run Module 6 (Differential Gene Expression Analysis)<br>
            Both analyses must be completed before running correlation analysis."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  req(mirna_no_zero_reactive(), mrna_no_zero_reactive())
  
  mrna_no_zero <- mrna_no_zero_reactive()
  mirna_no_zero <- mirna_no_zero_reactive()
  
  cor <- negative_cor(mrna_data = mrna_no_zero, 
                      mirna_data = mirna_no_zero, 
                      cut.off = 0, 
                      method = "spearman")
  cor <- as.data.frame(cor)
  cor$Correlation <- as.numeric(cor$Correlation)
  
  neg_cor <- cor
  neg_cor <- dplyr::rename(neg_cor,ensembl_gene_id = Gene)
  neg_cor$ensembl_gene_id <- gsub("\\..*","",neg_cor$ensembl_gene_id)
  
  hsmart <- mart()
  
  mapping <- getBM(
    attributes = c('ensembl_gene_id',org_string(),'description'), 
    filters = 'ensembl_gene_id',
    values = neg_cor$ensembl_gene_id,
    mart = hsmart
  )
  
  neg_cor <- merge(neg_cor, mapping, by="ensembl_gene_id")
  neg_cor_reactive(neg_cor)
  
  output$num_neg_correlations <- renderText({
    if (!is.null(neg_cor)) {
      paste("The number of miRNA-gene correlations is",nrow(neg_cor))
    } else {
      paste("No miRNA-gene correlations to show!")
      return()
    }
  })
  
  output$NegativeCorrelationTable <- renderTable({
    if (!is.null(neg_cor)) {
      head(neg_cor)
    } else {
      return()
    }
  })
  
})

observeEvent(input$show_database_support, {
  req(neg_cor_reactive(),predicted_target_reactive())
  
  predicted <- as.data.frame(predicted_target_reactive()) %>%
    as_tibble() %>%
    dplyr::select(ensembl_gene_id,!!sym(org_string()),miRNA)
  neg_cor_pairs <- as.data.frame(neg_cor_reactive()) %>% 
    as_tibble() %>%
    dplyr::select(ensembl_gene_id,!!sym(org_string()), miRNA)
  
  identical_rows <- inner_join(predicted,neg_cor_pairs, relationship = "many-to-many")
  supported_neg_cor_reactive(identical_rows)
  
  if ( nrow(identical_rows) == 0) {
    output$degMessage <- renderText("None of the DEGs are supported by database predictions")
  } else {
    output$degMessage <- renderText(paste("The number of supported miRNA-gene pairs is", nrow(identical_rows)))
  }
  
  output$SupportedNegativeCorrelationTable <- renderTable({
    if (!is.null(identical_rows)) {
      identical_rows <- identical_rows %>% filter(!!sym(org_string())!="")
      head(identical_rows)
    } else {
      return()
    }
  })
  
})

observeEvent(input$CorrelationPlotButton, {
  req(neg_cor_reactive(),supported_neg_cor_reactive())
  neg_cor <- neg_cor_reactive()
  supported_neg_cor <- supported_neg_cor_reactive()
  
  if (input$select_supported_only) { 
    neg_cor <- neg_cor %>%
      semi_join(supported_neg_cor, by = c("ensembl_gene_id", org_string(), "miRNA"))
  }
  
  miRNAs_with_high_correlation <- neg_cor %>%
    mutate(abs_Correlation = abs(Correlation)) %>%
    filter(abs_Correlation >= abs(input$correlationCutoff)) %>%
    mutate(abs_logratio_miRNA = abs(as.numeric(logratio_miRNA))) %>%
    filter(abs_logratio_miRNA >= input$logratio_miRNA_Cutoff) %>%
    pull(miRNA) %>%
    unique()
  
  genes_filtered <- neg_cor %>%
    mutate(abs_logratio_gene = abs(as.numeric(logratio_gene))) %>%
    filter(abs_logratio_gene >= input$logratio_gene_Cutoff) %>%
    pull(!!sym(org_string()))
  
  filtered_neg_cor <- neg_cor %>%
    filter(miRNA %in% miRNAs_with_high_correlation) %>% 
    filter(!!sym(org_string()) %in% genes_filtered)
  
  mean_correlations <- filtered_neg_cor %>%
    group_by(miRNA) %>%
    summarise(Mean_Correlation = mean(Correlation)) %>%
    arrange(Mean_Correlation)
  
  filtered_neg_cor$miRNA <- factor(filtered_neg_cor$miRNA, levels = mean_correlations$miRNA)
  
  plot <- ggplot(filtered_neg_cor, aes(x = !!sym(org_string()), y = miRNA, size = abs(Correlation), color = Correlation)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = "Negative Correlation between miRNA and Gene",
         x = "Gene",
         y = "miRNA",
         size = "Correlation (absolute value)",
         color = "Correlation Value") +
    theme_cowplot(12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  output$correlationPlot <- renderPlot({
    plot
  })
  
  correlation_plot_reactive(plot)
})

# Download handlers (outside observeEvents)
output$downloadNegCor <- downloadHandler(
  filename = function() {
    paste("neg_cor_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    req(neg_cor_reactive())
    write.csv(neg_cor_reactive(), file, row.names = FALSE)
  }
)

output$downloadsupportednegcor <- downloadHandler(
  filename = function() {
    paste("Supported_Negative_Correlations", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(supported_neg_cor_reactive())
    write.csv(supported_neg_cor_reactive(), file, row.names = FALSE, quote = FALSE)
  }
)