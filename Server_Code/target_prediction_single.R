# Module 5: Predict target genes for a single DE miRNA

observeEvent(input$runAnalysis, {
  if (is.null(res_significant()) || nrow(res_significant()) == 0) {
    showNotification(
      HTML("Please run Module 3 (Differential miRNA Expression Analysis) first.<br>
            Click on the 'Differential miRNA Expression Analysis' tab and run DESeq2 analysis."),
      type = "warning",
      duration = 10
    )
    return()
  }
  
  req(res_significant())
  
  gene_notification <- showNotification("Retrieving predicted targets...", type = "message", duration = NULL)
  
  mir <- input$DEMiRNA
  
  hsmart <- mart()
  
  predictions <- getPredictedTargets(mir, 
                                     species = org_abbrev(), 
                                     method = 'geom', 
                                     min_src = 2)
  
  if (is.null(predictions) || nrow(predictions) == 0) {
    mir_no_prefix <- sub(paste0("^",org_abbrev(),"-"), "", mir)
    message(paste("Trying again without prefix for miRNA:", mir_no_prefix))
    
    predictions <- getPredictedTargets(mir_no_prefix, 
                                       species = org_abbrev(), 
                                       method = 'geom', 
                                       min_src = 2)
    
    if (is.null(predictions) || nrow(predictions) == 0) {
      message(paste("No targets found for miRNA (after removing prefix):", mir_no_prefix))
      return()
    }
  }
  
  df <- as.data.frame(predictions)
  
  df$entrezgene_id <- row.names(df)
  mapping <- getBM(
    attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_string()), 
    filters = 'entrezgene_id',
    values = df$entrezgene_id,
    mart = hsmart
  )
  
  df <- merge(df, mapping, by = "entrezgene_id")
  
  df$miRNA <- mir
  
  finalgeneResults <- df[, c("entrezgene_id",
                             "rank_final",
                             "ensembl_gene_id",
                             "description",
                             org_string(),
                             "miRNA")]
  
  allPredictedGenes <- unique(finalgeneResults$entrezgene_id)
  Single_miRNA_PredictedGenes_reactive(allPredictedGenes)
  
  output$genetargets <- renderTable({
    if (!is.null(finalgeneResults)) {
      head(finalgeneResults)
    } else {
      return()
    }
  })
  
  removeNotification(gene_notification)
})

observeEvent(input$runORA_single_miRNA, {
  req(Single_miRNA_PredictedGenes_reactive())
  
  pathway_notification <- showNotification("Looking for enriched pathways...", type = "message", duration = NULL)
  
  Single_miRNA_Predicted_Genes <- Single_miRNA_PredictedGenes_reactive()
  
  background_genes <- keys(org_db(), keytype = "ENTREZID")
  
  if (input$single_miRNA_analysis_type == "GO") {
    ora_results <- enrichGO(gene = Single_miRNA_Predicted_Genes,
                            universe = background_genes,
                            OrgDb = org_db(),
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE)
  } else if (input$single_miRNA_analysis_type == "Reactome") {
    ora_results <- enrichPathway(gene = Single_miRNA_Predicted_Genes,
                                 organism = organism_reactive(),
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
  }
  
  Single_miRNA_Pathways_reactive(ora_results)
  
  output$resultsTable <- renderTable({
    if (!is.null(ora_results)) {
      head(ora_results)
    } else {
      return()
    }
  })
  
  removeNotification(pathway_notification)
})

observeEvent(input$plot_pathways_single_miRNA, {
  
  ora_results <- Single_miRNA_Pathways_reactive()
  
  output$ora_dotplot <- renderPlot({
    if (!is.null(ora_results)) {
      plot <- dotplot(ora_results, showCategory = 10) +
        ggtitle("Dotplot of Overrepresented Terms")
      dotplot_reactive(plot)
      plot
    }
  })
  
  output$ora_barplot <- renderPlot({
    if (!is.null(ora_results)) {
      plot <- barplot(ora_results, showCategory = 10) +
        ggtitle("Barplot of Overrepresented Terms")
      barplot_reactive(plot)
      plot
    }
  })
})

# Download handlers (outside observeEvents)  
output$downloadgeneResults <- downloadHandler(
  filename = function() {
    req(input$DEMiRNA)
    paste("miRNA", input$DEMiRNA, "predicted_gene_results", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(input$DEMiRNA)
    mir <- input$DEMiRNA
    
    hsmart <- mart()
    predictions <- getPredictedTargets(mir, species = org_abbrev(), method = 'geom', min_src = 2)
    
    if (is.null(predictions) || nrow(predictions) == 0) {
      mir_no_prefix <- sub(paste0("^",org_abbrev(),"-"), "", mir)
      predictions <- getPredictedTargets(mir_no_prefix, species = org_abbrev(), method = 'geom', min_src = 2)
    }
    
    if (!is.null(predictions) && nrow(predictions) > 0) {
      df <- as.data.frame(predictions)
      df$entrezgene_id <- row.names(df)
      mapping <- getBM(
        attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', org_string()), 
        filters = 'entrezgene_id',
        values = df$entrezgene_id,
        mart = hsmart
      )
      df <- merge(df, mapping, by = "entrezgene_id")
      df$miRNA <- mir
      finalgeneResults <- df[, c("entrezgene_id", "rank_final", "ensembl_gene_id", "description", org_string(), "miRNA")]
      write.csv(finalgeneResults, file, row.names = FALSE, quote = FALSE)
    }
  }
)

output$downloadFinalResults <- downloadHandler(
  filename = function() {
    req(input$DEMiRNA)
    paste("miRNA", input$DEMiRNA, "pathway_results", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(Single_miRNA_Pathways_reactive())
    ora_results <- Single_miRNA_Pathways_reactive()
    
    # Handle different types of results
    if (is.null(ora_results)) {
      # Create empty data frame with appropriate columns
      ora_df <- data.frame(
        ID = character(),
        Description = character(),
        GeneRatio = character(),
        BgRatio = character(),
        pvalue = numeric(),
        p.adjust = numeric(),
        qvalue = numeric(),
        geneID = character(),
        Count = integer()
      )
    } else if (inherits(ora_results, "enrichResult")) {
      # Extract result slot from S4 object
      ora_df <- ora_results@result
    } else {
      # Already a data frame
      ora_df <- as.data.frame(ora_results)
    }
    
    write.csv(ora_df, file, row.names = FALSE, quote = FALSE)
  }
)