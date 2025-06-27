# Module 4: Predict target genes for all DE miRNA

observeEvent(input$predict_all_genes, {
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
  mirVector <- res_significant()$miRNA
  
  ensembl_notification <- showNotification("Contacting Ensembl website...", type = "message", duration = NULL)
  
  allPredictedGenes <- c()
  geneResults <- list()
  
  hsmart <- mart()
  
  if (is.null(hsmart)) {
    removeNotification(ensembl_notification)
    return()
  }
  
  removeNotification(ensembl_notification)
  
  withProgress(message = 'Searching databases', value = 0, {
    n <- length(mirVector)
    
    for (i in 1:n) {
      mir <- mirVector[i]
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
          next
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
      
      df <- merge(df, mapping, by="entrezgene_id")
      
      df$miRNA <- mir
      
      geneResults[[mir]] <- df[,c("entrezgene_id",
                                  "rank_final",
                                  "ensembl_gene_id",
                                  "description",
                                  org_string(),
                                  "miRNA")]
      
      incProgress(1/n, detail = paste("Getting predicted targets for DE miRNA #", i))
    }
  })
  
  finalgeneResults <- do.call(rbind,geneResults)
  predicted_target_reactive(finalgeneResults)
  
  output$all_genetargets <- renderTable({
    if (!is.null(geneResults)) {
      head(finalgeneResults)
    } else {
      return()
    }
  })
  
  # Final results already stored in predicted_target_reactive
  
  # Determine which fold change column to use (DESeq2 vs limma)
  res_sig_df <- res_significant() %>% as.data.frame()
  fc_column <- if ("log2FoldChange" %in% colnames(res_sig_df)) "log2FoldChange" else "logFC"
  
  upregulated_miRNA <- res_sig_df %>% 
    filter(!!sym(fc_column) > 0 ) %>%
    dplyr::select(miRNA, !!fc_column)
  upregulated_miRNA <- merge(upregulated_miRNA, finalgeneResults, by= "miRNA")
  upregulated_miRNA_reactive(upregulated_miRNA)
  upregulated_miRNA_genes_reactive(unique(upregulated_miRNA$entrezgene_id))
  
  downregulated_miRNA <- res_sig_df %>% 
    filter(!!sym(fc_column) < 0 ) %>%
    dplyr::select(miRNA, !!fc_column)
  downregulated_miRNA <- merge(downregulated_miRNA, finalgeneResults, by= "miRNA")
  downregulated_miRNA_reactive(downregulated_miRNA)
  downregulated_miRNA_genes_reactive(unique(downregulated_miRNA$entrezgene_id))
})

observeEvent(input$runORA, {
  req((input$hsa || input$mmu), 
      upregulated_miRNA_reactive(), 
      downregulated_miRNA_reactive())
  
  pathway_notification <- showNotification("Looking for enriched pathways. This may take some time", type = "message", duration = NULL)
  
  if (input$select_up_or_down_miRNA == "Up") {
    allPredictedGenes <- upregulated_miRNA_genes_reactive()
  } else if (input$select_up_or_down_miRNA == "Down") {
    allPredictedGenes <- downregulated_miRNA_genes_reactive()
  }
  
  print(head(allPredictedGenes))
  
  background_genes <- keys(org_db(), keytype = "ENTREZID")
  
  if (input$analysis_type == "GO") {
    ora_results <- enrichGO(gene = allPredictedGenes,
                            universe = background_genes,
                            OrgDb = org_db(),
                            keyType = "ENTREZID",
                            ont = "BP",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE)
  } else if (input$analysis_type == "Reactome") {
    ora_results <- enrichPathway(gene = allPredictedGenes,
                                 organism = organism_reactive(),
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
  }
  
  All_miRNA_Pathways_reactive(ora_results %>% as.data.frame())
  
  target_genes <- ora_results@result$geneID
  pathways <- ora_results@result$ID
  pathway_description <- ora_results@result$Description
  fold_enrichment <- ora_results@result$FoldEnrichment
  q_value <- ora_results@result$qvalue
  
  target_genes_list <- strsplit(target_genes, "/")
  
  pathway_genes_df <- data.frame(
    pathway = rep(pathways, times = sapply(target_genes_list, length)),
    description = rep(pathway_description, times = sapply(target_genes_list, length)),
    enrichment = rep(fold_enrichment, times = sapply(target_genes_list, length)),
    q_value = rep(q_value, times = sapply(target_genes_list, length)),
    gene = unlist(target_genes_list)
  )
  pathway_genes_df <- pathway_genes_df %>% filter(q_value <0.05)
  gene_count <- pathway_genes_df %>%
    group_by(pathway) %>%
    summarise(total_genes = n_distinct(gene))
  
  if (input$select_up_or_down_miRNA == "Up") {
    mirna_gene_results <- upregulated_miRNA_reactive()
  } else if (input$select_up_or_down_miRNA == "Down") {
    mirna_gene_results <- downregulated_miRNA_reactive()
  }
  
  print(head(allPredictedGenes))
  
  miRNA_df <- as.data.frame(mirna_gene_results) %>% dplyr::select(org_string(),"miRNA")
  mapped_df <- miRNA_df %>%
    inner_join(pathway_genes_df, by = setNames("gene", org_string()), relationship = "many-to-many")
  mapped_df <- mapped_df %>%
    group_by(miRNA, pathway) %>%
    mutate(gene_count = n()) %>%
    ungroup()
  mapped_df <- mapped_df %>% 
    inner_join(gene_count, by = "pathway")
  mapped_df <- mapped_df %>%
    group_by(miRNA, pathway) %>%
    mutate(coverage = gene_count/total_genes) %>%
    ungroup()
  
  miRNA_mapped_pathways(mapped_df)
  
  output$all_resultsTable <- renderTable({
    if (!is.null(ora_results)) {
      head(ora_results)
    } else {
      return()
    }
  })
  
  removeNotification(pathway_notification)
  
})

# Download handlers (outside observeEvents)
output$all_downloadgeneResults <- downloadHandler(
  filename = function() {
    paste("gene_results", Sys.Date(), ".csv", sep = "_")
  },
  content = function(file) {
    req(predicted_target_reactive())
    write.csv(predicted_target_reactive(), file, row.names = FALSE, quote = FALSE)
  }
)

output$all_downloadFinalResults <- downloadHandler(
  filename = function() {
    paste("ora_results_all_predicted_genes_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    req(All_miRNA_Pathways_reactive())
    data <- All_miRNA_Pathways_reactive()
    
    # Handle NULL or empty results
    if (is.null(data) || (is.data.frame(data) && nrow(data) == 0)) {
      # Create empty data frame with appropriate columns
      data <- data.frame(
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
    } else if (!is.data.frame(data)) {
      data <- as.data.frame(data)
    }
    
    write.csv(data, file, row.names = FALSE, quote = FALSE)
  }
)