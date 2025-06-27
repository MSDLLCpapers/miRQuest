# Module 1: Exploratory Visualization

observeEvent(input$plotButton, {
  req(longData(),metadata())
  
  metadataData <- metadata()
  ld <- longData()
  print(names(ld))
  
  ld <- ld %>% 
    group_by(SampleID) %>%
    mutate(Percentage = Count / sum(Count) * 100)
  
  ld <- merge(metadataData, ld, by= "SampleID")
  print(dim(ld))
  
  ld <- ld %>% 
    group_by(SampleID) %>%
    filter(Percentage >= input$percentageFilter) %>%
    mutate(Percentage = Count / sum(Count) * 100)
  print(dim(ld))
  
  cols_large <- unlist(lapply(1:nrow(palettes_d_names), function(i) {
    package_name <- palettes_d_names$package[i]
    palette_name <- palettes_d_names$palette[i]
    num_colors <- palettes_d_names$length[i]
    
    palette_colors <- paletteer_d(paste0(package_name, "::", palette_name), num_colors)
    
    return(palette_colors)
  }))
  
  scrambled_cols <- unique(sample(cols_large))
  
  if (input$meanToggle) {
    req(input$groupBy)
    print(names(ld))
    
    ld_mean <- ld %>%
      group_by(!!sym(input$groupBy), X) %>%
      summarise(MeanPercentage = mean(Percentage)) %>% 
      ungroup() %>% 
      group_by(!!sym(input$groupBy)) %>% 
      mutate(Percentage=MeanPercentage/sum(MeanPercentage)*100)
    print(dim(ld_mean))
    
    output$miRNAPlot <- renderPlot({
      ggplot(ld_mean, aes(x = get(input$groupBy), y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1))
    })
    
    # Store plot for download
    stacked_plot_reactive(function() {
      ggplot(ld_mean, aes(x = get(input$groupBy), y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1))
    })
  } else {
    req(input$groupBy)
    
    output$miRNAPlot <- renderPlot({
      ggplot(ld, aes(x = SampleID, y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
    
    # Store plot for download
    stacked_plot_reactive(function() {
      ggplot(ld, aes(x = SampleID, y = Percentage, fill = X)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(as.formula(paste("~", input$groupBy)), scales = "free_x") +
        labs(x = "Sample", y = "miRNA (%)", title = "Relative miRNA Expression") +
        theme_cowplot(12) +
        scale_fill_manual(name="miRNA", values = scrambled_cols) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
  }
})

observeEvent(input$MDSButton, {
  
  df_wide <- wideData() 
  metadataData <- metadata()
  
  dataDGE<-DGEList(counts=df_wide,genes=rownames(df_wide))
  o <- order(rowSums(dataDGE$counts), decreasing=TRUE)
  dataDGE <- dataDGE[o,]
  
  dataNorm <- calcNormFactors(dataDGE)
  MDSdata <- plotMDS(dataNorm)
  
  MDSxy = data.frame(x=MDSdata$x, y=MDSdata$y)
  colnames(MDSxy) = c(paste(MDSdata$axislabel, '1'), paste(MDSdata$axislabel, '2'))
  mds_xy <- MDSxy
  mds_xy$Dim1 <- mds_xy$`Leading logFC dim 1`
  mds_xy$Dim2 <- mds_xy$`Leading logFC dim 2`
  mds_xy$SampleID <- gsub("X","",colnames(df_wide))
  mds_xy <- merge(mds_xy, metadataData,by="SampleID")
  
  output$MDSPlot <- renderPlot({
    req(input$groupBy)
    grouping_variable <- input$groupBy
    
    plot <- ggplot2::ggplot(mds_xy, aes(x=Dim1, y=Dim2, colour= get(grouping_variable))) + 
      geom_point(aes(fill= get(grouping_variable)), colour="black", pch=21, size=3) +
      scale_fill_viridis_d(name="") +
      ggtitle("Comparison by MDS") +
      xlab("Dim1") +
      ylab("Dim2") +
      cowplot::theme_cowplot(12) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position="top", legend.justification = "center") 
    
    if (input$showSampleID) {
      plot <- plot + geom_text(aes(label=SampleID), vjust=-1, size=3) +
        scale_fill_viridis_d(name="")
    }
    
    print(plot) 
    
    # Store in reactive for download handler
    mds_xy_reactive(mds_xy)
    
    # Store plot for download
    mds_plot_reactive(function() {
      ggplot2::ggplot(mds_xy, aes(x=Dim1, y=Dim2, colour= get(input$groupBy))) + 
        geom_point(aes(fill= get(input$groupBy)), colour="black", pch=21, size=3) +
        scale_fill_viridis_d(name="") +
        ggtitle("Comparison by MDS") +
        xlab("Dim1") +
        ylab("Dim2") +
        cowplot::theme_cowplot(16) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="top", legend.justification = "center") +
        (if (input$showSampleID) geom_text(aes(label=SampleID), vjust=-1, size=3) else NULL)
    })
  })
})

observeEvent(input$heatmapButton, { 
  countData <- wideData()
  metadataData <- metadata()
  
  metadataData <- metadataData %>% filter(SampleID %in% colnames(countData))
  countData <- countData[,metadataData$SampleID]
  metadataData$SampleID <- factor(metadataData$SampleID, levels =colnames(countData))
  metadataData <- metadataData[order(metadataData$SampleID), ]
  
  designFormula <- as.formula(paste("~", input$groupBy))
  dds <- DESeqDataSetFromMatrix(countData = countData, 
                                colData = metadataData, 
                                design = designFormula)
  
  dds<- DESeq(dds)
  norm_cts <- vst(dds, blind = TRUE, sum( rowMeans( counts(dds, normalized=TRUE)) > 5 )) %>% assay()
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:input$countThreshold]
  
  print(select)
  
  df <- as.data.frame(colData(dds)) %>% 
    dplyr::select(c(input$groupBy))
  
  if (input$rowScale && input$colScale) {
    showNotification("You can only check one scaling option (row or column).", type = "error")
    return()
  }
  
  pheatmap_params <- list(
    mat = norm_cts[select,],
    annotation_col = df,
    show_rownames = input$showRowLabels,
    show_colnames = input$showColLabels
  )
  
  if (input$rowCluster) {
    pheatmap_params$cluster_rows <- TRUE
  } else {
    pheatmap_params$cluster_rows <- FALSE
  }
  
  if (input$colCluster) {
    pheatmap_params$cluster_cols <- TRUE
  } else {
    pheatmap_params$cluster_cols <- FALSE
  }
  
  if (input$rowScale) {
    pheatmap_params$scale <- "row"
  } else if (input$colScale) {
    pheatmap_params$scale <- "column"
  } else {
    pheatmap_params$scale <- "none"
  }
  
  heatmap <- do.call(pheatmap, pheatmap_params)
  
  output$heatmap <- renderPlot({
    heatmap
  })
  
  # Download handler for heatmap
  output$downloadHeatmap <- downloadHandler(
    filename = function() {
      paste("exploratory_heatmap_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      png(file, width = 800, height = 600)
      print(heatmap)
      dev.off()
    }
  )
})

# Download handler for MDS data (outside observeEvent)
output$downloadData <- downloadHandler(
  filename = function() {
    paste("mds_xy_", Sys.Date(), ".csv", sep="")
  },
  content = function(file) {
    req(mds_xy_reactive())
    write.csv(mds_xy_reactive(), file, row.names = FALSE)
  }
)