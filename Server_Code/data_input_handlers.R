## Lines 15- 19 of this script are taken from https://github.com/nf-core/smrnaseq/blob/master/bin/edgeR_miRBase.r (MIT license).
## If you use this software, please also cite: 10.5281/zenodo.3456879 (nf-core/smrnaseq)
# Handler for "Use Default Data" button
observeEvent(input$useDefaultData, {
  # Load miRNA data
  mirna_path <- here::here("inst", "Data", "COREAD_Downsampled_Tumor_vs_Normal_Mature_miRNA.csv")
  data <- read.csv(mirna_path, header = TRUE)
  
  df_wide <- data
  df_wide[is.na(df_wide)] <- 0
  row.names(df_wide) <- df_wide$X
  df_wide <- df_wide %>% 
    dplyr::select(-c(X))

  row_sub = apply(df_wide, 1, function(row) all(row == 0))
  df_wide <- df_wide[!row_sub, , drop = FALSE]
  
  drop_colsum_zero <- (colSums(df_wide, na.rm = TRUE) != 0)
  df_wide <- df_wide[, drop_colsum_zero]
  
  wideData(df_wide)
  
  ld <- data %>%
    pivot_longer(cols = -1, names_to = "SampleID", values_to = "Count")
  
  longData(ld)
  output$mirna_status <- renderText("Demo data loaded (COREAD)")
  
  # Load metadata
  metadata_path <- here::here("inst", "Data", "COREAD_Tumor_vs_Normal_Metadata.csv")
  metadataData <- read.csv(metadata_path, header = TRUE)
  metadata(metadataData)
  updateSelectInput(session, "designColumn", choices = names(metadataData))
  updateSelectInput(session, "groupBy", choices = names(metadataData))
  updateSelectInput(session, "groupBy_Boxplot", choices = names(metadataData))
  updateSelectInput(session, "metadataColumn", choices = names(metadataData))
  output$metadata_status <- renderText("Demo data loaded (COREAD)")
  
  # Load mRNA data
  mrna_path <- here::here("inst", "Data", "COREAD_Downsampled_Tumor_vs_Normal_mRNA.csv")
  mrna_data <- read.csv(mrna_path, header = TRUE)
  mrna_wideData(mrna_data)
  output$mrna_status <- renderText("Demo data loaded (COREAD)")
  
  # Set species to human for COREAD data
  updateCheckboxInput(session, "hsa", value = TRUE)
  updateCheckboxInput(session, "mmu", value = FALSE)
  
  # Select tissue_type in all metadata dropdowns if it exists
  if ("tissue_type" %in% names(metadataData)) {
    updateSelectInput(session, "designColumn", selected = "tissue_type")
    updateSelectInput(session, "groupBy", selected = "tissue_type")
    updateSelectInput(session, "groupBy_Boxplot", selected = "tissue_type")
    updateSelectInput(session, "metadataColumn", selected = "tissue_type")
    
    # Click the plot button to show stacked bar chart
    shinyjs::delay(500, {
      shinyjs::click("plotButton")
    })
  }
  
  showNotification("Default COREAD data loaded successfully!", type = "message", duration = 3)
  
  # Show workflow guidance
  showNotification(
    HTML("Next steps:<br>
          1. Go to 'Differential miRNA Expression Analysis' tab and click 'Run DESeq2 Analysis'<br>
          2. For correlation analysis, also run 'Differential Gene Expression Analysis'"),
    type = "message",
    duration = 10
  )
})

# Module 0A: Read and Store MiRNA Counts Data
observeEvent(input$countTable, {
  req(input$countTable)
  
  data <- read.csv(input$countTable$datapath, header = TRUE)
  
  if (is.null(data)) return()
  
  df_wide <- data
  df_wide[is.na(df_wide)] <- 0
  row.names(df_wide) <- df_wide$X
  df_wide <- df_wide %>% 
    dplyr::select(-c(X))
  
  row_sub = apply(df_wide, 1, function(row) all(row == 0))
  df_wide <- df_wide[!row_sub, , drop = FALSE]
  
  drop_colsum_zero <- (colSums(df_wide, na.rm = TRUE) != 0)
  df_wide <- df_wide[, drop_colsum_zero]
  
  wideData(df_wide)
  
  ld <- data %>%
    pivot_longer(cols = -1, names_to = "SampleID", values_to = "Count")
  
  longData(ld)
  
  # Update status to show user data is loaded
  output$mirna_status <- renderText("User data loaded")
})

# Module 0B: Read and Store Metadata
observeEvent(input$metadataFile, {
  req(input$metadataFile)
  
  metadataData <- read.csv(input$metadataFile$datapath, header = TRUE)
  
  metadata(metadataData)
  
  updateSelectInput(session, "designColumn", choices = names(metadataData))
  updateSelectInput(session, "groupBy", choices = names(metadataData))
  updateSelectInput(session, "groupBy_Boxplot", choices = names(metadataData))
  updateSelectInput(session, "metadataColumn", choices = names(metadataData))
  
  # Update status to show user data is loaded
  output$metadata_status <- renderText("User data loaded")
})

# Module 0C: Read and Store mRNA Data
observeEvent(input$mrna_countTable, {
  req(input$mrna_countTable)
  
  data <- read.csv(input$mrna_countTable$datapath, header = TRUE)
  
  if (is.null(data)) return()
  
  mrna_wideData(data)
  
  # Update status to show user data is loaded
  output$mrna_status <- renderText("User data loaded")
})
