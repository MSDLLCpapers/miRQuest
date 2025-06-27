# Define server logic
server <- function(input, output, session) {
  
  ### Reactive values, organized by module ---
  
  # Module 0A: Read and Store MiRNA Counts Data
  longData <- reactiveVal() 
  wideData <- reactiveVal()  
  
  # Module 0B: Read and Store Metadata
  metadata <- reactiveVal()  
  
  
  # Module 0C: Assess species selection and build mart objects accordingly
  mart <- reactive({
    validate(
      need(input$hsa || input$mmu, "Please select at least one species."),
      need(!(input$hsa && input$mmu), "Please select only one species at a time.")
    )
    if (input$hsa) {
      tryCatch({
        useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
      }, warning = function(w) {
        showNotification(paste("Cannot contact Ensembl, please try again later: ", w$message), type = "error")
        return(NULL)
      }, error = function(e) {
        showNotification(paste("Cannot contact Ensembl, please try again later: ", e$message), type = "error")
        return(NULL)
      })
    } else if (input$mmu) {
      tryCatch({
        useEnsembl(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl", 
                   mirror = "useast")
      }, warning = function(w) {
        showNotification(paste("Cannot contact Ensembl, please try again later", w$message), type = "error")
        return(NULL)
      }, error = function(e) {
        showNotification(paste("Cannot contact Ensembl, please try again later", e$message), type = "error")
        return(NULL)
      })
    } else {
      return(NULL)
    }
  })
  
  org_string <- reactive({
    validate(
      need(input$hsa || input$mmu, "Please select at least one species."),
      need(!(input$hsa && input$mmu), "Please select only one species at a time.")
    )
    if (input$hsa) {
      return("hgnc_symbol")
    } else if (input$mmu) {
      return("external_gene_name")
    } else {
      return(NULL)
    }
  })
  
  org_abbrev <- reactive({
    validate(
      need(input$hsa || input$mmu, "Please select at least one species."),
      need(!(input$hsa && input$mmu), "Please select only one species at a time.")
    )
    if (input$hsa) {
      return("hsa")
    } else if (input$mmu) {
      return("mmu")
    } else {
      return(NULL)
    }
  })
  
  org_db <- reactive({
    validate(
      need(input$hsa || input$mmu, "Please select at least one species."),
      need(!(input$hsa && input$mmu), "Please select only one species at a time.")
    )
    if (input$hsa) {
      return(org.Hs.eg.db)
    } else if (input$mmu) {
      return(org.Mm.eg.db)
    } else {
      return(NULL)
    }
  })
  
  organism_reactive <- reactive({
    validate(
      need(input$hsa || input$mmu, "Please select at least one species.")
    )
    if (input$hsa) {
      return("human")  
    } else if (input$mmu) {
      return("mouse")  
    } else {
      return(NULL)
    }
  })
  
  # Module 1: Exploratory Visualization
  boxplot_reactive <- reactiveVal()
  mds_xy_reactive <- reactiveVal()
  stacked_plot_reactive <- reactiveVal()
  mds_plot_reactive <- reactiveVal()
  heatmap_reactive <- reactiveVal()
  volcano_plot_reactive <- reactiveVal()
  significant_heatmap_reactive <- reactiveVal()
  correlation_plot_reactive <- reactiveVal()
  chord_plot_reactive <- reactiveVal()
  network_plot_reactive <- reactiveVal()
  dotplot_reactive <- reactiveVal()
  barplot_reactive <- reactiveVal()
  
  # Module 2: Data Subsetting
  subset_data_reactive <- reactiveVal()
  
  # Module 3: Do Differential miRNA Expression Analysis with DESEQ2
  DESEQ_obj <- reactiveVal()
  DE_miRNA_results_reactive <- reactiveVal()
  res_significant <- reactiveVal()
  norm_cts_reactive <- reactiveVal()
  heatmap_annotation_reactive <- reactiveVal()
  mirna_data_reactive <- reactiveVal()
  
  # Module 4: Predict target genes for all DE miRNA
  predicted_target_reactive <- reactiveVal()
  upregulated_miRNA_genes_reactive <- reactiveVal()
  downregulated_miRNA_genes_reactive <- reactiveVal()
  upregulated_miRNA_reactive <- reactiveVal()
  downregulated_miRNA_reactive <- reactiveVal()
  All_miRNA_Pathways_reactive <- reactiveVal()
  miRNA_mapped_pathways <- reactiveVal()
  
  # Module 5: Predict target genes for a single DE miRNA
  Single_miRNA_PredictedGenes_reactive <- reactiveVal()
  Single_miRNA_Pathways_reactive <- reactiveVal()
  
  # Module 6: Do Differential Gene Expression Analysis with DESEQ2
  mrna_res_significant <- reactiveVal()
  mrna_data_reactive <- reactiveVal()
  mrna_wideData <- reactiveVal()
  
  # Module 7: miRNA-mRNA Correlation Analysis
  neg_cor_reactive <- reactiveVal()
  supported_neg_cor_reactive <- reactiveVal()
  
  # Species validation observer
  observeEvent(input$submit, {
    output$error_message <- renderText({
      if (input$hsa && input$mmu) {
        return("Please select only one species at a time.")
      } else if (input$hsa || input$mmu) {
        return("Success!")
      } else {
        return("Please select at least one species.")
      }
    })
  })
  
  # Source modular server files
  files_to_source <- list.files("Server_Code", pattern = "\\.R$", ignore.case = TRUE)
  
  cat("\n\n IN SERVER.R. Sourcing ... \n\n")
  
  for (file in files_to_source) {
    print(file)
    source(
      paste0("Server_Code/", file),
      local = TRUE
    )
  }
  
  # Auto-click the "Use Default Data" button when the app starts
  observe({
    shinyjs::delay(100, {
      shinyjs::click("useDefaultData")
    })
  }) %>% bindEvent(TRUE, once = TRUE)
}