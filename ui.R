# Define UI for application

ui <- fluidPage(
  useShinyjs(),
  #useBusyIndicators(spinners=FALSE,pulse=FALSE,fade=TRUE),
  titlePanel("miRQuest: Interactive Analysis of MicroRNA Sequencing Data"),
  theme=bs_theme(version=5, bootswatch = "cosmo"),
  
  sidebarLayout(
    sidebarPanel(
      actionButton("useDefaultData", "Use Default Data", class = "btn-primary", style = "width: 100%; margin-bottom: 20px;"),
      tags$hr(),
      h6("Upload raw miRNA counts x sample table"),
      fileInput("countTable", "Upload miRNA Count Table", 
                accept = c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
      textOutput("mirna_status"),
      tags$hr(),
      h6("Upload sample metadata; assumes first column is named 'SampleID'"),
      fileInput("metadataFile", "Upload Metadata File", 
                accept = c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
      textOutput("metadata_status"),
      tags$hr(),
      h6("Optionally, add RNA-Seq data for the same samples"),
      fileInput("mrna_countTable", "Upload RNA-Seq Count Table", 
                accept = c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
      textOutput("mrna_status"),
      tags$hr(),
      h6("Select species"),
      checkboxInput("hsa", "Homo sapiens", FALSE),
      checkboxInput("mmu", "Mus musculus", FALSE),
      tags$hr(),
      fluidRow(
        column(8, actionButton("submit", "Submit")),  # Submit button
        column(7, textOutput("error_message"))  # Error message next to the button
      ),
      width=3
    ),
    
    mainPanel(
      tabsetPanel(
        
        tabPanel("Exploratory Visualization",
                 selectInput("groupBy", "Select Metadata Column to Group Plots", choices = NULL),
                 #card(card_header("card1"),
                 actionButton("plotButton", "Show miRnome as Stacked Column Chart"),
                 checkboxInput("meanToggle", "Calculate Mean Relative Abundance", value = FALSE),  # New Checkbox
                 sliderInput("percentageFilter", 
                             "% miRNAs that comprise the sample", 
                             min = 0, 
                             max = 5, 
                             value = 1, 
                             step = 0.2),
                 plotOutput("miRNAPlot"),
                 downloadButton("downloadStackedColumn", "Download Plot"),
                 tags$hr(),
                 # ),
                 #card(card_header("card2"),
                 actionButton("MDSButton", "Show miRnome as a Single Point in 2D Space"),
                 checkboxInput("showSampleID", "Show Sample ID on MDS Plot", value = FALSE),
                 plotOutput("MDSPlot"),
                 downloadButton("downloadData", "Download Dataframe"),
                 downloadButton("downloadPlot", "Download Plot"),
                 tags$hr(),
                 #),
                 #card(card_header("card3"),
                 actionButton("heatmapButton", "Show miRnome as Heatmap"),
                 sliderInput("countThreshold", "Top n number of miRNA with the most counts", value = 20, min = 5, max=50, step = 5),
                 # Organizing checkboxes into three rows with two checkboxes per row
                 fluidRow(
                   column(6, checkboxInput("colCluster", "Cluster by Sample", value = TRUE)),
                   column(6, checkboxInput("rowCluster", "Cluster by miRNA", value = FALSE))
                 ),
                 fluidRow(
                   column(6, checkboxInput("rowScale", "Scale by miRNA", value = TRUE)),
                   column(6, checkboxInput("colScale", "Scale by Sample", value = FALSE))  
                 ),
                 fluidRow(
                   column(6, checkboxInput("showRowLabels", "Show Row Labels", value = FALSE)),  
                   column(6, checkboxInput("showColLabels", "Show Column Labels", value = FALSE))   
                 ),
                 plotOutput("heatmap"),
                 downloadButton("downloadHeatmap", "Download Heatmap"),
                 tags$hr()
                 
                 #)
        ),
        tabPanel("Data Subsetting",
                 selectInput("metadataColumn", "Select Metadata Column", choices = NULL),
                 uiOutput("filterValueInput"),  # Dynamic UI for filter values
                 actionButton("subsetDataButton", "Generate Subset"),
                 downloadButton("downloadSubset", "Download Subset as CSV"),  # Download button
                 tableOutput("subsetDataTable")  # Table to display the subsetted data
                 
        ),
        tabPanel("Differential miRNA Expression Analysis",
                 selectInput("designColumn", "Select Design Column", choices = NULL),
                 actionButton("runDESeq", "Run DESeq2 Analysis"),
                 numericInput("DEM_padj_filter", "Adjusted P-value Threshold", 
                              value = 0.05, min = 0, max = 1, step = 0.01),
                 #actionButton("showDESeqTable",  "Preview Differentially Expressed miRNA"), 
                 textOutput("dem_count"), 
                 downloadButton("download_significant_miRNA", "Download Significantly Differentially Expressed miRNA"),
                 tableOutput("resSignificantTable"), 
                 tags$hr(),
                 selectInput("selectedMiRNA", "Select Differentially Expressed miRNA", choices = NULL),
                 selectInput("groupBy_Boxplot", "Select Metadata Column For Plot Grouping", choices = NULL),
                 plotOutput("boxplot"),
                 downloadButton("downloadBoxplot", "Download Plot"),
                 tags$hr(),
                 numericInput("DEM_log2FoldChange_filter", "Absolute Log2 Fold Change Threshold", 
                              min = 0, 
                              value = 1, 
                              step = 0.5),
                 sliderInput("volcano_plot_label_size", "Label Size",
                             min = 1, 
                             max = 6, 
                             value=3,
                             step=1),
                 # future_expansion_code: limma::voom analysis button
                # actionButton("runVoom", "Run limma::voom Analysis"),  # New button for limma::voom
                 actionButton("Volcano_Plot_Button", "Show DE miRNA as Volcano Plot"),
                 plotOutput("volcanoPlot"),
                 downloadButton("downloadVolcano", "Download Volcano Plot"),
                 # future_expansion_code: limma volcano plot download
                 # downloadButton("downloadVolcanoLimma", "Download limma Volcano Plot"),
                 tags$hr(),
                 actionButton("significant_heatmap_button", "Show DE miRNA as Heatmap"),
                 checkboxInput("significant_heatmap_showRowLabels", "Show Row Labels", value = FALSE),
                 checkboxInput("significant_heatmap_showColLabels", "Show Column Labels", value = FALSE),
                 plotOutput("significant_miRNA_heatmap"),
                 downloadButton("downloadSignificantHeatmap", "Download Plot"),
                 tags$hr(),
                 # future_expansion_code: DE miRNA data table for debugging/correlation analysis
                 # h5("DE miRNA Data for Correlation Analysis:"),
                 # tableOutput("mirna_DE_table"),
                 # downloadButton("download_mirna_DE_table", "Download DE miRNA Data"),
                 # tags$hr()
                 # Plot Outputs
                 
                 
                 
        ),
        
        tabPanel("Single miRNA Target Prediction and Pathway Enrichment",
                 # Analysis related outputs
                 
                 selectInput("DEMiRNA", "Select Differentially Expressed miRNA", choices = NULL),
                 actionButton("runAnalysis", "Predict Gene Targets for a Single miRNA"),
                 downloadButton("downloadgeneResults", "Download Predicted Genes"),
                 tableOutput("genetargets"), 
                 tags$hr(),
                 selectInput("single_miRNA_analysis_type", "Select Pathway Database",
                             choices = c("GO" = "GO", "Reactome" = "Reactome"),
                             selected = "Reactome"),  
                 actionButton("runORA_single_miRNA", "Perform Pathway Overrepresentation Analysis"), 
                 downloadButton("downloadFinalResults", "Download Enriched Pathways"),
                 tableOutput("resultsTable"),
                 tags$hr(),
                 actionButton("plot_pathways_single_miRNA", "Plot Enriched Pathways"), 
              layout_columns(
                card(card_header("Dotplot"),
                 plotOutput("ora_dotplot"),
                 downloadButton("downloadDotplot", "Download Plot")
                ),
                card(card_header("Barplot"),
                 plotOutput("ora_barplot"),
                 downloadButton("downloadBarplot","Download Plot")
                )
              )
        ),
        
        tabPanel("Aggregated Pathway Enrichment of Predicted Genes",
                 # Analysis related outputs
                 actionButton("predict_all_genes", "Predict Targeted Genes for all DE miRNA"),  
                 downloadButton("all_downloadgeneResults", "Download Predicted Genes"),
                 tableOutput("all_genetargets"), 
                 tags$hr(),
                 selectInput("select_up_or_down_miRNA", "Select Up- or Down-regulated miRNA",
                             choices = c("Up" = "Up", "Down" = "Down"),
                             selected = "Up"),
                 selectInput("analysis_type", "Select Pathway Database",
                             choices = c("GO" = "GO", "Reactome" = "Reactome"),
                             selected = "Reactome"),  
                 actionButton("runORA", "Perform Pathway Overrepresentation Analysis"),
                 downloadButton("all_downloadFinalResults", "Download Enriched Pathways"),
                 tableOutput("all_resultsTable"),
                 tags$hr(),
                 
                 
                 # New inputs for user-defined parameters
                 h6("The number of pathways shown will be affected by coverage selection."),
                 numericInput("num_pathways", "Maximum Number of Pathways to Show:", value = 10, min = 2),
                 h6("Coverage indicates the minimum fraction of genes in a pathway an miRNA must target to be shown."),
                 numericInput("min_coverage", "Minimum Coverage (0 to 1):", value = 0.1, min = 0, max = 1, step = 0.01),
                 #absolutePanel(draggable=TRUE,
                 layout_columns(
                   card(card_header("Chord Plot"),
                        actionButton("showChordPlot", "Show Chord Plot"),
                        tags$div(
                          style = "overflow: visible; width: 100%; height: 600px;",
                          plotOutput("chordPlot", width = "100%", height = "600px")),
                        downloadButton("downloadChord", "Download Plot")
                        #    tags$div(
                        # style = "overflow: visible; width: 100%; height: 600px;", 
                        # plotOutput("chordPlot", width = "100%", height = "600px")),
                   ),
                   # New action buttons for plotting
                   card(card_header("Network Plot"),
                        actionButton("showNetworkPlot", "Show Network Plot"),
                        plotOutput("networkPlot"),
                        downloadButton("downloadNetwork", "Download Plot")
                   )
                 )
                 # )   
        ),
        tabPanel("Differential Gene Expression Analysis",
                 h6("Uses the same model design as is selected for miRNA"),
                 actionButton("mrna_runDESeq", "Run DESeq2 Analysis"),  
                 numericInput("padj_filter", "Adjusted P-value Threshold", 
                              value = 0.05, min = 0, max = 1, step = 0.01),
                 textOutput("deg_count"),  # Output for displaying the number of DEGs
                 downloadButton("download_DEG_Table", "Download Differentially Expressed Genes"),
                 tableOutput("mrna_resSignificantTable"), 
                 tags$hr()
        ),
        
        tabPanel("miRNA - Gene correlation analysis",
                 fluidRow(
                   column(6, numericInput("filter_DEG_correlation", "mRNA Log2 Fold Change Threshold", 
                                          value = 1, min = 0, step = 0.1)),
                   column(6, numericInput("filter_DEM_correlation", "miRNA Log2 Fold Change Threshold", 
                                          value = 1, min = 0, step = 0.1))
                 ),
                 textOutput("DEG_DEM_filter_message"),
                 tags$hr(),
                 actionButton("runCorrelationAnalysis", "Run Correlation Analysis"), 
                 textOutput("num_neg_correlations"),
                 downloadButton("downloadNegCor", "Download Negative Correlations"),
                 tableOutput("NegativeCorrelationTable"),
                 tags$hr(),
                 actionButton("show_database_support", "See Correlations with Database Support"), 
                 textOutput("degMessage"),
                 downloadButton("downloadsupportednegcor", "Download Supported Negative Correlations"),
                 tableOutput("SupportedNegativeCorrelationTable"),
                 tags$hr(),
                 actionButton("CorrelationPlotButton", "Show Correlation Plot"),
                 checkboxInput("select_supported_only", 
                               label = "Show Supported miRNA-Gene Pairs Only", 
                               value = FALSE),
                 fluidRow(
                   column(4, numericInput("correlationCutoff", "Correlation Cutoff:", 
                                          value = -0.5, min = -1, max = 0, step = 0.1)),  
                   column(4, numericInput("logratio_miRNA_Cutoff", "Absolute DE miRNA log2FoldChange Cutoff:", 
                                          value = 1, min = 0, step = 1)),  
                   column(4, numericInput("logratio_gene_Cutoff", "Absolute DE Gene log2FoldChange Cutoff", 
                                          value = 1, min = 0, step = 1))  
                 ),
                 plotOutput("correlationPlot"),
                 downloadButton("downloadCorrelationPlot", "Download Plot")
        )
      )
    )
  )
)