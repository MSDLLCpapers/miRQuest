# Future expansion code - commented out code from original server.R

# Reactive values for table visibility
# showDESeqTable <- reactiveVal(FALSE)
# genetable <- reactiveVal(FALSE)
# mrna_showDESeqTable <- reactiveVal(FALSE)
# show_database_support <- reactiveVal(FALSE)

# Observe the button to toggle visibility of tables
# observeEvent(input$show_database_support, {
#   show_database_support(!show_database_support())
# })
# observeEvent(input$showDESeqTable, {
#   showDESeqTable(!showDESeqTable())
# })
# observeEvent(input$mrna_showDESeqTable, {
#   mrna_showDESeqTable(!mrna_showDESeqTable())
# })
# observeEvent(input$genetable, {
#   genetable(!genetable())
# })

# Alternative data reading examples
# data <- read.csv("Tumor_vs_Normal_TCGA_COAD_miRNA_counts.csv")
# metadataData <- read.csv("Tumor_vs_Normal_TCGA_COAD_metadata.csv")
# countData <- read.csv("COREAD_Downsampled_Tumor_vs_Normal_Mature_miRNA.csv", row.names=1)
# countData <- read.csv("syn_miRNA_raw_counts.csv")
# metadataData <- read.csv("Metadata.csv")
# mRNAcountData <- read.csv("COREAD_Downsampled_Tumor_vs_Normal_mRNA.csv", row.names=1)
# metadataData <- read.csv("COREAD_Tumor_vs_Normal_Metadata.csv")

# Commented out filtering logic
# Identify factor levels 
# unique_groups <- unique(metadataData[[input$designColumn]])
# 
# Calculate the lengths of each group
# group1_length <- sum(metadataData[[input$designColumn]] == unique_groups[1])
# group2_length <- sum(metadataData[[input$designColumn]] == unique_groups[2])
# 
# Assign smallest_group_size to the smaller of the two lengths
# smallest_group_size <- min(group1_length, group2_length)
# print(smallest_group_size)
# 
# Filter for genes that have more than 5 counts in more than the smallest group size
# genes_to_keep <- rowSums((counts(dds) >= 5)) >= smallest_group_size
# dds <- dds[genes_to_keep, ]

# Alternative result filtering
# res <- results(dds, contrast=c("Outcome","Responder","Nonresponder")) 

# Alternative boxplot code
# req(input$selectedMiRNA, longData())
# 
# selected_miRNA <- input$selectedMiRNA
# ld <- longData()
# print(ld$X)
# 
# Filter the data for the selected miRNA
# ld_filtered <- ld %>% filter(X == selected_miRNA)
# print(ld_filtered)
# 
# Ensure that counts are greater than 0
# ld_filtered <- ld_filtered %>% filter(Count > 0)
# ld_filtered <- ld_filtered %>%
#   mutate(Count = ifelse(Count == 0, 1, Count))
# 
# Check if filtered data is available
# if (nrow(ld_filtered) == 0) {
#   return(NULL)  # Do not plot if there's no data
# }
# 
# Calculate log(Count)
# ld_filtered <- ld_filtered %>%
#   mutate(log_Count = log(Count))
# 
# ggplot(ld_filtered, aes(x = !!sym(grouping_variable), y = log_Count)) +
#   geom_boxplot(fill = "lightblue") +
#   geom_jitter(color = "darkblue", width = 0.2, alpha = 0.5) +
#   facet_wrap(as.formula(paste("~", grouping_variable)), scales = "free_x") +
#   labs(title = paste("Boxplot of", selected_miRNA),
#        x = "Sample",
#        y = "log(Count)") +
#   theme_cowplot(12) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Alternative volcano plot comments
# ylim(0, 15) +

# Parallel processing example
# Define the function to process each miRNA
# process_mir <- function(mir) {
#   Get predicted targets for the current miRNA
#   predictions <- getPredictedTargets(mir, species = 'hsa', method = 'geom', min_src = 2)
#   
#   Check if predictions are empty
#   if (is.null(predictions) || nrow(predictions) == 0) {
#     Remove the "hsa-" prefix and try again
#     mir_no_prefix <- sub("^hsa-", "", mir)
#     message(paste("Trying again without prefix for miRNA:", mir_no_prefix))
#     
#     predictions <- getPredictedTargets(mir_no_prefix, species = 'hsa', method = 'geom', min_src = 2)
#     
#     if (is.null(predictions) || nrow(predictions) == 0) {
#       message(paste("No targets found for miRNA (after removing prefix):", mir_no_prefix))
#       return(NULL)  # Return NULL if no predictions found
#     }
#   }
#   
#   Convert predictions to a data frame
#   df <- as.data.frame(predictions)
#   
#   Grab list of predicted genes 
#   df$entrezgene_id <- row.names(df)
#   mapping <- getBM(
#     attributes = c('ensembl_gene_id', 'description', 'entrezgene_id', 'hgnc_symbol'), 
#     filters = 'entrezgene_id',
#     values = df$entrezgene_id,
#     mart = hsmart
#   )
#   
#   df <- merge(df, mapping, by = "entrezgene_id")
#   df$miRNA <- mir
#   
#   return(df[, c("entrezgene_id", "rank_final", "ensembl_gene_id", "description", org_string(), "miRNA")])
# }
# 
# Set up a cluster
# numCores <- detectCores() - 1  # Use one less than the total number of cores
# cl <- makeCluster(numCores)
# 
# Export necessary variables and functions to the cluster
# clusterExport(cl, varlist = c("getPredictedTargets", "getBM", "hsmart"))
# 
# Parallel processing
# geneResults <- parLapply(cl, mirVector, process_mir)
# 
# Stop the cluster
# stopCluster(cl)
# 
# Remove NULL results (if any)
# geneResults <- geneResults[!sapply(geneResults, is.null)]
# 
# Combine results into a single data frame
# finalgeneResults <- do.call(rbind, geneResults)

# Store full list of predicted genes 
# allPredictedGenes <- unique(finalgeneResults$entrezgene_id)
# allPredictedGenes_reactive(allPredictedGenes)

# Alternative library loading
# library(parallel)
# library(stringr)
# library(shinycssloaders)

# Alternative testing code
# Test chord diagram
# 
# mat = matrix(1:9, 3)
# rownames(mat) = c("super extra duper long string name","pwy 2", "pwy 3")
# colnames(mat) = LETTERS[1:3]
# mat
# 
# df = data.frame(from = rep(rownames(mat), times = ncol(mat)),
#                 to = rep(colnames(mat), each = nrow(mat)),
#                 value = as.vector(mat),
#                 stringsAsFactors = FALSE)
# df$from <- stringr::str_wrap(df$from, width=20)
# df$from
# 
# chordDiagram(df)

# Modal dialog examples
# Network Plot Modal
# observeEvent(input$showNetworkPlot, {
#   showModal(modalDialog(
#     title = "Network Plot",
#     plotOutput("networkPlotModal"),
#     easyClose = TRUE,
#     footer = NULL,
#     size='xl'
#   ))
# })
# 
# Chord Plot Modal
# observeEvent(input$showChordPlot, {
#   showModal(modalDialog(
#     title = "Chord Plot",
#     plotOutput("chordPlotModal"),
#     easyClose = TRUE,
#     footer = NULL,
#     size='xl'
#   ))
# })

# UI theme alternative
# useBusyIndicators(spinners=FALSE,pulse=FALSE,fade=TRUE),

# Error checking examples
# Check if neg_cor_reactive or predicted_target_reactive is NULL
# if (nrow(as.data.frame(neg_cor_reactive())==0) || as.data.frame(nrow(predicted_target_reactive()==0))) {
#   output$degMessage <- renderText("Error: One or both of the required data sources are not available.")
#   
#   output$SupportedNegativeCorrelationTable <- renderTable({
#     return()  # Return nothing if there's an error
#   })
#   return()  # Exit the observeEvent early
# }

# Alternative plot titles
# title_text <- paste(levels_group[2], "vs", levels_group[1])
# subtitle_text <- paste("miRNA ~", grouping_variable)

# Alternative plot outputs
# output$chordPlotModal <- renderPlot({
# output$networkPlotModal <- renderPlot({

# Alternative positioning
# absolutePanel(draggable=TRUE,