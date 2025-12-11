# Module 2: Filter count table by metadata column value
output$filterValueInput <- renderUI({
  req(input$metadataColumn)
  metadataData <- metadata()
  
  unique_values <- unique(metadataData[[input$metadataColumn]])
  
  selectInput("filterValue", "Select Filter Value", choices = unique_values)
})

observeEvent(input$subsetDataButton, {
  req(wideData(), metadata(), input$metadataColumn, input$filterValue)
  
  wd <- wideData()
  metadataData <- metadata()
  
  selected_column <- input$metadataColumn
  filter_value <- input$filterValue
  
  samples <- metadataData %>%
    filter(!!sym(selected_column) == filter_value) %>%
    pull(SampleID)
  
  subsets_data <- wd[,colnames(wd) %in% samples]
  
  row.names(subsets_data)<-row.names(wd)
  print(row.names(subsets_data))
  
  # Store in reactive for download handler
  subset_data_reactive(subsets_data)
  
  output$subsetDataTable <- renderTable({
    head(subsets_data)
  })
})

# Download handler (outside observeEvent)
output$downloadSubset <- downloadHandler(
  filename = function() {
    paste("subset_data_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    req(subset_data_reactive())
    write.csv(subset_data_reactive(), file, row.names = TRUE)
  }
)