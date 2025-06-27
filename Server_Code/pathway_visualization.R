# Pathway Visualization - Chord Plot and Network Plot

observeEvent(input$showChordPlot,{
  req(miRNA_mapped_pathways())
  
  chord_data <- miRNA_mapped_pathways() %>% 
    dplyr::select(c("miRNA", "description", "enrichment", "coverage")) %>%
    unique()
  
  top_pathways <- All_miRNA_Pathways_reactive() %>% 
    as.data.frame() %>% 
    arrange(qvalue) %>% 
    slice_head(n = input$num_pathways) %>%
    pull(Description)
  
  filtered_data <- chord_data %>%
    filter(description %in% top_pathways) %>% 
    filter(coverage >= input$min_coverage)
  
  filtered_data$description <- stringr::str_wrap(filtered_data$description, width=30)
  filtered_data$enrichment <- as.numeric(filtered_data$enrichment)
  
  mat <- as.matrix(xtabs(enrichment ~ miRNA + description, data = filtered_data))
  
  chord_plot <- function() {
    circos.clear()
    chordDiagram(mat, 
                 transparency = 0.5, 
                 grid.col = NULL,
                 annotationTrack = "grid")
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, 
                  CELL_META$ylim[1], 
                  CELL_META$sector.index,
                  facing = "clockwise", 
                  niceFacing = TRUE, 
                  adj = c(0.8),
                  cex = 1)
    }, bg.border = NA)
  }
  
  output$chordPlot <- renderPlot({
    chord_plot()
  })
  
  chord_plot_reactive(chord_plot)
})

observeEvent(input$showNetworkPlot, {
  req(miRNA_mapped_pathways())
  chord_data <- miRNA_mapped_pathways() %>% 
    dplyr::select(c("miRNA", "description", "enrichment", "coverage")) %>%
    unique()
  
  top_pathways <- All_miRNA_Pathways_reactive() %>% 
    as.data.frame() %>% 
    arrange(qvalue) %>% 
    slice_head(n = input$num_pathways) %>%
    pull(Description)
  
  filtered_data <- chord_data %>%
    filter(description %in% top_pathways) %>% 
    filter(coverage >= input$min_coverage)
  
  edges <- filtered_data[, c("miRNA", "description")]
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  node_colors <- data.frame(name = V(g)$name,
                            type = ifelse(V(g)$name %in% filtered_data$miRNA, "miRNA", "description"))
  
  V(g)$color <- ifelse(node_colors$type == "miRNA", "lightblue", "purple")
  
  set.seed(123)
  plot <- ggraph(g, layout = "fr") + 
    geom_edge_link(color = "grey") +
    geom_node_point(aes(color = color), size = 5) +
    geom_node_text(aes(label = name), repel = TRUE) +
    scale_color_identity() +
    theme_graph(fg_text_colour = 'white')
  
  output$networkPlot <- renderPlot({
    plot
  })
  
  network_plot_reactive(plot)
})