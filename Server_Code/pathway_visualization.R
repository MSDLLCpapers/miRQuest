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

# observeEvent(input$showNetworkPlot, {
#   req(miRNA_mapped_pathways())
#   chord_data <- miRNA_mapped_pathways() %>% 
#     dplyr::select(c("miRNA", "description", "enrichment", "coverage")) %>%
#     unique()
#   
#   top_pathways <- All_miRNA_Pathways_reactive() %>% 
#     as.data.frame() %>% 
#     arrange(qvalue) %>% 
#     slice_head(n = input$num_pathways) %>%
#     pull(Description)
#   
#   filtered_data <- chord_data %>%
#     filter(description %in% top_pathways) %>% 
#     filter(coverage >= input$min_coverage)
#   
#   edges <- filtered_data[, c("miRNA", "description")]
#   
#   g <- graph_from_data_frame(edges, directed = FALSE)
#   
#   node_colors <- data.frame(name = V(g)$name,
#                             type = ifelse(V(g)$name %in% filtered_data$miRNA, "miRNA", "description"))
#   
#   V(g)$color <- ifelse(node_colors$type == "miRNA", "lightblue", "purple")
#   
#   set.seed(123)
#   plot <- ggraph(g, layout = "fr") + 
#     geom_edge_link(color = "grey") +
#     geom_node_point(aes(color = color), size = 5) +
#     geom_node_text(aes(label = name), repel = TRUE) +
#     scale_color_identity() +
#     theme_graph(fg_text_colour = 'white')
#   
#   output$networkPlot <- renderPlot({
#     plot
#   })
#   
#   network_plot_reactive(plot)
# })

observeEvent(input$showNetworkPlot, {
  req(miRNA_mapped_pathways(), All_miRNA_Pathways_reactive())
  
  # Prepare data
  chord_data <- miRNA_mapped_pathways() %>%
    dplyr::select(miRNA, description, enrichment, coverage) %>%
    dplyr::distinct()
  
  top_pathways <- All_miRNA_Pathways_reactive() %>%
    as.data.frame() %>%
    dplyr::arrange(qvalue) %>%
    dplyr::slice_head(n = input$num_pathways) %>%
    dplyr::pull(Description)
  
  filtered_data <- chord_data %>%
    dplyr::filter(description %in% top_pathways) %>%
    dplyr::filter(coverage >= input$min_coverage)
  
  # Build unique edges
  edges <- filtered_data %>%
    dplyr::transmute(
      from = miRNA,
      to   = description,
      title = paste0(
        "<b>miRNA:</b> ", miRNA,
        "<br><b>Pathway:</b> ", description,
        ifelse(!is.na(enrichment), paste0("<br><b>Enrichment:</b> ", enrichment), ""),
        ifelse(!is.na(coverage),   paste0("<br><b>Coverage:</b> ", coverage), "")
      )
    ) %>%
    dplyr::distinct(from, to, .keep_all = TRUE)
  
  # Nodes: miRNA and pathways
  miRNA_nodes <- edges %>%
    dplyr::distinct(from) %>%
    dplyr::transmute(
      id = from, label = from,
      group = "miRNA", color = "#6EC5FF"
    )
  
  pathway_nodes <- edges %>%
    dplyr::distinct(to) %>%
    dplyr::transmute(
      id = to, label = to,
      group = "pathway", color = "#B184F0"
    )
  
  nodes <- dplyr::bind_rows(miRNA_nodes, pathway_nodes)
  
  # Compute degree for tooltips and optional pruning
  deg <- dplyr::bind_rows(
    edges %>% dplyr::count(from, name = "deg") %>% dplyr::rename(id = from),
    edges %>% dplyr::count(to,   name = "deg") %>% dplyr::rename(id = to)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(degree = sum(deg), .groups = "drop")
  
  nodes <- nodes %>%
    dplyr::left_join(deg, by = "id") %>%
    dplyr::mutate(degree = ifelse(is.na(degree), 0L, degree)) %>%
    dplyr::mutate(
      title = ifelse(group == "miRNA",
                     paste0("<b>miRNA:</b> ", label, "<br><b>Degree:</b> ", degree),
                     paste0("<b>Pathway:</b> ", label, "<br><b>Degree:</b> ", degree))
    )
  
  # Optional: prune small nodes to keep it readable
  if (!is.null(input$path_min_degree) && input$path_min_degree > 0) {
    keep_ids <- nodes %>% dplyr::filter(degree >= input$path_min_degree) %>% dplyr::pull(id)
    edges <- edges %>% dplyr::filter(from %in% keep_ids, to %in% keep_ids)
    nodes <- nodes %>% dplyr::filter(id %in% keep_ids)
  }
  
  # Recompute degree after pruning (for accurate tooltips)
  final_deg <- dplyr::bind_rows(
    edges %>% dplyr::count(from, name = "deg") %>% dplyr::rename(id = from),
    edges %>% dplyr::count(to,   name = "deg") %>% dplyr::rename(id = to)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(degree = sum(deg), .groups = "drop")
  
  nodes <- nodes %>%
    dplyr::select(-degree) %>%
    dplyr::left_join(final_deg, by = "id") %>%
    dplyr::mutate(degree = ifelse(is.na(degree), 0L, degree)) %>%
    dplyr::mutate(
      title = ifelse(group == "miRNA",
                     paste0("<b>miRNA:</b> ", label, "<br><b>Degree (shown):</b> ", degree),
                     paste0("<b>Pathway:</b> ", label, "<br><b>Degree (shown):</b> ", degree))
    )
  
  # Render interactive network
  output$networkPlot <- visNetwork::renderVisNetwork({
    visNetwork::visNetwork(nodes, edges, height = "700px", width = "100%") %>%
      visNetwork::visNodes(shape = "dot", size = 18, font = list(size = 16)) %>%
      visNetwork::visEdges(smooth = FALSE, color = list(color = "#B3B3B3"), arrows = "none") %>%
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = TRUE
      ) %>%
      visNetwork::visLegend(
        addNodes = list(
          list(label = "miRNA",   shape = "dot", color = "#6EC5FF"),
          list(label = "Pathway", shape = "dot", color = "#B184F0")
        ),
        useGroups = FALSE
      ) %>%
      visNetwork::visInteraction(hover = TRUE, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
      visNetwork::visIgraphLayout(randomSeed = 123)  # FR-like layout
  })
  
  # to store the widget for reuse
  #network_plot_reactive(NULL)  # or store nodes/edges if you need them elsewhere
})

