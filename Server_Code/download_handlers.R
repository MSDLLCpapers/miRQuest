# All download handlers must be defined at the top level, not inside observeEvents

# Module 1: Exploratory Visualization Downloads
output$downloadStackedColumn <- downloadHandler(
  filename = function() {
    paste("Stacked_Column_Chart_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- stacked_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$downloadPlot <- downloadHandler(
  filename = function() {
    paste("MDS_Plot_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- mds_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# Module 3: Differential miRNA Expression Downloads  
output$downloadVolcano <- downloadHandler(
  filename = function() {
    paste("Volcano_Plot_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- volcano_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$downloadSignificantHeatmap <- downloadHandler(
  filename = function() {
    paste("Significant_Heatmap_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- significant_heatmap_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# Module 5: Single miRNA Target Prediction Downloads
output$downloadDotplot <- downloadHandler(
  filename = function() {
    paste("Pathway_Dotplot_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- dotplot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$downloadBarplot <- downloadHandler(
  filename = function() {
    paste("Pathway_Barplot_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- barplot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# Module 4: Aggregated Pathway Downloads
output$downloadChord <- downloadHandler(
  filename = function() {
    paste("Chord_Plot_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- chord_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

# output$downloadNetwork <- downloadHandler(
#   filename = function() {
#     paste("Network_Plot_", Sys.Date(), ".png", sep="")
#   },
#   content = function(file) {
#     png(file, width = 800, height = 600)
#     plot_func <- network_plot_reactive()
#     if (is.function(plot_func)) {
#       print(plot_func())
#     } else {
#       print(plot_func)
#     }
#     dev.off()
#   }
# )

# Module 7: Correlation Analysis Downloads
output$downloadCorrelationPlot <- downloadHandler(
  filename = function() {
    paste("Correlation_Plot_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- correlation_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)

output$corr_downloadChord <- downloadHandler(
  filename = function() {
    paste("Chord_Plot_", Sys.Date(), ".png", sep="")
  },
  content = function(file) {
    png(file, width = 800, height = 600)
    plot_func <- chord_plot_reactive()
    if (is.function(plot_func)) {
      print(plot_func())
    } else {
      print(plot_func)
    }
    dev.off()
  }
)