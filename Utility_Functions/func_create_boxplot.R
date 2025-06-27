# Define a function to create the boxplot
create_boxplot <- function(dds, selected_miRNA, grouping_variable) {
  df <- plotCounts(dds = dds, gene = selected_miRNA, intgroup = grouping_variable, returnData = TRUE)
  
  ggplot(df, aes(x = !!sym(grouping_variable), y = count, fill = !!sym(grouping_variable))) + 
    geom_boxplot(fill = "lightblue") +
    geom_jitter(color = "darkblue", width = 0.2, alpha = 0.5) +
    facet_wrap(as.formula(paste("~", grouping_variable)), scales = "free_x") +
    labs(title = paste("Boxplot of", selected_miRNA),
         x = "Sample",
         y = "Normalized Counts") +
    theme_cowplot(12)
}