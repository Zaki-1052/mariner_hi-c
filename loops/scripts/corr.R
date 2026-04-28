# correlation_analysis.R

library(pheatmap)
library(viridis)
counts <- readRDS("outputs/full/06_counts_matrix.rds")
# Sample correlation
cor_matrix <- cor(counts, use = "complete.obs")

# Simple correlation plot
pheatmap(cor_matrix, 
         display_numbers = TRUE,
         color = viridis(100),
         main = "Sample Correlation")

# --- Loop-level heatmap (top variable loops) ---
# Calculate coefficient of variation for each loop
cv <- apply(counts, 1, function(x) sd(x) / mean(x))
top_variable <- names(sort(cv, decreasing = TRUE)[1:50])

# Heatmap of top variable loops
counts_scaled <- t(scale(t(counts[top_variable,])))
pheatmap(counts_scaled,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 50 Variable Loops",
         show_rownames = FALSE,
         labels_col = c("Control", "Mutant"))
