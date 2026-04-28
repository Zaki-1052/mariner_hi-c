# qc_plots.R - Comprehensive QC visualization

library(ggplot2)
library(reshape2)
library(corrplot)
counts <- readRDS("outputs/full/06_counts_matrix.rds")
# --- Distribution Plots ---
# Convert to long format for ggplot
counts_long <- melt(counts)
colnames(counts_long) <- c("Loop", "Sample", "Count")

# Density plots of count distributions
p1 <- ggplot(counts_long, aes(x = Count + 1, fill = Sample)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Count Distribution by Sample",
       x = "Hi-C Contact Count + 1 (log10)",
       y = "Density") +
  theme(legend.position = "top")

# Box plots
p2 <- ggplot(counts_long, aes(x = Sample, y = Count + 1, fill = Sample)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Count Distribution Comparison",
       y = "Hi-C Contact Count + 1 (log10)")

# Combine plots
library(patchwork)
p1 + p2

# --- MA Plot ---
# Calculate M (log fold change) and A (average expression)
A <- 0.5 * (log2(counts[,1] + 1) + log2(counts[,2] + 1))
M <- log2(counts[,2] + 1) - log2(counts[,1] + 1)

ma_data <- data.frame(A = A, M = M)

ggplot(ma_data, aes(x = A, y = M)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), color = "blue", linetype = "dotted") +
  theme_minimal() +
  labs(title = "MA Plot: Mutant vs Control",
       x = "Average log2(Count + 1)",
       y = "log2 Fold Change (mut/ctrl)")

