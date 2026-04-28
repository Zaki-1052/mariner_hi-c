counts <- readRDS("outputs/full/06_counts_matrix.rds")
print(dim(counts))
print(summary(as.vector(counts)))
print(paste("Correlation:", round(cor(counts[,1], counts[,2]), 3)))

# Quick MA plot
plot(rowMeans(log2(counts + 1)), 
     log2(counts[,2] + 1) - log2(counts[,1] + 1),
     xlab = "Average Expression", ylab = "Log2 FC",
     main = "MA Plot", pch = 16, col = rgb(0,0,0,0.5))
abline(h = 0, col = "red")
