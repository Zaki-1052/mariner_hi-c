# /path/to/local/mariner/analysis_explore.R

library(tidyverse)
library(edgeR)
library(pheatmap)

# Load the count matrix
counts <- readRDS("outputs/full/06_counts_matrix.rds")
dim(counts)  # Should be 150 x 2

# Basic inspection
head(counts)
summary(counts)

# Load the edgeR object for full context
y <- readRDS("outputs/full/06_edgeR_input.rds")

# View loop coordinates
head(y$genes)

