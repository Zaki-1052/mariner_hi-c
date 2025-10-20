
library(tidyverse)
library(DT)
library(plotly)

counts <- readRDS("outputs/full/06_counts_matrix.rds")
y <- readRDS("outputs/full/06_edgeR_input.rds")
metadata <- readRDS("outputs/full/05_metadata.rds")
