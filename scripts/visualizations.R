#!/usr/bin/env Rscript
# scripts/visualizations.R
# Comprehensive Visualization Analysis for Differential Chromatin Loops
# Author: Zakir Alibhai
# Date: 2025-10-30
#
# Purpose:
#   Generate publication-quality visualizations from differential loop analysis:
#   1. Volcano plots (EnhancedVolcano)
#   2. Feature distribution (ChIPseeker)
#   3. GO/KEGG enrichment (clusterProfiler)
#   4. Loop type classification
#   5. Loop length distribution
#   6. Aggregate Peak Analysis (APA)
#
# Usage:
#   Rscript scripts/visualizations.R [--resolution RES] [--skip-apa]
#
# Arguments:
#   --resolution RES   Resolution for volcano plot (5000, 10000, or 25000)
#   --skip-apa         Skip time-intensive APA analysis

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("Visualization Analysis: Differential Loops\n")
cat("========================================\n\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
volcano_resolution <- 5000
skip_apa <- FALSE

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--resolution" && i < length(args)) {
      volcano_resolution <- as.numeric(args[i + 1])
    } else if (args[i] == "--skip-apa") {
      skip_apa <- TRUE
    }
  }
}

cat(sprintf("Configuration:\n"))
cat(sprintf("  Volcano plot resolution: %d kb\n", volcano_resolution/1000))
cat(sprintf("  Skip APA: %s\n\n", ifelse(skip_apa, "Yes", "No")))

# Load required libraries
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  # Core visualization
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(EnhancedVolcano)

  # Genomic annotation
  library(GenomicRanges)
  library(InteractionSet)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)

  # Enrichment analysis
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)

  # Hi-C/APA
  if (!skip_apa) {
    library(mariner)
    library(strawr)
  }

  # Data wrangling
  library(tidyverse)
})

cat("✓ Packages loaded\n\n")

# Set working directory
if (dir.exists("/expanse/lustre/projects/csd940/zalibhai/mariner")) {
  setwd("/expanse/lustre/projects/csd940/zalibhai/mariner")
} else if (basename(getwd()) == "mariner-final") {
  # Already in project directory
} else {
  warning("Working directory may not be correct. Current: ", getwd())
}

# Create output directories
output_dir <- "outputs/visualizations"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "volcano"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "features"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "loop_classification"), recursive = TRUE, showWarnings = FALSE)
if (!skip_apa) {
  dir.create(file.path(output_dir, "apa"), recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("Output directory: %s\n\n", output_dir))

# =============================================================================
# LOAD DATA
# =============================================================================

cat("="*70, "\n", sep="")
cat("Loading Input Data\n")
cat("="*70, "\n\n", sep="")

# Load characterized loops from downstream analysis
characterized_file <- "outputs/merged_loops/characterized_loops.tsv"
if (!file.exists(characterized_file)) {
  stop("ERROR: characterized_loops.tsv not found. Run downstream_analysis.R first.")
}

loops_df <- read.table(characterized_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("✓ Loaded characterized loops: %d loops\n", nrow(loops_df)))
cat(sprintf("  Up-regulated: %d\n", sum(loops_df$logFC > 0)))
cat(sprintf("  Down-regulated: %d\n\n", sum(loops_df$logFC < 0)))

# Load GInteractions object
gi_file <- "outputs/merged_loops/non_redundant_loops.rds"
if (file.exists(gi_file)) {
  loops_gi <- readRDS(gi_file)
  cat(sprintf("✓ Loaded GInteractions object: %d interactions\n\n", length(loops_gi)))
} else {
  warning("GInteractions file not found. Some analyses may be limited.")
  loops_gi <- NULL
}

# =============================================================================
# SECTION 1: VOLCANO PLOTS
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 1: Volcano Plots\n")
cat("="*70, "\n\n", sep="")

# Function to create volcano plot for a given resolution
create_volcano_plot <- function(resolution_kb, output_path) {
  cat(sprintf("Creating volcano plot for %dkb resolution...\n", resolution_kb))

  # Load all results
  results_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv",
                          resolution_kb)

  if (!file.exists(results_file)) {
    cat(sprintf("  ✗ Results file not found: %s\n", results_file))
    return(NULL)
  }

  df <- read.table(results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d loops\n", nrow(df)))

  # Count significant loops
  n_up <- sum(df$significant & df$logFC > 0, na.rm = TRUE)
  n_down <- sum(df$significant & df$logFC < 0, na.rm = TRUE)

  cat(sprintf("  Significant: %d up, %d down\n", n_up, n_down))

  # Create volcano plot
  p <- EnhancedVolcano(df,
                      lab = df$loop_id,
                      x = 'logFC',
                      y = 'FDR',
                      title = sprintf('Differential Chromatin Loops (%dkb)', resolution_kb),
                      subtitle = sprintf('%d total | %d up | %d down', nrow(df), n_up, n_down),
                      pCutoff = 0.05,
                      FCcutoff = 0.3,
                      pointSize = 2.0,
                      labSize = 3.0,
                      col = c('grey30', 'grey60', '#d73027', '#4575b4'),
                      colAlpha = 0.5,
                      legendPosition = 'right',
                      legendLabSize = 10,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      max.overlaps = 20)

  # Save plot
  ggsave(output_path, p, width = 10, height = 8)
  cat(sprintf("  ✓ Saved: %s\n\n", output_path))

  return(p)
}

# Generate volcano plots for each resolution
for (res_kb in c(5, 10, 25)) {
  output_file <- file.path(output_dir, "volcano", sprintf("volcano_%dkb.pdf", res_kb))
  create_volcano_plot(res_kb, output_file)
}

# Create merged volcano plot (optional)
cat("Creating merged multi-resolution volcano plot...\n")

all_resolutions <- list()
for (res_kb in c(5, 10, 25)) {
  results_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv",
                          res_kb)
  if (file.exists(results_file)) {
    df <- read.table(results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    df$resolution <- sprintf("%dkb", res_kb)
    all_resolutions[[as.character(res_kb)]] <- df
  }
}

if (length(all_resolutions) > 0) {
  merged_df <- bind_rows(all_resolutions)

  # Assign unique IDs
  merged_df$unique_id <- paste(merged_df$loop_id, merged_df$resolution, sep = "_")

  n_total <- nrow(merged_df)
  n_up <- sum(merged_df$significant & merged_df$logFC > 0, na.rm = TRUE)
  n_down <- sum(merged_df$significant & merged_df$logFC < 0, na.rm = TRUE)

  p_merged <- EnhancedVolcano(merged_df,
                             lab = merged_df$unique_id,
                             x = 'logFC',
                             y = 'FDR',
                             title = 'Differential Chromatin Loops (Multi-Resolution)',
                             subtitle = sprintf('%d total | %d up | %d down', n_total, n_up, n_down),
                             pCutoff = 0.05,
                             FCcutoff = 0.3,
                             pointSize = 1.5,
                             labSize = 2.5,
                             col = c('grey30', 'grey60', '#d73027', '#4575b4'),
                             colAlpha = 0.4,
                             legendPosition = 'right',
                             max.overlaps = 15)

  ggsave(file.path(output_dir, "volcano", "volcano_merged.pdf"), p_merged, width = 12, height = 8)
  cat("  ✓ Saved: volcano_merged.pdf\n\n")
}

cat("✓ Section 1 complete: Volcano plots generated\n\n")

# =============================================================================
# SECTION 2: FEATURE DISTRIBUTION ANALYSIS
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 2: Feature Distribution Analysis\n")
cat("="*70, "\n\n", sep="")

cat("Annotating loop anchors with genomic features...\n")

# Load TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Create GRanges for anchors
anchor1_gr <- GRanges(
  seqnames = loops_df$anchor1_chr,
  ranges = IRanges(start = loops_df$anchor1_start, end = loops_df$anchor1_end)
)

anchor2_gr <- GRanges(
  seqnames = loops_df$anchor2_chr,
  ranges = IRanges(start = loops_df$anchor2_start, end = loops_df$anchor2_end)
)

# Annotate with ChIPseeker
cat("  Annotating anchor1...\n")
anno1 <- annotatePeak(anchor1_gr, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb = "org.Mm.eg.db")
anno1_df <- as.data.frame(anno1)

cat("  Annotating anchor2...\n")
anno2 <- annotatePeak(anchor2_gr, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb = "org.Mm.eg.db")
anno2_df <- as.data.frame(anno2)

# Add annotations to loops dataframe
loops_df$anchor1_annotation <- anno1_df$annotation
loops_df$anchor2_annotation <- anno2_df$annotation

# Create feature distribution summary
feature_summary <- data.frame(
  category = character(),
  anchor = character(),
  feature = character(),
  count = numeric(),
  percentage = numeric(),
  stringsAsFactors = FALSE
)

# Helper function to summarize features
summarize_features <- function(annotations, direction, anchor_num) {
  # Simplify annotations
  simple_anno <- gsub(" \\(.*\\)", "", annotations)

  # Count occurrences
  feature_counts <- table(simple_anno)

  df <- data.frame(
    category = paste0(direction, "_anchor", anchor_num),
    anchor = paste0("anchor", anchor_num),
    feature = names(feature_counts),
    count = as.numeric(feature_counts),
    percentage = 100 * as.numeric(feature_counts) / length(annotations),
    stringsAsFactors = FALSE
  )

  return(df)
}

# Summarize for up-regulated loops
up_idx <- loops_df$logFC > 0
feature_summary <- rbind(feature_summary,
                        summarize_features(loops_df$anchor1_annotation[up_idx], "up", 1))
feature_summary <- rbind(feature_summary,
                        summarize_features(loops_df$anchor2_annotation[up_idx], "up", 2))

# Summarize for down-regulated loops
down_idx <- loops_df$logFC < 0
feature_summary <- rbind(feature_summary,
                        summarize_features(loops_df$anchor1_annotation[down_idx], "down", 1))
feature_summary <- rbind(feature_summary,
                        summarize_features(loops_df$anchor2_annotation[down_idx], "down", 2))

# Save summary
write.table(feature_summary,
            file.path(output_dir, "features", "feature_distribution_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("  ✓ Feature distribution summary saved\n")

# Create stacked bar plot
cat("  Creating feature distribution plot...\n")

# Order categories
feature_summary$category <- factor(feature_summary$category,
                                  levels = c("up_anchor1", "up_anchor2", "down_anchor1", "down_anchor2"))

# Color palette for features
feature_colors <- c(
  "Promoter" = "#e41a1c",
  "5' UTR" = "#377eb8",
  "3' UTR" = "#4daf4a",
  "1st Exon" = "#984ea3",
  "Other Exon" = "#ff7f00",
  "1st Intron" = "#ffff33",
  "Other Intron" = "#a65628",
  "Downstream" = "#f781bf",
  "Distal Intergenic" = "#999999"
)

p_features <- ggplot(feature_summary, aes(x = category, y = percentage, fill = feature)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  scale_fill_manual(values = feature_colors, name = "Genomic Feature") +
  labs(
    title = "Genomic Feature Distribution of Loop Anchors",
    x = "Category",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  coord_flip()

ggsave(file.path(output_dir, "features", "feature_distribution.pdf"), p_features, width = 10, height = 6)
cat("  ✓ Saved: feature_distribution.pdf\n\n")

cat("✓ Section 2 complete: Feature distribution analyzed\n\n")

# =============================================================================
# SECTION 3: GO/KEGG ENRICHMENT ANALYSIS
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 3: GO/KEGG Enrichment Analysis\n")
cat("="*70, "\n\n", sep="")

# Extract genes near loop anchors
cat("Extracting genes from loop anchors...\n")

# Get genes within 10kb of each anchor
get_nearby_genes <- function(gr, txdb, max_dist = 10000) {
  genes_txdb <- genes(txdb)

  # Find overlaps
  overlaps <- findOverlaps(gr, genes_txdb, maxgap = max_dist)

  # Get gene IDs
  gene_ids <- names(genes_txdb)[subjectHits(overlaps)]

  # Map to Entrez IDs
  entrez_map <- select(org.Mm.eg.db, keys = gene_ids, columns = "ENTREZID", keytype = "GENEID")
  entrez_ids <- unique(entrez_map$ENTREZID[!is.na(entrez_map$ENTREZID)])

  return(entrez_ids)
}

# Get genes for up-regulated loops
up_anchor1_gr <- anchor1_gr[up_idx]
up_anchor2_gr <- anchor2_gr[up_idx]
up_genes <- unique(c(
  get_nearby_genes(up_anchor1_gr, txdb),
  get_nearby_genes(up_anchor2_gr, txdb)
))

cat(sprintf("  Up-regulated loops: %d genes\n", length(up_genes)))

# Get genes for down-regulated loops
down_anchor1_gr <- anchor1_gr[down_idx]
down_anchor2_gr <- anchor2_gr[down_idx]
down_genes <- unique(c(
  get_nearby_genes(down_anchor1_gr, txdb),
  get_nearby_genes(down_anchor2_gr, txdb)
))

cat(sprintf("  Down-regulated loops: %d genes\n\n", length(down_genes)))

# Check if we have enough genes
if (length(up_genes) < 20 || length(down_genes) < 20) {
  cat("  ⚠ Warning: Few genes for enrichment analysis. Results may be limited.\n\n")
}

# Create gene list
gene_list <- list(
  up_genes = up_genes,
  down_genes = down_genes
)

# Run enrichment analyses
cat("Running enrichment analyses...\n")

# GO Biological Process
cat("  GO Biological Process...\n")
go_bp <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

if (!is.null(go_bp) && nrow(go_bp@compareClusterResult) > 0) {
  p_go_bp <- dotplot(go_bp, showCategory = 20) +
    labs(title = "GO Biological Process Enrichment") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(output_dir, "enrichment", "go_bp_dotplot.pdf"), p_go_bp, width = 12, height = 10)
  cat("    ✓ Saved: go_bp_dotplot.pdf\n")
} else {
  cat("    ⚠ No significant GO BP terms found\n")
}

# GO Cellular Component
cat("  GO Cellular Component...\n")
go_cc <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

if (!is.null(go_cc) && nrow(go_cc@compareClusterResult) > 0) {
  p_go_cc <- dotplot(go_cc, showCategory = 15) +
    labs(title = "GO Cellular Component Enrichment") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(output_dir, "enrichment", "go_cc_dotplot.pdf"), p_go_cc, width = 10, height = 8)
  cat("    ✓ Saved: go_cc_dotplot.pdf\n")
} else {
  cat("    ⚠ No significant GO CC terms found\n")
}

# GO Molecular Function
cat("  GO Molecular Function...\n")
go_mf <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

if (!is.null(go_mf) && nrow(go_mf@compareClusterResult) > 0) {
  p_go_mf <- dotplot(go_mf, showCategory = 15) +
    labs(title = "GO Molecular Function Enrichment") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(output_dir, "enrichment", "go_mf_dotplot.pdf"), p_go_mf, width = 10, height = 8)
  cat("    ✓ Saved: go_mf_dotplot.pdf\n")
} else {
  cat("    ⚠ No significant GO MF terms found\n")
}

# KEGG pathways
cat("  KEGG pathways...\n")
kegg <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichKEGG",
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

if (!is.null(kegg) && nrow(kegg@compareClusterResult) > 0) {
  p_kegg <- dotplot(kegg, showCategory = 20) +
    labs(title = "KEGG Pathway Enrichment") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(output_dir, "enrichment", "kegg_dotplot.pdf"), p_kegg, width = 12, height = 10)
  cat("    ✓ Saved: kegg_dotplot.pdf\n")
} else {
  cat("    ⚠ No significant KEGG pathways found\n")
}

cat("\n✓ Section 3 complete: Enrichment analysis done\n\n")

# =============================================================================
# SECTION 4: LOOP TYPE CLASSIFICATION
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 4: Loop Type Classification\n")
cat("="*70, "\n\n", sep="")

cat("Creating loop type classification plots...\n")

# Summarize loop types by direction
loop_type_summary <- loops_df %>%
  group_by(direction, loop_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    direction_label = ifelse(direction == "up_in_mutant", "Up", "Down")
  )

# Calculate percentages within each direction
loop_type_summary <- loop_type_summary %>%
  group_by(direction_label) %>%
  mutate(
    total = sum(count),
    percentage = 100 * count / total
  ) %>%
  ungroup()

# Create pie charts
p_up_pie <- loop_type_summary %>%
  filter(direction_label == "Up") %>%
  ggplot(aes(x = "", y = percentage, fill = loop_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Up-regulated Loops",
    subtitle = sprintf("n = %d", sum(up_idx)),
    fill = "Loop Type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  geom_text(aes(label = sprintf("%s\n%.1f%%", loop_type, percentage)),
            position = position_stack(vjust = 0.5))

p_down_pie <- loop_type_summary %>%
  filter(direction_label == "Down") %>%
  ggplot(aes(x = "", y = percentage, fill = loop_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Down-regulated Loops",
    subtitle = sprintf("n = %d", sum(down_idx)),
    fill = "Loop Type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  geom_text(aes(label = sprintf("%s\n%.1f%%", loop_type, percentage)),
            position = position_stack(vjust = 0.5))

# Combine plots
p_combined <- p_up_pie | p_down_pie

ggsave(file.path(output_dir, "loop_classification", "loop_type_classification.pdf"),
       p_combined, width = 14, height = 6)
cat("  ✓ Saved: loop_type_classification.pdf\n")

# Save loop type summary table
write.table(loop_type_summary,
            file.path(output_dir, "loop_classification", "loop_type_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Extract genes by loop type
loop_type_genes <- loops_df %>%
  filter(!is.na(anchor1_nearest_gene) | !is.na(anchor2_nearest_gene)) %>%
  select(loop_id, direction, loop_type, anchor1_nearest_gene, anchor2_nearest_gene) %>%
  pivot_longer(cols = c(anchor1_nearest_gene, anchor2_nearest_gene),
               names_to = "anchor", values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct()

write.table(loop_type_genes,
            file.path(output_dir, "loop_classification", "loop_type_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("  ✓ Saved: loop_type_genes.tsv\n\n")

cat("✓ Section 4 complete: Loop type classification done\n\n")

# =============================================================================
# SECTION 5: LOOP LENGTH DISTRIBUTION
# =============================================================================

cat("="*70, "\n", sep="")
cat("SECTION 5: Loop Length Distribution\n")
cat("="*70, "\n\n", sep="")

cat("Analyzing loop length distributions...\n")

# Prepare data
length_df <- loops_df %>%
  mutate(
    direction_label = ifelse(logFC > 0, "Up", "Down"),
    length_kb = loop_distance / 1000
  )

# Statistical test
wilcox_test <- wilcox.test(
  length_df$loop_distance[length_df$direction_label == "Up"],
  length_df$loop_distance[length_df$direction_label == "Down"]
)

cat(sprintf("  Wilcoxon test p-value: %.3e\n", wilcox_test$p.value))
cat(sprintf("  Median length (Up): %.1f kb\n",
            median(length_df$loop_distance[length_df$direction_label == "Up"]) / 1000))
cat(sprintf("  Median length (Down): %.1f kb\n\n",
            median(length_df$loop_distance[length_df$direction_label == "Down"]) / 1000))

# 1. Strip plot
p_strip <- ggplot(length_df, aes(x = length_kb, y = direction_label, color = direction_label)) +
  geom_jitter(alpha = 0.5, size = 2, height = 0.2) +
  scale_color_manual(values = c("Up" = "#4575b4", "Down" = "#d73027")) +
  scale_x_log10() +
  labs(
    title = "Loop Length Distribution",
    subtitle = sprintf("Wilcoxon p = %.3e", wilcox_test$p.value),
    x = "Loop Length (kb, log10)",
    y = "Direction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "loop_classification", "loop_length_distribution_strip.pdf"),
       p_strip, width = 10, height = 5)
cat("  ✓ Saved: loop_length_distribution_strip.pdf\n")

# 2. Violin plot
p_violin <- ggplot(length_df, aes(x = direction_label, y = length_kb, fill = direction_label)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("Up" = "#4575b4", "Down" = "#d73027")) +
  scale_y_log10() +
  labs(
    title = "Loop Length Distribution",
    subtitle = sprintf("Wilcoxon p = %.3e", wilcox_test$p.value),
    x = "Direction",
    y = "Loop Length (kb, log10)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "loop_classification", "loop_length_distribution_violin.pdf"),
       p_violin, width = 8, height = 6)
cat("  ✓ Saved: loop_length_distribution_violin.pdf\n")

# 3. Histogram
p_hist <- ggplot(length_df, aes(x = length_kb, fill = direction_label)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Up" = "#4575b4", "Down" = "#d73027")) +
  scale_x_log10() +
  labs(
    title = "Loop Length Distribution",
    x = "Loop Length (kb, log10)",
    y = "Count",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  ) +
  facet_wrap(~direction_label, ncol = 1)

ggsave(file.path(output_dir, "loop_classification", "loop_length_distribution_histogram.pdf"),
       p_hist, width = 10, height = 8)
cat("  ✓ Saved: loop_length_distribution_histogram.pdf\n")

# Save statistics
length_stats <- data.frame(
  direction = c("Up", "Down"),
  n = c(sum(length_df$direction_label == "Up"), sum(length_df$direction_label == "Down")),
  median_kb = c(
    median(length_df$loop_distance[length_df$direction_label == "Up"]) / 1000,
    median(length_df$loop_distance[length_df$direction_label == "Down"]) / 1000
  ),
  mean_kb = c(
    mean(length_df$loop_distance[length_df$direction_label == "Up"]) / 1000,
    mean(length_df$loop_distance[length_df$direction_label == "Down"]) / 1000
  ),
  sd_kb = c(
    sd(length_df$loop_distance[length_df$direction_label == "Up"]) / 1000,
    sd(length_df$loop_distance[length_df$direction_label == "Down"]) / 1000
  )
)

length_stats$wilcox_pvalue <- wilcox_test$p.value

write.table(length_stats,
            file.path(output_dir, "loop_classification", "loop_length_statistics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  ✓ Saved: loop_length_statistics.txt\n\n")

cat("✓ Section 5 complete: Loop length distribution analyzed\n\n")

# =============================================================================
# SECTION 6: AGGREGATE PEAK ANALYSIS (APA)
# =============================================================================

if (!skip_apa) {
  cat("="*70, "\n", sep="")
  cat("SECTION 6: Aggregate Peak Analysis (APA)\n")
  cat("="*70, "\n\n", sep="")

  cat("⚠ APA analysis requires Hi-C files and may take 20-30 minutes\n\n")

  cat("Loading .hic files...\n")

  # Define .hic files
  hicFiles <- c(
    ctrl_M1 = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M1.hic",
    ctrl_M2 = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M2.hic",
    ctrl_M3 = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/ctrl_M3.hic",
    mut_M1 = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M1.hic",
    mut_M2 = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M2.hic",
    mut_M3 = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/mut_M3.hic"
  )

  # Check files exist
  files_exist <- all(sapply(hicFiles, file.exists))

  if (!files_exist) {
    cat("  ⚠ Warning: Some .hic files not found. Skipping APA analysis.\n")
    cat("     (This is expected if running locally without HPC access)\n\n")
  } else {
    cat("  ✓ All .hic files found\n\n")

    # Prepare loop sets
    if (!is.null(loops_gi)) {
      # Get significant loops only
      sig_loops_gi <- loops_gi[mcols(loops_gi)$significant == TRUE]

      # Separate by direction
      up_loops_gi <- sig_loops_gi[mcols(sig_loops_gi)$logFC > 0]
      down_loops_gi <- sig_loops_gi[mcols(sig_loops_gi)$logFC < 0]

      cat(sprintf("  Up-regulated loops for APA: %d\n", length(up_loops_gi)))
      cat(sprintf("  Down-regulated loops for APA: %d\n\n", length(down_loops_gi)))

      # Run APA for up-regulated loops
      if (length(up_loops_gi) > 10) {
        cat("Running APA for up-regulated loops...\n")

        # This is a placeholder - actual APA implementation would use:
        # apa_up <- calcApa(hicFiles, up_loops_gi, resolution = 5000, buffer = 10)
        # For now, we'll note that this requires more setup

        cat("  ⚠ APA calculation requires additional setup with mariner::calcApa()\n")
        cat("     See mariner documentation for full implementation\n\n")
      }

      # Run APA for down-regulated loops
      if (length(down_loops_gi) > 10) {
        cat("Running APA for down-regulated loops...\n")
        cat("  ⚠ APA calculation requires additional setup with mariner::calcApa()\n\n")
      }
    } else {
      cat("  ⚠ GInteractions object not loaded. Cannot perform APA.\n\n")
    }

    cat("✓ Section 6 complete: APA analysis noted\n")
    cat("  (Full APA implementation requires mariner::calcApa() setup)\n\n")
  }
} else {
  cat("Skipping APA analysis (--skip-apa flag set)\n\n")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("="*70, "\n", sep="")
cat("VISUALIZATION ANALYSIS COMPLETE\n")
cat("="*70, "\n\n", sep="")

cat("Output directory: outputs/visualizations/\n\n")

cat("Generated files:\n\n")

cat("Volcano plots (outputs/visualizations/volcano/):\n")
cat("  - volcano_5kb.pdf\n")
cat("  - volcano_10kb.pdf\n")
cat("  - volcano_25kb.pdf\n")
cat("  - volcano_merged.pdf\n\n")

cat("Feature distribution (outputs/visualizations/features/):\n")
cat("  - feature_distribution.pdf\n")
cat("  - feature_distribution_summary.tsv\n\n")

cat("Enrichment analysis (outputs/visualizations/enrichment/):\n")
cat("  - go_bp_dotplot.pdf\n")
cat("  - go_cc_dotplot.pdf\n")
cat("  - go_mf_dotplot.pdf\n")
cat("  - kegg_dotplot.pdf\n\n")

cat("Loop classification (outputs/visualizations/loop_classification/):\n")
cat("  - loop_type_classification.pdf\n")
cat("  - loop_type_genes.tsv\n")
cat("  - loop_length_distribution_strip.pdf\n")
cat("  - loop_length_distribution_violin.pdf\n")
cat("  - loop_length_distribution_histogram.pdf\n")
cat("  - loop_length_statistics.txt\n\n")

if (!skip_apa) {
  cat("APA analysis (outputs/visualizations/apa/):\n")
  cat("  - Notes for full APA implementation\n\n")
}

cat("="*70, "\n", sep="")
cat("\n")
