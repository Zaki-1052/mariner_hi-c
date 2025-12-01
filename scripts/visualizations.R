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
  library(RColorBrewer)  # For color palettes with 10+ categories
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
  library(yaml)
})

cat("✓ Packages loaded\n\n")

# Load paths configuration and set working directory
config <- yaml::read_yaml("config/paths_config.yaml")
base_dir <- config$project$base_dir
if (dir.exists(base_dir)) {
  setwd(base_dir)
  cat(sprintf("  Working directory: %s\n", base_dir))
} else {
  warning("Configured base directory does not exist: ", base_dir)
  cat(sprintf("  Using current directory: %s\n", getwd()))
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

cat("\n========================================\n")
cat("Loading Input Data\n")
cat("\n========================================\n")

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

cat("\n========================================\n")
cat("SECTION 1: Volcano Plots\n")
cat("\n========================================\n")

# Function to create publication-quality volcano plot matching screenshot styling
create_publication_volcano <- function(results_file, output_path, title_suffix = "") {
  cat(sprintf("Creating enhanced volcano plot: %s\n", basename(output_path)))

  # Validate input file
  if (!file.exists(results_file)) {
    cat(sprintf("  ✗ Results file not found: %s\n", results_file))
    return(NULL)
  }

  # Load data
  df <- read.table(results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d loops\n", nrow(df)))

  # Calculate counts for annotations
  # Significant = FDR < 0.05 AND |logFC| > 0.3 (matches plot_color "red"/"blue")
  n_total <- nrow(df)
  n_up <- sum(df$significant & df$logFC > 0, na.rm = TRUE)
  n_down <- sum(df$significant & df$logFC < 0, na.rm = TRUE)

  cat(sprintf("  Significant: %d up, %d down\n", n_up, n_down))

  # Get data ranges for annotation positioning
  logfc_range <- range(df$logFC, na.rm = TRUE)
  fdr_range <- range(df$FDR[df$FDR > 0], na.rm = TRUE)
  neg_log10_fdr_max <- -log10(min(fdr_range))

  # Create EnhancedVolcano with publication settings
  p <- EnhancedVolcano(
    df,
    lab = NA,  # Hide individual labels for cleaner look
    x = 'logFC',
    y = 'FDR',
    title = paste0('WT vs KO ', title_suffix, ' Differential Loops'),
    subtitle = 'EnhancedVolcano',
    pCutoff = 0.05,
    FCcutoff = 0.3,
    pointSize = 2.5,
    labSize = 0,  # No point labels
    col = c('black', 'grey', 'red', 'darkred'),  # NS, FC-only, p-only, both
    colAlpha = 0.6,
    legendPosition = 'top',
    legendLabSize = 11,
    legendIconSize = 4.0,
    legendLabels = c('NS', bquote(Log[2]~FC), 'p-value', bquote(p-value~and~log[2]~FC)),
    drawConnectors = FALSE,  # Cleaner look
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 0.8,
    borderColour = 'black',
    xlab = 'LogFC Pixel Enrichment',
    ylab = 'FDR',
    xlim = c(logfc_range[1] - 0.2, logfc_range[2] + 0.2),
    caption = NULL  # Remove default caption
  )

  # Add custom text annotations
  p <- p +
    # Up-regulated count (top left, blue)
    annotate(
      "text",
      x = logfc_range[1] + 0.3,
      y = neg_log10_fdr_max * 0.98,
      label = as.character(n_up),
      color = "#3366CC",
      size = 6,
      fontface = "bold",
      hjust = 0
    ) +
    # Down-regulated count (top right, blue)
    annotate(
      "text",
      x = logfc_range[2] - 0.3,
      y = neg_log10_fdr_max * 0.98,
      label = as.character(n_down),
      color = "#3366CC",
      size = 6,
      fontface = "bold",
      hjust = 1
    ) +
    # Total variables (bottom right, black)
    annotate(
      "text",
      x = logfc_range[2] - 0.1,
      y = 0.15,
      label = sprintf("total = %d variables", n_total),
      color = "black",
      size = 3.5,
      hjust = 1
    ) +
    # Enhanced theme
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(b = 10)),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
      axis.text = element_text(size = 10),
      legend.position = "top",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    )

  # Save plot
  ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  cat(sprintf("  ✓ Saved: %s\n\n", output_path))

  return(p)
}

# Legacy function for backward compatibility
create_volcano_plot <- function(resolution_kb, output_path) {
  results_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv",
                          resolution_kb)
  title_suffix <- sprintf("%dkb", resolution_kb)
  return(create_publication_volcano(results_file, output_path, title_suffix))
}

# Generate volcano plots for each resolution
for (res_kb in c(5, 10, 25)) {
  output_file <- file.path(output_dir, "volcano", sprintf("volcano_%dkb.pdf", res_kb))
  create_volcano_plot(res_kb, output_file)
}

# Create merged multi-resolution volcano plot
cat("Creating merged multi-resolution volcano plot...\n")

# Check if merged_all_results exists (created by downstream_analysis.R)
merged_all_results_file <- "outputs/merged_loops/merged_all_results.tsv"
if (file.exists(merged_all_results_file)) {
  # Use the properly merged non-redundant dataset
  output_file <- file.path(output_dir, "volcano", "volcano_merged_multiresolution.pdf")
  create_publication_volcano(merged_all_results_file, output_file, "Multi-Resolution")
} else {
  # Fallback: create simple merged dataset
  cat("  Note: merged_all_results.tsv not found. Creating simple merged dataset...\n")

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

    # Save temporary file for function to read
    temp_file <- tempfile(fileext = ".tsv")
    write.table(merged_df, temp_file, sep = "\t", quote = FALSE, row.names = FALSE)

    output_file <- file.path(output_dir, "volcano", "volcano_merged_multiresolution.pdf")
    create_publication_volcano(temp_file, output_file, "Multi-Resolution")

    # Clean up
    unlink(temp_file)
  } else {
    cat("  ⚠ No resolution files found for merged volcano plot\n\n")
  }
}

# Create volcano plot for merged non-redundant loops
cat("Creating volcano plot for merged non-redundant loops...\n")

characterized_file <- "outputs/merged_loops/characterized_loops.tsv"
if (file.exists(characterized_file)) {
  output_file <- file.path(output_dir, "volcano", "volcano_nonredundant.pdf")
  create_publication_volcano(characterized_file, output_file, "Non-Redundant")
} else {
  cat("  ⚠ characterized_loops.tsv not found. Skipping non-redundant volcano plot.\n\n")
}

cat("✓ Section 1 complete: Enhanced volcano plots generated\n\n")

# =============================================================================
# SECTION 2: FEATURE DISTRIBUTION ANALYSIS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 2: Feature Distribution Analysis\n")
cat("\n========================================\n")

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

# Extract loops with promoter-annotated anchors
promoter_loops <- loops_df %>%
  filter(grepl("Promoter", anchor1_annotation) | grepl("Promoter", anchor2_annotation))

# Get genes from promoter anchors
promoter_genes <- data.frame(
  loop_id = promoter_loops$loop_id,
  direction = ifelse(promoter_loops$logFC > 0, "up", "down"),
  anchor1_annotation = promoter_loops$anchor1_annotation,
  anchor2_annotation = promoter_loops$anchor2_annotation,
  anchor1_gene = anno1_df$SYMBOL[match(promoter_loops$loop_id, loops_df$loop_id)],
  anchor2_gene = anno2_df$SYMBOL[match(promoter_loops$loop_id, loops_df$loop_id)]
)

write.table(promoter_genes,
            file.path(output_dir, "features", "promoter_annotated_loops.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("  Saved %d promoter-annotated loops to features/promoter_annotated_loops.tsv\n", nrow(promoter_genes)))

# === Get ALL overlapping genes per anchor (not just nearest) ===
cat("  Extracting all overlapping genes per anchor...\n")

# Get all genes from TxDb
all_genes_gr <- genes(txdb)

# Function to get all overlapping genes for a GRanges object
get_all_overlapping_genes <- function(anchor_gr, genes_gr, org_db) {
  # Find direct overlaps (no gap)
  overlaps <- findOverlaps(anchor_gr, genes_gr, maxgap = 0)

  # Build dataframe with anchor index and gene info
  result <- data.frame(
    anchor_idx = queryHits(overlaps),
    entrez_id = names(genes_gr)[subjectHits(overlaps)],
    stringsAsFactors = FALSE
  )

  # Map Entrez IDs to gene symbols
  if (nrow(result) > 0) {
    symbol_map <- AnnotationDbi::select(org_db,
                                         keys = unique(result$entrez_id),
                                         columns = "SYMBOL",
                                         keytype = "ENTREZID")
    result$symbol <- symbol_map$SYMBOL[match(result$entrez_id, symbol_map$ENTREZID)]
  } else {
    result$symbol <- character(0)
  }

  return(result)
}

# Get overlapping genes for both anchors
anchor1_all_genes_df <- get_all_overlapping_genes(anchor1_gr, all_genes_gr, org.Mm.eg.db)
anchor2_all_genes_df <- get_all_overlapping_genes(anchor2_gr, all_genes_gr, org.Mm.eg.db)

# Collapse to comma-separated gene lists per anchor
collapse_genes <- function(gene_df, n_anchors) {
  # Create empty vector for all anchors
  gene_lists <- rep(NA_character_, n_anchors)

  if (nrow(gene_df) > 0) {
    collapsed <- aggregate(symbol ~ anchor_idx, data = gene_df,
                           FUN = function(x) paste(unique(na.omit(x)), collapse = ","))
    gene_lists[collapsed$anchor_idx] <- collapsed$symbol
  }

  return(gene_lists)
}

# Add to loops_df
loops_df$anchor1_all_genes <- collapse_genes(anchor1_all_genes_df, nrow(loops_df))
loops_df$anchor2_all_genes <- collapse_genes(anchor2_all_genes_df, nrow(loops_df))

# Count genes per anchor
loops_df$anchor1_gene_count <- sapply(strsplit(as.character(loops_df$anchor1_all_genes), ","),
                                       function(x) sum(!is.na(x) & x != ""))
loops_df$anchor2_gene_count <- sapply(strsplit(as.character(loops_df$anchor2_all_genes), ","),
                                       function(x) sum(!is.na(x) & x != ""))

# Save complete annotation file
all_genes_annotation <- loops_df %>%
  mutate(direction = ifelse(logFC > 0, "up", "down")) %>%
  select(loop_id, logFC, FDR, direction,
         anchor1_chr, anchor1_start, anchor1_end,
         anchor1_annotation, anchor1_all_genes, anchor1_gene_count,
         anchor2_chr, anchor2_start, anchor2_end,
         anchor2_annotation, anchor2_all_genes, anchor2_gene_count)

write.table(all_genes_annotation,
            file.path(output_dir, "features", "anchors_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("  Saved complete gene annotations to features/anchors_all_genes.tsv\n"))
cat(sprintf("  Anchor1: median %d genes, max %d genes\n",
            median(loops_df$anchor1_gene_count, na.rm = TRUE),
            max(loops_df$anchor1_gene_count, na.rm = TRUE)))
cat(sprintf("  Anchor2: median %d genes, max %d genes\n",
            median(loops_df$anchor2_gene_count, na.rm = TRUE),
            max(loops_df$anchor2_gene_count, na.rm = TRUE)))

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
  "Exon" = "#e78ac3",
  "1st Intron" = "#ffff33",
  "Other Intron" = "#a65628",
  "Intron" = "#66c2a5",
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

cat("\n========================================\n")
cat("SECTION 3: GO/KEGG Enrichment Analysis\n")
cat("\n========================================\n")

# Extract genes near loop anchors
cat("Extracting genes from loop anchors...\n")

# Get genes within 10kb of each anchor
get_nearby_genes <- function(gr, txdb, max_dist = 10000) {
  genes_txdb <- genes(txdb)

  # Find overlaps
  overlaps <- findOverlaps(gr, genes_txdb, maxgap = max_dist)

  # Get gene IDs (these are already Entrez IDs from TxDb.Mmusculus.UCSC.mm10.knownGene)
  gene_ids <- names(genes_txdb)[subjectHits(overlaps)]

  # Return unique Entrez IDs (no mapping needed - IDs from TxDb are already Entrez format)
  entrez_ids <- unique(gene_ids[!is.na(gene_ids)])

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

cat("\n========================================\n")
cat("SECTION 4: Loop Type Classification\n")
cat("\n========================================\n")

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

# Create custom color palette for up to 10 loop types
# Using RColorBrewer Set3 + Paired for better distinction
n_loop_types <- length(unique(loops_df$loop_type))
if (n_loop_types <= 12) {
  loop_colors <- RColorBrewer::brewer.pal(max(3, n_loop_types), "Set3")
} else {
  # Combine Set3 (12) + Paired (12) if needed
  loop_colors <- c(RColorBrewer::brewer.pal(12, "Set3"),
                   RColorBrewer::brewer.pal(12, "Paired"))[1:n_loop_types]
}

# Create pie charts with improved labeling for 10 categories
p_up_pie <- loop_type_summary %>%
  filter(direction_label == "Up") %>%
  ggplot(aes(x = "", y = percentage, fill = loop_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = loop_colors) +
  labs(
    title = "Up-regulated Loops",
    subtitle = sprintf("n = %d", sum(up_idx)),
    fill = "Loop Type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.text = element_text(size = 8),  # Smaller legend for more categories
    legend.key.size = unit(0.4, "cm")
  ) +
  # Only show percentage labels for slices > 5% to reduce clutter
  geom_text(aes(label = ifelse(percentage > 5, sprintf("%.1f%%", percentage), "")),
            position = position_stack(vjust = 0.5), size = 3)

p_down_pie <- loop_type_summary %>%
  filter(direction_label == "Down") %>%
  ggplot(aes(x = "", y = percentage, fill = loop_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = loop_colors) +
  labs(
    title = "Down-regulated Loops",
    subtitle = sprintf("n = %d", sum(down_idx)),
    fill = "Loop Type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.text = element_text(size = 8),  # Smaller legend for more categories
    legend.key.size = unit(0.4, "cm")
  ) +
  # Only show percentage labels for slices > 5% to reduce clutter
  geom_text(aes(label = ifelse(percentage > 5, sprintf("%.1f%%", percentage), "")),
            position = position_stack(vjust = 0.5), size = 3)

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

cat("\n========================================\n")
cat("SECTION 5: Loop Length Distribution\n")
cat("\n========================================\n")

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
# APA analysis has been moved to a standalone script for comprehensive
# processing across all resolutions and loop sets.
#
# To run APA analysis:
#   Rscript scripts/apa_analysis.R
#
# Options:
#   --resolution RESOLUTION  Specific resolution (5000, 10000, 25000)
#   --loops LOOPS           Loop set ("merged", "resolution_specific")
#
# Examples:
#   # Run all resolutions and loop sets (default)
#   Rscript scripts/apa_analysis.R
#
#   # Run only 5kb resolution with merged loops
#   Rscript scripts/apa_analysis.R --resolution 5000 --loops merged
#
# Outputs: outputs/apa_results/
#
# Note: APA analysis requires .hic files and takes ~20-30 minutes per resolution

if (!skip_apa) {
  cat("\n========================================\n")
  cat("SECTION 6: Aggregate Peak Analysis (APA)\n")
  cat("========================================\n\n")
  cat("ℹ APA analysis is available as a standalone script\n")
  cat("  Run: Rscript scripts/apa_analysis.R\n")
  cat("  See script header for options and examples\n\n")
} else {
  cat("Skipping APA analysis (--skip-apa flag set)\n\n")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n========================================\n")
cat("VISUALIZATION ANALYSIS COMPLETE\n")
cat("\n========================================\n")

cat("Output directory: outputs/visualizations/\n\n")

cat("Generated files:\n\n")

cat("Volcano plots (outputs/visualizations/volcano/):\n")
cat("  - volcano_5kb.pdf\n")
cat("  - volcano_10kb.pdf\n")
cat("  - volcano_25kb.pdf\n")
cat("  - volcano_merged.pdf (all resolutions combined)\n")
cat("  - volcano_merged_nonredundant.pdf (final non-redundant loop set)\n\n")

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

cat("\n========================================\n")
cat("\n")
