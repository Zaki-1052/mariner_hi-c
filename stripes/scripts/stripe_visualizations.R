#!/usr/bin/env Rscript
# stripes/scripts/stripe_visualizations.R
# Stripe Visualization & Analysis Pipeline
# Author: Zakir Alibhai
# Date: 2025-12-30
#
# Purpose:
#   Generate publication-quality visualizations for differential stripe analysis:
#   1. Filter to medium/high confidence stripes
#   2. Volcano plots (relaxed thresholds + confidence colored)
#   3. Length distribution analysis
#   4. ChIP-seq-based anchor classification
#   5. GO/KEGG enrichment analysis
#   6. Summary statistics
#
# Usage:
#   Rscript scripts/stripe_visualizations.R [timepoint]
#
#   timepoint: "early" or "late" (default: runs both)
#
# Example:
#   Rscript scripts/stripe_visualizations.R early
#   Rscript scripts/stripe_visualizations.R late

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("Stripe Visualization & Analysis Pipeline\n")
cat("========================================\n\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  timepoints <- args[1]
  if (!timepoints %in% c("early", "late")) {
    stop("Invalid timepoint. Use 'early' or 'late'")
  }
  timepoints <- c(timepoints)
} else {
  timepoints <- c("early", "late")
}

cat(sprintf("Timepoints to process: %s\n\n", paste(timepoints, collapse = ", ")))

# Load required packages
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  # Core visualization
  library(ggplot2)
  library(patchwork)
  library(EnhancedVolcano)
  library(RColorBrewer)

  # Genomic annotation
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)

  # Enrichment analysis
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)

  # Data wrangling
  library(dplyr)
  library(tidyr)
  library(yaml)
})

cat("Packages loaded\n\n")

# Load configuration
config <- yaml::read_yaml("config/stripe_config.yaml")

# Set working directory
base_dir <- config$project$base_dir
if (dir.exists(base_dir)) {
  setwd(base_dir)
  cat(sprintf("Working directory: %s\n", base_dir))
} else {
  # Try stripes subdirectory
  if (dir.exists("stripes")) {
    setwd("stripes")
  }
  cat(sprintf("Using current directory: %s\n", getwd()))
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Create volcano plot with relaxed thresholds
#'
#' @param stripe_df Data frame with stripe data
#' @param output_path Output file path
#' @param title Plot title
create_stripe_volcano_relaxed <- function(stripe_df, output_path, title) {
  cat(sprintf("Creating relaxed threshold volcano plot: %s\n", basename(output_path)))

  # Calculate counts
  n_total <- nrow(stripe_df)
  n_lost <- sum(stripe_df$direction == "lost", na.rm = TRUE)
  n_gained <- sum(stripe_df$direction == "gained", na.rm = TRUE)

  # Get data ranges for annotation positioning
  logfc_range <- range(stripe_df$logFC, na.rm = TRUE)
  fdr_range <- range(stripe_df$FDR[stripe_df$FDR > 0], na.rm = TRUE)
  neg_log10_fdr_max <- -log10(min(fdr_range, na.rm = TRUE))

  # Create EnhancedVolcano with relaxed thresholds
  p <- EnhancedVolcano(
    stripe_df,
    lab = NA,
    x = 'logFC',
    y = 'FDR',
    title = title,
    subtitle = 'Relaxed thresholds: FDR < 0.10, |logFC| > 0.2',
    pCutoff = 0.10,
    FCcutoff = 0.2,
    pointSize = 2.5,
    labSize = 0,
    col = c('grey80', 'grey60', 'grey60', 'red3'),
    colAlpha = 0.7,
    legendPosition = 'top',
    legendLabSize = 10,
    legendIconSize = 3.0,
    legendLabels = c('NS', bquote(Log[2]~FC), 'p-value', bquote(p~and~log[2]~FC)),
    drawConnectors = FALSE,
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 0.8,
    borderColour = 'black',
    xlab = bquote(Log[2]~Fold~Change~(Mutant/Control)),
    ylab = bquote(-Log[10]~FDR),
    xlim = c(min(logfc_range) - 0.1, max(logfc_range) + 0.1),
    caption = NULL
  )

  # Add custom annotations
  p <- p +
    annotate(
      "text",
      x = logfc_range[1] + 0.1,
      y = neg_log10_fdr_max * 0.95,
      label = sprintf("Lost: %d", n_lost),
      color = "#3366CC",
      size = 5,
      fontface = "bold",
      hjust = 0
    ) +
    annotate(
      "text",
      x = logfc_range[2] - 0.1,
      y = neg_log10_fdr_max * 0.95,
      label = sprintf("Gained: %d", n_gained),
      color = "#CC3366",
      size = 5,
      fontface = "bold",
      hjust = 1
    ) +
    annotate(
      "text",
      x = logfc_range[2] - 0.05,
      y = 0.1,
      label = sprintf("Total: %d stripes", n_total),
      color = "black",
      size = 3.5,
      hjust = 1
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15)
    )

  ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  cat(sprintf("  Saved: %s\n", output_path))

  return(p)
}

#' Create volcano plot with confidence tier coloring
#'
#' @param stripe_df Data frame with stripe data
#' @param output_path Output file path
#' @param title Plot title
create_stripe_volcano_confidence <- function(stripe_df, output_path, title) {
  cat(sprintf("Creating confidence-colored volcano plot: %s\n", basename(output_path)))

  # Create combined color category
  stripe_df <- stripe_df %>%
    mutate(
      color_category = case_when(
        direction == "lost" & direction_confidence == "high" ~ "Lost (High)",
        direction == "lost" & direction_confidence == "medium" ~ "Lost (Medium)",
        direction == "lost" & direction_confidence == "low" ~ "Lost (Low)",
        direction == "gained" & direction_confidence == "high" ~ "Gained (High)",
        direction == "gained" & direction_confidence == "medium" ~ "Gained (Medium)",
        direction == "gained" & direction_confidence == "low" ~ "Gained (Low)",
        direction == "unchanged" ~ "Unchanged",
        TRUE ~ "Other"
      ),
      neg_log10_fdr = -log10(FDR)
    )

  # Define colors
  color_palette <- c(
    "Lost (High)" = "#0000FF",
    "Lost (Medium)" = "#6666FF",
    "Lost (Low)" = "#CCCCFF",
    "Gained (High)" = "#FF0000",
    "Gained (Medium)" = "#FF6666",
    "Gained (Low)" = "#FFCCCC",
    "Unchanged" = "#999999",
    "Other" = "#DDDDDD"
  )

  # Calculate counts for legend
  n_lost_high <- sum(stripe_df$color_category == "Lost (High)", na.rm = TRUE)
  n_lost_med <- sum(stripe_df$color_category == "Lost (Medium)", na.rm = TRUE)
  n_lost_low <- sum(stripe_df$color_category == "Lost (Low)", na.rm = TRUE)
  n_gained_high <- sum(stripe_df$color_category == "Gained (High)", na.rm = TRUE)
  n_gained_med <- sum(stripe_df$color_category == "Gained (Medium)", na.rm = TRUE)
  n_gained_low <- sum(stripe_df$color_category == "Gained (Low)", na.rm = TRUE)
  n_unchanged <- sum(stripe_df$color_category == "Unchanged", na.rm = TRUE)

  # Create plot
  p <- ggplot(stripe_df, aes(x = logFC, y = neg_log10_fdr, color = color_category)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(
      values = color_palette,
      name = "Direction (Confidence)",
      labels = c(
        "Lost (High)" = sprintf("Lost High (%d)", n_lost_high),
        "Lost (Medium)" = sprintf("Lost Medium (%d)", n_lost_med),
        "Lost (Low)" = sprintf("Lost Low (%d)", n_lost_low),
        "Gained (High)" = sprintf("Gained High (%d)", n_gained_high),
        "Gained (Medium)" = sprintf("Gained Medium (%d)", n_gained_med),
        "Gained (Low)" = sprintf("Gained Low (%d)", n_gained_low),
        "Unchanged" = sprintf("Unchanged (%d)", n_unchanged),
        "Other" = "Other"
      )
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.10), linetype = "dotted", color = "grey60", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "grey60", linewidth = 0.5) +
    labs(
      title = title,
      subtitle = "Colored by direction and confidence tier",
      x = bquote(Log[2]~Fold~Change~(Mutant/Control)),
      y = bquote(-Log[10]~FDR),
      caption = "Dashed: FDR=0.05, |logFC|=0.3 | Dotted: FDR=0.10, |logFC|=0.2"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8, color = "grey50"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 9),
      plot.margin = margin(15, 15, 15, 15)
    )

  ggsave(output_path, p, width = 12, height = 8, dpi = 300)
  cat(sprintf("  Saved: %s\n", output_path))

  return(p)
}

#' Annotate stripe anchors with ChIP-seq data
#'
#' @param stripe_df Data frame with stripe data
#' @param k27ac_path Path to H3K27ac peaks BED file
#' @param k4me1_path Path to H3K4me1 peaks BED file (can be NULL)
#' @param tss_threshold Distance to TSS for promoter classification
#' @return Data frame with annotation columns added
annotate_stripe_anchors <- function(stripe_df, k27ac_path, k4me1_path = NULL,
                                     tss_threshold = 2000) {
  cat("\nAnnotating stripe anchors with ChIP-seq data...\n")

  # Load TxDb for TSS
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes <- genes(txdb)
  tss <- resize(genes, width = 1, fix = "start")
  cat(sprintf("  Loaded %d gene TSS positions\n", length(tss)))

  # Create GRanges for stripe anchors
  anchor_gr <- GRanges(
    seqnames = stripe_df$chr,
    ranges = IRanges(start = stripe_df$anchor_x1, end = stripe_df$anchor_x2)
  )
  cat(sprintf("  Created GRanges for %d stripe anchors\n", length(anchor_gr)))

  # Calculate distance to nearest TSS
  cat("  Calculating distance to nearest TSS...\n")
  nearest_tss <- distanceToNearest(anchor_gr, tss)
  distance_to_tss <- rep(NA_real_, length(anchor_gr))
  nearest_gene_idx <- rep(NA_integer_, length(anchor_gr))

  if (length(nearest_tss) > 0) {
    distance_to_tss[queryHits(nearest_tss)] <- mcols(nearest_tss)$distance
    nearest_gene_idx[queryHits(nearest_tss)] <- subjectHits(nearest_tss)
  }

  # Get gene names
  nearest_gene <- rep(NA_character_, length(anchor_gr))
  valid_idx <- !is.na(nearest_gene_idx)
  nearest_gene[valid_idx] <- names(genes)[nearest_gene_idx[valid_idx]]

  # Load H3K27ac peaks
  cat(sprintf("  Loading H3K27ac peaks: %s\n", k27ac_path))
  if (file.exists(k27ac_path)) {
    k27ac_gr <- import(k27ac_path, format = "bed")
    h3k27ac_overlap <- countOverlaps(anchor_gr, k27ac_gr) > 0
    cat(sprintf("    Loaded %d peaks, %d anchors overlap\n",
                length(k27ac_gr), sum(h3k27ac_overlap)))
  } else {
    cat("    WARNING: File not found, setting all overlaps to FALSE\n")
    h3k27ac_overlap <- rep(FALSE, length(anchor_gr))
  }

  # Load H3K4me1 peaks (if available)
  if (!is.null(k4me1_path) && file.exists(k4me1_path)) {
    cat(sprintf("  Loading H3K4me1 peaks: %s\n", k4me1_path))
    k4me1_gr <- import(k4me1_path, format = "bed")
    h3k4me1_overlap <- countOverlaps(anchor_gr, k4me1_gr) > 0
    cat(sprintf("    Loaded %d peaks, %d anchors overlap\n",
                length(k4me1_gr), sum(h3k4me1_overlap)))
  } else {
    cat("  H3K4me1 peaks not available, setting all overlaps to FALSE\n")
    h3k4me1_overlap <- rep(FALSE, length(anchor_gr))
  }

  # Classify anchor types
  cat("  Classifying anchor types...\n")
  anchor_type <- rep("Other", length(anchor_gr))

  # Promoter: H3K27ac+ AND <=2kb from TSS
  is_promoter <- h3k27ac_overlap & !is.na(distance_to_tss) & distance_to_tss <= tss_threshold
  anchor_type[is_promoter] <- "Promoter"

  # Active Enhancer: H3K27ac+ AND >2kb from TSS
  is_active_enhancer <- h3k27ac_overlap & (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_active_enhancer] <- "Active_Enhancer"

  # Poised Enhancer: H3K4me1+ (no H3K27ac) AND >2kb from TSS
  is_poised_enhancer <- !h3k27ac_overlap & h3k4me1_overlap &
                        (is.na(distance_to_tss) | distance_to_tss > tss_threshold)
  anchor_type[is_poised_enhancer] <- "Poised_Enhancer"

  # Print summary
  cat("\n  Anchor type distribution:\n")
  for (type in c("Promoter", "Active_Enhancer", "Poised_Enhancer", "Other")) {
    count <- sum(anchor_type == type)
    pct <- 100 * count / length(anchor_type)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", type, count, pct))
  }

  # Add columns to data frame
  stripe_df$distance_to_tss <- distance_to_tss
  stripe_df$nearest_gene <- nearest_gene
  stripe_df$h3k27ac_overlap <- h3k27ac_overlap
  stripe_df$h3k4me1_overlap <- h3k4me1_overlap
  stripe_df$anchor_type <- anchor_type

  cat("\n  Annotation complete\n")

  return(stripe_df)
}

# =============================================================================
# MAIN PROCESSING LOOP
# =============================================================================

for (timepoint in timepoints) {

  cat("\n")
  cat("########################################\n")
  cat(sprintf("Processing: %s timepoint\n", toupper(timepoint)))
  cat("########################################\n\n")

  # Define paths
  input_file <- sprintf("outputs/%s/res_5kb/04_final_differential_stripes.tsv", timepoint)
  output_dir <- sprintf("outputs/%s/res_5kb", timepoint)
  plots_dir <- file.path(output_dir, "plots")

  # Create output directories
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Check input file exists
  if (!file.exists(input_file)) {
    cat(sprintf("WARNING: Input file not found: %s\n", input_file))
    cat("Skipping this timepoint.\n")
    next
  }

  # =============================================================================
  # SECTION 1: LOAD DATA AND FILTER
  # =============================================================================

  cat("\n========================================\n")
  cat("Section 1: Loading and Filtering Data\n")
  cat("========================================\n\n")

  # Load stripe data
  stripe_df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("Loaded %d stripes from %s\n", nrow(stripe_df), basename(input_file)))

  # Summary by direction
  cat("\nStripe counts by direction:\n")
  direction_counts <- table(stripe_df$direction)
  for (dir in names(direction_counts)) {
    cat(sprintf("  %-12s: %d\n", dir, direction_counts[dir]))
  }

  # Summary by confidence
  cat("\nStripe counts by direction_confidence:\n")
  conf_counts <- table(stripe_df$direction_confidence)
  for (conf in c("high", "medium", "low")) {
    if (conf %in% names(conf_counts)) {
      cat(sprintf("  %-12s: %d\n", conf, conf_counts[conf]))
    }
  }

  # Filter to medium/high confidence
  stripe_medium_high <- stripe_df %>%
    filter(direction_confidence %in% c("high", "medium"))

  cat(sprintf("\nFiltered to medium/high confidence: %d stripes\n", nrow(stripe_medium_high)))

  # Save filtered list
  filtered_file <- file.path(output_dir, "05_medium_high_confidence_stripes.tsv")
  write.table(stripe_medium_high, filtered_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved: %s\n", filtered_file))

  # =============================================================================
  # SECTION 2: VOLCANO PLOTS
  # =============================================================================

  cat("\n========================================\n")
  cat("Section 2: Volcano Plots\n")
  cat("========================================\n\n")

  # Version A: Relaxed thresholds
  volcano_relaxed_file <- file.path(plots_dir, "volcano_relaxed_thresholds.pdf")
  create_stripe_volcano_relaxed(
    stripe_df,
    volcano_relaxed_file,
    sprintf("Differential Stripes - %s Timepoint", tools::toTitleCase(timepoint))
  )

  # Version B: Confidence colored
  volcano_conf_file <- file.path(plots_dir, "volcano_confidence_colored.pdf")
  create_stripe_volcano_confidence(
    stripe_df,
    volcano_conf_file,
    sprintf("Differential Stripes - %s Timepoint", tools::toTitleCase(timepoint))
  )

  cat("Volcano plots complete\n")

  # =============================================================================
  # SECTION 3: LENGTH DISTRIBUTION ANALYSIS
  # =============================================================================

  cat("\n========================================\n")
  cat("Section 3: Length Distribution Analysis\n")
  cat("========================================\n\n")

  # Prepare data for length analysis (lost vs gained only)
  length_df <- stripe_df %>%
    filter(direction %in% c("lost", "gained")) %>%
    mutate(
      direction_label = ifelse(direction == "lost", "Lost", "Gained"),
      stripe_length_kb = stripe_length / 1000,
      anchor_width_kb = anchor_width / 1000
    )

  cat(sprintf("Analyzing %d stripes (lost + gained)\n", nrow(length_df)))

  # Statistical tests
  if (nrow(filter(length_df, direction == "lost")) > 1 &&
      nrow(filter(length_df, direction == "gained")) > 1) {

    wilcox_length <- wilcox.test(
      length_df$stripe_length[length_df$direction == "lost"],
      length_df$stripe_length[length_df$direction == "gained"]
    )

    wilcox_width <- wilcox.test(
      length_df$anchor_width[length_df$direction == "lost"],
      length_df$anchor_width[length_df$direction == "gained"]
    )

    cat(sprintf("Stripe length: Wilcoxon p = %.3e\n", wilcox_length$p.value))
    cat(sprintf("Anchor width: Wilcoxon p = %.3e\n", wilcox_width$p.value))

    # Summary statistics
    stats_df <- length_df %>%
      group_by(direction_label) %>%
      summarise(
        n = n(),
        stripe_length_median = median(stripe_length) / 1000,
        stripe_length_mean = mean(stripe_length) / 1000,
        stripe_length_min = min(stripe_length) / 1000,
        stripe_length_max = max(stripe_length) / 1000,
        anchor_width_median = median(anchor_width) / 1000,
        anchor_width_mean = mean(anchor_width) / 1000,
        .groups = "drop"
      )

    stats_df$stripe_length_wilcox_p <- wilcox_length$p.value
    stats_df$anchor_width_wilcox_p <- wilcox_width$p.value

    # Save statistics
    stats_file <- file.path(plots_dir, "length_statistics.tsv")
    write.table(stats_df, stats_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("Saved: %s\n", stats_file))

  } else {
    cat("Insufficient data for statistical tests\n")
    wilcox_length <- list(p.value = NA)
    wilcox_width <- list(p.value = NA)
  }

  # Color palette
  direction_colors <- c("Lost" = "#4575b4", "Gained" = "#d73027")

  # Plot 1: Strip plot
  p_strip <- ggplot(length_df, aes(x = stripe_length_kb, y = direction_label, color = direction_label)) +
    geom_jitter(alpha = 0.6, size = 2.5, height = 0.2) +
    scale_color_manual(values = direction_colors) +
    labs(
      title = sprintf("Stripe Length Distribution - %s", tools::toTitleCase(timepoint)),
      subtitle = sprintf("Wilcoxon p = %.3e", wilcox_length$p.value),
      x = "Stripe Length (kb)",
      y = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )

  ggsave(file.path(plots_dir, "length_distribution_strip.pdf"), p_strip, width = 10, height = 5, dpi = 300)
  cat("Saved: length_distribution_strip.pdf\n")

  # Plot 2: Violin plot
  p_violin <- ggplot(length_df, aes(x = direction_label, y = stripe_length_kb, fill = direction_label)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = direction_colors) +
    labs(
      title = sprintf("Stripe Length Distribution - %s", tools::toTitleCase(timepoint)),
      subtitle = sprintf("Wilcoxon p = %.3e", wilcox_length$p.value),
      x = "Direction",
      y = "Stripe Length (kb)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )

  ggsave(file.path(plots_dir, "length_distribution_violin.pdf"), p_violin, width = 8, height = 6, dpi = 300)
  cat("Saved: length_distribution_violin.pdf\n")

  # Plot 3: Histogram
  p_hist <- ggplot(length_df, aes(x = stripe_length_kb, fill = direction_label)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = direction_colors) +
    facet_wrap(~direction_label, ncol = 1) +
    labs(
      title = sprintf("Stripe Length Distribution - %s", tools::toTitleCase(timepoint)),
      x = "Stripe Length (kb)",
      y = "Count",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )

  ggsave(file.path(plots_dir, "length_distribution_histogram.pdf"), p_hist, width = 10, height = 8, dpi = 300)
  cat("Saved: length_distribution_histogram.pdf\n")

  cat("Length distribution analysis complete\n")

  # =============================================================================
  # SECTION 4: CHIP-SEQ ANCHOR CLASSIFICATION
  # =============================================================================

  cat("\n========================================\n")
  cat("Section 4: ChIP-seq Anchor Classification\n")
  cat("========================================\n\n")

  # Get ChIP-seq file paths based on timepoint
  if (timepoint == "early") {
    k27ac_path <- config$chipseq_peaks$h3k27ac_p12
    k4me1_path <- NULL  # Not available for early timepoint
  } else {
    k27ac_path <- config$chipseq_peaks$h3k27ac
    k4me1_path <- config$chipseq_peaks$h3k4me1
  }

  # Annotate stripes
  stripe_annotated <- annotate_stripe_anchors(
    stripe_df,
    k27ac_path = k27ac_path,
    k4me1_path = k4me1_path,
    tss_threshold = 2000
  )

  # Save annotated stripes
  annotated_file <- file.path(output_dir, "05_all_annotated_stripes.tsv")
  write.table(stripe_annotated, annotated_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved: %s\n", annotated_file))

  # Create anchor type distribution plot
  anchor_type_counts <- stripe_annotated %>%
    count(anchor_type) %>%
    mutate(percentage = 100 * n / sum(n))

  type_colors <- c(
    "Promoter" = "#e41a1c",
    "Active_Enhancer" = "#377eb8",
    "Poised_Enhancer" = "#4daf4a",
    "Other" = "#999999"
  )

  p_anchor_type <- ggplot(anchor_type_counts, aes(x = reorder(anchor_type, -n), y = n, fill = anchor_type)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_text(aes(label = sprintf("%d (%.1f%%)", n, percentage)), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = type_colors) +
    labs(
      title = sprintf("Stripe Anchor Type Distribution - %s", tools::toTitleCase(timepoint)),
      subtitle = sprintf("n = %d stripes", nrow(stripe_annotated)),
      x = "Anchor Type",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(file.path(plots_dir, "anchor_type_distribution.pdf"), p_anchor_type, width = 8, height = 6, dpi = 300)
  cat("Saved: anchor_type_distribution.pdf\n")

  # Create anchor type by direction plot
  anchor_direction <- stripe_annotated %>%
    filter(direction %in% c("lost", "gained")) %>%
    count(direction, anchor_type) %>%
    group_by(direction) %>%
    mutate(percentage = 100 * n / sum(n)) %>%
    ungroup()

  p_anchor_direction <- ggplot(anchor_direction, aes(x = direction, y = n, fill = anchor_type)) +
    geom_bar(stat = "identity", position = "stack", color = "black", alpha = 0.8) +
    scale_fill_manual(values = type_colors) +
    labs(
      title = sprintf("Anchor Type by Stripe Direction - %s", tools::toTitleCase(timepoint)),
      x = "Direction",
      y = "Count",
      fill = "Anchor Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plots_dir, "anchor_type_by_direction.pdf"), p_anchor_direction, width = 8, height = 6, dpi = 300)
  cat("Saved: anchor_type_by_direction.pdf\n")

  cat("ChIP-seq classification complete\n")

  # =============================================================================
  # SECTION 5: GO/KEGG ENRICHMENT ANALYSIS
  # =============================================================================

  cat("\n========================================\n")
  cat("Section 5: GO/KEGG Enrichment Analysis\n")
  cat("========================================\n\n")

  # Get genes near stripe anchors
  cat("Extracting genes near stripe anchors...\n")

  # Get unique genes from annotated stripes
  lost_genes <- stripe_annotated %>%
    filter(direction == "lost") %>%
    pull(nearest_gene) %>%
    unique() %>%
    na.omit()

  gained_genes <- stripe_annotated %>%
    filter(direction == "gained") %>%
    pull(nearest_gene) %>%
    unique() %>%
    na.omit()

  cat(sprintf("  Lost stripe genes: %d\n", length(lost_genes)))
  cat(sprintf("  Gained stripe genes: %d\n", length(gained_genes)))

  # Convert to Entrez IDs (genes from TxDb are already Entrez IDs)
  lost_entrez <- as.character(lost_genes)
  gained_entrez <- as.character(gained_genes)

  # Check if we have enough genes
  min_genes <- 5
  run_enrichment <- length(lost_entrez) >= min_genes && length(gained_entrez) >= min_genes

  if (run_enrichment) {
    cat("\nRunning enrichment analysis...\n")

    gene_list <- list(
      lost = lost_entrez,
      gained = gained_entrez
    )

    enrichment_results <- list()

    # GO Biological Process
    cat("  GO Biological Process...\n")
    tryCatch({
      go_bp <- compareCluster(
        geneCluster = gene_list,
        fun = "enrichGO",
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.2
      )

      if (!is.null(go_bp) && nrow(go_bp@compareClusterResult) > 0) {
        p_go_bp <- dotplot(go_bp, showCategory = 15) +
          labs(title = sprintf("GO Biological Process - %s", tools::toTitleCase(timepoint))) +
          theme(plot.title = element_text(face = "bold"))
        ggsave(file.path(plots_dir, "go_bp_dotplot.pdf"), p_go_bp, width = 12, height = 10, dpi = 300)
        cat("    Saved: go_bp_dotplot.pdf\n")
        enrichment_results$go_bp <- go_bp@compareClusterResult
      } else {
        cat("    No significant GO BP terms found\n")
      }
    }, error = function(e) {
      cat(sprintf("    Error in GO BP: %s\n", e$message))
    })

    # GO Cellular Component
    cat("  GO Cellular Component...\n")
    tryCatch({
      go_cc <- compareCluster(
        geneCluster = gene_list,
        fun = "enrichGO",
        OrgDb = org.Mm.eg.db,
        ont = "CC",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.2
      )

      if (!is.null(go_cc) && nrow(go_cc@compareClusterResult) > 0) {
        p_go_cc <- dotplot(go_cc, showCategory = 15) +
          labs(title = sprintf("GO Cellular Component - %s", tools::toTitleCase(timepoint))) +
          theme(plot.title = element_text(face = "bold"))
        ggsave(file.path(plots_dir, "go_cc_dotplot.pdf"), p_go_cc, width = 10, height = 8, dpi = 300)
        cat("    Saved: go_cc_dotplot.pdf\n")
        enrichment_results$go_cc <- go_cc@compareClusterResult
      } else {
        cat("    No significant GO CC terms found\n")
      }
    }, error = function(e) {
      cat(sprintf("    Error in GO CC: %s\n", e$message))
    })

    # GO Molecular Function
    cat("  GO Molecular Function...\n")
    tryCatch({
      go_mf <- compareCluster(
        geneCluster = gene_list,
        fun = "enrichGO",
        OrgDb = org.Mm.eg.db,
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.2
      )

      if (!is.null(go_mf) && nrow(go_mf@compareClusterResult) > 0) {
        p_go_mf <- dotplot(go_mf, showCategory = 15) +
          labs(title = sprintf("GO Molecular Function - %s", tools::toTitleCase(timepoint))) +
          theme(plot.title = element_text(face = "bold"))
        ggsave(file.path(plots_dir, "go_mf_dotplot.pdf"), p_go_mf, width = 10, height = 8, dpi = 300)
        cat("    Saved: go_mf_dotplot.pdf\n")
        enrichment_results$go_mf <- go_mf@compareClusterResult
      } else {
        cat("    No significant GO MF terms found\n")
      }
    }, error = function(e) {
      cat(sprintf("    Error in GO MF: %s\n", e$message))
    })

    # KEGG pathways
    cat("  KEGG pathways...\n")
    tryCatch({
      kegg <- compareCluster(
        geneCluster = gene_list,
        fun = "enrichKEGG",
        organism = "mmu",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.2
      )

      if (!is.null(kegg) && nrow(kegg@compareClusterResult) > 0) {
        p_kegg <- dotplot(kegg, showCategory = 15) +
          labs(title = sprintf("KEGG Pathways - %s", tools::toTitleCase(timepoint))) +
          theme(plot.title = element_text(face = "bold"))
        ggsave(file.path(plots_dir, "kegg_dotplot.pdf"), p_kegg, width = 12, height = 10, dpi = 300)
        cat("    Saved: kegg_dotplot.pdf\n")
        enrichment_results$kegg <- kegg@compareClusterResult
      } else {
        cat("    No significant KEGG pathways found\n")
      }
    }, error = function(e) {
      cat(sprintf("    Error in KEGG: %s\n", e$message))
    })

    # Save enrichment results
    if (length(enrichment_results) > 0) {
      enrichment_df <- bind_rows(enrichment_results, .id = "ontology")
      write.table(enrichment_df, file.path(output_dir, "enrichment_results.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat(sprintf("Saved: enrichment_results.tsv (%d terms)\n", nrow(enrichment_df)))
    }

  } else {
    cat(sprintf("Skipping enrichment analysis: insufficient genes (need >= %d per group)\n", min_genes))
  }

  cat("Enrichment analysis complete\n")

  # =============================================================================
  # SECTION 6: SUMMARY STATISTICS
  # =============================================================================

  cat("\n========================================\n")
  cat("Section 6: Summary Statistics\n")
  cat("========================================\n\n")

  # Generate summary text
  summary_lines <- c(
    "================================================================================",
    sprintf("STRIPE VISUALIZATION SUMMARY - %s TIMEPOINT", toupper(timepoint)),
    "================================================================================",
    sprintf("Date: %s", Sys.Date()),
    "",
    "STRIPE COUNTS",
    "-------------",
    sprintf("Total stripes: %d", nrow(stripe_df)),
    "",
    "By direction:",
    sprintf("  Lost (control_only):    %d", sum(stripe_df$direction == "lost")),
    sprintf("  Gained (mutant_only):   %d", sum(stripe_df$direction == "gained")),
    sprintf("  Unchanged (shared):     %d", sum(stripe_df$direction == "unchanged")),
    "",
    "By confidence tier:",
    sprintf("  High confidence:   %d", sum(stripe_df$direction_confidence == "high")),
    sprintf("  Medium confidence: %d", sum(stripe_df$direction_confidence == "medium")),
    sprintf("  Low confidence:    %d", sum(stripe_df$direction_confidence == "low")),
    "",
    sprintf("Medium/High confidence stripes: %d", nrow(stripe_medium_high)),
    "",
    "ANCHOR CLASSIFICATION",
    "---------------------"
  )

  for (type in c("Promoter", "Active_Enhancer", "Poised_Enhancer", "Other")) {
    count <- sum(stripe_annotated$anchor_type == type)
    pct <- 100 * count / nrow(stripe_annotated)
    summary_lines <- c(summary_lines, sprintf("  %-20s: %5d (%.1f%%)", type, count, pct))
  }

  summary_lines <- c(summary_lines,
    "",
    "LENGTH STATISTICS",
    "-----------------",
    sprintf("Stripe length (median): %.1f kb", median(stripe_df$stripe_length) / 1000),
    sprintf("Stripe length (range): %.1f - %.1f kb",
            min(stripe_df$stripe_length) / 1000, max(stripe_df$stripe_length) / 1000),
    sprintf("Anchor width (median): %.1f kb", median(stripe_df$anchor_width) / 1000),
    "",
    "DIFFERENTIAL STATISTICS",
    "-----------------------",
    sprintf("logFC range: %.3f to %.3f", min(stripe_df$logFC), max(stripe_df$logFC)),
    sprintf("logFC median: %.3f", median(stripe_df$logFC)),
    sprintf("FDR < 0.05: %d stripes", sum(stripe_df$FDR < 0.05)),
    sprintf("FDR < 0.10: %d stripes", sum(stripe_df$FDR < 0.10)),
    "",
    "DIRECTIONAL CONSISTENCY",
    "-----------------------",
    sprintf("Consistent (logFC matches direction): %d (%.1f%%)",
            sum(stripe_df$direction_consistent, na.rm = TRUE),
            100 * mean(stripe_df$direction_consistent, na.rm = TRUE)),
    "",
    "OUTPUT FILES",
    "------------",
    sprintf("05_medium_high_confidence_stripes.tsv: %d stripes", nrow(stripe_medium_high)),
    sprintf("05_all_annotated_stripes.tsv: %d stripes", nrow(stripe_annotated)),
    "Volcano plots: volcano_relaxed_thresholds.pdf, volcano_confidence_colored.pdf",
    "Length distribution: length_distribution_*.pdf, length_statistics.tsv",
    "Anchor classification: anchor_type_distribution.pdf, anchor_type_by_direction.pdf",
    "Enrichment: go_*_dotplot.pdf, kegg_dotplot.pdf, enrichment_results.tsv",
    "",
    "================================================================================"
  )

  # Write summary file
  summary_file <- file.path(output_dir, "05_visualization_summary.txt")
  writeLines(summary_lines, summary_file)
  cat(sprintf("Saved: %s\n", summary_file))

  # Print summary to console
  cat("\n")
  cat(paste(summary_lines, collapse = "\n"))
  cat("\n")

}

# =============================================================================
# FINAL MESSAGE
# =============================================================================

cat("\n")
cat("########################################\n")
cat("STRIPE VISUALIZATION PIPELINE COMPLETE\n")
cat("########################################\n\n")

cat("Processed timepoints:\n")
for (tp in timepoints) {
  cat(sprintf("  - %s\n", tp))
}

cat("\nOutput directories:\n")
for (tp in timepoints) {
  cat(sprintf("  - outputs/%s/res_5kb/\n", tp))
  cat(sprintf("    - plots/\n"))
}

cat("\n")
