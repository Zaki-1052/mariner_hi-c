#!/usr/bin/env Rscript
# scripts/compartment_volcano_plot.R
# Volcano Plot for Chromatin Compartment (PC1) Differential Analysis
# Author: Zakir Alibhai
# Date: 2025-11-24
#
# Purpose:
#   Generate publication-quality volcano plots from compartment differential
#   analysis performed using getDiffExpression.pl with PC1 values (HOMER/Hi-C).
#
#   This script creates EnhancedVolcano plots showing PC1 value differences
#   between control and mutant conditions, representing A/B compartment shifts.
#
#   Generates TWO versions of the volcano plot:
#   1. Relaxed thresholds (FDR < 0.15, |Diff| > 0.15) - exploratory analysis
#   2. Standard thresholds (FDR < 0.05, |Diff| > 0.30) - publication quality
#
# Biological Interpretation:
#   PC1 values indicate chromatin compartment identity:
#   - Positive PC1 = A compartment (active, euchromatin)
#   - Negative PC1 = B compartment (inactive, heterochromatin)
#
#   Difference interpretation (ctrl vs. mut):
#   - Positive Difference = shift toward A compartment in mutant (more active)
#   - Negative Difference = shift toward B compartment in mutant (more inactive)
#
# Input Format:
#   Tab-delimited file from getDiffExpression.pl with annotatePeaks output:
#   - PeakID (chr-position format)
#   - Chr, Start, End, Strand
#   - Gene annotations (Gene Name, Annotation, etc.)
#   - PC1 values per sample
#   - "ctrl vs. mut Difference" (effect size)
#   - "ctrl vs. mut p-value" (raw p-value)
#   - "ctrl vs. mut adj. p-value" (FDR-adjusted)
#
# Usage:
#   Rscript scripts/compartment_volcano_plot.R [INPUT_FILE] [OPTIONS]
#
# Arguments:
#   INPUT_FILE              Path to compartment differential file
#                           (default: tad_analysis/diffcompartments.txt)
#   --output DIR            Output directory (default: outputs/compartment_analysis/)
#   --title TEXT            Custom plot title (default: auto-generated)
#   --width WIDTH           Plot width in inches (default: 12)
#   --height HEIGHT         Plot height in inches (default: 10)
#   --label-genes           Label top genes on the plot (default: FALSE)
#   --n-labels N            Number of top genes to label (default: 10)
#
# Examples:
#   # Basic usage with defaults
#   Rscript scripts/compartment_volcano_plot.R
#
#   # Custom input file
#   Rscript scripts/compartment_volcano_plot.R tad_analysis/diffcompartments.txt
#
#   # With gene labeling
#   Rscript scripts/compartment_volcano_plot.R --label-genes --n-labels 20
#
# Output:
#   - compartment_volcano_relaxed.pdf - Volcano plot with relaxed thresholds
#   - compartment_volcano_standard.pdf - Volcano plot with standard thresholds
#   - compartment_significant_relaxed.tsv - Significant regions (relaxed)
#   - compartment_significant_standard.tsv - Significant regions (standard)
#   - compartment_volcano_summary.txt - Summary statistics
#   - compartment_all_annotated.tsv - Full dataset with annotations
#
# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("================================================\n")
cat("Compartment (PC1) Differential Analysis: Volcano Plot\n")
cat("================================================\n\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default parameters
input_file <- "tad_analysis/diffcompartments.txt"  # TODO: not in data/
# Original: output_dir <- "outputs/compartment_analysis"
tsv_dir <- "data/tsvs/figure_1_tads_boundaries_compartments"
plot_dir <- "data/plots/figure_1_tads_boundaries_compartments"
custom_title <- NULL
plot_width <- 12
plot_height <- 10
label_genes <- FALSE
n_labels <- 10

# Define threshold sets
thresholds <- list(
  relaxed = list(fdr = 0.15, fc = 0.15, name = "relaxed"),
  standard = list(fdr = 0.05, fc = 0.30, name = "standard")
)

# Parse arguments
i <- 1
while (i <= length(args)) {
  if (args[i] == "--output" && i < length(args)) {
    tsv_dir <- args[i + 1]  # Original: output_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--title" && i < length(args)) {
    custom_title <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--width" && i < length(args)) {
    plot_width <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--height" && i < length(args)) {
    plot_height <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--label-genes") {
    label_genes <- TRUE
    i <- i + 1
  } else if (args[i] == "--n-labels" && i < length(args)) {
    n_labels <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (!startsWith(args[i], "--")) {
    # Positional argument = input file
    input_file <- args[i]
    i <- i + 1
  } else {
    i <- i + 1
  }
}

# Validate input file exists
if (!file.exists(input_file)) {
  stop(sprintf("ERROR: Input file not found: %s", input_file))
}

cat(sprintf("Configuration:\n"))
cat(sprintf("  Input file: %s\n", input_file))
cat(sprintf("  Output directory (TSVs): %s\n", tsv_dir))  # Original: output_dir
cat(sprintf("  Output directory (plots): %s\n", plot_dir))
cat(sprintf("  Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))
cat(sprintf("  Label genes: %s\n", ifelse(label_genes, "Yes", "No")))
cat(sprintf("  Threshold sets: relaxed (FDR<0.15, |Diff|>0.15), "))
cat(sprintf("standard (FDR<0.05, |Diff|>0.30)\n\n"))

# Load required libraries
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(EnhancedVolcano)
})

# Load multi-format output utility for PDF + SVG + JPEG output
# Original: source("scripts/utils/multi_format_output.R")
source("data/scripts/_shared/multi_format_output.R")

# Create output directories
# Original: dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

cat("================================================\n")
cat("Loading Compartment Differential Expression Data\n")
cat("================================================\n\n")

# Read compartment differential expression file
cat(sprintf("Reading: %s\n", input_file))
comp_data <- read.table(input_file, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, comment.char = "",
                        check.names = FALSE, quote = "")

cat(sprintf("Loaded %d genomic regions\n\n", nrow(comp_data)))

# Display column names for debugging
cat("Column names in input file:\n")
print(names(comp_data)[1:min(10, length(names(comp_data)))])
cat("...\n\n")

# Identify key columns (handle variable naming from getDiffExpression.pl)
peak_id_col <- names(comp_data)[1]  # First column is PeakID
chr_col <- "Chr"
start_col <- "Start"
end_col <- "End"
gene_name_col <- "Gene Name"
annotation_col <- "Annotation"

difference_col <- grep("ctrl vs\\. mut Difference",
                       names(comp_data), value = TRUE, ignore.case = TRUE)[1]
pvalue_col <- grep("ctrl vs\\. mut p-value$",
                   names(comp_data), value = TRUE, ignore.case = TRUE)[1]
adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value",
                       names(comp_data), value = TRUE, ignore.case = TRUE)[1]

# Validate required columns exist
if (is.na(difference_col)) {
  stop("ERROR: Could not find 'ctrl vs. mut Difference' column")
}
if (is.na(adj_pvalue_col)) {
  stop("ERROR: Could not find 'ctrl vs. mut adj. p-value' column")
}

cat(sprintf("Identified key columns:\n"))
cat(sprintf("  Peak ID: %s\n", peak_id_col))
cat(sprintf("  Difference (effect size): %s\n", difference_col))
cat(sprintf("  Adjusted p-value: %s\n", adj_pvalue_col))

# Identify per-sample PC1 columns (bedGraph avg columns) for ctrl/mut averages
bedgraph_cols <- grep("bedGraph avg", names(comp_data), value = TRUE)
stopifnot("Expected 6 bedGraph avg columns (3 ctrl + 3 mut)" = length(bedgraph_cols) == 6)
ctrl_pc1_cols <- grep("ctrl", bedgraph_cols, value = TRUE)
mut_pc1_cols  <- grep("mut", bedgraph_cols, value = TRUE)
stopifnot("Expected 3 ctrl PC1 columns" = length(ctrl_pc1_cols) == 3)
stopifnot("Expected 3 mut PC1 columns"  = length(mut_pc1_cols) == 3)

ctrl_avg_PC1 <- rowMeans(sapply(ctrl_pc1_cols, function(col) as.numeric(comp_data[[col]])))
mut_avg_PC1  <- rowMeans(sapply(mut_pc1_cols, function(col) as.numeric(comp_data[[col]])))

cat(sprintf("  Ctrl PC1 columns: %d identified\n", length(ctrl_pc1_cols)))
cat(sprintf("  Mut PC1 columns:  %d identified\n", length(mut_pc1_cols)))
cat("\n")

# Create clean working dataframe with standardized column names
comp_df <- data.frame(
  Region_ID = comp_data[[peak_id_col]],
  Chr = if (chr_col %in% names(comp_data)) comp_data[[chr_col]] else NA,
  Start = if (start_col %in% names(comp_data)) comp_data[[start_col]] else NA,
  End = if (end_col %in% names(comp_data)) comp_data[[end_col]] else NA,
  Gene_Name = if (gene_name_col %in% names(comp_data)) {
    comp_data[[gene_name_col]]
  } else {
    NA
  },
  Annotation = if (annotation_col %in% names(comp_data)) {
    comp_data[[annotation_col]]
  } else {
    NA
  },
  Difference = as.numeric(comp_data[[difference_col]]),
  pvalue = if (!is.na(pvalue_col)) as.numeric(comp_data[[pvalue_col]]) else NA,
  adj_pvalue = as.numeric(comp_data[[adj_pvalue_col]]),
  ctrl_avg_PC1 = ctrl_avg_PC1,
  mut_avg_PC1 = mut_avg_PC1,
  stringsAsFactors = FALSE
)

# Calculate region size
if (!is.na(comp_df$Start[1]) && !is.na(comp_df$End[1])) {
  comp_df$Region_size <- comp_df$End - comp_df$Start
}

# Remove rows with NA values in critical columns
n_before <- nrow(comp_df)
comp_df <- comp_df[!is.na(comp_df$Difference) & !is.na(comp_df$adj_pvalue), ]
n_after <- nrow(comp_df)

if (n_before > n_after) {
  cat(sprintf("Removed %d regions with missing values\n", n_before - n_after))
}

# Handle zero or very small p-values (set floor to prevent Inf in -log10)
min_pval <- 1e-300
comp_df$adj_pvalue[comp_df$adj_pvalue == 0] <- min_pval
comp_df$adj_pvalue[comp_df$adj_pvalue < min_pval] <- min_pval

# Classify compartment shift direction
# Positive Difference = shift toward A compartment (more active) in mutant
# Negative Difference = shift toward B compartment (more inactive) in mutant
comp_df$direction <- ifelse(comp_df$Difference > 0,
                            "B_to_A_in_Mutant",
                            "A_to_B_in_Mutant")
comp_df$direction_label <- ifelse(comp_df$Difference > 0,
                                  "More Active (B->A)",
                                  "More Inactive (A->B)")

# Create label column for plotting (use gene name if available, otherwise region ID)
comp_df$plot_label <- ifelse(
  !is.na(comp_df$Gene_Name) & comp_df$Gene_Name != "" & comp_df$Gene_Name != "-",
  comp_df$Gene_Name,
  comp_df$Region_ID
)

n_total <- nrow(comp_df)

# Data range info
cat("\nData ranges:\n")
cat(sprintf("  Difference: [%.3f, %.3f]\n",
            min(comp_df$Difference), max(comp_df$Difference)))
cat(sprintf("  Adj. p-value: [%.2e, %.2f]\n",
            min(comp_df$adj_pvalue), max(comp_df$adj_pvalue)))
cat(sprintf("  -log10(adj. p-value): [%.2f, %.2f]\n\n",
            min(-log10(comp_df$adj_pvalue)), max(-log10(comp_df$adj_pvalue))))

# =============================================================================
# VOLCANO PLOT GENERATION FUNCTION
# =============================================================================

generate_volcano_plot <- function(df, fdr_threshold, fc_threshold, threshold_name,
                                  plot_title_base, tsv_dir, plot_dir, plot_width,
                                  plot_height, label_genes, n_labels) {
                                  # Original: output_dir parameter split into tsv_dir + plot_dir

  cat(sprintf("\n--- Generating %s threshold plot ", threshold_name))
  cat(sprintf("(FDR < %.2f, |Diff| > %.2f) ---\n", fdr_threshold, fc_threshold))

  # Calculate statistics for this threshold
  n_total <- nrow(df)
  n_sig_both <- sum((df$adj_pvalue < fdr_threshold) &
                      (abs(df$Difference) > fc_threshold))
  n_fc_up <- sum(df$Difference > fc_threshold)
  n_fc_down <- sum(df$Difference < -fc_threshold)
  n_fdr_sig <- sum(df$adj_pvalue < fdr_threshold)
  n_sig <- sum((df$adj_pvalue < fdr_threshold) |
                 (abs(df$Difference) > fc_threshold))
  n_ns <- n_total - n_sig

  cat(sprintf("  Total regions: %d\n", n_total))
  cat(sprintf("  Both criteria (FDR & FC): %d\n", n_sig_both))
  cat(sprintf("  FC significant: %d (%d B->A, %d A->B)\n",
              n_fc_up + n_fc_down, n_fc_up, n_fc_down))
  cat(sprintf("  FDR significant: %d\n", n_fdr_sig))

  # Get data ranges for annotation positioning
  diff_range <- range(df$Difference, na.rm = TRUE)
  neg_log10_fdr_max <- -log10(min(df$adj_pvalue[df$adj_pvalue > 0], na.rm = TRUE))

  # Calculate appropriate y-axis limits
  y_max <- -log10(min(df$adj_pvalue))
  y_lim <- c(0, max(2, y_max * 1.15))

  # Create plot title
  if (is.null(plot_title_base)) {
    plot_title <- sprintf("Compartment Switching: Control vs Mutant (%s)",
                          threshold_name)
  } else {
    plot_title <- sprintf("%s (%s)", plot_title_base, threshold_name)
  }

  subtitle <- sprintf(
    "PC1 Analysis | FDR < %.2f, |Diff| > %.2f | B->A = more active in mutant",
    fdr_threshold, fc_threshold
  )

  # Determine labeling
  if (label_genes) {
    lab_size <- 3
    select_lab <- df$plot_label
    # Only label top significant genes
    df$label_this <- FALSE
    sig_df <- df[(df$adj_pvalue < fdr_threshold) &
                   (abs(df$Difference) > fc_threshold), ]
    if (nrow(sig_df) > 0) {
      top_by_fc <- head(sig_df[order(-abs(sig_df$Difference)), ], n_labels)
      df$label_this[rownames(df) %in% rownames(top_by_fc)] <- TRUE
    }
    select_lab <- ifelse(df$label_this, df$plot_label, "")
  } else {
    lab_size <- 0
    select_lab <- df$plot_label
  }

  # Create EnhancedVolcano
  p <- EnhancedVolcano(
    df,
    lab = select_lab,
    x = "Difference",
    y = "adj_pvalue",
    title = plot_title,
    subtitle = subtitle,
    pCutoff = fdr_threshold,
    FCcutoff = fc_threshold,
    pointSize = 1.5,  # Smaller points for large dataset
    labSize = lab_size,
    col = c("grey80", "grey50", "steelblue", "firebrick3"),
    colAlpha = 0.5,
    legendPosition = "top",
    legendLabSize = 10,
    legendIconSize = 3.5,
    legendLabels = c("NS", "Difference", "FDR", "FDR and Difference"),
    drawConnectors = label_genes,
    widthConnectors = 0.5,
    colConnectors = "grey30",
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = "full",
    borderWidth = 0.8,
    borderColour = "black",
    xlab = bquote(bold("Ctrl vs Mut Difference (PC1)")),
    ylab = bquote(bold("-log"[10] * "(Adjusted p-value)")),
    xlim = c(diff_range[1] - 0.5, diff_range[2] + 0.5),
    ylim = y_lim,
    caption = NULL
  )

  # Add annotations
  p <- p +
    # B->A count (right side - more active in mutant)
    annotate("text", x = diff_range[2] - 0.2, y = neg_log10_fdr_max * 0.95,
             label = sprintf("B->A: %d", n_fc_up), color = "firebrick3",
             size = 5, fontface = "bold", hjust = 1) +
    # A->B count (left side - more inactive in mutant)
    annotate("text", x = diff_range[1] + 0.2, y = neg_log10_fdr_max * 0.95,
             label = sprintf("A->B: %d", n_fc_down), color = "steelblue",
             size = 5, fontface = "bold", hjust = 0) +
    # Total regions count
    annotate("text", x = diff_range[2] - 0.1, y = 0.3,
             label = sprintf("total = %d regions", n_total), color = "black",
             size = 3.5, hjust = 1) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0,
                                margin = margin(b = 2)),
      plot.subtitle = element_text(size = 10, hjust = 0,
                                   margin = margin(b = 10)),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title.x = element_text(size = 12, face = "bold",
                                  margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, face = "bold",
                                  margin = margin(r = 10)),
      axis.text = element_text(size = 10),
      legend.position = "top",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      plot.margin = margin(15, 15, 15, 15)
    )

  # Save plot in multiple formats (PDF, SVG, JPEG)
  # Original: output_base <- file.path(output_dir, sprintf("compartment_volcano_%s", threshold_name))
  output_base <- file.path(plot_dir, sprintf("1D_compartment_volcano_%s", threshold_name))
  save_multiformat_ggplot(p, output_base, width = plot_width, height = plot_height)

  # Save significant regions for this threshold
  df$significant <- (df$adj_pvalue < fdr_threshold) &
                    (abs(df$Difference) > fc_threshold)
  sig_regions <- df[df$significant, ]
  sig_regions <- sig_regions[order(sig_regions$adj_pvalue), ]

  # Original: output_sig <- file.path(output_dir, sprintf("compartment_significant_%s.tsv", threshold_name))
  output_sig <- file.path(tsv_dir,
                          sprintf("1D_compartment_significant_%s.tsv", threshold_name))
  write.table(sig_regions, output_sig, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s (%d regions)\n", output_sig, nrow(sig_regions)))

  # Return statistics for summary
  list(
    threshold_name = threshold_name,
    fdr = fdr_threshold,
    fc = fc_threshold,
    n_total = n_total,
    n_sig_both = n_sig_both,
    n_fdr_sig = n_fdr_sig,
    n_fc_up = n_fc_up,
    n_fc_down = n_fc_down,
    n_sig = n_sig,
    n_ns = n_ns
  )
}

# =============================================================================
# GENERATE VOLCANO PLOTS FOR BOTH THRESHOLD SETS
# =============================================================================

cat("================================================\n")
cat("Generating Volcano Plots\n")
cat("================================================\n")

results_list <- list()

for (thresh in thresholds) {
  result <- generate_volcano_plot(
    df = comp_df,
    fdr_threshold = thresh$fdr,
    fc_threshold = thresh$fc,
    threshold_name = thresh$name,
    plot_title_base = custom_title,
    tsv_dir = tsv_dir,  # Original: output_dir = output_dir
    plot_dir = plot_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    label_genes = label_genes,
    n_labels = n_labels
  )
  results_list[[thresh$name]] <- result
}

# =============================================================================
# SAVE COMBINED SUMMARY
# =============================================================================

cat("\n================================================\n")
cat("Saving Summary\n")
cat("================================================\n\n")

# Build summary text
summary_text <- c(
  "Compartment (PC1) Differential Analysis - Summary Statistics",
  "============================================================",
  "",
  sprintf("Analysis Date: %s", Sys.Date()),
  sprintf("Input File: %s", input_file),
  sprintf("Output Directory (TSVs): %s", tsv_dir),  # Original: output_dir
  sprintf("Output Directory (plots): %s", plot_dir),
  "",
  sprintf("Total regions analyzed: %d", n_total),
  "",
  "Biological Interpretation:",
  "  Positive Difference = shift toward A compartment (more active) in mutant",
  "  Negative Difference = shift toward B compartment (more inactive) in mutant",
  "",
  "Data Ranges:",
  sprintf("  Difference: [%.3f, %.3f]",
          min(comp_df$Difference), max(comp_df$Difference)),
  sprintf("  Adj. p-value: [%.2e, %.2f]",
          min(comp_df$adj_pvalue), max(comp_df$adj_pvalue)),
  ""
)

# Add results for each threshold
for (result in results_list) {
  summary_text <- c(summary_text,
    "------------------------------------------------------------",
    sprintf("%s Thresholds (FDR < %.2f, |Diff| > %.2f)",
            toupper(result$threshold_name), result$fdr, result$fc),
    "------------------------------------------------------------",
    sprintf("  Both criteria (FDR & FC): %d", result$n_sig_both),
    sprintf("  FDR significant only: %d", result$n_fdr_sig),
    sprintf("  FC significant: %d", result$n_fc_up + result$n_fc_down),
    sprintf("    - B->A (more active, Diff > %.2f): %d",
            result$fc, result$n_fc_up),
    sprintf("    - A->B (more inactive, Diff < -%.2f): %d",
            result$fc, result$n_fc_down),
    sprintf("  At least one criterion: %d (%.1f%%)",
            result$n_sig, 100 * result$n_sig / result$n_total),
    sprintf("  Non-significant: %d (%.1f%%)",
            result$n_ns, 100 * result$n_ns / result$n_total),
    ""
  )
}

# Add top regions by effect size
summary_text <- c(summary_text,
  "------------------------------------------------------------",
  "Top 15 Regions by Effect Size (largest |Difference|):",
  "------------------------------------------------------------"
)

comp_df_sorted <- comp_df[order(-abs(comp_df$Difference)), ]
top_regions_fc <- head(comp_df_sorted[, c("Region_ID", "Gene_Name", "Difference",
                                           "adj_pvalue", "direction_label")], 15)
for (i in seq_len(nrow(top_regions_fc))) {
  gene_info <- if (!is.na(top_regions_fc$Gene_Name[i]) &&
                   top_regions_fc$Gene_Name[i] != "" &&
                   top_regions_fc$Gene_Name[i] != "-") {
    sprintf(" (%s)", top_regions_fc$Gene_Name[i])
  } else {
    ""
  }
  summary_text <- c(summary_text,
    sprintf("%2d. %s%s | Diff=%.3f | FDR=%.2e | %s",
            i, top_regions_fc$Region_ID[i], gene_info,
            top_regions_fc$Difference[i], top_regions_fc$adj_pvalue[i],
            top_regions_fc$direction_label[i]))
}

# Add top regions by p-value
summary_text <- c(summary_text,
  "",
  "Top 15 Regions by FDR (smallest adj. p-value):",
  "-----------------------------------------------"
)
comp_df_sorted_p <- comp_df[order(comp_df$adj_pvalue), ]
top_regions_p <- head(comp_df_sorted_p[, c("Region_ID", "Gene_Name", "Difference",
                                            "adj_pvalue", "direction_label")], 15)
for (i in seq_len(nrow(top_regions_p))) {
  gene_info <- if (!is.na(top_regions_p$Gene_Name[i]) &&
                   top_regions_p$Gene_Name[i] != "" &&
                   top_regions_p$Gene_Name[i] != "-") {
    sprintf(" (%s)", top_regions_p$Gene_Name[i])
  } else {
    ""
  }
  summary_text <- c(summary_text,
    sprintf("%2d. %s%s | Diff=%.3f | FDR=%.2e | %s",
            i, top_regions_p$Region_ID[i], gene_info,
            top_regions_p$Difference[i], top_regions_p$adj_pvalue[i],
            top_regions_p$direction_label[i]))
}

# Gene enrichment summary (if gene names available)
if (any(!is.na(comp_df$Gene_Name) & comp_df$Gene_Name != "" &
        comp_df$Gene_Name != "-")) {

  # Count significant genes for standard threshold
  sig_std <- comp_df[(comp_df$adj_pvalue < 0.05) &
                      (abs(comp_df$Difference) > 0.30), ]
  sig_genes <- sig_std$Gene_Name[!is.na(sig_std$Gene_Name) &
                                  sig_std$Gene_Name != "" &
                                  sig_std$Gene_Name != "-"]
  unique_sig_genes <- unique(sig_genes)

  summary_text <- c(summary_text,
    "",
    "------------------------------------------------------------",
    "Gene Summary (Standard Thresholds):",
    "------------------------------------------------------------",
    sprintf("  Regions with gene annotations: %d",
            sum(!is.na(comp_df$Gene_Name) & comp_df$Gene_Name != "" &
                comp_df$Gene_Name != "-")),
    sprintf("  Significant regions with genes: %d", length(sig_genes)),
    sprintf("  Unique significant genes: %d", length(unique_sig_genes))
  )

  if (length(unique_sig_genes) > 0 && length(unique_sig_genes) <= 50) {
    summary_text <- c(summary_text,
      "",
      "Significant genes (alphabetical):",
      paste("  ", paste(sort(unique_sig_genes), collapse = ", "))
    )
  }
}

# Original: output_summary <- file.path(output_dir, "compartment_volcano_summary.txt")
output_summary <- file.path(tsv_dir, "1D_compartment_volcano_summary.txt")
writeLines(summary_text, output_summary)
cat(sprintf("Saved summary: %s\n", output_summary))

# Save full annotated dataset
# Original: output_all <- file.path(output_dir, "compartment_all_annotated.tsv")
output_all <- file.path(tsv_dir, "1D_compartment_all_annotated.tsv")
write.table(comp_df, output_all, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved full dataset: %s\n\n", output_all))

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================\n\n")

cat(sprintf("Output directory (TSVs): %s\n", tsv_dir))  # Original: output_dir
cat(sprintf("Output directory (plots): %s\n\n", plot_dir))

cat("Generated files:\n")
cat("  - compartment_volcano_relaxed.pdf (FDR<0.15, |Diff|>0.15)\n")
cat("  - compartment_volcano_standard.pdf (FDR<0.05, |Diff|>0.30)\n")
cat(sprintf("  - compartment_significant_relaxed.tsv (%d regions)\n",
            results_list$relaxed$n_sig_both))
cat(sprintf("  - compartment_significant_standard.tsv (%d regions)\n",
            results_list$standard$n_sig_both))
cat("  - compartment_volcano_summary.txt\n")
cat("  - compartment_all_annotated.tsv\n\n")

cat("Summary of compartment switching:\n")
cat(sprintf("  Relaxed thresholds:  %d regions (%d B->A, %d A->B)\n",
            results_list$relaxed$n_sig_both,
            results_list$relaxed$n_fc_up,
            results_list$relaxed$n_fc_down))
cat(sprintf("  Standard thresholds: %d regions (%d B->A, %d A->B)\n",
            results_list$standard$n_sig_both,
            results_list$standard$n_fc_up,
            results_list$standard$n_fc_down))

cat("\n================================================\n\n")
