#!/usr/bin/env Rscript
# scripts/tad_volcano_plot.R
# Volcano Plot for TAD Differential Expression Analysis
# Author: Zakir Alibhai
# Date: 2025-11-20
#
# Purpose:
#   Generate publication-quality volcano plots from TAD differential expression
#   analysis performed using getDiffExpression.pl (HOMER/Hi-C).
#
#   This script creates EnhancedVolcano plots showing TAD inclusion ratio
#   differences between control and mutant conditions, with significance
#   testing via adjusted p-values.
#
# Input Format:
#   Tab-delimited file from getDiffExpression.pl with columns:
#   - TAD name (chr:start-end)
#   - chr1, start1, end1, chr2, start2, end2
#   - InclusionRatio(IR) and replicate IR values
#   - "ctrl vs. mut Difference" (effect size, analogous to logFC)
#   - "ctrl vs. mut p-value" (raw p-value)
#   - "ctrl vs. mut adj. p-value" (FDR-adjusted, primary significance metric)
#
# Usage:
#   Rscript scripts/tad_volcano_plot.R [INPUT_FILE] [OPTIONS]
#
# Arguments:
#   INPUT_FILE              Path to TAD differential expression file (required)
#   --fdr THRESHOLD         FDR p-value cutoff (default: 0.15, lenient for TAD data)
#   --fc THRESHOLD          Difference cutoff (default: 0.15, moderate effect size)
#   --output DIR            Output directory (default: outputs/tad_analysis/)
#   --title TEXT            Custom plot title (default: auto-generated)
#   --width WIDTH           Plot width in inches (default: 10)
#   --height HEIGHT         Plot height in inches (default: 8)
#
# Examples:
#   # Basic usage with defaults
#   Rscript scripts/tad_volcano_plot.R Bap1.diff.tad.txt
#
#   # Custom thresholds
#   Rscript scripts/tad_volcano_plot.R Bap1.diff.tad.txt --fdr 0.1 --fc 0.5
#
#   # Custom output location
#   Rscript scripts/tad_volcano_plot.R Bap1.diff.tad.txt --output outputs/custom/
#
# Output:
#   - tad_volcano_plot.pdf - Publication-quality volcano plot
#   - tad_significant.tsv - Filtered significant TADs (FDR < threshold)
#   - tad_volcano_summary.txt - Summary statistics
#
# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("TAD Differential Expression: Volcano Plot\n")
cat("========================================\n\n")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default parameters
input_file <- NULL
fdr_threshold <- 0.15  # More lenient for TAD data (typical range: 0.1-0.4)
fc_threshold <- 0.15   # Moderate effect size for inclusion ratio differences
output_dir <- "outputs/tad_analysis"
custom_title <- NULL
plot_width <- 10
plot_height <- 8

# Parse arguments
if (length(args) == 0) {
  cat("ERROR: No input file specified\n\n")
  cat("Usage: Rscript scripts/tad_volcano_plot.R INPUT_FILE [OPTIONS]\n")
  cat("  --fdr THRESHOLD         FDR cutoff (default: 0.15)\n")
  cat("  --fc THRESHOLD          Difference cutoff (default: 0.15)\n")
  cat("  --output DIR            Output directory (default: outputs/tad_analysis/)\n")
  cat("  --title TEXT            Custom title\n")
  cat("  --width WIDTH           Plot width (default: 10)\n")
  cat("  --height HEIGHT         Plot height (default: 8)\n\n")
  quit(status = 1)
}

# First argument is always the input file
input_file <- args[1]

# Parse optional arguments
i <- 2
while (i <= length(args)) {
  if (args[i] == "--fdr" && i < length(args)) {
    fdr_threshold <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--fc" && i < length(args)) {
    fc_threshold <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--output" && i < length(args)) {
    output_dir <- args[i + 1]
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
cat(sprintf("  FDR threshold: %.3f\n", fdr_threshold))
cat(sprintf("  FC threshold: %.3f\n", fc_threshold))
cat(sprintf("  Output directory: %s\n", output_dir))
cat(sprintf("  Plot dimensions: %.1f × %.1f inches\n\n", plot_width, plot_height))

# Load required libraries
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(EnhancedVolcano)
})

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

cat("========================================\n")
cat("Loading TAD Differential Expression Data\n")
cat("========================================\n\n")

# Read TAD differential expression file
cat(sprintf("Reading: %s\n", input_file))
tad_data <- read.table(input_file, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE)

cat(sprintf("✓ Loaded %d TADs\n\n", nrow(tad_data)))

# Display column names for debugging
cat("Column names in input file:\n")
print(names(tad_data))
cat("\n")

# Identify key columns (handle variable naming from getDiffExpression.pl)
# The script produces columns with exact names we need to match
tad_name_col <- names(tad_data)[1]  # First column is TAD name
difference_col <- grep("ctrl vs\\. mut Difference", names(tad_data), value = TRUE, ignore.case = TRUE)[1]
pvalue_col <- grep("ctrl vs\\. mut p-value$", names(tad_data), value = TRUE, ignore.case = TRUE)[1]
adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value", names(tad_data), value = TRUE, ignore.case = TRUE)[1]

# Validate required columns exist
if (is.na(difference_col)) {
  stop("ERROR: Could not find 'ctrl vs. mut Difference' column")
}
if (is.na(adj_pvalue_col)) {
  stop("ERROR: Could not find 'ctrl vs. mut adj. p-value' column")
}

cat(sprintf("Identified key columns:\n"))
cat(sprintf("  TAD name: %s\n", tad_name_col))
cat(sprintf("  Difference (effect size): %s\n", difference_col))
cat(sprintf("  Adjusted p-value: %s\n", adj_pvalue_col))
cat("\n")

# Create clean working dataframe with standardized column names
tad_df <- data.frame(
  TAD_name = tad_data[[tad_name_col]],
  chr1 = tad_data$chr1,
  start1 = tad_data$start1,
  end1 = tad_data$end1,
  Difference = as.numeric(tad_data[[difference_col]]),
  pvalue = if (!is.na(pvalue_col)) as.numeric(tad_data[[pvalue_col]]) else NA,
  adj_pvalue = as.numeric(tad_data[[adj_pvalue_col]]),
  stringsAsFactors = FALSE
)

# Calculate TAD size
tad_df$TAD_size <- tad_df$end1 - tad_df$start1

# Remove rows with NA values in critical columns
n_before <- nrow(tad_df)
tad_df <- tad_df[!is.na(tad_df$Difference) & !is.na(tad_df$adj_pvalue), ]
n_after <- nrow(tad_df)

if (n_before > n_after) {
  cat(sprintf("⚠ Removed %d TADs with missing values\n", n_before - n_after))
}

# Handle zero or very small p-values (set floor to prevent Inf in -log10)
min_pval <- 1e-300
tad_df$adj_pvalue[tad_df$adj_pvalue == 0] <- min_pval
tad_df$adj_pvalue[tad_df$adj_pvalue < min_pval] <- min_pval

# Create significance flag
tad_df$significant <- (tad_df$adj_pvalue < fdr_threshold) & (abs(tad_df$Difference) > fc_threshold)

# Classify direction
tad_df$direction <- ifelse(tad_df$Difference > 0, "Increased_in_Mutant", "Decreased_in_Mutant")

# Calculate summary statistics
n_total <- nrow(tad_df)

# Counts for different categories (matching EnhancedVolcano logic)
# Category 1: Both FDR and FC significant (dark red)
n_sig_both <- sum(tad_df$significant)

# Category 2: FC significant (includes both grey and dark red points)
n_fc_up <- sum(tad_df$Difference > fc_threshold)  # Right side (increased in mutant)
n_fc_down <- sum(tad_df$Difference < -fc_threshold)  # Left side (decreased in mutant)

# Category 3: FDR significant (includes both red and dark red points)
n_fdr_sig <- sum(tad_df$adj_pvalue < fdr_threshold)

# For plot annotations, show counts of TADs exceeding FC threshold (visible colored points on each side)
n_up <- n_fc_up
n_down <- n_fc_down

# Overall significance (at least one criterion)
n_sig <- sum((tad_df$adj_pvalue < fdr_threshold) | (abs(tad_df$Difference) > fc_threshold))
n_ns <- n_total - n_sig

cat("\n========================================\n")
cat("Summary Statistics\n")
cat("========================================\n\n")
cat(sprintf("Total TADs analyzed: %d\n", n_total))
cat(sprintf("\nSignificance categories:\n"))
cat(sprintf("  Both criteria (FDR < %.3f AND |Diff| > %.3f): %d\n",
            fdr_threshold, fc_threshold, n_sig_both))
cat(sprintf("  FDR significant only: %d\n", n_fdr_sig))
cat(sprintf("  FC significant only: %d\n", n_fc_up + n_fc_down))
cat(sprintf("    - Increased in Mutant (Diff > %.3f): %d\n", fc_threshold, n_fc_up))
cat(sprintf("    - Decreased in Mutant (Diff < -%.3f): %d\n", fc_threshold, n_fc_down))
cat(sprintf("  At least one criterion: %d (%.1f%%)\n", n_sig, 100 * n_sig / n_total))
cat(sprintf("  Non-significant: %d (%.1f%%)\n\n", n_ns, 100 * n_ns / n_total))

# Data range info
cat("Data ranges:\n")
cat(sprintf("  Difference: [%.3f, %.3f]\n", min(tad_df$Difference), max(tad_df$Difference)))
cat(sprintf("  Adj. p-value: [%.2e, %.2f]\n", min(tad_df$adj_pvalue), max(tad_df$adj_pvalue)))
cat(sprintf("  -log10(adj. p-value): [%.2f, %.2f]\n\n",
            min(-log10(tad_df$adj_pvalue)), max(-log10(tad_df$adj_pvalue))))

# =============================================================================
# CREATE VOLCANO PLOT
# =============================================================================

cat("========================================\n")
cat("Generating Volcano Plot\n")
cat("========================================\n\n")

# Determine plot title
if (is.null(custom_title)) {
  plot_title <- "TAD Differential Expression: Control vs Mutant"
} else {
  plot_title <- custom_title
}

# Get data ranges for annotation positioning
diff_range <- range(tad_df$Difference, na.rm = TRUE)
fdr_range <- range(tad_df$adj_pvalue[tad_df$adj_pvalue > 0], na.rm = TRUE)
neg_log10_fdr_max <- -log10(min(fdr_range))

cat("Creating EnhancedVolcano plot...\n")

# Calculate appropriate y-axis limits based on actual data
y_max <- -log10(min(tad_df$adj_pvalue))
y_lim <- c(0, max(1.5, y_max * 1.15))  # At least 1.5, or 15% above max

# Create EnhancedVolcano with publication settings matching pipeline style
p <- EnhancedVolcano(
  tad_df,
  lab = tad_df$TAD_name,
  x = 'Difference',
  y = 'adj_pvalue',
  title = plot_title,
  subtitle = 'TAD Inclusion Ratio Analysis',
  pCutoff = fdr_threshold,
  FCcutoff = fc_threshold,
  pointSize = 2.5,
  labSize = 0,  # Hide individual labels for cleaner look
  col = c('black', 'grey50', 'red', 'darkred'),  # NS, FC-only, p-only, both
  colAlpha = 0.6,
  legendPosition = 'top',
  legendLabSize = 11,
  legendIconSize = 4.0,
  legendLabels = c('NS', 'Difference', 'FDR', 'FDR and Difference'),
  drawConnectors = FALSE,
  gridlines.major = TRUE,
  gridlines.minor = FALSE,
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black',
  xlab = bquote(bold('Ctrl vs Mut Difference')),
  ylab = bquote(bold('-log'[10]*'(Adjusted p-value)')),
  xlim = c(diff_range[1] - 0.2, diff_range[2] + 0.2),
  ylim = y_lim,  # Use calculated y-axis limits
  caption = NULL
)

# Add custom text annotations matching pipeline style
p <- p +
  # Increased in mutant count (top right, blue)
  annotate(
    "text",
    x = diff_range[2] - 0.1,
    y = neg_log10_fdr_max * 0.98,
    label = as.character(n_up),
    color = "#3366CC",
    size = 6,
    fontface = "bold",
    hjust = 1
  ) +
  # Decreased in mutant count (top left, blue)
  annotate(
    "text",
    x = diff_range[1] + 0.1,
    y = neg_log10_fdr_max * 0.98,
    label = as.character(n_down),
    color = "#3366CC",
    size = 6,
    fontface = "bold",
    hjust = 0
  ) +
  # Total TADs (bottom right, black)
  annotate(
    "text",
    x = diff_range[2] - 0.05,
    y = 0.15,
    label = sprintf("total = %d TADs", n_total),
    color = "black",
    size = 3.5,
    hjust = 1
  ) +
  # Enhanced theme matching pipeline
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
output_plot <- file.path(output_dir, "tad_volcano_plot.pdf")
ggsave(output_plot, p, width = plot_width, height = plot_height, dpi = 300)
cat(sprintf("✓ Saved: %s\n\n", output_plot))

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("========================================\n")
cat("Saving Results\n")
cat("========================================\n\n")

# Save significant TADs
sig_tads <- tad_df[tad_df$significant, ]
sig_tads <- sig_tads[order(sig_tads$adj_pvalue), ]  # Sort by significance

output_sig <- file.path(output_dir, "tad_significant.tsv")
write.table(sig_tads, output_sig, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("✓ Saved significant TADs: %s\n", output_sig))
cat(sprintf("  %d TADs (FDR < %.3f, |Diff| > %.3f)\n\n", nrow(sig_tads), fdr_threshold, fc_threshold))

# Save summary statistics
summary_text <- c(
  "TAD Differential Expression Volcano Plot - Summary Statistics",
  "==============================================================",
  "",
  sprintf("Analysis Date: %s", Sys.Date()),
  sprintf("Input File: %s", input_file),
  sprintf("Output Directory: %s", output_dir),
  "",
  "Thresholds:",
  sprintf("  FDR cutoff: %.3f", fdr_threshold),
  sprintf("  Difference cutoff: %.3f", fc_threshold),
  "",
  "Results:",
  sprintf("  Total TADs analyzed: %d", n_total),
  "",
  "Significance categories:",
  sprintf("  Both criteria (FDR < %.3f AND |Diff| > %.3f): %d", fdr_threshold, fc_threshold, n_sig_both),
  sprintf("  FDR significant only: %d", n_fdr_sig),
  sprintf("  FC significant (|Diff| > %.3f): %d", fc_threshold, n_fc_up + n_fc_down),
  sprintf("    - Increased in Mutant (Diff > %.3f): %d", fc_threshold, n_fc_up),
  sprintf("    - Decreased in Mutant (Diff < -%.3f): %d", fc_threshold, n_fc_down),
  sprintf("  At least one criterion: %d (%.1f%%)", n_sig, 100 * n_sig / n_total),
  sprintf("  Non-significant: %d (%.1f%%)", n_ns, 100 * n_ns / n_total),
  "",
  "Data Ranges:",
  sprintf("  Difference: [%.3f, %.3f]", min(tad_df$Difference), max(tad_df$Difference)),
  sprintf("  Adj. p-value: [%.2e, %.2f]", min(tad_df$adj_pvalue), max(tad_df$adj_pvalue)),
  "",
  "Top 10 TADs by Effect Size (largest |Difference|):",
  "---------------------------------------------------"
)

# Add top 10 TADs by absolute difference
tad_df_sorted <- tad_df[order(-abs(tad_df$Difference)), ]
top_tads_fc <- head(tad_df_sorted[, c("TAD_name", "Difference", "adj_pvalue", "direction")], 10)
for (i in 1:nrow(top_tads_fc)) {
  summary_text <- c(summary_text,
                   sprintf("%2d. %s | Diff=%.3f | FDR=%.2e | %s",
                          i, top_tads_fc$TAD_name[i], top_tads_fc$Difference[i],
                          top_tads_fc$adj_pvalue[i], top_tads_fc$direction[i]))
}

# Add top 10 by p-value if any are significant
if (n_fdr_sig > 0) {
  summary_text <- c(summary_text,
                   "",
                   "Top 10 TADs by FDR (smallest adj. p-value):",
                   "--------------------------------------------")
  tad_df_sorted_p <- tad_df[order(tad_df$adj_pvalue), ]
  top_tads_p <- head(tad_df_sorted_p[, c("TAD_name", "Difference", "adj_pvalue", "direction")], 10)
  for (i in 1:nrow(top_tads_p)) {
    summary_text <- c(summary_text,
                     sprintf("%2d. %s | Diff=%.3f | FDR=%.2e | %s",
                            i, top_tads_p$TAD_name[i], top_tads_p$Difference[i],
                            top_tads_p$adj_pvalue[i], top_tads_p$direction[i]))
  }
}

output_summary <- file.path(output_dir, "tad_volcano_summary.txt")
writeLines(summary_text, output_summary)
cat(sprintf("✓ Saved summary: %s\n\n", output_summary))

# Save full annotated dataset
output_all <- file.path(output_dir, "tad_all_annotated.tsv")
write.table(tad_df, output_all, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("✓ Saved full dataset: %s\n\n", output_all))

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat(sprintf("Output directory: %s\n\n", output_dir))

cat("Generated files:\n")
cat(sprintf("  - tad_volcano_plot.pdf (%d × %d inches)\n", plot_width, plot_height))
cat(sprintf("  - tad_significant.tsv (%d TADs)\n", nrow(sig_tads)))
cat("  - tad_volcano_summary.txt\n")
cat("  - tad_all_annotated.tsv\n\n")

if (n_sig_both == 0 && n_fdr_sig == 0 && (n_fc_up + n_fc_down) == 0) {
  cat("⚠ WARNING: No TADs meet either threshold\n")
  cat("  Consider relaxing thresholds with --fdr or --fc options\n\n")
} else {
  cat(sprintf("✓ Found TADs meeting criteria:\n"))
  cat(sprintf("  Both criteria (FDR & FC): %d\n", n_sig_both))
  cat(sprintf("  FDR significant: %d\n", n_fdr_sig))
  cat(sprintf("  FC significant: %d (%d increased, %d decreased in mutant)\n\n",
              n_fc_up + n_fc_down, n_fc_up, n_fc_down))
}

cat("========================================\n\n")
