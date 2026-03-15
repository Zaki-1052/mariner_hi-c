#!/usr/bin/env Rscript
# scripts/tad_volcano_plot.R
# Volcano Plot for TAD Differential Expression Analysis
# Author: Zakir Alibhai
# Date: 2025-11-20 (Updated: 2025-11-24)
#
# Purpose:
#   Generate publication-quality volcano plots from TAD differential analysis.
#   Supports two input formats:
#     1. HOMER getDiffExpression.pl output (inclusion ratio differences)
#     2. TADCompare output (Gap_Score Z-scores, p-values computed internally)
#
#   Format is auto-detected from column names. Generates TWO threshold versions:
#   - HOMER:      Relaxed (FDR<0.15, |Diff|>0.15), Standard (FDR<0.05, |Diff|>0.30)
#   - TADCompare: Relaxed (FDR<0.15, |Z|>2.0),     Standard (FDR<0.05, |Z|>3.0)
#
# Input Formats:
#   HOMER (getDiffExpression.pl):
#   - TAD name, chr1, start1, end1, "ctrl vs. mut Difference", p-value, adj. p-value
#
#   TADCompare (tadcompare_final_filtered.tsv):
#   - chr, Boundary, Gap_Score, TAD_Score1, TAD_Score2, Differential, Enriched_In, Type
#   - P-values are derived from Gap_Score as Z-scores: 2*pnorm(-|Z|), then BH FDR
#   - Gap_Score sign convention depends on matrix order in TADCompare call
#
# Usage:
#   Rscript scripts/tad_volcano_plot.R [INPUT_FILE] [OPTIONS]
#
# Arguments:
#   INPUT_FILE              Path to TAD differential expression file
#                           (default: tad_analysis/Bap1.diff.tad.txt)
#   --output DIR            Output directory (default: outputs/tad_analysis/)
#   --title TEXT            Custom plot title (default: auto-generated)
#   --width WIDTH           Plot width in inches (default: 10)
#   --height HEIGHT         Plot height in inches (default: 8)
#
# Examples:
#   # Basic usage with defaults (generates both threshold versions)
#   Rscript scripts/tad_volcano_plot.R
#
#   # Custom input file
#   Rscript scripts/tad_volcano_plot.R tad_analysis/Bap1.diff.tad.txt
#
#   # Custom output location
#   Rscript scripts/tad_volcano_plot.R --output outputs/custom/
#
# Output:
#   - tad_volcano_relaxed.pdf - Volcano plot with relaxed thresholds
#   - tad_volcano_standard.pdf - Volcano plot with standard thresholds
#   - tad_significant_relaxed.tsv - Significant TADs (FDR < 0.15, |Diff| > 0.15)
#   - tad_significant_standard.tsv - Significant TADs (FDR < 0.05, |Diff| > 0.30)
#   - tad_volcano_summary.txt - Summary statistics for both thresholds
#   - tad_all_annotated.tsv - Full dataset with annotations
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
input_file <- "tad_analysis/Bap1.diff.tad.txt"
output_dir <- "outputs/tad_analysis"
custom_title <- NULL
plot_width <- 10
plot_height <- 8

# Thresholds and labels are set after format detection (see below)

# Parse arguments
i <- 1
while (i <= length(args)) {
  if (args[i] == "--output" && i < length(args)) {
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
cat(sprintf("  Output directory: %s\n", output_dir))
cat(sprintf("  Plot dimensions: %.1f x %.1f inches\n\n", plot_width, plot_height))

# Load required libraries
cat("Loading required packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(EnhancedVolcano)
})

# Load multi-format output utility for PDF + SVG + JPEG output
source("scripts/utils/multi_format_output.R")

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

cat(sprintf("Loaded %d TADs\n\n", nrow(tad_data)))

# Display column names for debugging
cat("Column names in input file:\n")
print(names(tad_data))
cat("\n")

# Auto-detect input format from column names
col_names <- names(tad_data)
is_homer <- any(grepl("ctrl vs\\. mut Difference", col_names, ignore.case = TRUE))
is_tadcompare <- all(c("Gap_Score", "Boundary", "Differential") %in% col_names)

if (!is_homer && !is_tadcompare) {
  stop(sprintf("Unrecognized input format. Expected HOMER or TADCompare columns.\nFound: %s",
               paste(col_names, collapse = ", ")))
}

input_format <- if (is_homer) "homer" else "tadcompare"
cat(sprintf("Detected input format: %s\n\n", toupper(input_format)))

if (input_format == "homer") {
  # === HOMER getDiffExpression.pl format ===
  tad_name_col <- col_names[1]
  difference_col <- grep("ctrl vs\\. mut Difference", col_names, value = TRUE, ignore.case = TRUE)[1]
  pvalue_col <- grep("ctrl vs\\. mut p-value$", col_names, value = TRUE, ignore.case = TRUE)[1]
  adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value", col_names, value = TRUE, ignore.case = TRUE)[1]

  stopifnot("Could not find 'ctrl vs. mut Difference' column" = !is.na(difference_col))
  stopifnot("Could not find 'ctrl vs. mut adj. p-value' column" = !is.na(adj_pvalue_col))

  cat(sprintf("Identified key columns:\n"))
  cat(sprintf("  TAD name: %s\n", tad_name_col))
  cat(sprintf("  Difference (effect size): %s\n", difference_col))
  cat(sprintf("  Adjusted p-value: %s\n\n", adj_pvalue_col))

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
  tad_df$TAD_size <- tad_df$end1 - tad_df$start1

  # HOMER thresholds (inclusion ratio difference scale)
  thresholds <- list(
    relaxed = list(fdr = 0.15, fc = 0.15, name = "relaxed"),
    standard = list(fdr = 0.05, fc = 0.30, name = "standard")
  )
  x_axis_label <- bquote(bold('Ctrl vs Mut Difference'))
  analysis_type <- "TAD Inclusion Ratio Analysis"
  effect_label <- "Difference"

} else {
  # === TADCompare format ===
  # Gap_Score is a Z-score; compute p-values and FDR from it
  stopifnot("TADCompare file missing 'chr' column" = "chr" %in% col_names)

  cat("Computing p-values from Gap_Score Z-scores (BH FDR correction)...\n")
  pvals <- 2 * pnorm(-abs(tad_data$Gap_Score))
  adj_pvals <- p.adjust(pvals, method = "BH")
  cat(sprintf("  FDR range: [%.2e, %.4f]\n\n", min(adj_pvals), max(adj_pvals)))

  # Gap_Score used directly as effect size (Difference)
  # Sign convention: depends on matrix order in TADCompare call.
  # For early timepoint (matrices swapped): positive = mutant-enriched.
  # For standard matrix order: positive = control-enriched (negate with --negate if needed).
  tad_df <- data.frame(
    TAD_name = paste0(tad_data$chr, ":", tad_data$Boundary),
    chr1 = tad_data$chr,
    start1 = tad_data$Boundary,
    end1 = tad_data$Boundary,
    Difference = tad_data$Gap_Score,
    pvalue = pvals,
    adj_pvalue = adj_pvals,
    stringsAsFactors = FALSE
  )
  tad_df$TAD_size <- 0

  # TADCompare thresholds (Z-score scale)
  thresholds <- list(
    relaxed = list(fdr = 0.15, fc = 2.0, name = "relaxed"),
    standard = list(fdr = 0.05, fc = 3.0, name = "standard")
  )
  x_axis_label <- bquote(bold('Boundary Gap Score'))
  analysis_type <- "TADCompare Boundary Analysis"
  effect_label <- "Gap Score"
}

cat(sprintf("Thresholds: relaxed (FDR<%.2f, |%s|>%.2f), standard (FDR<%.2f, |%s|>%.2f)\n",
            thresholds$relaxed$fdr, effect_label, thresholds$relaxed$fc,
            thresholds$standard$fdr, effect_label, thresholds$standard$fc))

# Remove rows with NA values in critical columns
n_before <- nrow(tad_df)
tad_df <- tad_df[!is.na(tad_df$Difference) & !is.na(tad_df$adj_pvalue), ]
n_after <- nrow(tad_df)

if (n_before > n_after) {
  cat(sprintf("Removed %d TADs with missing values\n", n_before - n_after))
}

# Handle zero or very small p-values (set floor to prevent Inf in -log10)
min_pval <- 1e-300
tad_df$adj_pvalue[tad_df$adj_pvalue == 0] <- min_pval
tad_df$adj_pvalue[tad_df$adj_pvalue < min_pval] <- min_pval

# Classify direction
tad_df$direction <- ifelse(tad_df$Difference > 0, "Increased_in_Mutant", "Decreased_in_Mutant")

n_total <- nrow(tad_df)

# Data range info
cat("\nData ranges:\n")
cat(sprintf("  Difference: [%.3f, %.3f]\n", min(tad_df$Difference), max(tad_df$Difference)))
cat(sprintf("  Adj. p-value: [%.2e, %.2f]\n", min(tad_df$adj_pvalue), max(tad_df$adj_pvalue)))
cat(sprintf("  -log10(adj. p-value): [%.2f, %.2f]\n\n",
            min(-log10(tad_df$adj_pvalue)), max(-log10(tad_df$adj_pvalue))))

# =============================================================================
# VOLCANO PLOT GENERATION FUNCTION
# =============================================================================

generate_volcano_plot <- function(df, fdr_threshold, fc_threshold, threshold_name,
                                  plot_title_base, output_dir, plot_width, plot_height,
                                  x_label = bquote(bold('Ctrl vs Mut Difference')),
                                  analysis_type = "TAD Inclusion Ratio Analysis",
                                  effect_label = "Difference") {

  cat(sprintf("\n--- Generating %s threshold plot (FDR < %.2f, |Diff| > %.2f) ---\n",
              threshold_name, fdr_threshold, fc_threshold))

  # Calculate statistics for this threshold
  n_total <- nrow(df)
  n_sig_both <- sum((df$adj_pvalue < fdr_threshold) & (abs(df$Difference) > fc_threshold))
  n_fc_up <- sum(df$Difference > fc_threshold)
  n_fc_down <- sum(df$Difference < -fc_threshold)
  n_fdr_sig <- sum(df$adj_pvalue < fdr_threshold)
  n_sig <- sum((df$adj_pvalue < fdr_threshold) | (abs(df$Difference) > fc_threshold))
  n_ns <- n_total - n_sig

  cat(sprintf("  Total TADs: %d\n", n_total))
  cat(sprintf("  Both criteria (FDR & FC): %d\n", n_sig_both))
  cat(sprintf("  FC significant: %d (%d up, %d down)\n", n_fc_up + n_fc_down, n_fc_up, n_fc_down))
  cat(sprintf("  FDR significant: %d\n", n_fdr_sig))

  # Get data ranges for annotation positioning
  diff_range <- range(df$Difference, na.rm = TRUE)
  neg_log10_fdr_max <- -log10(min(df$adj_pvalue[df$adj_pvalue > 0], na.rm = TRUE))

  # Calculate appropriate y-axis limits
  y_max <- -log10(min(df$adj_pvalue))
  y_lim <- c(0, max(1.5, y_max * 1.15))

  # Create plot title
  if (is.null(plot_title_base)) {
    plot_title <- sprintf("TAD Differential Expression: Control vs Mutant (%s)", threshold_name)
  } else {
    plot_title <- sprintf("%s (%s)", plot_title_base, threshold_name)
  }

  subtitle <- sprintf("%s | FDR < %.2f, |%s| > %.2f",
                      analysis_type, fdr_threshold, effect_label, fc_threshold)

  # Create EnhancedVolcano
  p <- EnhancedVolcano(
    df,
    lab = df$TAD_name,
    x = 'Difference',
    y = 'adj_pvalue',
    title = plot_title,
    subtitle = subtitle,
    pCutoff = fdr_threshold,
    FCcutoff = fc_threshold,
    pointSize = 2.5,
    labSize = 0,
    col = c('black', 'grey50', 'red', 'darkred'),
    colAlpha = 0.6,
    legendPosition = 'top',
    legendLabSize = 11,
    legendIconSize = 4.0,
    legendLabels = c('NS', effect_label, 'FDR', paste('FDR and', effect_label)),
    drawConnectors = FALSE,
    gridlines.major = TRUE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 0.8,
    borderColour = 'black',
    xlab = x_label,
    ylab = bquote(bold('-log'[10]*'(Adjusted p-value)')),
    xlim = c(diff_range[1] - 0.2, diff_range[2] + 0.2),
    ylim = y_lim,
    caption = NULL
  )

  # Add annotations
  p <- p +
    annotate("text", x = diff_range[2] - 0.1, y = neg_log10_fdr_max * 0.98,
             label = as.character(n_fc_up), color = "#3366CC", size = 6,
             fontface = "bold", hjust = 1) +
    annotate("text", x = diff_range[1] + 0.1, y = neg_log10_fdr_max * 0.98,
             label = as.character(n_fc_down), color = "#3366CC", size = 6,
             fontface = "bold", hjust = 0) +
    annotate("text", x = diff_range[2] - 0.05, y = 0.15,
             label = sprintf("total = %d TADs", n_total), color = "black",
             size = 3.5, hjust = 1) +
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

  # Save plot in multiple formats (PDF, SVG, JPEG)
  output_base <- file.path(output_dir, sprintf("tad_volcano_%s", threshold_name))
  save_multiformat_ggplot(p, output_base, width = plot_width, height = plot_height)

  # Save significant TADs for this threshold
  df$significant <- (df$adj_pvalue < fdr_threshold) & (abs(df$Difference) > fc_threshold)
  sig_tads <- df[df$significant, ]
  sig_tads <- sig_tads[order(sig_tads$adj_pvalue), ]

  output_sig <- file.path(output_dir, sprintf("tad_significant_%s.tsv", threshold_name))
  write.table(sig_tads, output_sig, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved: %s (%d TADs)\n", output_sig, nrow(sig_tads)))

  # Return statistics for summary
  return(list(
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
  ))
}

# =============================================================================
# GENERATE VOLCANO PLOTS FOR BOTH THRESHOLD SETS
# =============================================================================

cat("========================================\n")
cat("Generating Volcano Plots\n")
cat("========================================\n")

results_list <- list()

for (thresh in thresholds) {
  result <- generate_volcano_plot(
    df = tad_df,
    fdr_threshold = thresh$fdr,
    fc_threshold = thresh$fc,
    threshold_name = thresh$name,
    plot_title_base = custom_title,
    output_dir = output_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    x_label = x_axis_label,
    analysis_type = analysis_type,
    effect_label = effect_label
  )
  results_list[[thresh$name]] <- result
}

# =============================================================================
# SAVE COMBINED SUMMARY
# =============================================================================

cat("\n========================================\n")
cat("Saving Summary\n")
cat("========================================\n\n")

# Build summary text
summary_text <- c(
  "TAD Differential Expression Volcano Plot - Summary Statistics",
  "==============================================================",
  "",
  sprintf("Analysis Date: %s", Sys.Date()),
  sprintf("Input File: %s", input_file),
  sprintf("Input Format: %s", toupper(input_format)),
  sprintf("Output Directory: %s", output_dir),
  "",
  sprintf("Total boundaries analyzed: %d", n_total),
  "",
  "Data Ranges:",
  sprintf("  %s: [%.3f, %.3f]", effect_label, min(tad_df$Difference), max(tad_df$Difference)),
  sprintf("  Adj. p-value: [%.2e, %.2f]", min(tad_df$adj_pvalue), max(tad_df$adj_pvalue)),
  ""
)

# Add results for each threshold
for (result in results_list) {
  summary_text <- c(summary_text,
    "--------------------------------------------------------------",
    sprintf("%s Thresholds (FDR < %.2f, |%s| > %.2f)",
            toupper(result$threshold_name), result$fdr, effect_label, result$fc),
    "--------------------------------------------------------------",
    sprintf("  Both criteria (FDR & %s): %d", effect_label, result$n_sig_both),
    sprintf("  FDR significant only: %d", result$n_fdr_sig),
    sprintf("  %s significant: %d", effect_label, result$n_fc_up + result$n_fc_down),
    sprintf("    - Increased in Mutant (%s > %.2f): %d", effect_label, result$fc, result$n_fc_up),
    sprintf("    - Decreased in Mutant (%s < -%.2f): %d", effect_label, result$fc, result$n_fc_down),
    sprintf("  At least one criterion: %d (%.1f%%)", result$n_sig, 100 * result$n_sig / result$n_total),
    sprintf("  Non-significant: %d (%.1f%%)", result$n_ns, 100 * result$n_ns / result$n_total),
    ""
  )
}

# Add top 10 TADs by effect size
summary_text <- c(summary_text,
  "--------------------------------------------------------------",
  sprintf("Top 10 Boundaries by Effect Size (largest |%s|):", effect_label),
  "--------------------------------------------------------------"
)

tad_df_sorted <- tad_df[order(-abs(tad_df$Difference)), ]
top_tads_fc <- head(tad_df_sorted[, c("TAD_name", "Difference", "adj_pvalue", "direction")], 10)
for (i in 1:nrow(top_tads_fc)) {
  summary_text <- c(summary_text,
    sprintf("%2d. %s | Diff=%.3f | FDR=%.2e | %s",
            i, top_tads_fc$TAD_name[i], top_tads_fc$Difference[i],
            top_tads_fc$adj_pvalue[i], top_tads_fc$direction[i]))
}

# Add top 10 by p-value
summary_text <- c(summary_text,
  "",
  "Top 10 TADs by FDR (smallest adj. p-value):",
  "--------------------------------------------"
)
tad_df_sorted_p <- tad_df[order(tad_df$adj_pvalue), ]
top_tads_p <- head(tad_df_sorted_p[, c("TAD_name", "Difference", "adj_pvalue", "direction")], 10)
for (i in 1:nrow(top_tads_p)) {
  summary_text <- c(summary_text,
    sprintf("%2d. %s | Diff=%.3f | FDR=%.2e | %s",
            i, top_tads_p$TAD_name[i], top_tads_p$Difference[i],
            top_tads_p$adj_pvalue[i], top_tads_p$direction[i]))
}

output_summary <- file.path(output_dir, "tad_volcano_summary.txt")
writeLines(summary_text, output_summary)
cat(sprintf("Saved summary: %s\n", output_summary))

# Save full annotated dataset
output_all <- file.path(output_dir, "tad_all_annotated.tsv")
write.table(tad_df, output_all, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Saved full dataset: %s\n\n", output_all))

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat(sprintf("Output directory: %s\n\n", output_dir))

cat("Generated files:\n")
cat(sprintf("  - tad_volcano_relaxed.pdf (FDR<%.2f, |%s|>%.2f)\n",
            thresholds$relaxed$fdr, effect_label, thresholds$relaxed$fc))
cat(sprintf("  - tad_volcano_standard.pdf (FDR<%.2f, |%s|>%.2f)\n",
            thresholds$standard$fdr, effect_label, thresholds$standard$fc))
cat(sprintf("  - tad_significant_relaxed.tsv (%d TADs)\n", results_list$relaxed$n_sig_both))
cat(sprintf("  - tad_significant_standard.tsv (%d TADs)\n", results_list$standard$n_sig_both))
cat("  - tad_volcano_summary.txt\n")
cat("  - tad_all_annotated.tsv\n\n")

cat("Summary of significant TADs:\n")
cat(sprintf("  Relaxed thresholds:  %d TADs (both criteria)\n", results_list$relaxed$n_sig_both))
cat(sprintf("  Standard thresholds: %d TADs (both criteria)\n", results_list$standard$n_sig_both))

cat("\n========================================\n\n")
