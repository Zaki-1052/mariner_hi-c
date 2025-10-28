#!/usr/bin/env Rscript
# scripts/compare_resolutions.R
# Meta-analysis: Compare differential loop results across resolutions
# Author: Zakir
# Date: 2025-10-20

# =============================================================================
# SETUP
# =============================================================================

cat("\n========================================\n")
cat("Multi-Resolution Comparison Analysis\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(InteractionSet)
  library(VennDiagram)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)
})

# Set working directory
setwd("/expanse/lustre/projects/csd940/zalibhai/mariner")

# Create output directory
output_dir <- "outputs/resolution_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SECTION 1: LOAD RESULTS FROM ALL RESOLUTIONS
# =============================================================================

cat("SECTION 1: Loading Results\n")
cat("---------------------------\n")

resolutions <- c(5000, 10000, 25000)
results_list <- list()
coords_list <- list()
final_results_list <- list()

for (res in resolutions) {
  res_kb <- res / 1000
  cat(sprintf("Loading %dkb results... ", res_kb))

  # Load differential results
  results_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.rds", res_kb)
  coords_file <- sprintf("outputs/res_%dkb/03_binned.rds", res_kb)
  final_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/final_results.tsv", res_kb)

  if (!file.exists(results_file) || !file.exists(coords_file)) {
    cat(sprintf("✗ Missing files for %dkb\n", res_kb))
    next
  }

  results_list[[as.character(res)]] <- readRDS(results_file)
  coords_list[[as.character(res)]] <- readRDS(coords_file)

  # Load final results (stringent thresholds: |logFC| > 0.3, FDR < 0.03)
  if (file.exists(final_file)) {
    final_results_list[[as.character(res)]] <- read.table(final_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    cat(sprintf("✓ %d loops (+ %d final)\n",
                nrow(results_list[[as.character(res)]]),
                nrow(final_results_list[[as.character(res)]])))
  } else {
    cat(sprintf("✓ %d loops (final results not found)\n", nrow(results_list[[as.character(res)]])))
  }
}

if (length(results_list) == 0) {
  stop("ERROR: No resolution results found. Run multi-resolution pipeline first.")
}

cat(sprintf("\n✓ Loaded results from %d resolutions\n", length(results_list)))
if (length(final_results_list) > 0) {
  cat(sprintf("✓ Loaded final results (stringent) from %d resolutions\n\n", length(final_results_list)))
} else {
  cat("  Note: No final_results.tsv files found. Run generate_final_results.py to create them.\n\n")
}

# =============================================================================
# SECTION 2: BASIC STATISTICS PER RESOLUTION
# =============================================================================

cat("SECTION 2: Basic Statistics\n")
cat("----------------------------\n\n")

summary_stats <- data.frame(
  resolution_kb = numeric(),
  total_loops = numeric(),
  significant_fdr05 = numeric(),
  pct_significant = numeric(),
  up_in_mutant = numeric(),
  down_in_mutant = numeric(),
  median_logFC_all = numeric(),
  median_logFC_sig = numeric(),
  final_loops = numeric(),
  final_up = numeric(),
  final_down = numeric(),
  pct_final = numeric(),
  median_logFC_final = numeric(),
  stringsAsFactors = FALSE
)

for (res in names(results_list)) {
  res_kb <- as.numeric(res) / 1000
  df <- results_list[[res]]

  n_sig <- sum(df$significant, na.rm = TRUE)
  n_up <- sum(df$significant & df$logFC > 0, na.rm = TRUE)
  n_down <- sum(df$significant & df$logFC < 0, na.rm = TRUE)

  # Calculate final results stats (stringent thresholds)
  if (res %in% names(final_results_list)) {
    final_df <- final_results_list[[res]]
    n_final <- nrow(final_df)
    n_final_up <- sum(final_df$logFC > 0, na.rm = TRUE)
    n_final_down <- sum(final_df$logFC < 0, na.rm = TRUE)
    median_fc_final <- median(abs(final_df$logFC), na.rm = TRUE)
  } else {
    n_final <- NA
    n_final_up <- NA
    n_final_down <- NA
    median_fc_final <- NA
  }

  summary_stats <- rbind(summary_stats, data.frame(
    resolution_kb = res_kb,
    total_loops = nrow(df),
    significant_fdr05 = n_sig,
    pct_significant = 100 * n_sig / nrow(df),
    up_in_mutant = n_up,
    down_in_mutant = n_down,
    median_logFC_all = median(abs(df$logFC), na.rm = TRUE),
    median_logFC_sig = median(abs(df$logFC[df$significant]), na.rm = TRUE),
    final_loops = n_final,
    final_up = n_final_up,
    final_down = n_final_down,
    pct_final = if (!is.na(n_final)) 100 * n_final / nrow(df) else NA,
    median_logFC_final = median_fc_final
  ))
}

cat("Standard thresholds (FDR < 0.05):\n")
print(summary_stats[, c("resolution_kb", "total_loops", "significant_fdr05", "pct_significant",
                        "up_in_mutant", "down_in_mutant")])

if (any(!is.na(summary_stats$final_loops))) {
  cat("\nStringent thresholds (|logFC| > 0.3, FDR < 0.03):\n")
  print(summary_stats[, c("resolution_kb", "final_loops", "pct_final",
                          "final_up", "final_down", "median_logFC_final")])
}

# Save summary
write.table(summary_stats,
            file.path(output_dir, "summary_statistics_by_resolution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n✓ Summary statistics saved\n\n")

# =============================================================================
# SECTION 3: COORDINATE-BASED OVERLAP ANALYSIS
# =============================================================================

cat("SECTION 3: Loop Overlap Analysis\n")
cat("---------------------------------\n\n")

# Function to match loops across resolutions
# Uses a tolerance based on the coarser resolution
match_loops <- function(coords1, coords2, tolerance_kb = 10) {
  # Convert to data.frame with coordinates
  df1 <- data.frame(
    chr1 = as.character(seqnames(anchors(coords1, "first"))),
    start1 = start(anchors(coords1, "first")),
    end1 = end(anchors(coords1, "first")),
    chr2 = as.character(seqnames(anchors(coords1, "second"))),
    start2 = start(anchors(coords1, "second")),
    end2 = end(anchors(coords1, "second")),
    index1 = 1:length(coords1)
  )

  df2 <- data.frame(
    chr1 = as.character(seqnames(anchors(coords2, "first"))),
    start1 = start(anchors(coords2, "first")),
    end1 = end(anchors(coords2, "first")),
    chr2 = as.character(seqnames(anchors(coords2, "second"))),
    start2 = start(anchors(coords2, "second")),
    end2 = end(anchors(coords2, "second")),
    index2 = 1:length(coords2)
  )

  tolerance_bp <- tolerance_kb * 1000

  # Find overlaps
  matches <- data.frame(index1 = integer(), index2 = integer())

  for (i in 1:nrow(df1)) {
    # Same chromosome check
    same_chr <- df2$chr1 == df1$chr1[i] & df2$chr2 == df1$chr2[i]

    # Distance check (within tolerance for both anchors)
    dist1 <- abs(df2$start1 - df1$start1[i])
    dist2 <- abs(df2$start2 - df1$start2[i])

    close_enough <- same_chr & dist1 <= tolerance_bp & dist2 <= tolerance_bp

    if (any(close_enough)) {
      # Take closest match
      distances <- dist1 + dist2
      best_match <- which(close_enough)[which.min(distances[close_enough])]
      matches <- rbind(matches, data.frame(index1 = i, index2 = best_match))
    }
  }

  return(matches)
}

cat("Matching loops across resolutions (10kb tolerance)...\n")

# Match loops between all pairs of resolutions
overlap_matrix <- matrix(0, nrow = length(resolutions), ncol = length(resolutions))
rownames(overlap_matrix) <- paste0(resolutions/1000, "kb")
colnames(overlap_matrix) <- paste0(resolutions/1000, "kb")

for (i in 1:length(resolutions)) {
  for (j in 1:length(resolutions)) {
    res_i <- as.character(resolutions[i])
    res_j <- as.character(resolutions[j])

    if (i == j) {
      overlap_matrix[i, j] <- nrow(results_list[[res_i]])
    } else if (res_i %in% names(coords_list) && res_j %in% names(coords_list)) {
      matches <- match_loops(coords_list[[res_i]], coords_list[[res_j]], tolerance_kb = 10)
      overlap_matrix[i, j] <- nrow(matches)
    }
  }
}

cat("\nLoop overlap matrix:\n")
print(overlap_matrix)

# Save overlap matrix
write.table(overlap_matrix,
            file.path(output_dir, "loop_overlap_matrix.tsv"),
            sep = "\t", quote = FALSE)

cat("\n✓ Overlap analysis complete\n\n")

# =============================================================================
# SECTION 4: DIFFERENTIAL CONCORDANCE ANALYSIS
# =============================================================================

cat("SECTION 4: Differential Concordance\n")
cat("------------------------------------\n\n")

# For matched loops, check if they're differential at the same direction
cat("Analyzing concordance of differential calls...\n\n")

# Compare 5kb vs 10kb
if ("5000" %in% names(coords_list) && "10000" %in% names(coords_list)) {
  matches_5_10 <- match_loops(coords_list[["5000"]], coords_list[["10000"]], tolerance_kb = 10)

  res_5kb <- results_list[["5000"]]
  res_10kb <- results_list[["10000"]]

  matched_data <- data.frame(
    sig_5kb = res_5kb$significant[matches_5_10$index1],
    sig_10kb = res_10kb$significant[matches_5_10$index2],
    fc_5kb = res_5kb$logFC[matches_5_10$index1],
    fc_10kb = res_10kb$logFC[matches_5_10$index2]
  )

  # Concordance statistics
  both_sig <- sum(matched_data$sig_5kb & matched_data$sig_10kb, na.rm = TRUE)
  both_sig_same_dir <- sum(matched_data$sig_5kb & matched_data$sig_10kb &
                            sign(matched_data$fc_5kb) == sign(matched_data$fc_10kb), na.rm = TRUE)

  cat("5kb vs 10kb:\n")
  cat(sprintf("  Matched loops: %d\n", nrow(matched_data)))
  cat(sprintf("  Both significant: %d (%.1f%%)\n", both_sig, 100*both_sig/nrow(matched_data)))
  cat(sprintf("  Concordant direction: %d (%.1f%% of both-sig)\n",
              both_sig_same_dir, 100*both_sig_same_dir/both_sig))
  cat(sprintf("  Fold-change correlation: %.3f\n\n",
              cor(matched_data$fc_5kb, matched_data$fc_10kb, use = "complete.obs")))
}

# Compare 5kb vs 25kb
if ("5000" %in% names(coords_list) && "25000" %in% names(coords_list)) {
  matches_5_25 <- match_loops(coords_list[["5000"]], coords_list[["25000"]], tolerance_kb = 10)

  res_5kb <- results_list[["5000"]]
  res_25kb <- results_list[["25000"]]

  matched_data <- data.frame(
    sig_5kb = res_5kb$significant[matches_5_25$index1],
    sig_25kb = res_25kb$significant[matches_5_25$index2],
    fc_5kb = res_5kb$logFC[matches_5_25$index1],
    fc_25kb = res_25kb$logFC[matches_5_25$index2]
  )

  both_sig <- sum(matched_data$sig_5kb & matched_data$sig_25kb, na.rm = TRUE)
  both_sig_same_dir <- sum(matched_data$sig_5kb & matched_data$sig_25kb &
                            sign(matched_data$fc_5kb) == sign(matched_data$fc_25kb), na.rm = TRUE)

  cat("5kb vs 25kb:\n")
  cat(sprintf("  Matched loops: %d\n", nrow(matched_data)))
  cat(sprintf("  Both significant: %d (%.1f%%)\n", both_sig, 100*both_sig/nrow(matched_data)))
  cat(sprintf("  Concordant direction: %d (%.1f%% of both-sig)\n",
              both_sig_same_dir, 100*both_sig_same_dir/both_sig))
  cat(sprintf("  Fold-change correlation: %.3f\n\n",
              cor(matched_data$fc_5kb, matched_data$fc_25kb, use = "complete.obs")))
}

# Compare 10kb vs 25kb
if ("10000" %in% names(coords_list) && "25000" %in% names(coords_list)) {
  matches_10_25 <- match_loops(coords_list[["10000"]], coords_list[["25000"]], tolerance_kb = 10)

  res_10kb <- results_list[["10000"]]
  res_25kb <- results_list[["25000"]]

  matched_data <- data.frame(
    sig_10kb = res_10kb$significant[matches_10_25$index1],
    sig_25kb = res_25kb$significant[matches_10_25$index2],
    fc_10kb = res_10kb$logFC[matches_10_25$index1],
    fc_25kb = res_25kb$logFC[matches_10_25$index2]
  )

  both_sig <- sum(matched_data$sig_10kb & matched_data$sig_25kb, na.rm = TRUE)
  both_sig_same_dir <- sum(matched_data$sig_10kb & matched_data$sig_25kb &
                            sign(matched_data$fc_10kb) == sign(matched_data$fc_25kb), na.rm = TRUE)

  cat("10kb vs 25kb:\n")
  cat(sprintf("  Matched loops: %d\n", nrow(matched_data)))
  cat(sprintf("  Both significant: %d (%.1f%%)\n", both_sig, 100*both_sig/nrow(matched_data)))
  cat(sprintf("  Concordant direction: %d (%.1f%% of both-sig)\n",
              both_sig_same_dir, 100*both_sig_same_dir/both_sig))
  cat(sprintf("  Fold-change correlation: %.3f\n\n",
              cor(matched_data$fc_10kb, matched_data$fc_25kb, use = "complete.obs")))
}

cat("✓ Concordance analysis complete\n\n")

# =============================================================================
# SECTION 5: VISUALIZATIONS
# =============================================================================

cat("SECTION 5: Generating Visualizations\n")
cat("-------------------------------------\n\n")

# 1. Bar plot of differential loop counts by resolution (standard vs stringent)
cat("Creating differential loop count plot... ")
pdf(file.path(output_dir, "differential_loops_by_resolution.pdf"), width = 10, height = 6)

# Prepare data for standard thresholds
summary_standard <- summary_stats %>%
  select(resolution_kb, up_in_mutant, down_in_mutant) %>%
  pivot_longer(cols = c(up_in_mutant, down_in_mutant),
               names_to = "direction", values_to = "count") %>%
  mutate(
    direction = ifelse(direction == "up_in_mutant", "Up", "Down"),
    threshold = "Standard (FDR < 0.05)"
  )

# Prepare data for stringent thresholds (if available)
if (any(!is.na(summary_stats$final_loops))) {
  summary_final <- summary_stats %>%
    select(resolution_kb, final_up, final_down) %>%
    filter(!is.na(final_up)) %>%
    pivot_longer(cols = c(final_up, final_down),
                 names_to = "direction", values_to = "count") %>%
    mutate(
      direction = ifelse(direction == "final_up", "Up", "Down"),
      threshold = "Stringent (|logFC| > 0.3, FDR < 0.03)"
    )

  # Combine
  plot_data <- rbind(summary_standard, summary_final)

  ggplot(plot_data, aes(x = factor(resolution_kb), y = count, fill = interaction(direction, threshold))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    scale_fill_manual(
      name = "Threshold & Direction",
      values = c(
        "Up.Standard (FDR < 0.05)" = "#d73027",
        "Down.Standard (FDR < 0.05)" = "#4575b4",
        "Up.Stringent (|logFC| > 0.3, FDR < 0.03)" = "#a50026",
        "Down.Stringent (|logFC| > 0.3, FDR < 0.03)" = "#313695"
      ),
      labels = c(
        "Up (Standard)",
        "Down (Standard)",
        "Up (Stringent)",
        "Down (Stringent)"
      )
    ) +
    labs(
      title = "Differential Loops by Resolution and Threshold",
      x = "Resolution (kb)",
      y = "Number of Differential Loops"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
} else {
  # Only standard thresholds available
  ggplot(summary_standard, aes(x = factor(resolution_kb), y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#4575b4")) +
    labs(
      title = "Differential Loops by Resolution (Standard Thresholds)",
      x = "Resolution (kb)",
      y = "Number of Significant Loops (FDR < 0.05)",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
}

dev.off()
cat("✓\n")

# 1b. Filtering cascade plot (if final results available)
if (any(!is.na(summary_stats$final_loops))) {
  cat("Creating filtering cascade plot... ")
  pdf(file.path(output_dir, "filtering_cascade_by_resolution.pdf"), width = 10, height = 6)

  cascade_data <- summary_stats %>%
    filter(!is.na(final_loops)) %>%
    select(resolution_kb, total_loops, significant_fdr05, final_loops) %>%
    pivot_longer(cols = c(total_loops, significant_fdr05, final_loops),
                 names_to = "stage", values_to = "count") %>%
    mutate(
      stage = factor(stage,
                    levels = c("total_loops", "significant_fdr05", "final_loops"),
                    labels = c("Total Tested", "Significant\n(FDR < 0.05)", "Final\n(|logFC| > 0.3, FDR < 0.03)"))
    )

  ggplot(cascade_data, aes(x = stage, y = count, fill = factor(resolution_kb), group = factor(resolution_kb))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_line(aes(color = factor(resolution_kb)), position = position_dodge(width = 0.9), size = 1) +
    geom_point(aes(color = factor(resolution_kb)), position = position_dodge(width = 0.9), size = 3) +
    scale_fill_brewer(palette = "Set1", name = "Resolution") +
    scale_color_brewer(palette = "Set1", name = "Resolution") +
    labs(
      title = "Filtering Cascade: From Total Loops to Final Results",
      x = "Analysis Stage",
      y = "Number of Loops"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      axis.text.x = element_text(size = 10)
    )

  dev.off()
  cat("✓\n")
}

# 2. Fold-change correlation scatter plots
if ("5000" %in% names(results_list) && "10000" %in% names(results_list)) {
  cat("Creating fold-change correlation plots... ")

  # 5kb vs 10kb
  matches_5_10 <- match_loops(coords_list[["5000"]], coords_list[["10000"]], tolerance_kb = 10)

  plot_data <- data.frame(
    fc_5kb = results_list[["5000"]]$logFC[matches_5_10$index1],
    fc_10kb = results_list[["10000"]]$logFC[matches_5_10$index2],
    sig_5kb = results_list[["5000"]]$significant[matches_5_10$index1],
    sig_10kb = results_list[["10000"]]$significant[matches_5_10$index2]
  ) %>%
    mutate(
      status = case_when(
        sig_5kb & sig_10kb ~ "Both Significant",
        sig_5kb ~ "5kb Only",
        sig_10kb ~ "10kb Only",
        TRUE ~ "Neither"
      )
    )

  pdf(file.path(output_dir, "foldchange_correlation_5kb_vs_10kb.pdf"), width = 8, height = 7)

  ggplot(plot_data, aes(x = fc_5kb, y = fc_10kb, color = status)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c(
      "Both Significant" = "#1b7837",
      "5kb Only" = "#d73027",
      "10kb Only" = "#4575b4",
      "Neither" = "#999999"
    )) +
    labs(
      title = "Fold-Change Correlation: 5kb vs 10kb",
      subtitle = sprintf("Pearson r = %.3f (n = %d matched loops)",
                        cor(plot_data$fc_5kb, plot_data$fc_10kb, use = "complete.obs"),
                        nrow(plot_data)),
      x = "log2 Fold Change (5kb)",
      y = "log2 Fold Change (10kb)",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "right"
    )

  dev.off()
  cat("✓\n")
}

# 3. Venn diagram (if 3 resolutions available)
if (length(results_list) == 3) {
  cat("Creating Venn diagram... ")

  sig_5kb <- which(results_list[["5000"]]$significant)
  sig_10kb <- which(results_list[["10000"]]$significant)
  sig_25kb <- which(results_list[["25000"]]$significant)

  # Map to matched loops
  matches_5_10 <- match_loops(coords_list[["5000"]], coords_list[["10000"]], tolerance_kb = 10)
  matches_5_25 <- match_loops(coords_list[["5000"]], coords_list[["25000"]], tolerance_kb = 10)

  # Create sets based on 5kb as reference
  set_5kb <- sig_5kb
  set_10kb <- matches_5_10$index1[matches_5_10$index2 %in% sig_10kb]
  set_25kb <- matches_5_25$index1[matches_5_25$index2 %in% sig_25kb]

  venn_list <- list(
    "5kb" = set_5kb,
    "10kb" = set_10kb,
    "25kb" = set_25kb
  )

  pdf(file.path(output_dir, "venn_diagram_differential_loops.pdf"), width = 8, height = 8)

  venn.plot <- venn.diagram(
    x = venn_list,
    filename = NULL,
    category.names = c("5kb", "10kb", "25kb"),
    col = "transparent",
    fill = c("#d73027", "#4575b4", "#1b7837"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    cat.fontface = "bold",
    main = "Differential Loops Across Resolutions",
    main.cex = 1.8
  )

  grid.draw(venn.plot)

  dev.off()
  cat("✓\n")
}

cat("\n✓ All visualizations generated\n\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("========================================\n")
cat("MULTI-RESOLUTION ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Output directory: outputs/resolution_comparison/\n\n")

cat("Files generated:\n")
cat("  1. summary_statistics_by_resolution.tsv\n")
cat("     → Loop counts and significance rates (standard + stringent)\n\n")
cat("  2. loop_overlap_matrix.tsv\n")
cat("     → Overlap counts between resolutions\n\n")
cat("  3. differential_loops_by_resolution.pdf\n")
cat("     → Bar plot comparing standard vs stringent thresholds\n\n")
if (any(!is.na(summary_stats$final_loops))) {
  cat("  4. filtering_cascade_by_resolution.pdf\n")
  cat("     → Filtering progression from total to final results\n\n")
  cat("  5. foldchange_correlation_*.pdf\n")
  cat("     → Scatter plots of fold-change concordance\n\n")
  cat("  6. venn_diagram_differential_loops.pdf\n")
  cat("     → Venn diagram of shared differential loops\n\n")
} else {
  cat("  4. foldchange_correlation_*.pdf\n")
  cat("     → Scatter plots of fold-change concordance\n\n")
  cat("  5. venn_diagram_differential_loops.pdf\n")
  cat("     → Venn diagram of shared differential loops\n\n")
}

cat("Key findings:\n")
cat("\nStandard thresholds (FDR < 0.05):\n")
for (i in 1:nrow(summary_stats)) {
  res_kb <- summary_stats$resolution_kb[i]
  n_sig <- summary_stats$significant_fdr05[i]
  pct <- summary_stats$pct_significant[i]
  cat(sprintf("  %dkb: %d differential loops (%.1f%%)\n", res_kb, n_sig, pct))
}

if (any(!is.na(summary_stats$final_loops))) {
  cat("\nStringent thresholds (|logFC| > 0.3, FDR < 0.03):\n")
  for (i in 1:nrow(summary_stats)) {
    res_kb <- summary_stats$resolution_kb[i]
    if (!is.na(summary_stats$final_loops[i])) {
      n_final <- summary_stats$final_loops[i]
      pct_final <- summary_stats$pct_final[i]
      pct_of_sig <- 100 * n_final / summary_stats$significant_fdr05[i]
      cat(sprintf("  %dkb: %d loops (%.1f%% of total, %.1f%% of standard significant)\n",
                  res_kb, n_final, pct_final, pct_of_sig))
    }
  }
}

cat("\n========================================\n")
