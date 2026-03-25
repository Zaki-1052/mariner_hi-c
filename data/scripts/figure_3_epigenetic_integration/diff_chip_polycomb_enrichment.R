# scripts/diff_chip_polycomb_enrichment.R
# Task 3f: Differential H3K27me3/H2AK119ub Overlap with Loop Categories
# Author: Zakir Alibhai
# Date: 2026-02-01
#
# Purpose:
#   Assess whether differential H3K27me3 and H2AK119ub peaks (from DiffBind)
#   overlap with specific loop categories:
#   (a) Long-range lost loops (down_in_mutant, distance > 500kb)
#   (b) Short-range gained loops (up_in_mutant, distance < 500kb)
#   (c) Unchanged loops (FDR >= 0.1)
#
# Biological hypothesis (BAP1-KO):
#   - K27me3_down enriched at long-range lost loops (losing repression)
#   - H2AK119ub_down enriched at long-range lost loops (losing ubiquitination)
#
# Caveat: 400bp summit width may underestimate broader domains -> expect false negatives
#
# Usage:
#   Rscript scripts/diff_chip_polycomb_enrichment.R
#
# Output:
#   output/diff_chip_polycomb_enrichment/
#     ├── tables/
#     ├── plots/
#     └── enrichment_analysis_summary.txt

# ==============================================================================
# LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(pheatmap)
})

# Load multi-format output utility
source("data/scripts/_shared/multi_format_output.R") # Original: source("scripts/utils/multi_format_output.R")

cat("=",'%>%' |> rep(79) |> paste(collapse=""), "\n")
cat("Task 3f: Differential H3K27me3/H2AK119ub Overlap with Loop Categories\n")
cat("=",'%>%' |> rep(79) |> paste(collapse=""), "\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Output directories
TSV_DIR  <- "data/tsvs/figure_3_epigenetic_integration"   # Original: OUTPUT_DIR <- "output/diff_chip_polycomb_enrichment"
PLOT_DIR <- "data/plots/figure_3_epigenetic_integration"   # Original: (plots went under OUTPUT_DIR/plots)
OUTPUT_DIR <- TSV_DIR  # kept for summary text file writes

# Input files
# Use merged_all_results.tsv which contains ALL loops (including unchanged) across resolutions
LOOP_FILES <- list(
  # All loops merged across resolutions (includes unchanged)
  all_loops = "outputs/250402-late_outputs/merged_loops/merged_all_results.tsv",  # Note: repo-relative path, not bundled in data/
  # Polycomb shared from task 3c (differential only)
  polycomb_shared = "data/tsvs/supplemental/polycomb_shared_loops.tsv"  # Original: output/shared_anchor_analysis/late/polycomb_specific/tables/polycomb_shared_loops.tsv
)

# Differential peak files (late timepoint for K27me3, H2AK119ub is not timepoint-specific)
DIFF_PEAK_FILES <- list(
  K27me3_down = "peaks/new/adult_K27me3_down.bed",       # Note: repo-relative path, not bundled in data/
  K27me3_up = "peaks/new/adult_K27me3_up.bed",           # Note: repo-relative path, not bundled in data/
  H2AK119ub_down = "peaks/new/H2AK119ub_down.bed",      # Note: repo-relative path, not bundled in data/
  H2AK119ub_up = "peaks/new/H2AK119ub_up.bed"           # Note: repo-relative path, not bundled in data/
)

# Analysis parameters
DISTANCE_THRESHOLD <- 500000  # 500kb threshold for long/short range
FDR_UNCHANGED_THRESHOLD <- 0.1  # FDR >= 0.1 for unchanged loops

# Color palette
COLORS <- list(
  long_range_lost = "#2166ac",      # Blue
  short_range_gained = "#b2182b",   # Red
  unchanged = "#969696",            # Gray
  other_differential = "#756bb1",   # Purple
  K27me3_down = "#1b9e77",
  K27me3_up = "#d95f02",
  H2AK119ub_down = "#7570b3",
  H2AK119ub_up = "#e7298a"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Load ChIP-seq peaks from BED file
#' Adapted from annotate_loops_extended.R
#'
#' @param bed_path Path to BED file
#' @param mark_name Name of the mark (for logging)
#' @return GRanges object
load_chip_peaks <- function(bed_path, mark_name = "ChIP") {
  if (!file.exists(bed_path)) {
    stop(sprintf("%s bed file not found: %s", mark_name, bed_path))
  }
  cat(sprintf("  Loading %s peaks from: %s\n", mark_name, bed_path))

  df <- read.table(bed_path, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE)
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2, end = df$V3)
  )

  cat(sprintf("    Loaded %d peaks\n", length(gr)))
  gr
}

#' Load all differential peak files
#'
#' @param peak_files Named list of file paths
#' @return Named list of GRanges objects
load_all_diff_peaks <- function(peak_files) {
  cat("\n[1] Loading differential peak files...\n")
  peaks <- list()
  for (name in names(peak_files)) {
    peaks[[name]] <- load_chip_peaks(peak_files[[name]], name)
  }
  peaks
}

#' Categorize loops into analysis categories
#'
#' Categories:
#'   - long_range_lost: down_in_mutant AND distance > 500kb
#'   - short_range_gained: up_in_mutant AND distance < 500kb
#'   - unchanged: FDR >= 0.1
#'   - other_differential: everything else (FDR < 0.1 but not fitting above)
#'
#' @param loops_df Data frame with loops
#' @param distance_threshold Distance threshold (default 500kb)
#' @param fdr_threshold FDR threshold for unchanged (default 0.1)
#' @return Data frame with loop_category column
categorize_loops <- function(loops_df, distance_threshold = 500000, fdr_threshold = 0.1) {
  # Calculate loop_distance if not present
  if (!"loop_distance" %in% colnames(loops_df)) {
    # Calculate midpoints and distance
    loops_df <- loops_df %>%
      mutate(
        anchor1_mid = (start1 + end1) / 2,
        anchor2_mid = (start2 + end2) / 2,
        loop_distance = abs(anchor2_mid - anchor1_mid)
      )
    cat("  Calculated loop distances from coordinates\n")
  }

  loops_df %>%
    mutate(
      loop_category = case_when(
        direction == "down_in_mutant" & loop_distance > distance_threshold ~ "long_range_lost",
        direction == "up_in_mutant" & loop_distance < distance_threshold ~ "short_range_gained",
        FDR >= fdr_threshold ~ "unchanged",
        TRUE ~ "other_differential"
      ),
      loop_category = factor(loop_category,
                             levels = c("long_range_lost", "short_range_gained",
                                       "unchanged", "other_differential"))
    )
}

#' Compute overlaps between loop anchors and differential peaks
#'
#' @param loops_df Data frame with loops (must have anchor coordinates)
#' @param diff_peaks_list Named list of GRanges with differential peaks
#' @return Data frame with overlap columns for each mark
compute_diff_chip_overlaps <- function(loops_df, diff_peaks_list) {
  cat("\n[3] Computing overlaps between anchors and differential peaks...\n")

  # Create GRanges for anchors
  anchor1_gr <- GRanges(
    seqnames = loops_df$chr1,
    ranges = IRanges(start = loops_df$start1, end = loops_df$end1)
  )
  anchor2_gr <- GRanges(
    seqnames = loops_df$chr2,
    ranges = IRanges(start = loops_df$start2, end = loops_df$end2)
  )

  # Compute overlaps for each differential peak set
  for (mark in names(diff_peaks_list)) {
    peaks_gr <- diff_peaks_list[[mark]]

    anchor1_overlap <- countOverlaps(anchor1_gr, peaks_gr) > 0
    anchor2_overlap <- countOverlaps(anchor2_gr, peaks_gr) > 0

    # Either anchor overlapping
    col_name <- paste0(mark, "_overlap")
    loops_df[[col_name]] <- anchor1_overlap | anchor2_overlap

    # Also track individual anchors
    loops_df[[paste0(mark, "_anchor1")]] <- anchor1_overlap
    loops_df[[paste0(mark, "_anchor2")]] <- anchor2_overlap

    n_overlap <- sum(loops_df[[col_name]])
    cat(sprintf("  %s: %d loops overlap (%0.1f%%)\n",
                mark, n_overlap, 100 * n_overlap / nrow(loops_df)))
  }

  loops_df
}

#' Summarize overlap counts by loop category
#'
#' @param loops_df Data frame with categorized loops and overlap columns
#' @param mark_names Names of marks to summarize
#' @return Summary data frame
summarize_overlaps_by_category <- function(loops_df, mark_names) {
  summary_list <- list()

  for (cat in levels(loops_df$loop_category)) {
    cat_loops <- loops_df %>% filter(loop_category == cat)
    n_total <- nrow(cat_loops)

    for (mark in mark_names) {
      col_name <- paste0(mark, "_overlap")
      n_overlap <- sum(cat_loops[[col_name]])
      pct_overlap <- 100 * n_overlap / n_total

      summary_list[[length(summary_list) + 1]] <- tibble(
        loop_category = cat,
        mark = mark,
        n_total = n_total,
        n_overlap = n_overlap,
        pct_overlap = pct_overlap
      )
    }
  }

  bind_rows(summary_list)
}

#' Run Fisher's exact test for enrichment
#'
#' @param loops_df Data frame with categorized loops
#' @param mark Name of the mark
#' @param test_category Category to test (e.g., "long_range_lost")
#' @param reference_category Reference category (default "unchanged")
#' @return List with test results
run_fisher_test <- function(loops_df, mark, test_category, reference_category = "unchanged") {
  col_name <- paste0(mark, "_overlap")

  # Get counts for test category
  test_loops <- loops_df %>% filter(loop_category == test_category)
  test_overlap <- sum(test_loops[[col_name]])
  test_no_overlap <- nrow(test_loops) - test_overlap

  # Get counts for reference category
  ref_loops <- loops_df %>% filter(loop_category == reference_category)
  ref_overlap <- sum(ref_loops[[col_name]])
  ref_no_overlap <- nrow(ref_loops) - ref_overlap

  # 2x2 contingency table
  #                    Test      Reference
  # Overlap            a         b
  # No_Overlap         c         d

  mat <- matrix(c(test_overlap, ref_overlap,
                  test_no_overlap, ref_no_overlap),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("Overlap", "No_Overlap"),
                               c(test_category, reference_category)))

  # Fisher's exact test
  fisher_result <- tryCatch({
    fisher.test(mat)
  }, error = function(e) {
    list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
  })

  list(
    mark = mark,
    test_category = test_category,
    reference_category = reference_category,
    contingency_table = mat,
    test_n = nrow(test_loops),
    test_overlap = test_overlap,
    test_pct = 100 * test_overlap / nrow(test_loops),
    ref_n = nrow(ref_loops),
    ref_overlap = ref_overlap,
    ref_pct = 100 * ref_overlap / nrow(ref_loops),
    odds_ratio = as.numeric(fisher_result$estimate),
    ci_lower = fisher_result$conf.int[1],
    ci_upper = fisher_result$conf.int[2],
    p_value = fisher_result$p.value
  )
}

#' Run all Fisher's tests for enrichment analysis
#'
#' Tests:
#'   - Each mark (4) x Each comparison (3) = 12 tests
#'   - Comparisons: long_range_lost vs unchanged, short_range_gained vs unchanged,
#'                  long_range_lost vs short_range_gained
#'
#' @param loops_df Data frame with categorized loops
#' @param mark_names Names of marks to test
#' @return Data frame with all test results
run_all_fisher_tests <- function(loops_df, mark_names) {
  cat("\n[4] Running Fisher's exact tests for enrichment...\n")

  comparisons <- list(
    c("long_range_lost", "unchanged"),
    c("short_range_gained", "unchanged"),
    c("long_range_lost", "short_range_gained")
  )

  results_list <- list()

  for (mark in mark_names) {
    for (comp in comparisons) {
      test_cat <- comp[1]
      ref_cat <- comp[2]

      result <- run_fisher_test(loops_df, mark, test_cat, ref_cat)
      results_list[[length(results_list) + 1]] <- as_tibble(result[!names(result) %in% "contingency_table"])
    }
  }

  results_df <- bind_rows(results_list)

  # Apply BH-FDR correction across all tests
  results_df$fdr <- p.adjust(results_df$p_value, method = "BH")

  # Add significance indicators
  results_df <- results_df %>%
    mutate(
      significant = fdr < 0.05,
      enrichment = case_when(
        odds_ratio > 1 & significant ~ "enriched",
        odds_ratio < 1 & significant ~ "depleted",
        TRUE ~ "ns"
      )
    )

  cat(sprintf("  Completed %d tests, %d significant at FDR < 0.05\n",
              nrow(results_df), sum(results_df$significant, na.rm = TRUE)))

  results_df
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create bar plot of overlap percentages by loop category
#'
#' @param summary_df Summary data frame from summarize_overlaps_by_category
#' @param title Plot title
#' @return ggplot object
create_overlap_barplot <- function(summary_df, title = "Differential ChIP-seq Peak Overlap by Loop Category") {
  # Order marks for display
  summary_df <- summary_df %>%
    mutate(
      mark = factor(mark, levels = c("K27me3_down", "K27me3_up",
                                     "H2AK119ub_down", "H2AK119ub_up")),
      mark_label = case_when(
        mark == "K27me3_down" ~ "K27me3\ndown",
        mark == "K27me3_up" ~ "K27me3\nup",
        mark == "H2AK119ub_down" ~ "H2AK119ub\ndown",
        mark == "H2AK119ub_up" ~ "H2AK119ub\nup"
      ),
      category_label = case_when(
        loop_category == "long_range_lost" ~ "Long-range\nLost (>500kb)",
        loop_category == "short_range_gained" ~ "Short-range\nGained (<500kb)",
        loop_category == "unchanged" ~ "Unchanged",
        loop_category == "other_differential" ~ "Other\nDifferential"
      ),
      category_label = factor(category_label, levels = c(
        "Long-range\nLost (>500kb)", "Short-range\nGained (<500kb)",
        "Unchanged", "Other\nDifferential"
      ))
    )

  ggplot(summary_df, aes(x = category_label, y = pct_overlap, fill = mark)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%d", n_overlap)),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 2.5) +
    scale_fill_manual(
      values = c(
        "K27me3_down" = COLORS$K27me3_down,
        "K27me3_up" = COLORS$K27me3_up,
        "H2AK119ub_down" = COLORS$H2AK119ub_down,
        "H2AK119ub_up" = COLORS$H2AK119ub_up
      ),
      labels = c("K27me3 down", "K27me3 up", "H2AK119ub down", "H2AK119ub up")
    ) +
    labs(
      title = title,
      subtitle = sprintf("Overlap with differential peaks (400bp summits)"),
      x = "Loop Category",
      y = "% Loops with Overlap",
      fill = "Differential\nPeak Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 10),
      legend.position = "right"
    ) +
    ylim(0, max(summary_df$pct_overlap) * 1.15)
}

#' Create enrichment dot plot
#'
#' @param results_df Results from run_all_fisher_tests
#' @param title Plot title
#' @return ggplot object
create_enrichment_dotplot <- function(results_df, title = "Enrichment Analysis: Differential Peaks at Loop Anchors") {
  # Filter to main comparisons (vs unchanged)
  plot_df <- results_df %>%
    filter(reference_category == "unchanged") %>%
    mutate(
      mark = factor(mark, levels = c("K27me3_down", "K27me3_up",
                                     "H2AK119ub_down", "H2AK119ub_up")),
      mark_label = case_when(
        mark == "K27me3_down" ~ "K27me3 down",
        mark == "K27me3_up" ~ "K27me3 up",
        mark == "H2AK119ub_down" ~ "H2AK119ub down",
        mark == "H2AK119ub_up" ~ "H2AK119ub up"
      ),
      category_label = case_when(
        test_category == "long_range_lost" ~ "Long-range Lost",
        test_category == "short_range_gained" ~ "Short-range Gained"
      ),
      log2_or = log2(odds_ratio),
      neg_log10_fdr = -log10(fdr + 1e-10)
    )

  ggplot(plot_df, aes(x = log2_or, y = mark_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = neg_log10_fdr, color = significant), alpha = 0.8) +
    geom_errorbarh(aes(xmin = log2(ci_lower), xmax = log2(ci_upper)),
                   height = 0.2, alpha = 0.5) +
    facet_wrap(~category_label, ncol = 1) +
    scale_color_manual(values = c("TRUE" = "#d62728", "FALSE" = "gray60"),
                       labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS")) +
    scale_size_continuous(
      name = "-log10(FDR)",
      range = c(2, 8),
      breaks = c(1, 2, 3, 5)
    ) +
    labs(
      title = title,
      subtitle = "Comparison to unchanged loops (FDR >= 0.1)",
      x = "log2(Odds Ratio)",
      y = "",
      color = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

#' Create heatmap of overlap percentages
#'
#' @param summary_df Summary data frame
#' @param title Plot title
#' @return pheatmap object
create_overlap_heatmap <- function(summary_df, title = "Overlap Heatmap") {
  # Map category names for display
  category_labels <- c(
    "long_range_lost" = "Long-range Lost",
    "short_range_gained" = "Short-range Gained",
    "unchanged" = "Unchanged",
    "other_differential" = "Other Differential"
  )

  mark_labels <- c(
    "K27me3_down" = "K27me3 down",
    "K27me3_up" = "K27me3 up",
    "H2AK119ub_down" = "H2AK119ub down",
    "H2AK119ub_up" = "H2AK119ub up"
  )

  # Pivot to matrix format
  mat <- summary_df %>%
    select(loop_category, mark, pct_overlap) %>%
    pivot_wider(names_from = mark, values_from = pct_overlap) %>%
    column_to_rownames("loop_category") %>%
    as.matrix()

  # Clean row and column names dynamically
  rownames(mat) <- category_labels[rownames(mat)]
  colnames(mat) <- mark_labels[colnames(mat)]

  pheatmap(
    mat,
    main = title,
    color = colorRampPalette(c("white", "#fee8c8", "#e34a33"))(50),
    display_numbers = TRUE,
    number_format = "%.1f",
    number_color = "black",
    fontsize_number = 10,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "gray80",
    cellwidth = 60,
    cellheight = 30,
    angle_col = 45
  )
}

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

run_analysis <- function() {

  # Create output directories
  dirs <- c(TSV_DIR, PLOT_DIR)  # Original: file.path(OUTPUT_DIR, "tables"), file.path(OUTPUT_DIR, "plots")
  for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # --------------------------------------------------------------------------
  # Load differential peaks
  # --------------------------------------------------------------------------
  diff_peaks <- load_all_diff_peaks(DIFF_PEAK_FILES)
  mark_names <- names(diff_peaks)

  # --------------------------------------------------------------------------
  # Load and categorize loops - ALL LOOPS
  # --------------------------------------------------------------------------
  cat("\n[2] Loading and categorizing loops...\n")

  all_loops <- read.table(LOOP_FILES$all_loops, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE)
  cat(sprintf("  All loops loaded: %d\n", nrow(all_loops)))

  all_loops <- categorize_loops(all_loops, DISTANCE_THRESHOLD, FDR_UNCHANGED_THRESHOLD)

  cat("\n  Loop category distribution:\n")
  cat_counts <- table(all_loops$loop_category)
  for (cat in names(cat_counts)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", cat, cat_counts[cat],
                100 * cat_counts[cat] / nrow(all_loops)))
  }

  # --------------------------------------------------------------------------
  # Load and categorize loops - POLYCOMB SHARED
  # --------------------------------------------------------------------------
  polycomb_loops <- read.table(LOOP_FILES$polycomb_shared, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)
  cat(sprintf("\n  Polycomb shared loops loaded: %d\n", nrow(polycomb_loops)))

  polycomb_loops <- categorize_loops(polycomb_loops, DISTANCE_THRESHOLD, FDR_UNCHANGED_THRESHOLD)

  cat("\n  Polycomb loop category distribution:\n")
  pc_counts <- table(polycomb_loops$loop_category)
  for (cat in names(pc_counts)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", cat, pc_counts[cat],
                100 * pc_counts[cat] / nrow(polycomb_loops)))
  }

  # --------------------------------------------------------------------------
  # Compute overlaps - ALL LOOPS
  # --------------------------------------------------------------------------
  all_loops <- compute_diff_chip_overlaps(all_loops, diff_peaks)

  # --------------------------------------------------------------------------
  # Compute overlaps - POLYCOMB SHARED
  # --------------------------------------------------------------------------
  cat("\n  Computing overlaps for Polycomb shared loops...\n")
  polycomb_loops <- compute_diff_chip_overlaps(polycomb_loops, diff_peaks)

  # --------------------------------------------------------------------------
  # Summarize overlaps by category
  # --------------------------------------------------------------------------
  cat("\n[5] Summarizing overlaps by loop category...\n")

  summary_all <- summarize_overlaps_by_category(all_loops, mark_names)
  summary_polycomb <- summarize_overlaps_by_category(polycomb_loops, mark_names)

  # Save summary tables
  write.table(summary_all,
              file.path(TSV_DIR, "3C_overlap_summary_all_loops.tsv"),  # Original: file.path(OUTPUT_DIR, "tables", "overlap_summary_all_loops.tsv")
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_polycomb,
              file.path(TSV_DIR, "3C_overlap_summary_polycomb.tsv"),  # Original: file.path(OUTPUT_DIR, "tables", "overlap_summary_polycomb_shared.tsv")
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("  Saved overlap summary tables\n")

  # --------------------------------------------------------------------------
  # Run Fisher's exact tests
  # --------------------------------------------------------------------------
  enrichment_all <- run_all_fisher_tests(all_loops, mark_names)
  enrichment_polycomb <- run_all_fisher_tests(polycomb_loops, mark_names)

  # Save enrichment test results
  write.table(enrichment_all,
              file.path(TSV_DIR, "3C_enrichment_tests_all_loops.tsv"),  # Original: file.path(OUTPUT_DIR, "tables", "enrichment_tests_all_loops.tsv")
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(enrichment_polycomb,
              file.path(TSV_DIR, "3C_enrichment_tests_polycomb.tsv"),  # Original: file.path(OUTPUT_DIR, "tables", "enrichment_tests_polycomb_shared.tsv")
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("  Saved enrichment test results\n")

  # Save loops with overlap annotations
  write.table(all_loops,
              file.path(TSV_DIR, "3C_all_loops_with_diff_chip_overlap.tsv"),  # Original: file.path(OUTPUT_DIR, "tables", "all_loops_with_diff_chip_overlap.tsv")
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(polycomb_loops,
              file.path(TSV_DIR, "3C_polycomb_loops_with_diff_chip_overlap.tsv"),  # Original: file.path(OUTPUT_DIR, "tables", "polycomb_loops_with_diff_chip_overlap.tsv")
              sep = "\t", quote = FALSE, row.names = FALSE)

  # --------------------------------------------------------------------------
  # Generate visualizations
  # --------------------------------------------------------------------------
  cat("\n[6] Generating visualizations...\n")

  # Bar plots - All loops
  p_bar_all <- create_overlap_barplot(summary_all,
                                       "Differential ChIP-seq Peak Overlap - All Loops")
  save_multiformat_ggplot(p_bar_all,
                          file.path(PLOT_DIR, "3C_01_overlap_barplot_all_loops"),  # Original: file.path(OUTPUT_DIR, "plots", "01_overlap_barplot_all_loops")
                          width = 10, height = 7)

  # Bar plots - Polycomb shared
  p_bar_polycomb <- create_overlap_barplot(summary_polycomb,
                                            "Differential ChIP-seq Peak Overlap - Polycomb Shared Loops")
  save_multiformat_ggplot(p_bar_polycomb,
                          file.path(PLOT_DIR, "3C_01_overlap_barplot_polycomb_shared"),  # Original: file.path(OUTPUT_DIR, "plots", "01_overlap_barplot_polycomb_shared")
                          width = 10, height = 7)

  # Enrichment dot plots - All loops
  p_enrich_all <- create_enrichment_dotplot(enrichment_all,
                                             "Enrichment Analysis - All Loops")
  save_multiformat_ggplot(p_enrich_all,
                          file.path(PLOT_DIR, "3C_02_enrichment_dotplot_all_loops"),  # Original: file.path(OUTPUT_DIR, "plots", "02_enrichment_dotplot_all_loops")
                          width = 9, height = 7)

  # Enrichment dot plots - Polycomb shared
  p_enrich_polycomb <- create_enrichment_dotplot(enrichment_polycomb,
                                                  "Enrichment Analysis - Polycomb Shared Loops")
  save_multiformat_ggplot(p_enrich_polycomb,
                          file.path(PLOT_DIR, "3C_02_enrichment_dotplot_polycomb_shared"),  # Original: file.path(OUTPUT_DIR, "plots", "02_enrichment_dotplot_polycomb_shared")
                          width = 9, height = 7)

  # Heatmaps - create directly to avoid scoping issues with quote()
  # All loops heatmap
  heatmap_dir_all <- file.path(PLOT_DIR, "3C_03_heatmap_overlap_all_loops")  # Original: file.path(OUTPUT_DIR, "plots", "03_heatmap_overlap_all_loops")
  dir.create(heatmap_dir_all, recursive = TRUE, showWarnings = FALSE)

  pdf(file.path(heatmap_dir_all, "3C_03_heatmap_overlap_all_loops.pdf"), width = 8, height = 5)
  create_overlap_heatmap(summary_all, "Overlap % - All Loops")
  dev.off()

  svglite::svglite(file.path(heatmap_dir_all, "3C_03_heatmap_overlap_all_loops.svg"), width = 8, height = 5)
  create_overlap_heatmap(summary_all, "Overlap % - All Loops")
  dev.off()

  jpeg(file.path(heatmap_dir_all, "3C_03_heatmap_overlap_all_loops.jpg"), width = 8*300, height = 5*300, res = 300, quality = 95)
  create_overlap_heatmap(summary_all, "Overlap % - All Loops")
  dev.off()
  cat("  Saved: 03_heatmap_overlap_all_loops/{pdf,svg,jpg}\n")

  # Polycomb shared heatmap (only if there are non-zero categories)
  if (sum(pc_counts) > 0 && length(unique(summary_polycomb$loop_category[summary_polycomb$n_total > 0])) >= 2) {
    heatmap_dir_pc <- file.path(PLOT_DIR, "3C_03_heatmap_overlap_polycomb_shared")  # Original: file.path(OUTPUT_DIR, "plots", "03_heatmap_overlap_polycomb_shared")
    dir.create(heatmap_dir_pc, recursive = TRUE, showWarnings = FALSE)

    # Filter to categories with data
    summary_polycomb_filtered <- summary_polycomb %>%
      filter(n_total > 0)

    pdf(file.path(heatmap_dir_pc, "3C_03_heatmap_overlap_polycomb_shared.pdf"), width = 8, height = 5)
    create_overlap_heatmap(summary_polycomb_filtered, "Overlap % - Polycomb Shared")
    dev.off()

    svglite::svglite(file.path(heatmap_dir_pc, "3C_03_heatmap_overlap_polycomb_shared.svg"), width = 8, height = 5)
    create_overlap_heatmap(summary_polycomb_filtered, "Overlap % - Polycomb Shared")
    dev.off()

    jpeg(file.path(heatmap_dir_pc, "3C_03_heatmap_overlap_polycomb_shared.jpg"), width = 8*300, height = 5*300, res = 300, quality = 95)
    create_overlap_heatmap(summary_polycomb_filtered, "Overlap % - Polycomb Shared")
    dev.off()
    cat("  Saved: 03_heatmap_overlap_polycomb_shared/{pdf,svg,jpg}\n")
  } else {
    cat("  Skipped polycomb heatmap (insufficient categories with data)\n")
  }

  # --------------------------------------------------------------------------
  # Generate summary report
  # --------------------------------------------------------------------------
  cat("\n[7] Generating summary report...\n")

  report_lines <- c(
    "=",' ' |> rep(79) |> paste(collapse=""),
    "Task 3f: Differential H3K27me3/H2AK119ub Overlap with Loop Categories",
    "=",' ' |> rep(79) |> paste(collapse=""),
    "",
    sprintf("Analysis Date: %s", Sys.time()),
    "",
    "INPUT FILES:",
    "-------------",
    sprintf("All loops: %s (%d loops)", LOOP_FILES$all_loops, nrow(all_loops)),
    sprintf("Polycomb shared: %s (%d loops)", LOOP_FILES$polycomb_shared, nrow(polycomb_loops)),
    "",
    "Differential Peak Files:",
    sprintf("  K27me3_down: %d peaks", length(diff_peaks$K27me3_down)),
    sprintf("  K27me3_up: %d peaks", length(diff_peaks$K27me3_up)),
    sprintf("  H2AK119ub_down: %d peaks", length(diff_peaks$H2AK119ub_down)),
    sprintf("  H2AK119ub_up: %d peaks", length(diff_peaks$H2AK119ub_up)),
    "",
    "ANALYSIS PARAMETERS:",
    "--------------------",
    sprintf("Distance threshold for long/short range: %d bp", DISTANCE_THRESHOLD),
    sprintf("FDR threshold for unchanged loops: %.2f", FDR_UNCHANGED_THRESHOLD),
    "",
    "LOOP CATEGORY DISTRIBUTION - ALL LOOPS:",
    "---------------------------------------"
  )

  for (cat in names(cat_counts)) {
    report_lines <- c(report_lines,
                      sprintf("  %s: %d (%.1f%%)", cat, cat_counts[cat],
                              100 * cat_counts[cat] / nrow(all_loops)))
  }

  report_lines <- c(report_lines, "",
    "LOOP CATEGORY DISTRIBUTION - POLYCOMB SHARED:",
    "---------------------------------------------"
  )

  for (cat in names(pc_counts)) {
    report_lines <- c(report_lines,
                      sprintf("  %s: %d (%.1f%%)", cat, pc_counts[cat],
                              100 * pc_counts[cat] / nrow(polycomb_loops)))
  }

  # Significant enrichments
  report_lines <- c(report_lines, "",
    "SIGNIFICANT ENRICHMENTS (FDR < 0.05) - ALL LOOPS:",
    "-------------------------------------------------"
  )

  sig_all <- enrichment_all %>%
    filter(significant) %>%
    arrange(fdr)

  if (nrow(sig_all) > 0) {
    for (i in 1:nrow(sig_all)) {
      row <- sig_all[i, ]
      report_lines <- c(report_lines,
        sprintf("  %s at %s vs %s: OR=%.2f, FDR=%.2e (%s)",
                row$mark, row$test_category, row$reference_category,
                row$odds_ratio, row$fdr, row$enrichment))
    }
  } else {
    report_lines <- c(report_lines, "  No significant enrichments")
  }

  report_lines <- c(report_lines, "",
    "SIGNIFICANT ENRICHMENTS (FDR < 0.05) - POLYCOMB SHARED:",
    "-------------------------------------------------------"
  )

  sig_polycomb <- enrichment_polycomb %>%
    filter(significant) %>%
    arrange(fdr)

  if (nrow(sig_polycomb) > 0) {
    for (i in 1:nrow(sig_polycomb)) {
      row <- sig_polycomb[i, ]
      report_lines <- c(report_lines,
        sprintf("  %s at %s vs %s: OR=%.2f, FDR=%.2e (%s)",
                row$mark, row$test_category, row$reference_category,
                row$odds_ratio, row$fdr, row$enrichment))
    }
  } else {
    report_lines <- c(report_lines, "  No significant enrichments")
  }

  # Data patterns summary
  report_lines <- c(report_lines, "",
    "KEY DATA PATTERNS:",
    "------------------",
    "See README.md for detailed interpretation of results.",
    "",
    "CAVEATS:",
    "--------",
    "  - 400bp summit width from DiffBind may underestimate broader H3K27me3/H2AK119ub domains",
    "  - Expected false negatives due to narrow peak calls",
    "  - H2AK119ub files are not timepoint-specific",
    "  - Follow-up bigWig heatmap analysis recommended for broader domain visualization",
    "",
    "OUTPUT FILES:",
    "-------------",
    "  tables/overlap_summary_all_loops.tsv",
    "  tables/overlap_summary_polycomb_shared.tsv",
    "  tables/enrichment_tests_all_loops.tsv",
    "  tables/enrichment_tests_polycomb_shared.tsv",
    "  tables/all_loops_with_diff_chip_overlap.tsv",
    "  tables/polycomb_loops_with_diff_chip_overlap.tsv",
    "  plots/01_overlap_barplot_all_loops/{pdf,svg,jpg}",
    "  plots/01_overlap_barplot_polycomb_shared/{pdf,svg,jpg}",
    "  plots/02_enrichment_dotplot_all_loops/{pdf,svg,jpg}",
    "  plots/02_enrichment_dotplot_polycomb_shared/{pdf,svg,jpg}",
    "  plots/03_heatmap_overlap_all_loops/{pdf,svg,jpg}",
    "  plots/03_heatmap_overlap_polycomb_shared/{pdf,svg,jpg}"
  )

  writeLines(report_lines, file.path(TSV_DIR, "3C_enrichment_analysis_summary.txt"))  # Original: file.path(OUTPUT_DIR, "enrichment_analysis_summary.txt")
  cat("  Saved summary report\n")

  cat("\n", "=",'%>%' |> rep(50) |> paste(collapse=""), "\n")
  cat("Analysis complete!\n")
  cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
  cat("=",'%>%' |> rep(50) |> paste(collapse=""), "\n")

  # Return results for programmatic use
  invisible(list(
    all_loops = all_loops,
    polycomb_loops = polycomb_loops,
    summary_all = summary_all,
    summary_polycomb = summary_polycomb,
    enrichment_all = enrichment_all,
    enrichment_polycomb = enrichment_polycomb,
    diff_peaks = diff_peaks
  ))
}

# ==============================================================================
# EXECUTE
# ==============================================================================

if (!interactive()) {
  results <- run_analysis()
} else {
  cat("Script loaded. Run run_analysis() to execute.\n")
}
