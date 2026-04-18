# scripts/h2ak119ub_loop_integration.R
# H2AK119ub Loop Integration: Tasks 8b, 8c, and Section 3 Connection
# Author: Zakir Alibhai
# Date: 2026-02-10
#
# Purpose:
#   Integrate H2AK119ub (ubiquitination) data with differential chromatin loop
#   analysis. BAP1 is an H2AK119ub deubiquitinase; BAP1-KO accumulates K119ub
#   at Polycomb targets, altering chromatin contacts.
#
# Hypothesis: "Ubiquitination is a buffer to stop K27ac contact" -- once K119ub
# threshold is reached, long-range Polycomb contacts are destabilized.
#
# Sections:
#   A: Anchor-level K119ub annotation (peak overlaps)
#   B: CDF/Density distance distributions by K119ub status
#   C: Enrichment testing by loop category x chromatin state
#   D: Continuous signal correlation (loop logFC vs K119ub change)
#   E: Shared anchor K119ub analysis (loop switching sites)
#
# Usage:
#   Rscript scripts/h2ak119ub_loop_integration.R
#
# Output:
#   output/h2ak119ub_loop_integration/late/
#     tables/  (8 TSV files)
#     plots/   (19 multi-format plots)
#     analysis_summary.txt

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

source("scripts/utils/multi_format_output.R")

cat("================================================================================\n")
cat("H2AK119ub Loop Integration Analysis\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

OUTPUT_DIR <- "output/h2ak119ub_loop_integration/late"

INPUT_FILES <- list(
  all_loops        = "outputs/250402-late_outputs/merged_loops/merged_all_results.tsv",
  diff_loops       = "outputs/250402-late_outputs/merged_loops/characterized_loops.tsv",
  shared_anchors   = "output/shared_anchor_analysis/late/tables/shared_anchors.tsv",
  shared_loops     = "output/shared_anchor_analysis/late/tables/shared_anchor_loops.tsv",
  polycomb_shared  = "output/shared_anchor_analysis/late/polycomb_specific/tables/polycomb_shared_loops.tsv",
  signal           = "peaks/k119ub_anchor_signal.tsv"
)

PEAK_FILES <- list(
  K119ub_up   = "peaks/new/H2AK119ub_up.bed",
  K119ub_down = "peaks/new/H2AK119ub_down.bed",
  K119ub_ctrl = "peaks/intersect/P51_K119ub_ctrl_intersect.bed",
  K119ub_mut  = "peaks/intersect/P51_K119ub_mut_intersect.bed"
)

DISTANCE_THRESHOLD <- 500000
FDR_UNCHANGED_THRESHOLD <- 0.1
MIN_SAMPLES <- 10

DISTANCE_BINS <- list(
  "< 100kb"      = c(0, 100000),
  "100kb-500kb"   = c(100000, 500000),
  "500kb-1Mb"     = c(500000, 1000000),
  "> 1Mb"         = c(1000000, Inf)
)

COLORS <- list(
  lost     = "#d73027",
  gained   = "#4575b4",
  unchanged = "#969696",
  K119ub_up   = "#e7298a",
  K119ub_down = "#7570b3",
  enriched = "#d62728",
  depleted = "#1f77b4",
  ns       = "gray60"
)

DIRECTION_LABELS <- c(
  "down_in_mutant" = "Lost in BAP1-KO",
  "up_in_mutant"   = "Gained in BAP1-KO"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Load BED file into GRanges (handles 3-col and 6-col BED)
load_bed_peaks <- function(bed_path, mark_name = "ChIP") {
  if (!file.exists(bed_path)) stop(sprintf("%s bed file not found: %s", mark_name, bed_path))
  cat(sprintf("  Loading %s from: %s\n", mark_name, bed_path))

  df <- read.table(bed_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gr <- GRanges(seqnames = df$V1, ranges = IRanges(start = df$V2, end = df$V3))

  cat(sprintf("    %d peaks loaded\n", length(gr)))
  gr
}

#' Annotate loop anchors with K119ub peak overlaps
annotate_k119ub_overlaps <- function(loops_df, peaks_list) {
  # Build anchor GRanges
  anchor1_gr <- GRanges(seqnames = loops_df$chr1,
                         ranges = IRanges(start = loops_df$start1, end = loops_df$end1))
  anchor2_gr <- GRanges(seqnames = loops_df$chr2,
                         ranges = IRanges(start = loops_df$start2, end = loops_df$end2))

  for (mark in names(peaks_list)) {
    peaks_gr <- peaks_list[[mark]]
    a1_hit <- countOverlaps(anchor1_gr, peaks_gr) > 0
    a2_hit <- countOverlaps(anchor2_gr, peaks_gr) > 0

    loops_df[[paste0("anchor1_", mark)]] <- a1_hit
    loops_df[[paste0("anchor2_", mark)]] <- a2_hit
    loops_df[[paste0(mark, "_either")]] <- a1_hit | a2_hit

    n <- sum(a1_hit | a2_hit)
    cat(sprintf("    %s: %d loops (%0.1f%%) have overlap at either anchor\n",
                mark, n, 100 * n / nrow(loops_df)))
  }

  # Derive K119ub differential status (either-anchor level)
  loops_df$K119ub_diff_status <- case_when(
    loops_df$K119ub_up_either & loops_df$K119ub_down_either ~ "both",
    loops_df$K119ub_up_either ~ "gained_ub",
    loops_df$K119ub_down_either ~ "lost_ub",
    TRUE ~ "neither"
  )

  # Derive K119ub condition status
  loops_df$K119ub_condition_status <- case_when(
    loops_df$K119ub_ctrl_either & loops_df$K119ub_mut_either ~ "both",
    loops_df$K119ub_ctrl_either ~ "ctrl_only",
    loops_df$K119ub_mut_either ~ "mut_only",
    TRUE ~ "neither"
  )

  loops_df
}

#' Categorize loops into analysis groups
categorize_loops <- function(loops_df) {
  if (!"loop_distance" %in% colnames(loops_df)) {
    loops_df <- loops_df %>%
      mutate(loop_distance = abs((start2 + end2) / 2 - (start1 + end1) / 2))
  }

  loops_df %>%
    mutate(
      loop_category = factor(case_when(
        direction == "down_in_mutant" & loop_distance > DISTANCE_THRESHOLD ~ "long_range_lost",
        direction == "up_in_mutant" & loop_distance < DISTANCE_THRESHOLD ~ "short_range_gained",
        FDR >= FDR_UNCHANGED_THRESHOLD ~ "unchanged",
        TRUE ~ "other_differential"
      ), levels = c("long_range_lost", "short_range_gained", "unchanged", "other_differential"))
    )
}

#' Run Fisher's exact test for enrichment (2x2)
run_fisher_test <- function(loops_df, overlap_col, test_category, ref_category) {
  test_loops <- loops_df %>% filter(loop_category == test_category)
  ref_loops  <- loops_df %>% filter(loop_category == ref_category)

  a <- sum(test_loops[[overlap_col]])
  b <- sum(ref_loops[[overlap_col]])
  c <- nrow(test_loops) - a
  d <- nrow(ref_loops) - b

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(c("Overlap", "No_Overlap"), c(test_category, ref_category)))

  ft <- tryCatch(fisher.test(mat), error = function(e) {
    list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
  })

  tibble(
    overlap_col      = overlap_col,
    test_category    = test_category,
    ref_category     = ref_category,
    test_n           = nrow(test_loops),
    test_overlap     = a,
    test_pct         = 100 * a / max(nrow(test_loops), 1),
    ref_n            = nrow(ref_loops),
    ref_overlap      = b,
    ref_pct          = 100 * b / max(nrow(ref_loops), 1),
    odds_ratio       = as.numeric(ft$estimate),
    ci_lower         = ft$conf.int[1],
    ci_upper         = ft$conf.int[2],
    p_value          = ft$p.value
  )
}

#' Generate CDF plot comparing lost vs gained within a K119ub-filtered subset
generate_k119ub_cdf <- function(loops_df, subset_name, output_path, n_unchanged = NA) {
  lost   <- loops_df %>% filter(direction == "down_in_mutant")
  gained <- loops_df %>% filter(direction == "up_in_mutant")

  if (nrow(lost) < 5 || nrow(gained) < 5) {
    cat(sprintf("    Skipping CDF for %s: insufficient data (Lost: %d, Gained: %d)\n",
                subset_name, nrow(lost), nrow(gained)))
    return(invisible(NULL))
  }

  ks <- ks.test(lost$loop_distance, gained$loop_distance)
  med_lost <- median(lost$loop_distance)
  med_gained <- median(gained$loop_distance)

  plot_df <- loops_df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      direction_label = factor(ifelse(direction == "down_in_mutant",
                                       "Lost in BAP1-KO", "Gained in BAP1-KO"),
                                levels = c("Lost in BAP1-KO", "Gained in BAP1-KO")),
      loop_distance_kb = loop_distance / 1000
    )

  n_label <- sprintf("n = %d (Lost: %d, Gained: %d)", nrow(plot_df), nrow(lost), nrow(gained))
  if (!is.na(n_unchanged)) n_label <- paste0(n_label, sprintf(", Unchanged: %d", n_unchanged))

  p <- ggplot(plot_df, aes(x = loop_distance_kb, color = direction_label)) +
    stat_ecdf(geom = "step", linewidth = 1.2) +
    scale_color_manual(values = c("Lost in BAP1-KO" = COLORS$lost,
                                   "Gained in BAP1-KO" = COLORS$gained),
                        name = "Direction") +
    scale_x_log10(labels = comma, breaks = c(10, 100, 1000, 10000), limits = c(10, 100000)) +
    geom_vline(xintercept = med_lost / 1000, color = COLORS$lost, linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = med_gained / 1000, color = COLORS$gained, linetype = "dashed", alpha = 0.7) +
    annotate("text", x = 50000, y = 0.15,
             label = sprintf("KS p = %.2e", ks$p.value), hjust = 1, size = 4) +
    annotate("text", x = 50000, y = 0.08,
             label = sprintf("Median: Lost = %d kb, Gained = %d kb",
                             round(med_lost / 1000), round(med_gained / 1000)),
             hjust = 1, size = 3.5) +
    annotate("text", x = 50000, y = 0.22, label = n_label, hjust = 1, size = 3.5, fontface = "italic") +
    labs(title = sprintf("Loop Distance CDF: %s", subset_name),
         subtitle = ifelse(med_lost > med_gained,
                           "BAP1-KO preferentially loses long-range loops",
                           "Similar or reversed distance pattern"),
         x = "Loop Distance (kb, log scale)", y = "Cumulative Proportion") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
          legend.position = c(0.15, 0.85),
          legend.background = element_rect(fill = "white", color = NA),
          panel.grid.minor = element_blank())

  save_multiformat_ggplot(p, output_path, width = 8, height = 6)
  invisible(p)
}

#' Generate density plot comparing lost vs gained within a K119ub-filtered subset
generate_k119ub_density <- function(loops_df, subset_name, output_path) {
  lost   <- loops_df %>% filter(direction == "down_in_mutant")
  gained <- loops_df %>% filter(direction == "up_in_mutant")

  if (nrow(lost) < 5 || nrow(gained) < 5) {
    cat(sprintf("    Skipping density for %s: insufficient data (Lost: %d, Gained: %d)\n",
                subset_name, nrow(lost), nrow(gained)))
    return(invisible(NULL))
  }

  wx <- wilcox.test(lost$loop_distance, gained$loop_distance)
  med_lost <- median(lost$loop_distance)
  med_gained <- median(gained$loop_distance)

  plot_df <- loops_df %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(
      direction_label = factor(ifelse(direction == "down_in_mutant",
                                       "Lost in BAP1-KO", "Gained in BAP1-KO"),
                                levels = c("Lost in BAP1-KO", "Gained in BAP1-KO")),
      loop_distance_kb = loop_distance / 1000
    )

  p <- ggplot(plot_df, aes(x = loop_distance_kb, fill = direction_label)) +
    geom_density(alpha = 0.5, color = "black", linewidth = 0.5) +
    geom_vline(xintercept = med_lost / 1000, color = COLORS$lost, linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = med_gained / 1000, color = COLORS$gained, linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = c("Lost in BAP1-KO" = COLORS$lost,
                                  "Gained in BAP1-KO" = COLORS$gained),
                       name = "") +
    scale_x_log10(labels = comma, breaks = c(10, 100, 1000, 10000), limits = c(10, 100000)) +
    annotate("text", x = 20, y = 0,
             label = sprintf("Wilcoxon p = %.2e | n = %d", wx$p.value, nrow(plot_df)),
             hjust = 0, vjust = -0.5, size = 4) +
    labs(title = sprintf("Loop Distance Density: %s", subset_name),
         subtitle = ifelse(med_lost > med_gained,
                           "Lost loops are systematically longer | H2AK119ub anchored",
                           "Similar distance distributions | H2AK119ub anchored"),
         x = "Loop Distance (kb, log scale)", y = "Density") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
          legend.position = "top",
          panel.grid.minor = element_blank())

  save_multiformat_ggplot(p, output_path, width = 8, height = 6)
  invisible(p)
}

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

run_analysis <- function() {

  # Create output directories
  tables_dir <- file.path(OUTPUT_DIR, "tables")
  plots_dir  <- file.path(OUTPUT_DIR, "plots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

  # ==========================================================================
  # LOAD ALL DATA
  # ==========================================================================
  cat("\n[0] Loading input data...\n")

  # Peak files
  cat("  Loading K119ub peak files...\n")
  peaks <- list()
  for (name in names(PEAK_FILES)) {
    peaks[[name]] <- load_bed_peaks(PEAK_FILES[[name]], name)
  }

  # All loops (39,344)
  cat("\n  Loading all loops...\n")
  all_loops <- read.table(INPUT_FILES$all_loops, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("    All loops: %d\n", nrow(all_loops)))

  # Differential loops (2,910)
  cat("  Loading characterized (differential) loops...\n")
  diff_loops <- read.table(INPUT_FILES$diff_loops, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("    Differential loops: %d\n", nrow(diff_loops)))

  # Shared anchor data
  cat("  Loading shared anchor data...\n")
  shared_anchors <- read.table(INPUT_FILES$shared_anchors, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  shared_loops   <- read.table(INPUT_FILES$shared_loops, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  polycomb_shared <- read.table(INPUT_FILES$polycomb_shared, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("    Shared anchors: %d, Shared loops: %d, Polycomb shared: %d\n",
              nrow(shared_anchors), nrow(shared_loops), nrow(polycomb_shared)))

  # Signal data for Section D
  cat("  Loading K119ub signal data...\n")
  if (!file.exists(INPUT_FILES$signal)) {
    stop(sprintf("K119ub signal file not found: %s", INPUT_FILES$signal))
  }
  signal_data <- read.table(INPUT_FILES$signal, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("    Signal data: %d rows from %s\n", nrow(signal_data), INPUT_FILES$signal))

  # ==========================================================================
  # SECTION A: ANCHOR-LEVEL K119UB ANNOTATION
  # ==========================================================================
  cat("\n============================================================\n")
  cat("SECTION A: Anchor-Level K119ub Annotation\n")
  cat("============================================================\n")

  cat("\n  Annotating all loops with K119ub overlaps...\n")
  all_loops <- annotate_k119ub_overlaps(all_loops, peaks)
  all_loops <- categorize_loops(all_loops)

  cat("\n  Annotating differential loops with K119ub overlaps...\n")
  diff_loops <- annotate_k119ub_overlaps(diff_loops, peaks)

  # Compute loop_distance for diff_loops if missing
  if (!"loop_distance" %in% colnames(diff_loops)) {
    diff_loops$loop_distance <- abs(diff_loops$anchor2_midpoint - diff_loops$anchor1_midpoint)
  }
  if (!"loop_distance_kb" %in% colnames(diff_loops)) {
    diff_loops$loop_distance_kb <- diff_loops$loop_distance / 1000
  }

  # Report K119ub status distributions
  cat("\n  K119ub diff status (all loops):\n")
  status_table <- table(all_loops$K119ub_diff_status)
  for (s in names(status_table)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", s, status_table[s], 100 * status_table[s] / nrow(all_loops)))
  }

  cat("\n  K119ub diff status (differential loops):\n")
  status_table_diff <- table(diff_loops$K119ub_diff_status)
  for (s in names(status_table_diff)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", s, status_table_diff[s],
                100 * status_table_diff[s] / nrow(diff_loops)))
  }

  # Save annotated all loops
  write.table(all_loops, file.path(tables_dir, "loops_with_k119ub_annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n  Saved: loops_with_k119ub_annotation.tsv\n")

  # ==========================================================================
  # SECTION B: CDF/DENSITY BY K119UB STATUS
  # ==========================================================================
  cat("\n============================================================\n")
  cat("SECTION B: CDF/Density by K119ub Status\n")
  cat("============================================================\n")

  # Use diff_loops for directional CDF/density analysis
  directional <- diff_loops %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant"))
  cat(sprintf("\n  Directional loops: %d (Lost: %d, Gained: %d)\n",
              nrow(directional),
              sum(directional$direction == "down_in_mutant"),
              sum(directional$direction == "up_in_mutant")))

  # K119ub_up filtered analyses
  cat("\n  --- K119ub_up (gained ubiquitination) ---\n")

  # One-anchor
  one_up <- directional %>% filter(anchor1_K119ub_up | anchor2_K119ub_up)
  cat(sprintf("  K119ub_up one+ anchor: %d loops\n", nrow(one_up)))
  generate_k119ub_cdf(one_up, "K119ub-Gained Anchors (One+)",
                       file.path(plots_dir, "01_cdf_k119ub_up_one_anchor"))
  generate_k119ub_density(one_up, "K119ub-Gained Anchors (One+)",
                           file.path(plots_dir, "05_density_k119ub_up_one_anchor"))

  # Both-anchors
  both_up <- directional %>% filter(anchor1_K119ub_up & anchor2_K119ub_up)
  cat(sprintf("  K119ub_up both anchors: %d loops\n", nrow(both_up)))
  generate_k119ub_cdf(both_up, "K119ub-Gained Anchors (Both)",
                       file.path(plots_dir, "02_cdf_k119ub_up_both_anchors"))
  generate_k119ub_density(both_up, "K119ub-Gained Anchors (Both)",
                           file.path(plots_dir, "06_density_k119ub_up_both_anchors"))

  # K119ub_down filtered analyses
  cat("\n  --- K119ub_down (lost ubiquitination) ---\n")

  one_down <- directional %>% filter(anchor1_K119ub_down | anchor2_K119ub_down)
  cat(sprintf("  K119ub_down one+ anchor: %d loops\n", nrow(one_down)))
  generate_k119ub_cdf(one_down, "K119ub-Lost Anchors (One+)",
                       file.path(plots_dir, "03_cdf_k119ub_down_one_anchor"))
  generate_k119ub_density(one_down, "K119ub-Lost Anchors (One+)",
                           file.path(plots_dir, "07_density_k119ub_down_one_anchor"))

  both_down <- directional %>% filter(anchor1_K119ub_down & anchor2_K119ub_down)
  cat(sprintf("  K119ub_down both anchors: %d loops\n", nrow(both_down)))
  generate_k119ub_cdf(both_down, "K119ub-Lost Anchors (Both)",
                       file.path(plots_dir, "04_cdf_k119ub_down_both_anchors"))
  generate_k119ub_density(both_down, "K119ub-Lost Anchors (Both)",
                           file.path(plots_dir, "08_density_k119ub_down_both_anchors"))

  # Collect distance statistics
  subsets <- list(
    list(name = "K119ub_up_one",  df = one_up),
    list(name = "K119ub_up_both", df = both_up),
    list(name = "K119ub_down_one",  df = one_down),
    list(name = "K119ub_down_both", df = both_down)
  )

  dist_stats <- bind_rows(lapply(subsets, function(s) {
    lost <- s$df %>% filter(direction == "down_in_mutant")
    gained <- s$df %>% filter(direction == "up_in_mutant")

    ks_p <- if (nrow(lost) >= 2 && nrow(gained) >= 2) {
      tryCatch(ks.test(lost$loop_distance, gained$loop_distance)$p.value, error = function(e) NA)
    } else NA
    wx_p <- if (nrow(lost) >= 2 && nrow(gained) >= 2) {
      tryCatch(wilcox.test(lost$loop_distance, gained$loop_distance)$p.value, error = function(e) NA)
    } else NA

    tibble(
      subset = s$name, total = nrow(s$df),
      n_lost = nrow(lost), n_gained = nrow(gained),
      median_lost_kb = ifelse(nrow(lost) > 0, round(median(lost$loop_distance) / 1000, 1), NA),
      median_gained_kb = ifelse(nrow(gained) > 0, round(median(gained$loop_distance) / 1000, 1), NA),
      ks_pvalue = ks_p, wilcox_pvalue = wx_p
    )
  }))

  write.table(dist_stats, file.path(tables_dir, "distance_stats_by_k119ub.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n  Saved: distance_stats_by_k119ub.tsv\n")

  # ==========================================================================
  # SECTION C: ENRICHMENT TESTING
  # ==========================================================================
  cat("\n============================================================\n")
  cat("SECTION C: Enrichment Testing by Loop Category x Chromatin State\n")
  cat("============================================================\n")

  # --- C1: Global enrichment (all loops, vs unchanged) ---
  cat("\n  --- C1: Global K119ub enrichment by loop category ---\n")

  overlap_cols <- c("K119ub_up_either", "K119ub_down_either",
                     "K119ub_ctrl_either", "K119ub_mut_either")
  comparisons <- list(
    c("long_range_lost", "unchanged"),
    c("short_range_gained", "unchanged"),
    c("long_range_lost", "short_range_gained")
  )

  global_results <- bind_rows(lapply(overlap_cols, function(col) {
    bind_rows(lapply(comparisons, function(comp) {
      run_fisher_test(all_loops, col, comp[1], comp[2])
    }))
  }))
  global_results$fdr <- p.adjust(global_results$p_value, method = "BH")
  global_results$significant <- global_results$fdr < 0.05
  global_results$enrichment <- case_when(
    global_results$odds_ratio > 1 & global_results$significant ~ "enriched",
    global_results$odds_ratio < 1 & global_results$significant ~ "depleted",
    TRUE ~ "ns"
  )

  cat(sprintf("  Global: %d tests, %d significant\n",
              nrow(global_results), sum(global_results$significant, na.rm = TRUE)))

  # --- C2: Chromatin-state stratified enrichment (diff loops only) ---
  cat("\n  --- C2: K119ub enrichment stratified by chromatin state ---\n")

  # For chromatin-state stratification, use diff_loops which have anchor types
  anchor_types <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                     "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Other")

  # Assign directional category to diff_loops for testing
  diff_loops$dir_binary <- ifelse(diff_loops$direction == "down_in_mutant", "lost", "gained")
  diff_directional <- diff_loops %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant"))

  chromstate_results <- bind_rows(lapply(anchor_types, function(atype) {
    # Loops where either anchor has this chromatin type
    subset <- diff_directional %>%
      filter(anchor1_type == atype | anchor2_type == atype)

    if (nrow(subset) < MIN_SAMPLES) {
      return(tibble(chromatin_state = atype, overlap_col = NA_character_,
                     n_total = nrow(subset), note = "insufficient_data"))
    }

    bind_rows(lapply(c("K119ub_up_either", "K119ub_down_either"), function(col) {
      # 2x2: lost vs gained × K119ub overlap
      lost <- subset %>% filter(dir_binary == "lost")
      gained <- subset %>% filter(dir_binary == "gained")

      a <- sum(lost[[col]])
      b <- sum(gained[[col]])
      c <- nrow(lost) - a
      d <- nrow(gained) - b

      mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      ft <- tryCatch(fisher.test(mat), error = function(e) {
        list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
      })

      tibble(
        chromatin_state = atype,
        overlap_col     = col,
        n_total         = nrow(subset),
        n_lost          = nrow(lost),
        n_gained        = nrow(gained),
        lost_overlap    = a,
        lost_pct        = 100 * a / max(nrow(lost), 1),
        gained_overlap  = b,
        gained_pct      = 100 * b / max(nrow(gained), 1),
        odds_ratio      = as.numeric(ft$estimate),
        ci_lower        = ft$conf.int[1],
        ci_upper        = ft$conf.int[2],
        p_value         = ft$p.value
      )
    }))
  }))

  # Filter out insufficient data rows and apply FDR
  chromstate_results <- chromstate_results %>% filter(!is.na(overlap_col))
  if (nrow(chromstate_results) > 0) {
    chromstate_results$fdr <- p.adjust(chromstate_results$p_value, method = "BH")
    chromstate_results$significant <- chromstate_results$fdr < 0.05
  }

  cat(sprintf("  Chromatin-state: %d tests, %d significant\n",
              nrow(chromstate_results),
              sum(chromstate_results$significant, na.rm = TRUE)))

  # --- C3: Distance-stratified enrichment (all loops) ---
  cat("\n  --- C3: K119ub enrichment stratified by distance ---\n")

  distance_results <- bind_rows(lapply(names(DISTANCE_BINS), function(bin_name) {
    bounds <- DISTANCE_BINS[[bin_name]]
    subset <- all_loops %>%
      filter(loop_distance >= bounds[1], loop_distance < bounds[2])

    if (nrow(subset) < MIN_SAMPLES) {
      return(tibble(distance_bin = bin_name, note = "insufficient_data"))
    }

    bind_rows(lapply(c("K119ub_up_either", "K119ub_down_either"), function(col) {
      bind_rows(lapply(list(c("long_range_lost", "unchanged"),
                             c("short_range_gained", "unchanged")), function(comp) {
        test_sub <- subset %>% filter(loop_category == comp[1])
        ref_sub  <- subset %>% filter(loop_category == comp[2])

        if (nrow(test_sub) < 5 || nrow(ref_sub) < 5) {
          return(tibble(distance_bin = bin_name, overlap_col = col,
                         test_category = comp[1], ref_category = comp[2],
                         note = "insufficient_data"))
        }

        result <- run_fisher_test(subset, col, comp[1], comp[2])
        result$distance_bin <- bin_name
        result
      }))
    }))
  }))

  distance_results <- distance_results %>% filter(is.na(note) | note != "insufficient_data")
  if (nrow(distance_results) > 0 && "p_value" %in% colnames(distance_results)) {
    valid <- !is.na(distance_results$p_value)
    distance_results$fdr <- NA_real_
    distance_results$fdr[valid] <- p.adjust(distance_results$p_value[valid], method = "BH")
    distance_results$significant <- distance_results$fdr < 0.05
  }

  cat(sprintf("  Distance-stratified: %d tests\n", nrow(distance_results)))

  # Save enrichment tables
  write.table(global_results, file.path(tables_dir, "enrichment_global.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(chromstate_results, file.path(tables_dir, "enrichment_by_chromatin_state.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(distance_results, file.path(tables_dir, "enrichment_by_distance.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # --- C Plots ---
  cat("\n  Generating enrichment plots...\n")

  # Plot 09: Chromatin-state enrichment dotplot
  if (nrow(chromstate_results) > 0 && any(!is.na(chromstate_results$odds_ratio))) {
    plot_c1 <- chromstate_results %>%
      filter(!is.na(odds_ratio)) %>%
      mutate(
        mark_label = ifelse(grepl("up", overlap_col), "K119ub gained", "K119ub lost"),
        log2_or = log2(pmax(odds_ratio, 0.01)),
        chromatin_state = factor(chromatin_state, levels = rev(anchor_types))
      )

    p09 <- ggplot(plot_c1, aes(x = log2_or, y = chromatin_state)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(aes(size = -log10(fdr + 1e-10), color = significant), alpha = 0.8) +
      geom_errorbarh(aes(xmin = log2(pmax(ci_lower, 0.001)), xmax = log2(pmin(ci_upper, 1000))),
                     height = 0.2, alpha = 0.5) +
      facet_wrap(~mark_label, ncol = 1) +
      scale_color_manual(values = c("TRUE" = COLORS$enriched, "FALSE" = COLORS$ns),
                          labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS")) +
      scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
      labs(title = "K119ub Enrichment at Lost vs Gained Loops",
           subtitle = "Stratified by anchor chromatin state",
           x = "log2(Odds Ratio)", y = "", color = "Significance") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14),
            strip.background = element_rect(fill = "gray90"),
            strip.text = element_text(face = "bold"))

    save_multiformat_ggplot(p09, file.path(plots_dir, "09_enrichment_dotplot_by_chromatin_state"),
                             width = 9, height = 8)
  }

  # Plot 10: K119ub overlap % heatmap (loop direction x chromatin state)
  if (nrow(chromstate_results) > 0) {
    heatmap_data <- chromstate_results %>%
      filter(!is.na(odds_ratio), grepl("up", overlap_col)) %>%
      select(chromatin_state, lost_pct, gained_pct) %>%
      pivot_longer(cols = c(lost_pct, gained_pct),
                    names_to = "direction", values_to = "pct") %>%
      mutate(direction = ifelse(direction == "lost_pct", "Lost", "Gained"))

    if (nrow(heatmap_data) > 0) {
      mat <- heatmap_data %>%
        pivot_wider(names_from = direction, values_from = pct) %>%
        column_to_rownames("chromatin_state") %>%
        as.matrix()

      heatmap_dir <- file.path(plots_dir, "10_enrichment_heatmap_category_x_chromstate")
      dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

      for (fmt in c("pdf", "svg", "jpg")) {
        fpath <- file.path(heatmap_dir, paste0("10_enrichment_heatmap_category_x_chromstate.", fmt))
        if (fmt == "pdf") pdf(fpath, width = 6, height = 6)
        else if (fmt == "svg") svglite::svglite(fpath, width = 6, height = 6)
        else jpeg(fpath, width = 6 * 300, height = 6 * 300, res = 300, quality = 95)
        pheatmap(mat, main = "% K119ub-Gained Overlap by Anchor Type",
                 color = colorRampPalette(c("white", "#fee8c8", "#e34a33"))(50),
                 display_numbers = TRUE, number_format = "%.1f", number_color = "black",
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 border_color = "gray80", cellwidth = 60, cellheight = 30, angle_col = 45)
        dev.off()
      }
      cat("  Saved: 10_enrichment_heatmap_category_x_chromstate/{pdf,svg,jpg}\n")
    }
  }

  # Plot 11: Distance-stratified enrichment dotplot
  if (nrow(distance_results) > 0 && "odds_ratio" %in% colnames(distance_results)) {
    plot_c3 <- distance_results %>%
      filter(!is.na(odds_ratio)) %>%
      mutate(
        mark_label = ifelse(grepl("up", overlap_col), "K119ub gained", "K119ub lost"),
        log2_or = log2(pmax(odds_ratio, 0.01)),
        comparison = paste(test_category, "vs", ref_category),
        distance_bin = factor(distance_bin, levels = names(DISTANCE_BINS))
      )

    if (nrow(plot_c3) > 0) {
      p11 <- ggplot(plot_c3, aes(x = log2_or, y = comparison)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(aes(size = -log10(fdr + 1e-10), color = significant), alpha = 0.8) +
        geom_errorbarh(aes(xmin = log2(pmax(ci_lower, 0.001)), xmax = log2(pmin(ci_upper, 1000))),
                       height = 0.2, alpha = 0.5) +
        facet_grid(distance_bin ~ mark_label) +
        scale_color_manual(values = c("TRUE" = COLORS$enriched, "FALSE" = COLORS$ns)) +
        scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
        labs(title = "K119ub Enrichment by Distance Bin",
             x = "log2(Odds Ratio)", y = "", color = "Significance") +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 14),
              strip.background = element_rect(fill = "gray90"),
              strip.text = element_text(face = "bold", size = 9))

      save_multiformat_ggplot(p11, file.path(plots_dir, "11_enrichment_dotplot_by_distance"),
                               width = 10, height = 8)
    }
  }

  # ==========================================================================
  # SECTION D: CONTINUOUS SIGNAL CORRELATION
  # ==========================================================================
  cat("\n============================================================\n")
  cat("SECTION D: Continuous Signal Correlation\n")
  cat("============================================================\n")

  # Join K119ub anchor signal by coordinate match (chr:start-end)
  anchor_fc <- signal_data %>%
    filter(signal_class == "quantifiable") %>%
    select(anchor_id, log2fc) %>%
    distinct(anchor_id, .keep_all = TRUE)

  diff_with_signal <- diff_directional %>%
    mutate(
      anchor1_id = paste0(anchor1_chr, ":", anchor1_start, "-", anchor1_end),
      anchor2_id = paste0(anchor2_chr, ":", anchor2_start, "-", anchor2_end)
    ) %>%
    left_join(anchor_fc %>% rename(anchor1_k119ub_fc = log2fc),
              by = c("anchor1_id" = "anchor_id")) %>%
    left_join(anchor_fc %>% rename(anchor2_k119ub_fc = log2fc),
              by = c("anchor2_id" = "anchor_id"))

  diff_with_signal$mean_anchor_k119ub_fc <- rowMeans(
    cbind(diff_with_signal$anchor1_k119ub_fc, diff_with_signal$anchor2_k119ub_fc),
    na.rm = TRUE
  )

  corr_results <- tibble()
  logistic_results <- tibble()

  {
    # Filter to loops with quantifiable signal
    sig_loops <- diff_with_signal %>%
      filter(!is.na(mean_anchor_k119ub_fc) & is.finite(mean_anchor_k119ub_fc))
    cat(sprintf("  Loops with K119ub signal: %d / %d\n", nrow(sig_loops), nrow(diff_with_signal)))

    if (nrow(sig_loops) >= MIN_SAMPLES) {

      # Plot 12: Main scatter (loop logFC vs K119ub FC)
      cat("\n  Generating correlation scatter...\n")

      # Determine dominant anchor type for coloring
      sig_loops$dominant_anchor_type <- ifelse(
        sig_loops$anchor1_type %in% c("Polycomb", "Repressed_Promoter", "Bivalent_Promoter"),
        sig_loops$anchor1_type,
        ifelse(sig_loops$anchor2_type %in% c("Polycomb", "Repressed_Promoter", "Bivalent_Promoter"),
               sig_loops$anchor2_type,
               ifelse(sig_loops$anchor1_type %in% c("Active_Enhancer", "Active_Promoter"),
                      sig_loops$anchor1_type, sig_loops$anchor2_type)))

      # Simplify to 3 groups for coloring
      sig_loops$anchor_group <- case_when(
        sig_loops$dominant_anchor_type %in% c("Polycomb", "Repressed_Promoter", "Bivalent_Promoter") ~ "Polycomb/Repressive",
        sig_loops$dominant_anchor_type %in% c("Active_Enhancer", "Active_Promoter") ~ "Active",
        TRUE ~ "Other"
      )

      rho <- cor.test(sig_loops$mean_anchor_k119ub_fc, sig_loops$logFC, method = "spearman")

      p12 <- ggplot(sig_loops, aes(x = mean_anchor_k119ub_fc, y = logFC)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
        geom_point(aes(color = anchor_group), alpha = 0.4, size = 1.5) +
        geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
        scale_color_manual(values = c("Polycomb/Repressive" = "#756bb1",
                                       "Active" = "#e6550d",
                                       "Other" = "gray50"),
                            name = "Anchor Type") +
        annotate("text", x = Inf, y = Inf,
                 label = sprintf("Spearman rho = %.3f\np = %.2e\nn = %d (%s)",
                                 rho$estimate, rho$p.value, nrow(sig_loops), INPUT_FILES$signal),
                 hjust = 1.1, vjust = 1.5, size = 4) +
        labs(title = "Loop Contact Change vs K119ub Signal Change",
             subtitle = "BAP1-KO: H2AK119ub accumulation at anchors vs loop strength",
             x = "Mean Anchor K119ub log2FC (mut/ctrl)",
             y = "Loop logFC (positive = gained in BAP1-KO)") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14),
              plot.subtitle = element_text(size = 10, color = "grey40"),
              legend.position = c(0.15, 0.15),
              legend.background = element_rect(fill = alpha("white", 0.8)))

      save_multiformat_ggplot(p12, file.path(plots_dir, "12_scatter_loopFC_vs_k119ub_FC"),
                               width = 9, height = 8)

      # Plot 13: Faceted by chromatin state
      if (length(unique(sig_loops$anchor_group)) > 1) {
        facet_corrs <- sig_loops %>%
          group_by(anchor_group) %>%
          filter(n() >= 10) %>%
          ungroup()

        if (nrow(facet_corrs) > 0) {
          p13 <- ggplot(facet_corrs, aes(x = mean_anchor_k119ub_fc, y = logFC)) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
            geom_point(aes(color = anchor_group), alpha = 0.4, size = 1.5) +
            geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
            facet_wrap(~anchor_group, scales = "free") +
            scale_color_manual(values = c("Polycomb/Repressive" = "#756bb1",
                                           "Active" = "#e6550d",
                                           "Other" = "gray50")) +
            labs(title = "Loop logFC vs K119ub Change by Anchor Type",
                 x = "Mean Anchor K119ub log2FC", y = "Loop logFC") +
            theme_bw(base_size = 11) +
            theme(plot.title = element_text(face = "bold", size = 14),
                  strip.background = element_rect(fill = "gray90"),
                  strip.text = element_text(face = "bold"),
                  legend.position = "none")

          save_multiformat_ggplot(p13, file.path(plots_dir, "13_scatter_by_chromatin_state"),
                                   width = 12, height = 5)
        }
      }

      # Plot 14: Boxplot of H2AK119ub signal by loop direction
      n_lost   <- sum(sig_loops$direction == "down_in_mutant")
      n_gained <- sum(sig_loops$direction == "up_in_mutant")
      direction_labels <- c(
        "down_in_mutant" = sprintf("Lost\n(n = %s)",   format(n_lost,   big.mark = ",")),
        "up_in_mutant"   = sprintf("Gained\n(n = %s)", format(n_gained, big.mark = ","))
      )

      p14 <- ggplot(sig_loops, aes(x = factor(direction,
                                               levels = c("down_in_mutant", "up_in_mutant"),
                                               labels = direction_labels),
                                    y = mean_anchor_k119ub_fc, fill = direction)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = c("down_in_mutant" = COLORS$lost,
                                      "up_in_mutant" = COLORS$gained)) +
        stat_summary(fun = median, geom = "text",
                     aes(label = sprintf("%.3f", after_stat(y))),
                     vjust = -0.5, size = 3.5) +
        labs(title = "H2AK119ub Signal Change by Loop Direction",
             subtitle = sprintf("Wilcoxon p = %.2e   |   total n = %s loops",
                                tryCatch(wilcox.test(mean_anchor_k119ub_fc ~ direction,
                                                      data = sig_loops)$p.value,
                                         error = function(e) NA),
                                format(n_lost + n_gained, big.mark = ",")),
             x = "Loop Direction", y = "Mean Anchor H2AK119ub log2FC") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 14),
              legend.position = "none")

      save_multiformat_ggplot(p14, file.path(plots_dir, "14_boxplot_k119ub_by_loop_direction"),
                               width = 6, height = 7)

      # Stratified correlations for summary
      corr_results <- sig_loops %>%
        group_by(anchor_group) %>%
        filter(n() >= 10) %>%
        summarize(
          n = n(),
          spearman_rho = cor(mean_anchor_k119ub_fc, logFC, method = "spearman", use = "complete.obs"),
          spearman_p = tryCatch(
            cor.test(mean_anchor_k119ub_fc, logFC, method = "spearman")$p.value,
            error = function(e) NA
          ),
          .groups = "drop"
        )

      # Plot 15: Correlation summary heatmap
      if (nrow(corr_results) > 1) {
        # Add distance-stratified correlations
        sig_loops$dist_bin <- cut(sig_loops$loop_distance,
                                   breaks = c(0, 100000, 500000, 1000000, Inf),
                                   labels = names(DISTANCE_BINS))

        dist_corrs <- sig_loops %>%
          group_by(dist_bin, anchor_group) %>%
          filter(n() >= 10) %>%
          summarize(
            spearman_rho = cor(mean_anchor_k119ub_fc, logFC, method = "spearman", use = "complete.obs"),
            n = n(),
            .groups = "drop"
          )

        if (nrow(dist_corrs) > 0) {
          corr_mat <- dist_corrs %>%
            select(-n) %>%
            pivot_wider(names_from = anchor_group, values_from = spearman_rho) %>%
            column_to_rownames("dist_bin") %>%
            as.matrix()

          heatmap_dir_15 <- file.path(plots_dir, "15_correlation_summary_heatmap")
          dir.create(heatmap_dir_15, recursive = TRUE, showWarnings = FALSE)

          for (fmt in c("pdf", "svg", "jpg")) {
            fpath <- file.path(heatmap_dir_15, paste0("15_correlation_summary_heatmap.", fmt))
            if (fmt == "pdf") pdf(fpath, width = 7, height = 5)
            else if (fmt == "svg") svglite::svglite(fpath, width = 7, height = 5)
            else jpeg(fpath, width = 7 * 300, height = 5 * 300, res = 300, quality = 95)
            pheatmap(corr_mat,
                     main = "Spearman Correlation: Loop logFC vs K119ub Change",
                     color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(50),
                     breaks = seq(-0.5, 0.5, length.out = 51),
                     display_numbers = TRUE, number_format = "%.2f",
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     border_color = "gray80", cellwidth = 60, cellheight = 30)
            dev.off()
          }
          cat("  Saved: 15_correlation_summary_heatmap/{pdf,svg,jpg}\n")
        }
      }

      write.table(corr_results, file.path(tables_dir, "correlation_results.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      # "Ubiquitination buffer" logistic regression
      cat("\n  Testing ubiquitination buffer hypothesis...\n")
      sig_loops$is_lost <- as.integer(sig_loops$direction == "down_in_mutant")
      sig_loops$log_distance <- log10(sig_loops$loop_distance)

      logistic_fit <- tryCatch({
        glm(is_lost ~ mean_anchor_k119ub_fc + log_distance, data = sig_loops, family = binomial)
      }, error = function(e) {
        cat(sprintf("    Logistic regression failed: %s\n", e$message))
        NULL
      })

      if (!is.null(logistic_fit)) {
        logistic_summary <- summary(logistic_fit)
        logistic_results <- as_tibble(logistic_summary$coefficients, rownames = "term")
        colnames(logistic_results) <- c("term", "estimate", "std_error", "z_value", "p_value")
        logistic_results$odds_ratio <- exp(logistic_results$estimate)
        logistic_results$signal_file <- INPUT_FILES$signal

        write.table(logistic_results, file.path(tables_dir, "logistic_regression_results.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        cat(sprintf("    K119ub coefficient: %.4f (p = %.4e, OR = %.3f)\n",
                    logistic_results$estimate[logistic_results$term == "mean_anchor_k119ub_fc"],
                    logistic_results$p_value[logistic_results$term == "mean_anchor_k119ub_fc"],
                    logistic_results$odds_ratio[logistic_results$term == "mean_anchor_k119ub_fc"]))
      }
    }
  }

  # ==========================================================================
  # SECTION E: SHARED ANCHOR K119UB ANALYSIS
  # ==========================================================================
  cat("\n============================================================\n")
  cat("SECTION E: Shared Anchor K119ub Analysis\n")
  cat("============================================================\n")

  # Annotate shared anchor loops with K119ub
  shared_loops <- annotate_k119ub_overlaps(shared_loops, peaks)
  polycomb_shared <- annotate_k119ub_overlaps(polycomb_shared, peaks)

  # Build GRanges for shared anchors
  shared_anchor_gr <- GRanges(
    seqnames = shared_anchors$chr,
    ranges = IRanges(start = shared_anchors$start, end = shared_anchors$end)
  )

  # K119ub overlap at shared anchors
  for (mark in names(peaks)) {
    shared_anchors[[paste0(mark, "_overlap")]] <- countOverlaps(shared_anchor_gr, peaks[[mark]]) > 0
  }

  # Build non-shared anchor set from diff loops for comparison
  all_anchor1 <- diff_loops %>%
    select(chr = anchor1_chr, start = anchor1_start, end = anchor1_end) %>% distinct()
  all_anchor2 <- diff_loops %>%
    select(chr = anchor2_chr, start = anchor2_start, end = anchor2_end) %>% distinct()
  all_anchors <- bind_rows(all_anchor1, all_anchor2) %>% distinct()

  all_anchor_gr <- GRanges(seqnames = all_anchors$chr,
                            ranges = IRanges(start = all_anchors$start, end = all_anchors$end))

  # Identify non-shared anchors (no overlap with shared_anchor_gr)
  is_shared <- countOverlaps(all_anchor_gr, shared_anchor_gr) > 0
  nonshared_anchor_gr <- all_anchor_gr[!is_shared]

  for (mark in names(peaks)) {
    all_anchors[[paste0(mark, "_overlap")]] <- countOverlaps(all_anchor_gr, peaks[[mark]]) > 0
  }
  all_anchors$is_shared <- is_shared

  cat(sprintf("\n  Total unique anchors: %d (Shared: %d, Non-shared: %d)\n",
              nrow(all_anchors), sum(is_shared), sum(!is_shared)))

  # E1: Compare K119ub overlap rates at shared vs non-shared
  cat("\n  --- E1: Shared vs non-shared K119ub overlap ---\n")

  shared_comparison <- tibble()
  for (mark in c("K119ub_up", "K119ub_down", "K119ub_ctrl", "K119ub_mut")) {
    col <- paste0(mark, "_overlap")
    shared_n <- sum(all_anchors[[col]][all_anchors$is_shared])
    nonshared_n <- sum(all_anchors[[col]][!all_anchors$is_shared])
    shared_tot <- sum(all_anchors$is_shared)
    nonshared_tot <- sum(!all_anchors$is_shared)

    mat <- matrix(c(shared_n, nonshared_n, shared_tot - shared_n, nonshared_tot - nonshared_n),
                  nrow = 2, byrow = TRUE)
    ft <- tryCatch(fisher.test(mat), error = function(e) {
      list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
    })

    shared_comparison <- bind_rows(shared_comparison, tibble(
      mark = mark,
      shared_n = shared_n, shared_pct = 100 * shared_n / shared_tot,
      nonshared_n = nonshared_n, nonshared_pct = 100 * nonshared_n / nonshared_tot,
      odds_ratio = as.numeric(ft$estimate),
      p_value = ft$p.value
    ))
  }
  shared_comparison$fdr <- p.adjust(shared_comparison$p_value, method = "BH")

  write.table(shared_comparison, file.path(tables_dir, "shared_anchor_k119ub_comparison.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("  Shared vs non-shared K119ub enrichment:\n")
  for (i in seq_len(nrow(shared_comparison))) {
    r <- shared_comparison[i, ]
    cat(sprintf("    %s: shared=%.1f%%, non-shared=%.1f%%, OR=%.2f, FDR=%.3e\n",
                r$mark, r$shared_pct, r$nonshared_pct, r$odds_ratio, r$fdr))
  }

  # Plot 16: K119ub overlap at shared vs non-shared (grouped bar)
  plot_e1 <- all_anchors %>%
    mutate(anchor_type = ifelse(is_shared, "Shared", "Non-shared")) %>%
    select(anchor_type, K119ub_up_overlap, K119ub_down_overlap) %>%
    pivot_longer(-anchor_type, names_to = "mark", values_to = "has_overlap") %>%
    mutate(mark = gsub("_overlap", "", mark)) %>%
    group_by(anchor_type, mark) %>%
    summarize(pct = 100 * mean(has_overlap), n = sum(has_overlap), .groups = "drop")

  p16 <- ggplot(plot_e1, aes(x = anchor_type, y = pct, fill = mark)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = n), position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("K119ub_up" = COLORS$K119ub_up, "K119ub_down" = COLORS$K119ub_down),
                       labels = c("K119ub gained", "K119ub lost")) +
    labs(title = "K119ub Overlap: Shared vs Non-shared Anchors",
         subtitle = "Shared = loop-switching sites (lost + gained at same anchor)",
         x = "", y = "% Anchors with K119ub Overlap", fill = "DiffBind Peak") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14))

  save_multiformat_ggplot(p16, file.path(plots_dir, "16_k119ub_shared_vs_nonshared"),
                           width = 7, height = 6)

  # E2: K119ub at shared anchor loops — lost-end vs gained-end
  cat("\n  --- E2: K119ub at lost-end vs gained-end of shared loops ---\n")

  shared_loop_summary <- shared_loops %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(loop_direction = ifelse(direction == "down_in_mutant", "Lost", "Gained")) %>%
    group_by(loop_direction) %>%
    summarize(
      n = n(),
      pct_K119ub_up_anchor1 = 100 * mean(anchor1_K119ub_up),
      pct_K119ub_up_anchor2 = 100 * mean(anchor2_K119ub_up),
      pct_K119ub_up_either  = 100 * mean(K119ub_up_either),
      pct_K119ub_down_either = 100 * mean(K119ub_down_either),
      .groups = "drop"
    )

  # Plot 17: Paired comparison at shared loops
  paired_data <- shared_loops %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    select(loop_id, direction, K119ub_up_either, K119ub_down_either, K119ub_ctrl_either, K119ub_mut_either) %>%
    pivot_longer(cols = c(K119ub_up_either, K119ub_down_either),
                  names_to = "mark", values_to = "has_overlap") %>%
    mutate(
      direction_label = ifelse(direction == "down_in_mutant", "Lost Loop", "Gained Loop"),
      mark_label = ifelse(mark == "K119ub_up_either", "K119ub gained", "K119ub lost")
    ) %>%
    group_by(direction_label, mark_label) %>%
    summarize(pct = 100 * mean(has_overlap), n = sum(has_overlap), total = n(), .groups = "drop")

  if (nrow(paired_data) > 0) {
    p17 <- ggplot(paired_data, aes(x = direction_label, y = pct, fill = mark_label)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_text(aes(label = sprintf("%d/%d", n, total)),
                position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
      scale_fill_manual(values = c("K119ub gained" = COLORS$K119ub_up,
                                    "K119ub lost" = COLORS$K119ub_down)) +
      labs(title = "K119ub Overlap at Shared Anchor Loop Ends",
           subtitle = sprintf("n = %d shared anchor loops", nrow(shared_loops)),
           x = "Loop Direction", y = "% Loops with K119ub Overlap", fill = "") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14))

    save_multiformat_ggplot(p17, file.path(plots_dir, "17_paired_k119ub_at_shared_anchors"),
                             width = 7, height = 6)
  }

  # E3: Polycomb shared subset
  cat("\n  --- E3: Polycomb-specific shared anchor K119ub ---\n")

  polycomb_k119ub <- polycomb_shared %>%
    filter(direction %in% c("down_in_mutant", "up_in_mutant")) %>%
    mutate(direction_label = ifelse(direction == "down_in_mutant", "Lost", "Gained"))

  if (nrow(polycomb_k119ub) >= MIN_SAMPLES) {
    pc_summary <- polycomb_k119ub %>%
      group_by(direction_label) %>%
      summarize(
        n = n(),
        pct_K119ub_up = 100 * mean(K119ub_up_either),
        pct_K119ub_down = 100 * mean(K119ub_down_either),
        .groups = "drop"
      ) %>%
      pivot_longer(cols = c(pct_K119ub_up, pct_K119ub_down),
                    names_to = "mark", values_to = "pct") %>%
      mutate(mark_label = ifelse(mark == "pct_K119ub_up", "K119ub gained", "K119ub lost"))

    p18 <- ggplot(pc_summary, aes(x = direction_label, y = pct, fill = mark_label)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      scale_fill_manual(values = c("K119ub gained" = COLORS$K119ub_up,
                                    "K119ub lost" = COLORS$K119ub_down)) +
      labs(title = "K119ub at Polycomb Shared Anchor Loops",
           subtitle = sprintf("n = %d Polycomb shared loops", nrow(polycomb_k119ub)),
           x = "Loop Direction", y = "% Loops with K119ub Overlap", fill = "") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14))

    save_multiformat_ggplot(p18, file.path(plots_dir, "18_polycomb_shared_k119ub"),
                             width = 7, height = 6)
  } else {
    cat("    Insufficient Polycomb shared loops for plotting\n")
  }

  # E4: Shared anchor K119ub direction × switching pattern
  cat("\n  --- E4: K119ub status × switching pattern ---\n")

  shared_anchors$k119ub_direction <- case_when(
    shared_anchors$K119ub_up_overlap & shared_anchors$K119ub_down_overlap ~ "Both",
    shared_anchors$K119ub_up_overlap ~ "Gaining K119ub",
    shared_anchors$K119ub_down_overlap ~ "Losing K119ub",
    TRUE ~ "No Change"
  )

  shared_anchors$switching_bias <- case_when(
    shared_anchors$n_lost_loops > shared_anchors$n_gained_loops ~ "More Lost",
    shared_anchors$n_gained_loops > shared_anchors$n_lost_loops ~ "More Gained",
    TRUE ~ "Balanced"
  )

  switch_data <- shared_anchors %>%
    count(k119ub_direction, switching_bias) %>%
    group_by(k119ub_direction) %>%
    mutate(pct = 100 * n / sum(n))

  if (nrow(switch_data) > 0) {
    p19 <- ggplot(switch_data, aes(x = k119ub_direction, y = n, fill = switching_bias)) +
      geom_col(position = "fill") +
      scale_fill_manual(values = c("More Lost" = COLORS$lost,
                                    "More Gained" = COLORS$gained,
                                    "Balanced" = COLORS$unchanged)) +
      scale_y_continuous(labels = percent) +
      labs(title = "Shared Anchor Switching Pattern by K119ub Direction",
           subtitle = sprintf("n = %d shared anchors", nrow(shared_anchors)),
           x = "K119ub Status at Anchor", y = "Proportion", fill = "Switching Bias") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 14))

    save_multiformat_ggplot(p19, file.path(plots_dir, "19_shared_anchor_k119ub_direction"),
                             width = 8, height = 6)
  }

  # ==========================================================================
  # WRITE ANALYSIS SUMMARY
  # ==========================================================================
  cat("\n============================================================\n")
  cat("Writing analysis summary...\n")
  cat("============================================================\n")

  summary_lines <- c(
    "================================================================================",
    "H2AK119ub Loop Integration Analysis Summary",
    "================================================================================",
    sprintf("Date: %s", Sys.time()),
    sprintf("Signal source: %s", INPUT_FILES$signal),
    "",
    "INPUT FILES:",
    sprintf("  All loops: %s (%d)", INPUT_FILES$all_loops, nrow(all_loops)),
    sprintf("  Differential loops: %s (%d)", INPUT_FILES$diff_loops, nrow(diff_loops)),
    sprintf("  Shared anchors: %d, Shared loops: %d, Polycomb shared: %d",
            nrow(shared_anchors), nrow(shared_loops), nrow(polycomb_shared)),
    "",
    "PEAK FILES:",
    sprintf("  K119ub_up: %d peaks", length(peaks$K119ub_up)),
    sprintf("  K119ub_down: %d peaks", length(peaks$K119ub_down)),
    sprintf("  K119ub_ctrl: %d peaks", length(peaks$K119ub_ctrl)),
    sprintf("  K119ub_mut: %d peaks", length(peaks$K119ub_mut)),
    "",
    "SECTION A: K119ub Annotation",
    "----------------------------"
  )

  for (s in names(table(all_loops$K119ub_diff_status))) {
    n <- sum(all_loops$K119ub_diff_status == s)
    summary_lines <- c(summary_lines,
                        sprintf("  %s: %d (%.1f%%)", s, n, 100 * n / nrow(all_loops)))
  }

  summary_lines <- c(summary_lines, "",
    "SECTION B: Distance Analysis",
    "----------------------------"
  )
  for (i in seq_len(nrow(dist_stats))) {
    r <- dist_stats[i, ]
    summary_lines <- c(summary_lines,
      sprintf("  %s: n=%d (Lost:%d, Gained:%d), Med Lost=%.0f kb, Med Gained=%.0f kb, KS p=%.2e",
              r$subset, r$total, r$n_lost, r$n_gained,
              ifelse(is.na(r$median_lost_kb), 0, r$median_lost_kb),
              ifelse(is.na(r$median_gained_kb), 0, r$median_gained_kb),
              ifelse(is.na(r$ks_pvalue), 1, r$ks_pvalue)))
  }

  summary_lines <- c(summary_lines, "",
    "SECTION C: Enrichment Testing",
    "-----------------------------",
    sprintf("  Global tests: %d (%d significant)",
            nrow(global_results), sum(global_results$significant, na.rm = TRUE)),
    sprintf("  Chromatin-state tests: %d (%d significant)",
            nrow(chromstate_results), sum(chromstate_results$significant, na.rm = TRUE)),
    sprintf("  Distance-stratified tests: %d", nrow(distance_results))
  )

  if (nrow(global_results) > 0 && any(global_results$significant, na.rm = TRUE)) {
    summary_lines <- c(summary_lines, "", "  Significant global enrichments:")
    sig <- global_results %>% filter(significant) %>% arrange(fdr)
    for (i in seq_len(min(nrow(sig), 10))) {
      r <- sig[i, ]
      summary_lines <- c(summary_lines,
        sprintf("    %s: %s vs %s, OR=%.2f, FDR=%.2e",
                r$overlap_col, r$test_category, r$ref_category, r$odds_ratio, r$fdr))
    }
  }

  summary_lines <- c(summary_lines, "",
    "SECTION D: Continuous Signal Correlation",
    "----------------------------------------",
    sprintf("  Signal source: %s", INPUT_FILES$signal)
  )

  if (nrow(corr_results) > 0) {
    for (i in seq_len(nrow(corr_results))) {
      r <- corr_results[i, ]
      summary_lines <- c(summary_lines,
        sprintf("  %s: rho=%.3f, p=%.2e, n=%d", r$anchor_group, r$spearman_rho, r$spearman_p, r$n))
    }
  }

  if (nrow(logistic_results) > 0) {
    k119ub_row <- logistic_results %>% filter(term == "mean_anchor_k119ub_fc")
    if (nrow(k119ub_row) > 0) {
      summary_lines <- c(summary_lines, "",
        "  Logistic regression: P(lost) ~ K119ub_fc + log(distance)",
        sprintf("    K119ub_fc coefficient: %.4f (p = %.4e, OR = %.3f)",
                k119ub_row$estimate, k119ub_row$p_value, k119ub_row$odds_ratio),
        sprintf("    Interpretation: %s",
                ifelse(k119ub_row$estimate > 0,
                       "Higher K119ub change -> more likely LOST (supports buffer hypothesis)",
                       "Higher K119ub change -> more likely GAINED (contradicts buffer hypothesis)")))
    }
  }

  summary_lines <- c(summary_lines, "",
    "SECTION E: Shared Anchor K119ub",
    "-------------------------------",
    sprintf("  Shared anchors: %d", nrow(shared_anchors))
  )

  if (nrow(shared_comparison) > 0) {
    for (i in seq_len(nrow(shared_comparison))) {
      r <- shared_comparison[i, ]
      summary_lines <- c(summary_lines,
        sprintf("  %s: shared=%.1f%% vs non-shared=%.1f%%, OR=%.2f, FDR=%.3e",
                r$mark, r$shared_pct, r$nonshared_pct, r$odds_ratio, r$fdr))
    }
  }

  summary_lines <- c(summary_lines, "",
    "================================================================================",
    "OUTPUT FILES:",
    "  tables/loops_with_k119ub_annotation.tsv",
    "  tables/distance_stats_by_k119ub.tsv",
    "  tables/enrichment_global.tsv",
    "  tables/enrichment_by_chromatin_state.tsv",
    "  tables/enrichment_by_distance.tsv",
    "  tables/correlation_results.tsv",
    "  tables/logistic_regression_results.tsv",
    "  tables/shared_anchor_k119ub_comparison.tsv",
    "  plots/01-08: CDF/density plots",
    "  plots/09-11: Enrichment plots",
    "  plots/12-15: Correlation plots",
    "  plots/16-19: Shared anchor plots",
    "================================================================================"
  )

  writeLines(summary_lines, file.path(OUTPUT_DIR, "analysis_summary.txt"))
  cat("  Saved: analysis_summary.txt\n")

  cat("\n================================================================================\n")
  cat("Analysis complete!\n")
  cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
  cat("================================================================================\n")

  invisible(list(
    all_loops = all_loops,
    diff_loops = diff_loops,
    shared_anchors = shared_anchors,
    global_results = global_results,
    chromstate_results = chromstate_results,
    corr_results = corr_results,
    logistic_results = logistic_results,
    shared_comparison = shared_comparison,
    signal_file = INPUT_FILES$signal
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
