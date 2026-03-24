# scripts/shared_anchor_analysis.R
#
# Shared Anchor / Loop Switching Analysis
#
# Biological Question: Do anchors that participate in lost long-range loops
# also participate in gained short-range loops? Prove "loop switching" at
# the same genomic sites.
#
# Usage:
#   Rscript scripts/shared_anchor_analysis.R                      # Late timepoint (default)
#   Rscript scripts/shared_anchor_analysis.R --timepoint early    # Early timepoint
#   Rscript scripts/shared_anchor_analysis.R --timepoint both     # Both timepoints
#   Rscript scripts/shared_anchor_analysis.R --tolerance 20000    # Custom anchor tolerance (bp)
#
# Output:
#   output/shared_anchor_analysis/{timepoint}/
#     tables/     - TSV files with statistics
#     plots/      - Multi-format plots (PDF, SVG, JPEG)
#     apa_subsets/- RDS files for downstream APA analysis

# ==============================================================================
# 1. PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(readxl)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(InteractionSet)
})

# Load multi-format output utility
source("data/scripts/_shared/multi_format_output.R")  # Original: source("scripts/utils/multi_format_output.R")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

# Base directory
BASE_DIR <- getwd()

# Default parameters
DEFAULT_TOLERANCE <- 10000  # 10kb for anchor overlap matching

# DEG thresholds
DEG_PADJ_THRESHOLD <- 0.05
DEG_LFC_THRESHOLD <- 0.3

# GREAT-style gene association parameters
GREAT_UPSTREAM <- 5000
GREAT_DOWNSTREAM <- 1000
GREAT_MAX_EXTENSION <- 100000

# Timepoint-specific file mappings
TIMEPOINT_CONFIG <- list(
  late = list(
    loops_file = file.path(BASE_DIR, "data/upstream/loop_calls/late_characterized_loops.tsv"),  # Original: 25042-late_outputs/merged_loops/characterized_loops.tsv
    rna_file = file.path(BASE_DIR, "data/upstream/rna_seq/adult_rnaseq_results.xlsx"),  # Original: tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx
    output_dir_tsvs = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_analysis/late (tables)
    output_dir_plots = file.path(BASE_DIR, "data/plots/supplemental"),  # Original: output/shared_anchor_analysis/late (plots)
    output_dir = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_analysis/late
    label = "Late Timepoint"
  ),
  early = list(
    loops_file = file.path(BASE_DIR, "250831-early_outputs/merged_loops/characterized_loops.tsv"),  # TODO: not in data/
    rna_file = file.path(BASE_DIR, "tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx"),  # TODO: not in data/
    output_dir_tsvs = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_analysis/early (tables)
    output_dir_plots = file.path(BASE_DIR, "data/plots/supplemental"),  # Original: output/shared_anchor_analysis/early (plots)
    output_dir = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_analysis/early
    label = "Early Timepoint"
  )
)

# Chromatin state colors (matching existing pipeline)
CHROMATIN_COLORS <- c(
  "Active_Promoter" = "#E41A1C",
  "Repressed_Promoter" = "#377EB8",
  "Bivalent_Promoter" = "#984EA3",
  "Polycomb" = "#4DAF4A",
  "Active_Enhancer" = "#FF7F00",
  "Poised_Enhancer" = "#FFFF33",
  "Other" = "#A65628"
)

# Direction colors
DIRECTION_COLORS <- c(
  "down_in_mutant" = "#4575b4",  # Lost (blue)
  "up_in_mutant" = "#d73027"     # Gained (red)
)

# ==============================================================================
# 3. HELPER FUNCTIONS - Data Loading
# ==============================================================================

#' Load characterized loops from TSV
#' @param file_path Path to characterized_loops.tsv
#' @return tibble of loops with required columns
load_characterized_loops <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("Loops file not found: %s", file_path))
  }

  loops <- read_tsv(file_path, show_col_types = FALSE)

  # Verify required columns
  required_cols <- c("direction", "anchor1_chr", "anchor1_start", "anchor1_end",
                     "anchor2_chr", "anchor2_start", "anchor2_end", "loop_distance",
                     "anchor1_type", "anchor2_type")
  missing <- setdiff(required_cols, colnames(loops))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }

  cat(sprintf("  Loaded %d loops from %s\n", nrow(loops), basename(file_path)))
  cat(sprintf("    - down_in_mutant (lost): %d\n", sum(loops$direction == "down_in_mutant")))
  cat(sprintf("    - up_in_mutant (gained): %d\n", sum(loops$direction == "up_in_mutant")))

  return(loops)
}

#' Extract all anchors from loops as GRanges
#' @param loops tibble of loops
#' @param direction_filter Optional filter for direction
#' @return GRanges of all anchors with metadata
extract_anchors_as_granges <- function(loops, direction_filter = NULL) {
  if (!is.null(direction_filter)) {
    loops <- loops %>% filter(direction == direction_filter)
  }

  # Extract anchor1
  anchor1_gr <- GRanges(
    seqnames = loops$anchor1_chr,
    ranges = IRanges(start = loops$anchor1_start, end = loops$anchor1_end),
    loop_id = loops$loop_id,
    anchor_num = 1,
    direction = loops$direction,
    loop_distance = loops$loop_distance,
    anchor_type = loops$anchor1_type
  )

  # Extract anchor2
  anchor2_gr <- GRanges(
    seqnames = loops$anchor2_chr,
    ranges = IRanges(start = loops$anchor2_start, end = loops$anchor2_end),
    loop_id = loops$loop_id,
    anchor_num = 2,
    direction = loops$direction,
    loop_distance = loops$loop_distance,
    anchor_type = loops$anchor2_type
  )

  # Combine
  all_anchors <- c(anchor1_gr, anchor2_gr)
  return(all_anchors)
}

# ==============================================================================
# 4. TASK 1a: IDENTIFY SHARED ANCHORS
# ==============================================================================

#' Identify shared anchors between lost and gained loops
#' @param loops tibble of loops
#' @param tolerance Maximum distance (bp) for anchors to be considered overlapping
#' @return list with shared_anchors GRanges and summary statistics
identify_shared_anchors <- function(loops, tolerance = 10000) {
  cat("\n[Task 1a] Identifying shared anchors...\n")

  # Extract anchors by direction
  lost_anchors <- extract_anchors_as_granges(loops, "down_in_mutant")
  gained_anchors <- extract_anchors_as_granges(loops, "up_in_mutant")

  cat(sprintf("  Lost loop anchors: %d\n", length(lost_anchors)))
  cat(sprintf("  Gained loop anchors: %d\n", length(gained_anchors)))

  # Find overlapping anchors with tolerance (maxgap)
  overlaps <- findOverlaps(lost_anchors, gained_anchors, maxgap = tolerance)

  if (length(overlaps) == 0) {
    cat("  WARNING: No overlapping anchors found!\n")
    return(list(
      shared_anchors = GRanges(),
      shared_anchor_df = tibble(),
      lost_at_shared = tibble(),
      gained_at_shared = tibble(),
      stats = list(n_shared = 0, n_lost_only = length(unique(lost_anchors)), n_gained_only = length(unique(gained_anchors)))
    ))
  }

  # Get unique shared anchor positions (from lost anchors)
  shared_lost_idx <- unique(queryHits(overlaps))
  shared_gained_idx <- unique(subjectHits(overlaps))

  # Create merged shared anchor regions
  shared_from_lost <- lost_anchors[shared_lost_idx]
  shared_from_gained <- gained_anchors[shared_gained_idx]

  # Merge overlapping regions to get consensus shared anchors
  # Use GenomicRanges::reduce explicitly to avoid purrr conflict
  all_shared <- c(shared_from_lost, shared_from_gained)
  shared_reduced <- GenomicRanges::reduce(all_shared, ignore.strand = TRUE)

  # Build summary of shared anchors with loop counts
  # Pre-extract metadata for efficient access (convert to plain vectors to avoid Rle issues)
  lost_loop_ids <- as.character(mcols(lost_anchors)$loop_id)
  gained_loop_ids <- as.character(mcols(gained_anchors)$loop_id)

  shared_summary <- tibble()
  for (i in seq_along(shared_reduced)) {
    anchor_region <- shared_reduced[i]

    # Count overlapping lost loops
    lost_hits <- findOverlaps(anchor_region, lost_anchors, maxgap = tolerance)
    n_lost <- length(unique(lost_loop_ids[subjectHits(lost_hits)]))

    # Count overlapping gained loops
    gained_hits <- findOverlaps(anchor_region, gained_anchors, maxgap = tolerance)
    n_gained <- length(unique(gained_loop_ids[subjectHits(gained_hits)]))

    shared_summary <- bind_rows(shared_summary, tibble(
      anchor_id = paste0("shared_", i),
      chr = as.character(seqnames(anchor_region)),
      start = start(anchor_region),
      end = end(anchor_region),
      n_lost_loops = n_lost,
      n_gained_loops = n_gained,
      n_total_loops = n_lost + n_gained
    ))
  }

  # Get loops touching shared anchors
  lost_at_shared_idx <- unique(lost_loop_ids[queryHits(overlaps)])
  gained_at_shared_idx <- unique(gained_loop_ids[subjectHits(overlaps)])

  lost_at_shared <- loops %>% filter(loop_id %in% lost_at_shared_idx)
  gained_at_shared <- loops %>% filter(loop_id %in% gained_at_shared_idx)

  # Calculate non-shared anchor counts
  all_lost_loop_ids <- unique(lost_loop_ids)
  all_gained_loop_ids <- unique(gained_loop_ids)

  lost_only_count <- length(setdiff(all_lost_loop_ids, lost_at_shared_idx))
  gained_only_count <- length(setdiff(all_gained_loop_ids, gained_at_shared_idx))

  cat(sprintf("  Shared anchor regions (tolerance %d bp): %d\n", tolerance, length(shared_reduced)))
  cat(sprintf("  Lost loops touching shared anchors: %d\n", nrow(lost_at_shared)))
  cat(sprintf("  Gained loops touching shared anchors: %d\n", nrow(gained_at_shared)))
  cat(sprintf("  Lost loops with unique anchors: %d\n", lost_only_count))
  cat(sprintf("  Gained loops with unique anchors: %d\n", gained_only_count))

  return(list(
    shared_anchors = shared_reduced,
    shared_anchor_df = shared_summary,
    lost_at_shared = lost_at_shared,
    gained_at_shared = gained_at_shared,
    overlaps = overlaps,
    lost_anchors = lost_anchors,
    gained_anchors = gained_anchors,
    stats = list(
      n_shared = length(shared_reduced),
      n_lost_at_shared = nrow(lost_at_shared),
      n_gained_at_shared = nrow(gained_at_shared),
      n_lost_only = lost_only_count,
      n_gained_only = gained_only_count,
      tolerance = tolerance
    )
  ))
}

#' Create anchor overlap summary plot
#' @param shared_result Result from identify_shared_anchors
#' @return ggplot object
create_anchor_overlap_plot <- function(shared_result) {
  stats <- shared_result$stats

  plot_data <- tibble(
    category = c("Lost loops at\nshared anchors", "Gained loops at\nshared anchors",
                 "Lost loops at\nunique anchors", "Gained loops at\nunique anchors"),
    count = c(stats$n_lost_at_shared, stats$n_gained_at_shared,
              stats$n_lost_only, stats$n_gained_only),
    anchor_type = c("Shared", "Shared", "Unique", "Unique"),
    direction = c("Lost", "Gained", "Lost", "Gained")
  )

  p <- ggplot(plot_data, aes(x = category, y = count, fill = interaction(direction, anchor_type))) +
    geom_col(color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c(
      "Lost.Shared" = "#4575b4",
      "Gained.Shared" = "#d73027",
      "Lost.Unique" = "#91bfdb",
      "Gained.Unique" = "#fc8d59"
    ), name = "Category") +
    geom_text(aes(label = count), vjust = -0.3, size = 4, fontface = "bold") +
    labs(
      title = "Loop Distribution by Anchor Sharing",
      subtitle = sprintf("Shared anchor regions: %d (tolerance: %d bp)",
                        stats$n_shared, stats$tolerance),
      x = NULL,
      y = "Number of Loops"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    coord_cartesian(ylim = c(0, max(plot_data$count) * 1.15))

  return(p)
}

# ==============================================================================
# 5. TASK 1b: CHARACTERIZE SHARED VS NON-SHARED ANCHORS
# ==============================================================================

#' Characterize shared vs non-shared anchors
#' @param loops All loops
#' @param shared_result Result from identify_shared_anchors
#' @return list with comparison statistics and data
characterize_anchors <- function(loops, shared_result) {
  cat("\n[Task 1b] Characterizing shared vs non-shared anchors...\n")

  # Get loops at shared vs unique anchors
  shared_loop_ids <- c(shared_result$lost_at_shared$loop_id,
                       shared_result$gained_at_shared$loop_id)

  shared_loops <- loops %>% filter(loop_id %in% shared_loop_ids)
  unique_loops <- loops %>% filter(!loop_id %in% shared_loop_ids)

  # Extract anchor types for comparison
  get_anchor_types <- function(loop_df, label) {
    types1 <- loop_df %>% dplyr::select(anchor_type = anchor1_type) %>% mutate(anchor_class = label)
    types2 <- loop_df %>% dplyr::select(anchor_type = anchor2_type) %>% mutate(anchor_class = label)
    bind_rows(types1, types2)
  }

  anchor_type_data <- bind_rows(
    get_anchor_types(shared_loops, "Shared"),
    get_anchor_types(unique_loops, "Unique")
  )

  # Chi-square test for chromatin state distribution
  type_table <- table(anchor_type_data$anchor_type, anchor_type_data$anchor_class)
  chisq_result <- tryCatch({
    chisq.test(type_table)
  }, error = function(e) {
    warning("Chi-square test failed: ", e$message)
    list(p.value = NA, statistic = NA)
  })

  cat(sprintf("  Chromatin state Chi-square p-value: %.2e\n", chisq_result$p.value))

  # Compare distance to TSS
  shared_tss_dist <- c(shared_loops$anchor1_distance_to_tss, shared_loops$anchor2_distance_to_tss)
  unique_tss_dist <- c(unique_loops$anchor1_distance_to_tss, unique_loops$anchor2_distance_to_tss)

  tss_wilcox <- wilcox.test(shared_tss_dist, unique_tss_dist)
  cat(sprintf("  Distance to TSS Wilcoxon p-value: %.2e\n", tss_wilcox$p.value))

  # ChIP-seq mark comparison with Fisher's exact test
  chip_marks <- c("H3K27ac", "H3K4me1", "H3K27me3", "H3K4me3", "Bivalent_Promoter")
  chip_results <- list()

  for (mark in chip_marks) {
    col1 <- paste0("anchor1_", mark, "_overlap")
    col2 <- paste0("anchor2_", mark, "_overlap")

    if (col1 %in% colnames(shared_loops) && col2 %in% colnames(shared_loops)) {
      shared_pos <- sum(shared_loops[[col1]], na.rm = TRUE) + sum(shared_loops[[col2]], na.rm = TRUE)
      shared_neg <- 2 * nrow(shared_loops) - shared_pos
      unique_pos <- sum(unique_loops[[col1]], na.rm = TRUE) + sum(unique_loops[[col2]], na.rm = TRUE)
      unique_neg <- 2 * nrow(unique_loops) - unique_pos

      if (shared_pos > 0 || unique_pos > 0) {
        fisher_mat <- matrix(c(shared_pos, shared_neg, unique_pos, unique_neg), nrow = 2)
        fisher_test <- fisher.test(fisher_mat)

        chip_results[[mark]] <- tibble(
          mark = mark,
          shared_positive = shared_pos,
          shared_total = 2 * nrow(shared_loops),
          unique_positive = unique_pos,
          unique_total = 2 * nrow(unique_loops),
          odds_ratio = fisher_test$estimate,
          p_value = fisher_test$p.value,
          conf_low = fisher_test$conf.int[1],
          conf_high = fisher_test$conf.int[2]
        )
      }
    }
  }

  chip_df <- bind_rows(chip_results)

  return(list(
    anchor_type_data = anchor_type_data,
    type_table = type_table,
    chisq_result = chisq_result,
    tss_wilcox = tss_wilcox,
    shared_tss_dist = shared_tss_dist,
    unique_tss_dist = unique_tss_dist,
    chip_results = chip_df
  ))
}

#' Create chromatin state comparison plot
#' @param char_result Result from characterize_anchors
#' @return ggplot object
create_chromatin_state_plot <- function(char_result) {
  data <- char_result$anchor_type_data

  # Calculate proportions
  prop_data <- data %>%
    group_by(anchor_class, anchor_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(anchor_class) %>%
    mutate(prop = n / sum(n), total = sum(n)) %>%
    ungroup()

  p <- ggplot(prop_data, aes(x = anchor_class, y = prop, fill = anchor_type)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(values = CHROMATIN_COLORS, name = "Chromatin State") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Chromatin State Distribution",
      subtitle = sprintf("Chi-square p = %.2e", char_result$chisq_result$p.value),
      x = "Anchor Class",
      y = "Proportion of Anchors"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "right"
    )

  return(p)
}

#' Create ChIP enrichment dotplot
#' @param char_result Result from characterize_anchors
#' @return ggplot object
create_chip_enrichment_plot <- function(char_result) {
  if (nrow(char_result$chip_results) == 0) {
    return(NULL)
  }

  data <- char_result$chip_results %>%
    mutate(
      log2OR = log2(odds_ratio),
      sig = ifelse(p_value < 0.05, "Significant", "NS"),
      mark = factor(mark, levels = c("H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "Bivalent_Promoter"))
    )

  p <- ggplot(data, aes(x = log2OR, y = mark, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = log2(conf_low), xmax = log2(conf_high)), height = 0.2) +
    geom_point(size = 4) +
    scale_color_manual(values = c("Significant" = "#E41A1C", "NS" = "gray50"), name = "") +
    labs(
      title = "ChIP-seq Mark Enrichment at Shared vs Unique Anchors",
      subtitle = "Odds ratio > 1 = enriched at shared anchors",
      x = "log2(Odds Ratio)",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "bottom"
    )

  return(p)
}

# ==============================================================================
# 6. TASK 1c: PAIRED DISTANCE ANALYSIS
# ==============================================================================

#' Perform paired distance analysis at shared anchors
#' @param shared_result Result from identify_shared_anchors
#' @param tolerance bp tolerance for matching
#' @return list with paired statistics
paired_distance_analysis <- function(shared_result, tolerance = 10000) {
  cat("\n[Task 1c] Performing paired distance analysis...\n")

  lost_anchors <- shared_result$lost_anchors
  gained_anchors <- shared_result$gained_anchors
  overlaps <- shared_result$overlaps

  if (length(overlaps) == 0) {
    return(list(paired_data = tibble(), test_result = NULL))
  }

  # Pre-extract metadata for efficient access (convert to plain vectors)
  lost_loop_distances <- as.numeric(mcols(lost_anchors)$loop_distance)
  gained_loop_distances <- as.numeric(mcols(gained_anchors)$loop_distance)

  # Build paired data at shared anchor positions
  # For each overlap, record the anchor position and distances from both directions
  paired_records <- list()

  # Group by shared anchor region
  shared_anchors <- shared_result$shared_anchors

  for (i in seq_along(shared_anchors)) {
    anchor_region <- shared_anchors[i]

    # Find lost loop distances at this anchor
    lost_hits <- findOverlaps(anchor_region, lost_anchors, maxgap = tolerance)
    lost_distances <- unique(lost_loop_distances[subjectHits(lost_hits)])

    # Find gained loop distances at this anchor
    gained_hits <- findOverlaps(anchor_region, gained_anchors, maxgap = tolerance)
    gained_distances <- unique(gained_loop_distances[subjectHits(gained_hits)])

    if (length(lost_distances) > 0 && length(gained_distances) > 0) {
      paired_records[[i]] <- tibble(
        anchor_id = paste0("shared_", i),
        anchor_chr = as.character(seqnames(anchor_region)),
        anchor_start = start(anchor_region),
        anchor_end = end(anchor_region),
        median_lost_distance = median(lost_distances),
        median_gained_distance = median(gained_distances),
        mean_lost_distance = mean(lost_distances),
        mean_gained_distance = mean(gained_distances),
        n_lost = length(lost_distances),
        n_gained = length(gained_distances)
      )
    }
  }

  paired_data <- bind_rows(paired_records)

  if (nrow(paired_data) == 0) {
    cat("  WARNING: No paired data available for analysis\n")
    return(list(paired_data = tibble(), test_result = NULL))
  }

  cat(sprintf("  Shared anchors with paired data: %d\n", nrow(paired_data)))

  # Paired Wilcoxon signed-rank test
  paired_test <- wilcox.test(
    paired_data$median_lost_distance,
    paired_data$median_gained_distance,
    paired = TRUE,
    alternative = "greater"  # Test: lost > gained
  )

  # Effect size: paired Cohen's d
  diff <- paired_data$median_lost_distance - paired_data$median_gained_distance
  cohens_d <- mean(diff) / sd(diff)

  # Bootstrap 95% CI for median difference
  set.seed(42)
  n_boot <- 1000
  boot_medians <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(nrow(paired_data), replace = TRUE)
    boot_diff <- paired_data$median_lost_distance[idx] - paired_data$median_gained_distance[idx]
    boot_medians[b] <- median(boot_diff)
  }
  ci_95 <- quantile(boot_medians, c(0.025, 0.975))

  cat(sprintf("  Paired Wilcoxon p-value (lost > gained): %.2e\n", paired_test$p.value))
  cat(sprintf("  Median difference: %.0f bp (95%% CI: %.0f to %.0f)\n",
              median(diff), ci_95[1], ci_95[2]))
  cat(sprintf("  Cohen's d: %.2f\n", cohens_d))

  return(list(
    paired_data = paired_data,
    test_result = paired_test,
    cohens_d = cohens_d,
    median_diff = median(diff),
    ci_95 = ci_95,
    stats = list(
      median_lost = median(paired_data$median_lost_distance),
      median_gained = median(paired_data$median_gained_distance),
      mean_lost = mean(paired_data$median_lost_distance),
      mean_gained = mean(paired_data$median_gained_distance)
    )
  ))
}

#' Create paired distance scatter plot
#' @param distance_result Result from paired_distance_analysis
#' @return ggplot object
create_paired_scatter_plot <- function(distance_result) {
  if (nrow(distance_result$paired_data) == 0) {
    return(NULL)
  }

  data <- distance_result$paired_data

  # Calculate axis limits
  max_dist <- max(c(data$median_lost_distance, data$median_gained_distance), na.rm = TRUE)

  p <- ggplot(data, aes(x = median_gained_distance / 1e6, y = median_lost_distance / 1e6)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = n_lost + n_gained), alpha = 0.6, color = "#984EA3") +
    scale_size_continuous(name = "Total loops", range = c(2, 8)) +
    labs(
      title = "Loop Distance at Shared Anchors",
      subtitle = sprintf("Paired Wilcoxon p = %.2e\nMedian diff = %.0f kb",
                        distance_result$test_result$p.value,
                        distance_result$median_diff / 1000),
      x = "Gained Loop Distance (Mb)",
      y = "Lost Loop Distance (Mb)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "right"
    ) +
    coord_fixed(xlim = c(0, max_dist / 1e6 * 1.1), ylim = c(0, max_dist / 1e6 * 1.1))

  return(p)
}

#' Create distance violin plot for shared anchors
#' @param shared_result Result from identify_shared_anchors
#' @return ggplot object
create_distance_violin_plot <- function(shared_result) {
  # Combine lost and gained at shared anchors
  lost_data <- shared_result$lost_at_shared %>%
    dplyr::select(loop_id, loop_distance, direction) %>%
    mutate(direction_label = "Lost")

  gained_data <- shared_result$gained_at_shared %>%
    dplyr::select(loop_id, loop_distance, direction) %>%
    mutate(direction_label = "Gained")

  plot_data <- bind_rows(lost_data, gained_data)

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  # Mann-Whitney U test
  test_result <- wilcox.test(
    lost_data$loop_distance,
    gained_data$loop_distance,
    alternative = "greater"
  )

  p <- ggplot(plot_data, aes(x = direction_label, y = loop_distance / 1e6, fill = direction_label)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Lost" = "#4575b4", "Gained" = "#d73027")) +
    labs(
      title = "Loop Distance at Shared Anchors",
      subtitle = sprintf("Mann-Whitney p = %.2e (lost > gained)", test_result$p.value),
      x = NULL,
      y = "Loop Distance (Mb)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "none"
    )

  return(p)
}

# ==============================================================================
# 7. TASK 1d: APA SUBSETS EXPORT
# ==============================================================================

#' Export loop subsets for APA analysis
#' @param shared_result Result from identify_shared_anchors
#' @param output_dir Output directory
#' @return list with file paths
export_apa_subsets <- function(shared_result, output_dir) {
  cat("\n[Task 1d] Exporting APA subsets...\n")

  apa_dir <- file.path(output_dir, "apa_subsets")
  dir.create(apa_dir, recursive = TRUE, showWarnings = FALSE)

  # Get median distances
  lost_median <- median(shared_result$lost_at_shared$loop_distance, na.rm = TRUE)
  gained_median <- median(shared_result$gained_at_shared$loop_distance, na.rm = TRUE)

  # Lost long-range (above median of lost)
  lost_longrange <- shared_result$lost_at_shared %>%
    filter(loop_distance > lost_median)

  # Gained short-range (below median of gained)
  gained_shortrange <- shared_result$gained_at_shared %>%
    filter(loop_distance < gained_median)

  cat(sprintf("  Lost long-range loops (> %.0f kb): %d\n", lost_median / 1000, nrow(lost_longrange)))
  cat(sprintf("  Gained short-range loops (< %.0f kb): %d\n", gained_median / 1000, nrow(gained_shortrange)))

  # Convert to GInteractions for APA
  loops_to_gi <- function(loop_df) {
    if (nrow(loop_df) == 0) return(GInteractions())

    anchor1 <- GRanges(
      seqnames = loop_df$anchor1_chr,
      ranges = IRanges(start = loop_df$anchor1_start, end = loop_df$anchor1_end)
    )
    anchor2 <- GRanges(
      seqnames = loop_df$anchor2_chr,
      ranges = IRanges(start = loop_df$anchor2_start, end = loop_df$anchor2_end)
    )

    gi <- GInteractions(anchor1, anchor2)
    mcols(gi) <- loop_df %>% dplyr::select(-starts_with("anchor1_"), -starts_with("anchor2_"))
    return(gi)
  }

  lost_gi <- loops_to_gi(lost_longrange)
  gained_gi <- loops_to_gi(gained_shortrange)

  # Save RDS files
  lost_path <- file.path(apa_dir, "shared_lost_longrange.rds")
  gained_path <- file.path(apa_dir, "shared_gained_shortrange.rds")

  saveRDS(lost_gi, lost_path)
  saveRDS(gained_gi, gained_path)

  cat(sprintf("  Saved: %s\n", basename(lost_path)))
  cat(sprintf("  Saved: %s\n", basename(gained_path)))

  # Write summary
  summary_lines <- c(
    "=== APA Loop Subset Summary ===",
    "",
    sprintf("Date: %s", Sys.time()),
    "",
    "--- Lost Long-Range Loops ---",
    sprintf("Count: %d", nrow(lost_longrange)),
    sprintf("Distance threshold: > %.0f bp (median of lost at shared)", lost_median),
    sprintf("Distance range: %.0f - %.0f kb",
            min(lost_longrange$loop_distance, na.rm = TRUE) / 1000,
            max(lost_longrange$loop_distance, na.rm = TRUE) / 1000),
    sprintf("File: %s", basename(lost_path)),
    "",
    "--- Gained Short-Range Loops ---",
    sprintf("Count: %d", nrow(gained_shortrange)),
    sprintf("Distance threshold: < %.0f bp (median of gained at shared)", gained_median),
    sprintf("Distance range: %.0f - %.0f kb",
            min(gained_shortrange$loop_distance, na.rm = TRUE) / 1000,
            max(gained_shortrange$loop_distance, na.rm = TRUE) / 1000),
    sprintf("File: %s", basename(gained_path)),
    "",
    "=== End Summary ==="
  )

  summary_path <- file.path(apa_dir, "loop_subset_summary.txt")
  writeLines(summary_lines, summary_path)
  cat(sprintf("  Saved: %s\n", basename(summary_path)))

  # Also export as TSV for inspection
  if (nrow(lost_longrange) > 0) {
    write_tsv(lost_longrange, file.path(apa_dir, "shared_lost_longrange.tsv"))
  }
  if (nrow(gained_shortrange) > 0) {
    write_tsv(gained_shortrange, file.path(apa_dir, "shared_gained_shortrange.tsv"))
  }

  return(list(
    lost_path = lost_path,
    gained_path = gained_path,
    lost_longrange = lost_longrange,
    gained_shortrange = gained_shortrange,
    lost_median = lost_median,
    gained_median = gained_median
  ))
}

# ==============================================================================
# 8. TASK 1e: GENE EXPRESSION VIOLIN
# ==============================================================================

#' Load RNA-seq results
#' @param rna_file Path to RNA-seq Excel file
#' @return tibble of significant DEGs
load_rnaseq_degs <- function(rna_file) {
  if (!file.exists(rna_file)) {
    stop(sprintf("RNA-seq file not found: %s", rna_file))
  }

  rna_df <- read_excel(rna_file)

  # The column is named ensembl_gene_id but actually contains gene symbols
  if (!"ensembl_gene_id" %in% colnames(rna_df)) {
    stop("Expected column 'ensembl_gene_id' not found in RNA-seq file")
  }

  # Filter for significant DEGs
  deg_df <- rna_df %>%
    filter(!is.na(padj) & padj < DEG_PADJ_THRESHOLD) %>%
    filter(abs(log2FoldChange) > DEG_LFC_THRESHOLD) %>%
    dplyr::select(gene_symbol = ensembl_gene_id, log2FoldChange, padj, baseMean)

  cat(sprintf("  Loaded %d significant DEGs (padj < %g, |log2FC| > %g)\n",
              nrow(deg_df), DEG_PADJ_THRESHOLD, DEG_LFC_THRESHOLD))

  return(deg_df)
}

#' Convert gene symbols to Entrez IDs
#' @param gene_symbols Vector of gene symbols
#' @return tibble mapping symbols to Entrez IDs
convert_symbols_to_entrez <- function(gene_symbols) {
  mapping <- tryCatch({
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unique(gene_symbols),
      columns = "ENTREZID",
      keytype = "SYMBOL"
    )
  }, error = function(e) {
    warning(sprintf("Error converting gene symbols: %s", e$message))
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  })

  mapping <- mapping %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  cat(sprintf("  Converted %d / %d gene symbols to Entrez\n",
              nrow(mapping), length(unique(gene_symbols))))

  return(mapping)
}

#' Find genes at shared anchors using GREAT-style regulatory domains
#' @param shared_result Result from identify_shared_anchors
#' @return tibble of gene-anchor associations
find_anchor_genes <- function(shared_result) {
  # Get gene information from TxDb
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes_gr <- genes(txdb)

  # Get TSS positions
  tss_pos <- ifelse(as.character(strand(genes_gr)) == "-",
                    end(genes_gr), start(genes_gr))
  gene_strand <- as.character(strand(genes_gr))
  gene_chr <- as.character(seqnames(genes_gr))
  gene_ids <- names(genes_gr)

  gene_info <- data.frame(
    entrez_id = gene_ids,
    chr = gene_chr,
    tss = tss_pos,
    strand = gene_strand,
    stringsAsFactors = FALSE
  )

  # Calculate GREAT-style regulatory domains
  gene_info <- gene_info %>%
    arrange(chr, tss) %>%
    mutate(
      basal_start = ifelse(strand == "+", tss - GREAT_UPSTREAM, tss - GREAT_DOWNSTREAM),
      basal_end = ifelse(strand == "+", tss + GREAT_DOWNSTREAM, tss + GREAT_UPSTREAM),
      max_start = tss - GREAT_MAX_EXTENSION,
      max_end = tss + GREAT_MAX_EXTENSION
    ) %>%
    group_by(chr) %>%
    mutate(
      prev_basal_end = lag(basal_end, default = -Inf),
      next_basal_start = lead(basal_start, default = Inf),
      reg_start = pmax(max_start, prev_basal_end, 1),
      reg_end = pmin(max_end, next_basal_start)
    ) %>%
    ungroup() %>%
    mutate(
      reg_start = ifelse(reg_end < reg_start, basal_start, reg_start),
      reg_end = ifelse(reg_end < reg_start, basal_end, reg_end),
      reg_end = pmax(reg_end, reg_start)
    ) %>%
    dplyr::select(entrez_id, chr, tss, strand, reg_start, reg_end)

  # Create GRanges for gene regulatory domains
  gene_domains_gr <- GRanges(
    seqnames = gene_info$chr,
    ranges = IRanges(start = gene_info$reg_start, end = gene_info$reg_end),
    entrez_id = gene_info$entrez_id
  )

  # Find overlaps with shared anchors
  shared_anchors <- shared_result$shared_anchors

  if (length(shared_anchors) == 0) {
    return(tibble(entrez_id = character(), anchor_chr = character()))
  }

  overlaps <- findOverlaps(shared_anchors, gene_domains_gr, ignore.strand = TRUE)

  if (length(overlaps) == 0) {
    warning("No gene-anchor overlaps found")
    return(tibble(entrez_id = character(), anchor_chr = character()))
  }

  gene_anchor_df <- tibble(
    entrez_id = gene_domains_gr$entrez_id[subjectHits(overlaps)],
    anchor_chr = as.character(seqnames(shared_anchors)[queryHits(overlaps)]),
    anchor_start = start(shared_anchors)[queryHits(overlaps)],
    anchor_end = end(shared_anchors)[queryHits(overlaps)]
  ) %>%
    distinct()

  cat(sprintf("  Found %d gene-anchor associations at shared anchors\n", nrow(gene_anchor_df)))

  return(gene_anchor_df)
}

#' Create gene expression analysis
#' @param shared_result Result from identify_shared_anchors
#' @param deg_df DEG tibble
#' @return list with plot and data
gene_expression_analysis <- function(shared_result, deg_df) {
  cat("\n[Task 1e] Performing gene expression analysis...\n")

  # Get genes at shared anchors
  gene_anchor_df <- find_anchor_genes(shared_result)

  if (nrow(gene_anchor_df) == 0) {
    return(list(plot = NULL, data = tibble()))
  }

  # Convert DEG symbols to Entrez
  id_mapping <- convert_symbols_to_entrez(deg_df$gene_symbol)

  # Merge
  deg_with_entrez <- deg_df %>%
    inner_join(id_mapping, by = c("gene_symbol" = "SYMBOL"))

  shared_anchor_genes <- gene_anchor_df %>%
    inner_join(deg_with_entrez, by = c("entrez_id" = "ENTREZID")) %>%
    distinct(entrez_id, .keep_all = TRUE)

  cat(sprintf("  DEGs at shared anchors: %d\n", nrow(shared_anchor_genes)))

  if (nrow(shared_anchor_genes) < 5) {
    cat("  WARNING: Too few genes for meaningful violin plot\n")
    return(list(plot = NULL, data = shared_anchor_genes))
  }

  # One-sample Wilcoxon test (is median != 0?)
  wilcox_result <- wilcox.test(shared_anchor_genes$log2FoldChange, mu = 0)

  cat(sprintf("  One-sample Wilcoxon p-value (median != 0): %.2e\n", wilcox_result$p.value))
  cat(sprintf("  Median log2FC: %.3f\n", median(shared_anchor_genes$log2FoldChange)))

  # Create violin plot
  p <- ggplot(shared_anchor_genes, aes(x = "Shared Anchors", y = log2FoldChange)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_violin(fill = "#984EA3", alpha = 0.7, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.5) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1, color = "black") +
    labs(
      title = "DEG Expression at Shared Anchors",
      subtitle = sprintf("n = %d genes, Wilcoxon p = %.2e\nMedian log2FC = %.2f",
                        nrow(shared_anchor_genes), wilcox_result$p.value,
                        median(shared_anchor_genes$log2FoldChange)),
      x = NULL,
      y = expression(log[2]*"FC (BAP1-KO / WT)")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )

  return(list(
    plot = p,
    data = shared_anchor_genes,
    test = wilcox_result
  ))
}

# ==============================================================================
# 9. MAIN PIPELINE FUNCTION
# ==============================================================================

#' Process a single timepoint
#' @param timepoint "late" or "early"
#' @param tolerance Anchor overlap tolerance in bp
process_timepoint <- function(timepoint, tolerance = 10000) {
  cat("\n")
  cat("=======================================================\n")
  cat(sprintf("Processing %s timepoint\n", toupper(timepoint)))
  cat("=======================================================\n")

  config <- TIMEPOINT_CONFIG[[timepoint]]
  if (is.null(config)) {
    stop(sprintf("Unknown timepoint: %s", timepoint))
  }

  # Create output directories
  tables_dir <- config$output_dir_tsvs  # Original: file.path(config$output_dir, "tables")
  plots_dir <- config$output_dir_plots  # Original: file.path(config$output_dir, "plots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Step 1: Load loops
  cat("\n[Loading] Characterized loops...\n")
  loops <- load_characterized_loops(config$loops_file)

  # Task 1a: Identify shared anchors
  shared_result <- identify_shared_anchors(loops, tolerance)

  # Save shared anchor tables
  write_tsv(shared_result$shared_anchor_df, file.path(tables_dir, "shared_anchors.tsv"))  # Original: tables/shared_anchors.tsv
  cat(sprintf("  Saved: shared_anchors.tsv\n"))

  shared_loop_df <- bind_rows(
    shared_result$lost_at_shared %>% mutate(anchor_status = "shared"),
    shared_result$gained_at_shared %>% mutate(anchor_status = "shared")
  )
  write_tsv(shared_loop_df, file.path(tables_dir, "shared_anchor_loops.tsv"))  # Original: tables/shared_anchor_loops.tsv
  cat(sprintf("  Saved: shared_anchor_loops.tsv\n"))

  # Create and save anchor overlap plot
  overlap_plot <- create_anchor_overlap_plot(shared_result)
  save_multiformat_ggplot(overlap_plot, file.path(plots_dir, "anchor_overlap_summary"),
                          width = 8, height = 6)

  # Task 1b: Characterize anchors
  char_result <- characterize_anchors(loops, shared_result)

  # Save characterization table
  char_summary <- tibble(
    comparison = c("chromatin_state_chisq", "tss_distance_wilcox"),
    statistic = c(char_result$chisq_result$statistic, char_result$tss_wilcox$statistic),
    p_value = c(char_result$chisq_result$p.value, char_result$tss_wilcox$p.value)
  )
  write_tsv(char_summary, file.path(tables_dir, "anchor_characterization.tsv"))  # Original: tables/anchor_characterization.tsv
  if (nrow(char_result$chip_results) > 0) {
    write_tsv(char_result$chip_results, file.path(tables_dir, "shared_anchor_chip_enrichment.tsv"))  # Original: tables/chip_enrichment_results.tsv
  }

  # Create and save characterization plots
  chrom_plot <- create_chromatin_state_plot(char_result)
  save_multiformat_ggplot(chrom_plot, file.path(plots_dir, "chromatin_state_comparison"),
                          width = 8, height = 6)

  chip_plot <- create_chip_enrichment_plot(char_result)
  if (!is.null(chip_plot)) {
    save_multiformat_ggplot(chip_plot, file.path(plots_dir, "chip_enrichment_dotplot"),
                            width = 8, height = 5)
  }

  # Task 1c: Paired distance analysis
  distance_result <- paired_distance_analysis(shared_result, tolerance)

  if (nrow(distance_result$paired_data) > 0) {
    write_tsv(distance_result$paired_data, file.path(tables_dir, "shared_anchor_paired_distance.tsv"))  # Original: tables/paired_distance_stats.tsv

    # Create and save distance plots
    scatter_plot <- create_paired_scatter_plot(distance_result)
    if (!is.null(scatter_plot)) {
      save_multiformat_ggplot(scatter_plot, file.path(plots_dir, "paired_distance_scatter"),
                              width = 7, height = 7)
    }
  }

  violin_plot <- create_distance_violin_plot(shared_result)
  if (!is.null(violin_plot)) {
    save_multiformat_ggplot(violin_plot, file.path(plots_dir, "distance_violin_shared"),
                            width = 5, height = 6)
  }

  # Task 1d: APA subsets
  apa_result <- export_apa_subsets(shared_result, config$output_dir_tsvs)  # Original: config$output_dir

  # Task 1e: Gene expression
  deg_df <- load_rnaseq_degs(config$rna_file)
  expr_result <- gene_expression_analysis(shared_result, deg_df)

  if (!is.null(expr_result$plot)) {
    save_multiformat_ggplot(expr_result$plot, file.path(plots_dir, "expression_violin"),
                            width = 5, height = 6)
  }
  if (nrow(expr_result$data) > 0) {
    write_tsv(expr_result$data, file.path(tables_dir, "shared_anchor_genes.tsv"))  # Original: tables/shared_anchor_genes.tsv
  }

  # Generate summary report
  generate_summary_report(shared_result, char_result, distance_result,
                         apa_result, expr_result, config, tolerance)

  cat("\n")
  cat(sprintf("Completed %s timepoint\n", timepoint))
  cat(sprintf("  Output directory: %s\n", config$output_dir))

  return(list(
    shared_result = shared_result,
    char_result = char_result,
    distance_result = distance_result,
    apa_result = apa_result,
    expr_result = expr_result
  ))
}

#' Generate summary report
generate_summary_report <- function(shared_result, char_result, distance_result,
                                    apa_result, expr_result, config, tolerance) {
  report_lines <- c(
    "=================================================================",
    sprintf("SHARED ANCHOR ANALYSIS SUMMARY - %s", config$label),
    "=================================================================",
    "",
    sprintf("Date: %s", Sys.time()),
    sprintf("Input: %s", config$loops_file),
    sprintf("Tolerance: %d bp", tolerance),
    "",
    "--- SHARED ANCHOR IDENTIFICATION ---",
    sprintf("Shared anchor regions: %d", shared_result$stats$n_shared),
    sprintf("Lost loops at shared anchors: %d", shared_result$stats$n_lost_at_shared),
    sprintf("Gained loops at shared anchors: %d", shared_result$stats$n_gained_at_shared),
    sprintf("Lost loops at unique anchors: %d", shared_result$stats$n_lost_only),
    sprintf("Gained loops at unique anchors: %d", shared_result$stats$n_gained_only),
    "",
    "--- ANCHOR CHARACTERIZATION ---",
    sprintf("Chromatin state Chi-square p: %.2e", char_result$chisq_result$p.value),
    sprintf("TSS distance Wilcoxon p: %.2e", char_result$tss_wilcox$p.value),
    ""
  )

  if (!is.null(distance_result$test_result)) {
    report_lines <- c(report_lines,
      "--- PAIRED DISTANCE ANALYSIS ---",
      sprintf("Paired anchors tested: %d", nrow(distance_result$paired_data)),
      sprintf("Paired Wilcoxon p (lost > gained): %.2e", distance_result$test_result$p.value),
      sprintf("Median lost distance: %.0f bp", distance_result$stats$median_lost),
      sprintf("Median gained distance: %.0f bp", distance_result$stats$median_gained),
      sprintf("Median difference: %.0f bp", distance_result$median_diff),
      sprintf("Cohen's d: %.2f", distance_result$cohens_d),
      ""
    )
  }

  report_lines <- c(report_lines,
    "--- APA SUBSETS ---",
    sprintf("Lost long-range loops: %d", nrow(apa_result$lost_longrange)),
    sprintf("Gained short-range loops: %d", nrow(apa_result$gained_shortrange)),
    ""
  )

  if (!is.null(expr_result$test)) {
    report_lines <- c(report_lines,
      "--- GENE EXPRESSION ---",
      sprintf("DEGs at shared anchors: %d", nrow(expr_result$data)),
      sprintf("One-sample Wilcoxon p: %.2e", expr_result$test$p.value),
      sprintf("Median log2FC: %.3f", median(expr_result$data$log2FoldChange)),
      ""
    )
  }

  report_lines <- c(report_lines,
    "=================================================================",
    "END REPORT",
    "================================================================="
  )

  report_path <- file.path(config$output_dir_tsvs, "summary_report.txt")  # Original: file.path(config$output_dir, "summary_report.txt")
  writeLines(report_lines, report_path)
  cat(sprintf("\n  Saved: summary_report.txt\n"))
}

# ==============================================================================
# 10. COMMAND LINE INTERFACE
# ==============================================================================

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Default values
  timepoint_arg <- "late"
  tolerance_arg <- DEFAULT_TOLERANCE

  # Parse arguments
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint_arg <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--tolerance" && i < length(args)) {
      tolerance_arg <- as.numeric(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  cat("=================================================================\n")
  cat("SHARED ANCHOR / LOOP SWITCHING ANALYSIS\n")
  cat("=================================================================\n")
  cat(sprintf("Timepoint: %s\n", timepoint_arg))
  cat(sprintf("Tolerance: %d bp\n", tolerance_arg))
  cat(sprintf("Working directory: %s\n", getwd()))

  # Process timepoints
  if (tolower(timepoint_arg) == "both") {
    timepoints <- c("late", "early")
  } else if (tolower(timepoint_arg) %in% c("late", "early")) {
    timepoints <- tolower(timepoint_arg)
  } else {
    stop(sprintf("Invalid timepoint: %s. Use 'late', 'early', or 'both'", timepoint_arg))
  }

  results <- list()
  for (tp in timepoints) {
    results[[tp]] <- tryCatch({
      process_timepoint(tp, tolerance_arg)
    }, error = function(e) {
      cat(sprintf("\nERROR processing %s: %s\n", tp, e$message))
      return(NULL)
    })
  }

  cat("\n")
  cat("=================================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("=================================================================\n")

  for (tp in names(results)) {
    if (!is.null(results[[tp]])) {
      cat(sprintf("\n%s timepoint TSVs: %s\n", toupper(tp),
                  TIMEPOINT_CONFIG[[tp]]$output_dir_tsvs))
      cat(sprintf("%s timepoint plots: %s\n", toupper(tp),
                  TIMEPOINT_CONFIG[[tp]]$output_dir_plots))
    }
  }
}

# Run if called directly
if (!interactive()) {
  main()
}
