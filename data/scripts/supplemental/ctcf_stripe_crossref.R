# scripts/ctcf_stripe_crossref.R
#
# Cross-reference CTCF-anchored Lost Loops with Stripe Defects
#
# Biological Question: Do regions with lost CTCF-anchored loops show stripe defects?
# This tests whether disruption of CTCF-mediated loops and stripe architecture are
# coordinated in BAP1-KO, suggesting a common mechanism (e.g., cohesin/extrusion dysfunction).
#
# Usage:
#   Rscript scripts/ctcf_stripe_crossref.R                       # Late timepoint (default)
#   Rscript scripts/ctcf_stripe_crossref.R --timepoint early     # Early timepoint
#   Rscript scripts/ctcf_stripe_crossref.R --timepoint both      # Both timepoints
#   Rscript scripts/ctcf_stripe_crossref.R --tolerance 25000     # Custom anchor tolerance (bp)
#
# Output:
#   output/ctcf_stripe_crossref/{timepoint}/
#     tables/     - TSV files with overlap statistics
#     plots/      - Multi-format plots (PDF, SVG, JPEG)
#     summary_report.txt

# ==============================================================================
# 1. PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(futile.logger)
})

# Load multi-format output utility
source("data/scripts/_shared/multi_format_output.R")  # Original: source("scripts/utils/multi_format_output.R")

# Suppress VennDiagram logging
flog.threshold(ERROR)

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

BASE_DIR <- getwd()

DEFAULT_TOLERANCE <- 25000  # 25kb for anchor-stripe overlap matching

# Timepoint-specific file mappings
TIMEPOINT_CONFIG <- list(
  late = list(
    loops_file = file.path(BASE_DIR, "peaks/loop_annotation_extended/late/extended_characterized_loops.tsv"),  # Note: repo-relative path, not bundled in data/
    stripes_dir = file.path(BASE_DIR, "stripes/outputs/late"),  # Note: repo-relative path, not bundled in data/
    output_dir = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/ctcf_stripe_crossref/late
    output_dir_tsvs = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/ctcf_stripe_crossref/late (tables)
    output_dir_plots = file.path(BASE_DIR, "data/plots/supplemental"),  # Original: output/ctcf_stripe_crossref/late (plots)
    label = "Late Timepoint"
  ),
  early = list(
    loops_file = file.path(BASE_DIR, "peaks/loop_annotation_extended/early/extended_characterized_loops.tsv"),  # Note: repo-relative path, not bundled in data/
    stripes_dir = file.path(BASE_DIR, "stripes/outputs/early"),  # Note: repo-relative path, not bundled in data/
    output_dir = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/ctcf_stripe_crossref/early
    output_dir_tsvs = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/ctcf_stripe_crossref/early (tables)
    output_dir_plots = file.path(BASE_DIR, "data/plots/supplemental"),  # Original: output/ctcf_stripe_crossref/early (plots)
    label = "Early Timepoint"
  )
)

# Color palette
COLORS <- list(
  ctcf_with_stripe = "#E41A1C",      # Red - CTCF loops with stripe overlap
  ctcf_without_stripe = "#377EB8",   # Blue - CTCF loops without stripe overlap
  non_ctcf = "#4DAF4A",              # Green - Non-CTCF loops
  lost_stripe = "#984EA3",           # Purple - Lost stripes
  gained_stripe = "#FF7F00",         # Orange - Gained stripes
  overlap = "#FFFF33"                # Yellow - Overlapping regions
)

# ==============================================================================
# 3. HELPER FUNCTIONS - Data Loading
# ==============================================================================

#' Load extended characterized loops with CTCF filtering
#' @param file_path Path to extended_characterized_loops.tsv
#' @return tibble of loops
load_loops <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("Loops file not found: %s", file_path))
  }

  loops <- read_tsv(file_path, show_col_types = FALSE)

  # Verify required columns
  required_cols <- c("loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
                     "direction", "anchor1_CTCF_overlap", "anchor2_CTCF_overlap")
  missing <- setdiff(required_cols, colnames(loops))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }

  cat(sprintf("  Loaded %d loops from %s\n", nrow(loops), basename(file_path)))
  cat(sprintf("    - down_in_mutant (lost): %d\n", sum(loops$direction == "down_in_mutant")))
  cat(sprintf("    - up_in_mutant (gained): %d\n", sum(loops$direction == "up_in_mutant")))

  return(loops)
}

#' Filter loops for lost CTCF-anchored loops
#' @param loops tibble of loops
#' @return tibble of lost CTCF loops
filter_ctcf_lost_loops <- function(loops) {
  ctcf_lost <- loops %>%
    filter(direction == "down_in_mutant") %>%
    filter(anchor1_CTCF_overlap | anchor2_CTCF_overlap)

  cat(sprintf("  Lost loops with CTCF at any anchor: %d\n", nrow(ctcf_lost)))
  return(ctcf_lost)
}

#' Filter loops for lost non-CTCF loops (control comparison)
#' @param loops tibble of loops
#' @return tibble of lost non-CTCF loops
filter_non_ctcf_lost_loops <- function(loops) {
  non_ctcf_lost <- loops %>%
    filter(direction == "down_in_mutant") %>%
    filter(!anchor1_CTCF_overlap & !anchor2_CTCF_overlap)

  cat(sprintf("  Lost loops without CTCF at either anchor: %d\n", nrow(non_ctcf_lost)))
  return(non_ctcf_lost)
}

#' Load stripe data (lost and gained)
#' @param stripes_dir Path to stripes output directory
#' @return list with lost and gained stripes
load_stripes <- function(stripes_dir) {
  lost_file <- file.path(stripes_dir, "04_stripes_lost.tsv")
  gained_file <- file.path(stripes_dir, "04_stripes_gained.tsv")

  if (!file.exists(lost_file)) {
    stop(sprintf("Lost stripes file not found: %s", lost_file))
  }
  if (!file.exists(gained_file)) {
    stop(sprintf("Gained stripes file not found: %s", gained_file))
  }

  lost <- read_tsv(lost_file, show_col_types = FALSE)
  gained <- read_tsv(gained_file, show_col_types = FALSE)

  cat(sprintf("  Loaded stripes from %s\n", basename(stripes_dir)))
  cat(sprintf("    - Lost stripes: %d\n", nrow(lost)))
  cat(sprintf("    - Gained stripes: %d\n", nrow(gained)))

  return(list(lost = lost, gained = gained))
}

# ==============================================================================
# 4. GRANGES CONVERSION FUNCTIONS
# ==============================================================================

#' Convert loop anchors to GRanges (both anchors extracted)
#' @param loops tibble of loops
#' @return GRanges with all anchors
loop_anchors_to_granges <- function(loops) {
  if (nrow(loops) == 0) {
    return(GRanges())
  }

  # Extract anchor1
  anchor1_gr <- GRanges(
    seqnames = loops$chr1,
    ranges = IRanges(start = loops$start1, end = loops$end1),
    loop_id = loops$loop_id,
    anchor_num = 1L
  )

  # Extract anchor2
  anchor2_gr <- GRanges(
    seqnames = loops$chr2,
    ranges = IRanges(start = loops$start2, end = loops$end2),
    loop_id = loops$loop_id,
    anchor_num = 2L
  )

  # Combine
  all_anchors <- c(anchor1_gr, anchor2_gr)
  return(all_anchors)
}

#' Convert stripe anchors to GRanges
#' @param stripes tibble of stripes
#' @return GRanges of stripe anchors
stripe_anchors_to_granges <- function(stripes) {
  if (nrow(stripes) == 0) {
    return(GRanges())
  }

  GRanges(
    seqnames = stripes$chr,
    ranges = IRanges(start = stripes$anchor_x1, end = stripes$anchor_x2),
    stripe_id = stripes$stripe_id,
    confidence = stripes$confidence
  )
}

# ==============================================================================
# 5. OVERLAP ANALYSIS FUNCTIONS
# ==============================================================================

#' Find overlaps between loop anchors and stripe anchors
#' @param loop_anchors_gr GRanges of loop anchors
#' @param stripe_anchors_gr GRanges of stripe anchors
#' @param tolerance Maximum gap for overlap detection (bp)
#' @return list with overlap statistics
find_loop_stripe_overlaps <- function(loop_anchors_gr, stripe_anchors_gr, tolerance = 25000) {
  if (length(loop_anchors_gr) == 0 || length(stripe_anchors_gr) == 0) {
    return(list(
      overlaps = NULL,
      loop_ids_with_stripe = character(0),
      stripe_ids_with_loop = character(0),
      n_loop_anchors_overlapping = 0,
      n_stripes_overlapping = 0
    ))
  }

  # Find overlaps with tolerance window
  overlaps <- findOverlaps(loop_anchors_gr, stripe_anchors_gr, maxgap = tolerance)

  if (length(overlaps) == 0) {
    return(list(
      overlaps = overlaps,
      loop_ids_with_stripe = character(0),
      stripe_ids_with_loop = character(0),
      n_loop_anchors_overlapping = 0,
      n_stripes_overlapping = 0
    ))
  }

  # Get unique loop and stripe IDs that overlap
  loop_ids_with_stripe <- unique(as.character(loop_anchors_gr$loop_id[queryHits(overlaps)]))
  stripe_ids_with_loop <- unique(as.character(stripe_anchors_gr$stripe_id[subjectHits(overlaps)]))

  return(list(
    overlaps = overlaps,
    loop_ids_with_stripe = loop_ids_with_stripe,
    stripe_ids_with_loop = stripe_ids_with_loop,
    n_loop_anchors_overlapping = length(unique(queryHits(overlaps))),
    n_stripes_overlapping = length(stripe_ids_with_loop)
  ))
}

#' Compute Fisher's exact test for CTCF enrichment at stripe-defect regions
#' @param ctcf_loops tibble of CTCF lost loops
#' @param non_ctcf_loops tibble of non-CTCF lost loops
#' @param ctcf_with_stripe loop_ids of CTCF loops overlapping stripes
#' @param non_ctcf_with_stripe loop_ids of non-CTCF loops overlapping stripes
#' @return Fisher's test result
compute_fisher_enrichment <- function(ctcf_loops, non_ctcf_loops,
                                       ctcf_with_stripe, non_ctcf_with_stripe) {
  # Build 2x2 contingency table
  # Rows: CTCF status (CTCF, non-CTCF)
  # Columns: Stripe overlap (Yes, No)

  ctcf_yes <- sum(ctcf_loops$loop_id %in% ctcf_with_stripe)
  ctcf_no <- nrow(ctcf_loops) - ctcf_yes
  non_ctcf_yes <- sum(non_ctcf_loops$loop_id %in% non_ctcf_with_stripe)
  non_ctcf_no <- nrow(non_ctcf_loops) - non_ctcf_yes

  contingency <- matrix(
    c(ctcf_yes, ctcf_no, non_ctcf_yes, non_ctcf_no),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      c("CTCF_loops", "non_CTCF_loops"),
      c("Stripe_overlap", "No_overlap")
    )
  )

  cat("\n  Contingency table:\n")
  print(contingency)

  # Fisher's exact test
  fisher_result <- fisher.test(contingency)

  return(list(
    contingency = contingency,
    fisher_test = fisher_result,
    ctcf_pct = 100 * ctcf_yes / nrow(ctcf_loops),
    non_ctcf_pct = 100 * non_ctcf_yes / nrow(non_ctcf_loops)
  ))
}

#' Run permutation test for overlap enrichment
#' @param loop_anchors_gr GRanges of loop anchors
#' @param stripe_anchors_gr GRanges of stripe anchors
#' @param tolerance overlap tolerance in bp
#' @param n_perm number of permutations
#' @return list with permutation test results
run_permutation_test <- function(loop_anchors_gr, stripe_anchors_gr,
                                  tolerance = 25000, n_perm = 1000) {
  if (length(loop_anchors_gr) == 0 || length(stripe_anchors_gr) == 0) {
    return(list(observed = 0, null_mean = NA, null_sd = NA, p_value = NA))
  }

  # Observed overlap count (unique loops)
  obs_overlaps <- findOverlaps(loop_anchors_gr, stripe_anchors_gr, maxgap = tolerance)
  observed <- length(unique(loop_anchors_gr$loop_id[queryHits(obs_overlaps)]))

  cat(sprintf("  Running permutation test (n=%d)...\n", n_perm))

  # Get chromosome sizes for bounded shuffling
  chr_sizes <- seqlengths(loop_anchors_gr)
  if (all(is.na(chr_sizes))) {
    # Estimate from data
    chr_info <- data.frame(
      chr = as.character(seqnames(loop_anchors_gr)),
      max_pos = end(loop_anchors_gr)
    ) %>%
      group_by(chr) %>%
      summarize(size = max(max_pos) + 1000000, .groups = "drop")
    chr_sizes <- setNames(chr_info$size, chr_info$chr)
  }

  # Permutation: shuffle stripe positions within chromosomes
  set.seed(42)
  null_dist <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    # Shuffle stripes by random shift within chromosome bounds
    shifts <- sample(-500000:500000, length(stripe_anchors_gr), replace = TRUE)
    shuffled_start <- pmax(1, start(stripe_anchors_gr) + shifts)
    shuffled_end <- shuffled_start + width(stripe_anchors_gr) - 1

    shuffled_gr <- GRanges(
      seqnames = seqnames(stripe_anchors_gr),
      ranges = IRanges(start = shuffled_start, end = shuffled_end)
    )

    perm_overlaps <- findOverlaps(loop_anchors_gr, shuffled_gr, maxgap = tolerance)
    null_dist[i] <- length(unique(loop_anchors_gr$loop_id[queryHits(perm_overlaps)]))
  }

  # Calculate p-value (one-tailed: observed >= null)
  p_value <- (sum(null_dist >= observed) + 1) / (n_perm + 1)

  return(list(
    observed = observed,
    null_mean = mean(null_dist),
    null_sd = sd(null_dist),
    null_dist = null_dist,
    p_value = p_value,
    z_score = (observed - mean(null_dist)) / sd(null_dist)
  ))
}

#' Calculate distance to nearest stripe for non-overlapping loops
#' @param loop_anchors_gr GRanges of loop anchors
#' @param stripe_anchors_gr GRanges of stripe anchors
#' @param overlapping_loop_ids loop_ids that already overlap
#' @return tibble with distance information
calculate_distance_to_nearest_stripe <- function(loop_anchors_gr, stripe_anchors_gr,
                                                   overlapping_loop_ids) {
  if (length(loop_anchors_gr) == 0 || length(stripe_anchors_gr) == 0) {
    return(tibble(loop_id = character(), min_distance = numeric()))
  }

  # Get non-overlapping loop anchors
  non_overlap_idx <- !(loop_anchors_gr$loop_id %in% overlapping_loop_ids)
  non_overlap_gr <- loop_anchors_gr[non_overlap_idx]

  if (length(non_overlap_gr) == 0) {
    return(tibble(loop_id = character(), min_distance = numeric()))
  }

  # Find nearest stripe for each non-overlapping anchor
  nearest_idx <- nearest(non_overlap_gr, stripe_anchors_gr)

  # Handle NAs in nearest_idx (no hit on same chromosome)
  valid_idx <- !is.na(nearest_idx)

  if (sum(valid_idx) == 0) {
    return(tibble(loop_id = character(), min_distance = numeric()))
  }

  # Calculate distances only for valid pairs
  distances <- rep(NA_real_, length(non_overlap_gr))
  distances[valid_idx] <- distance(non_overlap_gr[valid_idx],
                                    stripe_anchors_gr[nearest_idx[valid_idx]])

  result <- tibble(
    loop_id = as.character(non_overlap_gr$loop_id),
    anchor_num = non_overlap_gr$anchor_num,
    min_distance = distances
  ) %>%
    filter(!is.na(min_distance)) %>%
    group_by(loop_id) %>%
    summarize(min_distance = min(min_distance, na.rm = TRUE), .groups = "drop")

  return(result)
}

# ==============================================================================
# 6. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create Venn diagram showing overlap
#' @param ctcf_loops CTCF lost loops
#' @param loops_with_stripe loop_ids overlapping stripes
#' @param output_path base path for output
create_venn_diagram <- function(ctcf_loops, loops_with_stripe, output_path) {
  # Calculate sets
  ctcf_set <- unique(ctcf_loops$loop_id)
  stripe_overlap_set <- loops_with_stripe

  n_ctcf_only <- length(setdiff(ctcf_set, stripe_overlap_set))
  n_overlap <- length(intersect(ctcf_set, stripe_overlap_set))
  n_stripe_only <- length(setdiff(stripe_overlap_set, ctcf_set))  # Should be 0 for our analysis

  # Create Venn diagram
  venn <- draw.pairwise.venn(
    area1 = length(ctcf_set),
    area2 = length(stripe_overlap_set),
    cross.area = n_overlap,
    category = c("Lost CTCF\nLoops", "Loops at\nLost Stripes"),
    fill = c(COLORS$ctcf_without_stripe, COLORS$lost_stripe),
    alpha = 0.5,
    cat.pos = c(-30, 30),
    cat.dist = c(0.05, 0.05),
    cat.cex = 1.2,
    cex = 1.5,
    fontfamily = "sans",
    cat.fontfamily = "sans",
    margin = 0.1,
    ind = FALSE
  )

  # Save in multiple formats
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  figure_name <- basename(output_path)
  output_dir <- file.path(dirname(output_path), figure_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # PDF
  pdf(file.path(output_dir, paste0(figure_name, ".pdf")), width = 8, height = 6)
  grid.draw(venn)
  grid.text("Lost CTCF Loops vs Lost Stripe Regions",
            y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()

  # SVG
  svglite::svglite(file.path(output_dir, paste0(figure_name, ".svg")), width = 8, height = 6)
  grid.draw(venn)
  grid.text("Lost CTCF Loops vs Lost Stripe Regions",
            y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()

  # JPEG
  jpeg(file.path(output_dir, paste0(figure_name, ".jpg")),
       width = 8 * 300, height = 6 * 300, res = 300, quality = 95)
  grid.draw(venn)
  grid.text("Lost CTCF Loops vs Lost Stripe Regions",
            y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()

  cat(sprintf("  Saved: %s/{pdf,svg,jpg}\n", figure_name))

  return(list(n_ctcf = length(ctcf_set), n_overlap = n_overlap))
}

#' Create enrichment bar plot comparing CTCF vs non-CTCF loops
#' @param enrichment_result Result from compute_fisher_enrichment
#' @param gained_pct Percentage of loops overlapping gained stripes (control)
#' @return ggplot object
create_enrichment_barplot <- function(enrichment_result, gained_ctcf_pct = 0, gained_non_ctcf_pct = 0) {
  plot_data <- tibble(
    category = c("CTCF Loops", "CTCF Loops", "non-CTCF Loops", "non-CTCF Loops"),
    stripe_type = c("Lost Stripes", "Gained Stripes", "Lost Stripes", "Gained Stripes"),
    percentage = c(enrichment_result$ctcf_pct, gained_ctcf_pct,
                   enrichment_result$non_ctcf_pct, gained_non_ctcf_pct)
  ) %>%
    mutate(
      category = factor(category, levels = c("CTCF Loops", "non-CTCF Loops")),
      stripe_type = factor(stripe_type, levels = c("Lost Stripes", "Gained Stripes"))
    )

  p <- ggplot(plot_data, aes(x = category, y = percentage, fill = stripe_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("Lost Stripes" = COLORS$lost_stripe,
                                  "Gained Stripes" = COLORS$gained_stripe),
                      name = "Stripe Type") +
    labs(
      title = "Overlap with Stripe Defect Regions",
      subtitle = sprintf("Fisher's exact p = %.2e, OR = %.2f",
                        enrichment_result$fisher_test$p.value,
                        enrichment_result$fisher_test$estimate),
      x = "Lost Loop Category",
      y = "% with Stripe Overlap"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "right"
    ) +
    coord_cartesian(ylim = c(0, max(plot_data$percentage) * 1.2))

  return(p)
}

#' Create distance distribution histogram
#' @param distance_data tibble with distance information
#' @return ggplot object
create_distance_histogram <- function(distance_data) {
  if (nrow(distance_data) == 0 || all(is.na(distance_data$min_distance))) {
    return(NULL)
  }

  # Filter out infinite values
  plot_data <- distance_data %>%
    filter(!is.na(min_distance) & is.finite(min_distance))

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  median_dist <- median(plot_data$min_distance, na.rm = TRUE)

  p <- ggplot(plot_data, aes(x = min_distance / 1e6)) +
    geom_histogram(bins = 50, fill = COLORS$ctcf_without_stripe,
                   color = "black", linewidth = 0.3, alpha = 0.7) +
    geom_vline(xintercept = median_dist / 1e6, linetype = "dashed",
               color = "red", linewidth = 1) +
    annotate("text", x = median_dist / 1e6, y = Inf,
             label = sprintf("Median: %.1f Mb", median_dist / 1e6),
             vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
    labs(
      title = "Distance to Nearest Lost Stripe",
      subtitle = "For non-overlapping CTCF lost loops",
      x = "Distance to Nearest Lost Stripe (Mb)",
      y = "Count"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )

  return(p)
}

#' Create permutation test null distribution plot
#' @param perm_result Result from run_permutation_test
#' @return ggplot object
create_permutation_plot <- function(perm_result) {
  if (is.null(perm_result$null_dist) || length(perm_result$null_dist) == 0) {
    return(NULL)
  }

  plot_data <- tibble(null_count = perm_result$null_dist)

  p <- ggplot(plot_data, aes(x = null_count)) +
    geom_histogram(bins = 30, fill = "gray70", color = "black", linewidth = 0.3) +
    geom_vline(xintercept = perm_result$observed, color = "red",
               linewidth = 1.5, linetype = "solid") +
    geom_vline(xintercept = perm_result$null_mean, color = "blue",
               linewidth = 1, linetype = "dashed") +
    annotate("text", x = perm_result$observed, y = Inf,
             label = sprintf("Observed: %d", perm_result$observed),
             vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
    labs(
      title = "Permutation Test: Loop-Stripe Overlap",
      subtitle = sprintf("p = %.4f, Z = %.2f", perm_result$p_value, perm_result$z_score),
      x = "Number of Overlapping Loops (permuted)",
      y = "Count"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )

  return(p)
}

# ==============================================================================
# 7. MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Process a single timepoint
#' @param timepoint "late" or "early"
#' @param tolerance Anchor overlap tolerance in bp
process_timepoint <- function(timepoint, tolerance = 25000) {
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

  # ===========================================================================
  # Step 1: Load Data
  # ===========================================================================
  cat("\n[Step 1] Loading data...\n")

  loops <- load_loops(config$loops_file)
  stripes <- load_stripes(config$stripes_dir)

  # Filter loops
  ctcf_lost_loops <- filter_ctcf_lost_loops(loops)
  non_ctcf_lost_loops <- filter_non_ctcf_lost_loops(loops)

  # Also get gained loops for comparison
  ctcf_gained_loops <- loops %>%
    filter(direction == "up_in_mutant") %>%
    filter(anchor1_CTCF_overlap | anchor2_CTCF_overlap)
  non_ctcf_gained_loops <- loops %>%
    filter(direction == "up_in_mutant") %>%
    filter(!anchor1_CTCF_overlap & !anchor2_CTCF_overlap)

  cat(sprintf("  Gained loops with CTCF: %d\n", nrow(ctcf_gained_loops)))
  cat(sprintf("  Gained loops without CTCF: %d\n", nrow(non_ctcf_gained_loops)))

  # ===========================================================================
  # Step 2: Convert to GRanges
  # ===========================================================================
  cat("\n[Step 2] Converting to GRanges...\n")

  ctcf_anchors_gr <- loop_anchors_to_granges(ctcf_lost_loops)
  non_ctcf_anchors_gr <- loop_anchors_to_granges(non_ctcf_lost_loops)
  lost_stripe_gr <- stripe_anchors_to_granges(stripes$lost)
  gained_stripe_gr <- stripe_anchors_to_granges(stripes$gained)

  cat(sprintf("  CTCF lost loop anchors: %d\n", length(ctcf_anchors_gr)))
  cat(sprintf("  Non-CTCF lost loop anchors: %d\n", length(non_ctcf_anchors_gr)))
  cat(sprintf("  Lost stripe anchors: %d\n", length(lost_stripe_gr)))
  cat(sprintf("  Gained stripe anchors: %d\n", length(gained_stripe_gr)))

  # ===========================================================================
  # Step 3: Find Overlaps with Lost Stripes
  # ===========================================================================
  cat("\n[Step 3] Finding overlaps with lost stripes...\n")
  cat(sprintf("  Using tolerance: %d bp\n", tolerance))

  ctcf_lost_overlap <- find_loop_stripe_overlaps(ctcf_anchors_gr, lost_stripe_gr, tolerance)
  non_ctcf_lost_overlap <- find_loop_stripe_overlaps(non_ctcf_anchors_gr, lost_stripe_gr, tolerance)

  cat(sprintf("  CTCF loops overlapping lost stripes: %d / %d (%.1f%%)\n",
              length(ctcf_lost_overlap$loop_ids_with_stripe),
              nrow(ctcf_lost_loops),
              100 * length(ctcf_lost_overlap$loop_ids_with_stripe) / nrow(ctcf_lost_loops)))
  cat(sprintf("  Non-CTCF loops overlapping lost stripes: %d / %d (%.1f%%)\n",
              length(non_ctcf_lost_overlap$loop_ids_with_stripe),
              nrow(non_ctcf_lost_loops),
              100 * length(non_ctcf_lost_overlap$loop_ids_with_stripe) / nrow(non_ctcf_lost_loops)))
  cat(sprintf("  Lost stripes with CTCF loop overlap: %d / %d\n",
              length(ctcf_lost_overlap$stripe_ids_with_loop),
              nrow(stripes$lost)))

  # ===========================================================================
  # Step 4: Find Overlaps with Gained Stripes (Control)
  # ===========================================================================
  cat("\n[Step 4] Finding overlaps with gained stripes (control)...\n")

  ctcf_gained_overlap <- find_loop_stripe_overlaps(ctcf_anchors_gr, gained_stripe_gr, tolerance)
  non_ctcf_gained_overlap <- find_loop_stripe_overlaps(non_ctcf_anchors_gr, gained_stripe_gr, tolerance)

  cat(sprintf("  CTCF lost loops overlapping gained stripes: %d (%.1f%%)\n",
              length(ctcf_gained_overlap$loop_ids_with_stripe),
              100 * length(ctcf_gained_overlap$loop_ids_with_stripe) / max(1, nrow(ctcf_lost_loops))))
  cat(sprintf("  Non-CTCF lost loops overlapping gained stripes: %d (%.1f%%)\n",
              length(non_ctcf_gained_overlap$loop_ids_with_stripe),
              100 * length(non_ctcf_gained_overlap$loop_ids_with_stripe) / max(1, nrow(non_ctcf_lost_loops))))

  # ===========================================================================
  # Step 5: Statistical Testing
  # ===========================================================================
  cat("\n[Step 5] Statistical testing...\n")

  # Fisher's exact test: CTCF vs non-CTCF enrichment at lost stripes
  cat("\n  --- Fisher's Exact Test: CTCF enrichment at lost stripes ---\n")
  fisher_result <- compute_fisher_enrichment(
    ctcf_lost_loops, non_ctcf_lost_loops,
    ctcf_lost_overlap$loop_ids_with_stripe,
    non_ctcf_lost_overlap$loop_ids_with_stripe
  )

  cat(sprintf("\n  Fisher's exact test p-value: %.2e\n", fisher_result$fisher_test$p.value))
  cat(sprintf("  Odds ratio: %.2f\n", fisher_result$fisher_test$estimate))
  cat(sprintf("  95%% CI: [%.2f, %.2f]\n",
              fisher_result$fisher_test$conf.int[1],
              fisher_result$fisher_test$conf.int[2]))

  # Permutation test
  cat("\n  --- Permutation Test: CTCF loop enrichment at lost stripes ---\n")
  perm_result <- run_permutation_test(ctcf_anchors_gr, lost_stripe_gr, tolerance, n_perm = 1000)

  cat(sprintf("  Observed overlapping loops: %d\n", perm_result$observed))
  cat(sprintf("  Null distribution mean: %.1f (SD: %.1f)\n", perm_result$null_mean, perm_result$null_sd))
  cat(sprintf("  Permutation p-value: %.4f\n", perm_result$p_value))
  cat(sprintf("  Z-score: %.2f\n", perm_result$z_score))

  # ===========================================================================
  # Step 6: Distance Analysis for Non-Overlapping Loops
  # ===========================================================================
  cat("\n[Step 6] Distance analysis for non-overlapping loops...\n")

  distance_data <- calculate_distance_to_nearest_stripe(
    ctcf_anchors_gr, lost_stripe_gr,
    ctcf_lost_overlap$loop_ids_with_stripe
  )

  if (nrow(distance_data) > 0 && any(!is.na(distance_data$min_distance))) {
    cat(sprintf("  Median distance to nearest stripe: %.0f bp\n",
                median(distance_data$min_distance, na.rm = TRUE)))
    cat(sprintf("  Loops analyzed: %d\n", nrow(distance_data)))
  }

  # ===========================================================================
  # Step 7: Save Tables
  # ===========================================================================
  cat("\n[Step 7] Saving tables...\n")

  # CTCF loops with stripe overlap
  ctcf_with_stripe <- ctcf_lost_loops %>%
    filter(loop_id %in% ctcf_lost_overlap$loop_ids_with_stripe)
  write_tsv(ctcf_with_stripe, file.path(tables_dir, "ctcf_stripe_ctcf_loops_with_stripe_overlap.tsv"))  # Original: tables/ctcf_loops_with_stripe_overlap.tsv
  cat(sprintf("  Saved: ctcf_stripe_ctcf_loops_with_stripe_overlap.tsv (%d rows)\n", nrow(ctcf_with_stripe)))

  # CTCF loops without stripe overlap
  ctcf_without_stripe <- ctcf_lost_loops %>%
    filter(!loop_id %in% ctcf_lost_overlap$loop_ids_with_stripe)
  write_tsv(ctcf_without_stripe, file.path(tables_dir, "ctcf_stripe_ctcf_loops_without_stripe_overlap.tsv"))  # Original: tables/ctcf_loops_without_stripe_overlap.tsv
  cat(sprintf("  Saved: ctcf_stripe_ctcf_loops_without_stripe_overlap.tsv (%d rows)\n", nrow(ctcf_without_stripe)))

  # Stripes with loop overlap
  stripes_with_loop <- stripes$lost %>%
    filter(stripe_id %in% ctcf_lost_overlap$stripe_ids_with_loop)
  write_tsv(stripes_with_loop, file.path(tables_dir, "ctcf_stripe_stripes_with_loop_overlap.tsv"))  # Original: tables/stripes_with_loop_overlap.tsv
  cat(sprintf("  Saved: ctcf_stripe_stripes_with_loop_overlap.tsv (%d rows)\n", nrow(stripes_with_loop)))

  # Enrichment statistics
  enrichment_stats <- tibble(
    test = c("fisher_exact", "permutation"),
    p_value = c(fisher_result$fisher_test$p.value, perm_result$p_value),
    statistic = c(fisher_result$fisher_test$estimate, perm_result$z_score),
    statistic_name = c("odds_ratio", "z_score"),
    ctcf_pct_overlap = c(fisher_result$ctcf_pct, NA),
    non_ctcf_pct_overlap = c(fisher_result$non_ctcf_pct, NA),
    observed_overlap = c(NA, perm_result$observed),
    null_mean = c(NA, perm_result$null_mean),
    null_sd = c(NA, perm_result$null_sd)
  )
  write_tsv(enrichment_stats, file.path(tables_dir, "ctcf_stripe_enrichment_statistics.tsv"))  # Original: tables/enrichment_statistics.tsv
  cat("  Saved: ctcf_stripe_enrichment_statistics.tsv\n")

  # Distance distribution
  if (nrow(distance_data) > 0) {
    write_tsv(distance_data, file.path(tables_dir, "ctcf_stripe_distance_to_nearest_stripe.tsv"))  # Original: tables/distance_to_nearest_stripe.tsv
    cat(sprintf("  Saved: ctcf_stripe_distance_to_nearest_stripe.tsv (%d rows)\n", nrow(distance_data)))
  }

  # Summary by loop type
  loop_type_summary <- ctcf_lost_loops %>%
    mutate(has_stripe_overlap = loop_id %in% ctcf_lost_overlap$loop_ids_with_stripe) %>%
    group_by(loop_type_extended, has_stripe_overlap) %>%
    summarize(n_loops = n(), .groups = "drop") %>%
    pivot_wider(names_from = has_stripe_overlap, values_from = n_loops,
                names_prefix = "stripe_overlap_", values_fill = 0)
  write_tsv(loop_type_summary, file.path(tables_dir, "ctcf_stripe_summary_by_loop_type.tsv"))  # Original: tables/summary_by_loop_type.tsv
  cat("  Saved: ctcf_stripe_summary_by_loop_type.tsv\n")

  # ===========================================================================
  # Step 8: Create Visualizations
  # ===========================================================================
  cat("\n[Step 8] Creating visualizations...\n")

  # Venn diagram
  venn_result <- create_venn_diagram(
    ctcf_lost_loops,
    ctcf_lost_overlap$loop_ids_with_stripe,
    file.path(plots_dir, "venn_overlap")
  )

  # Enrichment bar plot
  gained_ctcf_pct <- 100 * length(ctcf_gained_overlap$loop_ids_with_stripe) / max(1, nrow(ctcf_lost_loops))
  gained_non_ctcf_pct <- 100 * length(non_ctcf_gained_overlap$loop_ids_with_stripe) / max(1, nrow(non_ctcf_lost_loops))

  enrichment_plot <- create_enrichment_barplot(fisher_result, gained_ctcf_pct, gained_non_ctcf_pct)
  save_multiformat_ggplot(enrichment_plot, file.path(plots_dir, "enrichment_barplot"),
                          width = 8, height = 6)

  # Distance distribution
  dist_plot <- create_distance_histogram(distance_data)
  if (!is.null(dist_plot)) {
    save_multiformat_ggplot(dist_plot, file.path(plots_dir, "distance_distribution"),
                            width = 8, height = 6)
  }

  # Permutation test plot
  perm_plot <- create_permutation_plot(perm_result)
  if (!is.null(perm_plot)) {
    save_multiformat_ggplot(perm_plot, file.path(plots_dir, "permutation_test"),
                            width = 8, height = 6)
  }

  # ===========================================================================
  # Step 9: Generate Summary Report
  # ===========================================================================
  cat("\n[Step 9] Generating summary report...\n")

  report_lines <- c(
    "=================================================================",
    sprintf("CTCF-STRIPE CROSS-REFERENCE ANALYSIS - %s", config$label),
    "=================================================================",
    "",
    sprintf("Date: %s", Sys.time()),
    sprintf("Input loops: %s", config$loops_file),
    sprintf("Input stripes: %s", config$stripes_dir),
    sprintf("Overlap tolerance: %d bp", tolerance),
    "",
    "--- INPUT DATA SUMMARY ---",
    sprintf("Total loops: %d", nrow(loops)),
    sprintf("Lost loops (down_in_mutant): %d", sum(loops$direction == "down_in_mutant")),
    sprintf("  - With CTCF at any anchor: %d", nrow(ctcf_lost_loops)),
    sprintf("  - Without CTCF: %d", nrow(non_ctcf_lost_loops)),
    sprintf("Lost stripes: %d", nrow(stripes$lost)),
    sprintf("Gained stripes: %d", nrow(stripes$gained)),
    "",
    "--- OVERLAP RESULTS ---",
    sprintf("CTCF lost loops overlapping lost stripes: %d / %d (%.1f%%)",
            length(ctcf_lost_overlap$loop_ids_with_stripe),
            nrow(ctcf_lost_loops), fisher_result$ctcf_pct),
    sprintf("Non-CTCF lost loops overlapping lost stripes: %d / %d (%.1f%%)",
            length(non_ctcf_lost_overlap$loop_ids_with_stripe),
            nrow(non_ctcf_lost_loops), fisher_result$non_ctcf_pct),
    sprintf("Lost stripes with CTCF loop overlap: %d / %d",
            length(ctcf_lost_overlap$stripe_ids_with_loop), nrow(stripes$lost)),
    "",
    "--- STATISTICAL TESTS ---",
    "",
    "Fisher's Exact Test (CTCF vs non-CTCF enrichment):",
    sprintf("  p-value: %.2e", fisher_result$fisher_test$p.value),
    sprintf("  Odds ratio: %.2f", fisher_result$fisher_test$estimate),
    sprintf("  95%% CI: [%.2f, %.2f]",
            fisher_result$fisher_test$conf.int[1],
            fisher_result$fisher_test$conf.int[2]),
    "",
    "Permutation Test (n=1000):",
    sprintf("  Observed overlapping CTCF loops: %d", perm_result$observed),
    sprintf("  Null mean: %.1f (SD: %.1f)", perm_result$null_mean, perm_result$null_sd),
    sprintf("  p-value: %.4f", perm_result$p_value),
    sprintf("  Z-score: %.2f", perm_result$z_score),
    "",
    "--- INTERPRETATION ---"
  )

  # Add interpretation
  if (fisher_result$fisher_test$p.value < 0.05 && fisher_result$fisher_test$estimate > 1) {
    report_lines <- c(report_lines,
      "CTCF-anchored lost loops are SIGNIFICANTLY ENRICHED at lost stripe regions.",
      "This suggests coordinated disruption of CTCF/cohesin-mediated architecture in BAP1-KO.",
      sprintf("Enrichment: %.1fx higher in CTCF loops vs non-CTCF loops.", fisher_result$fisher_test$estimate)
    )
  } else if (fisher_result$fisher_test$p.value < 0.05 && fisher_result$fisher_test$estimate < 1) {
    report_lines <- c(report_lines,
      "CTCF-anchored lost loops are SIGNIFICANTLY DEPLETED at lost stripe regions.",
      "This suggests independent mechanisms for loop loss and stripe defects."
    )
  } else {
    report_lines <- c(report_lines,
      "No significant enrichment or depletion detected.",
      "CTCF loop loss and stripe defects may be independent phenomena."
    )
  }

  report_lines <- c(report_lines,
    "",
    "--- OUTPUT FILES ---",
    "tables/ctcf_loops_with_stripe_overlap.tsv",
    "tables/ctcf_loops_without_stripe_overlap.tsv",
    "tables/stripes_with_loop_overlap.tsv",
    "tables/enrichment_statistics.tsv",
    "tables/summary_by_loop_type.tsv",
    "plots/venn_overlap/",
    "plots/enrichment_barplot/",
    "plots/distance_distribution/",
    "plots/permutation_test/",
    "",
    "=================================================================",
    "END REPORT",
    "================================================================="
  )

  report_path <- file.path(config$output_dir_tsvs, "ctcf_stripe_summary_report.txt")  # Original: file.path(config$output_dir, "summary_report.txt")
  writeLines(report_lines, report_path)
  cat(sprintf("  Saved: ctcf_stripe_summary_report.txt\n"))

  cat("\n")
  cat(sprintf("Completed %s timepoint\n", timepoint))
  cat(sprintf("  Output directory: %s\n", config$output_dir))

  return(list(
    ctcf_lost_loops = ctcf_lost_loops,
    non_ctcf_lost_loops = non_ctcf_lost_loops,
    stripes = stripes,
    ctcf_overlap = ctcf_lost_overlap,
    non_ctcf_overlap = non_ctcf_lost_overlap,
    fisher_result = fisher_result,
    perm_result = perm_result,
    distance_data = distance_data
  ))
}

# ==============================================================================
# 8. COMMAND LINE INTERFACE
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
  cat("CTCF-STRIPE CROSS-REFERENCE ANALYSIS\n")
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
      cat(sprintf("\n%s timepoint output: %s\n", toupper(tp),
                  TIMEPOINT_CONFIG[[tp]]$output_dir))
    }
  }
}

# Run if called directly
if (!interactive()) {
  main()
}
