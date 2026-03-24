# scripts/shared_anchor_boundary_analysis.R
#
# Shared Anchor x TAD Boundary Cross-Reference Analysis
#
# Biological Question: Do the 212 shared anchors (loop-switching loci) colocalize
# with differential TAD boundaries more than non-shared differential loop anchors?
# Does the TAD structural change type (Merge) match the loop switching pattern?
#
# Four questions addressed:
#   Q1: Are shared anchors closer to differential TAD boundaries than non-shared?
#   Q2: Among shared anchors near boundaries, is Merge type enriched?
#   Q3: Are lost-loop partner anchors closer to boundaries than gained-loop partners?
#   Q4: At shared anchors near boundaries, does boundary enrichment direction match loop direction?
#
# Usage:
#   Rscript scripts/shared_anchor_boundary_analysis.R                     # Late (default)
#   Rscript scripts/shared_anchor_boundary_analysis.R --timepoint early   # Early timepoint
#   Rscript scripts/shared_anchor_boundary_analysis.R --timepoint both    # Both timepoints
#
# Output:
#   output/shared_anchor_boundary_analysis/{timepoint}/
#     tables/   - TSV files with per-question statistics
#     plots/    - Multi-format plots (PDF, SVG, JPEG)
#     summary_report.txt

# ==============================================================================
# 1. PACKAGE LOADING
# ==============================================================================

cat("=== Shared Anchor x TAD Boundary Cross-Reference Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Load multi-format output utility
source("data/scripts/_shared/multi_format_output.R")  # Original: source("scripts/utils/multi_format_output.R")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

BASE_DIR <- getwd()

# Distance thresholds for "near boundary" classification (bp)
DISTANCE_THRESHOLDS <- c(10000, 25000, 50000, 100000)

# Statistical parameters
N_PERMUTATIONS <- 1000
set.seed(42)

# Anchor overlap tolerance for control set exclusion
CONTROL_EXCLUSION_TOLERANCE <- 10000

# Timepoint-specific file mappings
TIMEPOINT_CONFIG <- list(
  late = list(
    shared_anchors_file = file.path(BASE_DIR, "data/tsvs/supplemental/shared_anchors.tsv"),  # Original: output/shared_anchor_analysis/late/tables/shared_anchors.tsv
    shared_loops_file = file.path(BASE_DIR, "data/tsvs/supplemental/shared_anchor_loops.tsv"),  # Original: output/shared_anchor_analysis/late/tables/shared_anchor_loops.tsv
    boundaries_file = file.path(BASE_DIR, "tads/results/late/final/tadcompare_final_filtered.tsv"),  # TODO: not in data/
    all_loops_file = file.path(BASE_DIR, "data/upstream/loop_calls/late_characterized_loops.tsv"),  # Original: outputs/250402-late_outputs/merged_loops/characterized_loops.tsv
    output_dir = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_boundary_analysis/late
    output_dir_tsvs = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_boundary_analysis/late (tables)
    output_dir_plots = file.path(BASE_DIR, "data/plots/supplemental"),  # Original: output/shared_anchor_boundary_analysis/late (plots)
    label = "Late Timepoint"
  ),
  early = list(
    shared_anchors_file = file.path(BASE_DIR, "output/shared_anchor_analysis/early/tables/shared_anchors.tsv"),  # TODO: not in data/
    shared_loops_file = file.path(BASE_DIR, "output/shared_anchor_analysis/early/tables/shared_anchor_loops.tsv"),  # TODO: not in data/
    boundaries_file = file.path(BASE_DIR, "tads/results/early/final/tadcompare_final_filtered.tsv"),  # TODO: not in data/
    all_loops_file = file.path(BASE_DIR, "outputs/250831-early_outputs/merged_loops/characterized_loops.tsv"),  # TODO: not in data/
    output_dir = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_boundary_analysis/early
    output_dir_tsvs = file.path(BASE_DIR, "data/tsvs/supplemental"),  # Original: output/shared_anchor_boundary_analysis/early (tables)
    output_dir_plots = file.path(BASE_DIR, "data/plots/supplemental"),  # Original: output/shared_anchor_boundary_analysis/early (plots)
    label = "Early Timepoint"
  )
)

# Color schemes (matching existing pipeline)
DIRECTION_COLORS <- c(
  "down_in_mutant" = "#4575b4",
  "up_in_mutant" = "#d73027"
)

BOUNDARY_TYPE_COLORS <- c(
  "Complex" = "#984EA3",
  "Merge" = "#4DAF4A",
  "Shifted" = "#FF7F00",
  "Split" = "#E41A1C",
  "Strength Change" = "#377EB8"
)

ENRICHED_COLORS <- c(
  "Matrix 1" = "#4575b4",
  "Matrix 2" = "#d73027"
)

# ==============================================================================
# 3. DATA LOADING FUNCTIONS
# ==============================================================================

#' Load shared anchor regions from TSV
#' @param file_path Path to shared_anchors.tsv
#' @return tibble with anchor_id, chr, start, end, n_lost_loops, n_gained_loops
load_shared_anchors <- function(file_path) {
  stopifnot(file.exists(file_path))
  anchors <- read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("anchor_id", "chr", "start", "end", "n_lost_loops", "n_gained_loops")
  missing <- setdiff(required_cols, colnames(anchors))
  if (length(missing) > 0) {
    stop(sprintf("shared_anchors.tsv missing required columns: %s", paste(missing, collapse = ", ")))
  }
  cat(sprintf("  Loaded %d shared anchors from %s\n", nrow(anchors), basename(file_path)))
  return(anchors)
}

#' Load shared anchor loops from TSV
#' @param file_path Path to shared_anchor_loops.tsv
#' @return tibble with loop details and anchor coordinates
load_shared_anchor_loops <- function(file_path) {
  stopifnot(file.exists(file_path))
  loops <- read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("loop_id", "direction", "anchor1_chr", "anchor1_start", "anchor1_end",
                     "anchor2_chr", "anchor2_start", "anchor2_end", "loop_distance")
  missing <- setdiff(required_cols, colnames(loops))
  if (length(missing) > 0) {
    stop(sprintf("shared_anchor_loops.tsv missing required columns: %s", paste(missing, collapse = ", ")))
  }
  cat(sprintf("  Loaded %d shared anchor loops (%d lost, %d gained)\n",
              nrow(loops),
              sum(loops$direction == "down_in_mutant"),
              sum(loops$direction == "up_in_mutant")))
  return(loops)
}

#' Convert shared anchors tibble to GRanges
#' @param anchors tibble from load_shared_anchors
#' @return GRanges with anchor metadata
shared_anchors_to_granges <- function(anchors) {
  GRanges(
    seqnames = anchors$chr,
    ranges = IRanges(start = anchors$start, end = anchors$end),
    anchor_id = anchors$anchor_id,
    n_lost_loops = anchors$n_lost_loops,
    n_gained_loops = anchors$n_gained_loops
  )
}

#' Load TAD boundaries and filter to differential only
#' @param file_path Path to tadcompare_final_filtered.tsv
#' @return tibble of differential boundaries
load_boundaries <- function(file_path) {
  stopifnot(file.exists(file_path))
  boundaries <- read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "Boundary", "Type", "Differential", "Enriched_In")
  missing <- setdiff(required_cols, colnames(boundaries))
  if (length(missing) > 0) {
    stop(sprintf("Boundaries file missing required columns: %s", paste(missing, collapse = ", ")))
  }

  # Filter to differential boundaries only
  diff_boundaries <- boundaries %>%
    filter(Type != "Non-Differential" & !is.na(Type))
  cat(sprintf("  Loaded %d total boundaries, %d differential\n", nrow(boundaries), nrow(diff_boundaries)))
  cat(sprintf("    By type: %s\n",
              paste(names(table(diff_boundaries$Type)), "=", table(diff_boundaries$Type), collapse = ", ")))
  return(diff_boundaries)
}

#' Convert boundaries to GRanges with buffer
#' @param boundaries tibble of differential boundaries
#' @param buffer Buffer around boundary position (bp, applied symmetrically)
#' @return GRanges with boundary metadata
boundaries_to_granges <- function(boundaries, buffer = 5000) {
  GRanges(
    seqnames = boundaries$chr,
    ranges = IRanges(
      start = pmax(1, boundaries$Boundary - buffer),
      end = boundaries$Boundary + buffer
    ),
    boundary_pos = boundaries$Boundary,
    type = boundaries$Type,
    enriched_in = boundaries$Enriched_In,
    gap_score = if ("Gap_Score" %in% colnames(boundaries)) boundaries$Gap_Score else NA_real_
  )
}

#' Load all differential loops for control set construction
#' @param file_path Path to characterized_loops.tsv
#' @return tibble of all differential loops
load_all_loops <- function(file_path) {
  stopifnot(file.exists(file_path))
  loops <- read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("anchor1_chr", "anchor1_start", "anchor1_end",
                     "anchor2_chr", "anchor2_start", "anchor2_end")
  missing <- setdiff(required_cols, colnames(loops))
  if (length(missing) > 0) {
    stop(sprintf("characterized_loops.tsv missing required columns: %s", paste(missing, collapse = ", ")))
  }
  cat(sprintf("  Loaded %d differential loops from %s\n", nrow(loops), basename(file_path)))
  return(loops)
}

# ==============================================================================
# 4. CONTROL SET CONSTRUCTION
# ==============================================================================

#' Build non-shared control anchors from all differential loops
#' Extract all unique anchors from 2,910 differential loops, reduce overlapping,
#' then remove any that overlap shared anchor regions (with tolerance).
#'
#' @param all_loops tibble of all differential loops
#' @param shared_anchors_gr GRanges of shared anchors
#' @param tolerance bp tolerance for overlap exclusion
#' @return GRanges of non-shared control anchors
build_nonshared_control_anchors <- function(all_loops, shared_anchors_gr, tolerance = CONTROL_EXCLUSION_TOLERANCE) {
  cat("\n[Building non-shared control anchor set]\n")

  # Extract anchor1 and anchor2 as GRanges
  a1 <- GRanges(
    seqnames = all_loops$anchor1_chr,
    ranges = IRanges(start = all_loops$anchor1_start, end = all_loops$anchor1_end)
  )
  a2 <- GRanges(
    seqnames = all_loops$anchor2_chr,
    ranges = IRanges(start = all_loops$anchor2_start, end = all_loops$anchor2_end)
  )

  # Combine and merge overlapping anchors
  all_anchors <- c(a1, a2)
  reduced <- GenomicRanges::reduce(all_anchors, ignore.strand = TRUE)
  cat(sprintf("  All differential loop anchors (reduced): %d\n", length(reduced)))

  # Remove any that overlap shared anchors (with tolerance)
  overlaps <- findOverlaps(reduced, shared_anchors_gr, maxgap = tolerance)
  shared_idx <- unique(queryHits(overlaps))

  if (length(shared_idx) > 0) {
    control_anchors <- reduced[-shared_idx]
  } else {
    control_anchors <- reduced
  }

  cat(sprintf("  Removed %d anchors overlapping shared anchors (tolerance %d bp)\n",
              length(shared_idx), tolerance))
  cat(sprintf("  Non-shared control anchors: %d\n", length(control_anchors)))

  return(control_anchors)
}

# ==============================================================================
# 5. Q1: SHARED vs NON-SHARED ANCHOR PROXIMITY TO BOUNDARIES
# ==============================================================================

#' Compute distance from each region in a GRanges to nearest differential boundary
#' @param regions_gr GRanges of query regions
#' @param boundaries_gr GRanges of differential boundaries (point-like with buffer)
#' @return numeric vector of distances (same length as regions_gr)
compute_distances_to_boundaries <- function(regions_gr, boundaries_gr) {
  dist_hits <- distanceToNearest(regions_gr, boundaries_gr, ignore.strand = TRUE)
  distances <- rep(NA_real_, length(regions_gr))
  if (length(dist_hits) > 0) {
    distances[queryHits(dist_hits)] <- mcols(dist_hits)$distance
  }
  return(distances)
}

#' Get nearest boundary metadata for each region
#' @param regions_gr GRanges of query regions
#' @param boundaries_gr GRanges of differential boundaries
#' @return tibble with nearest boundary type, enriched_in, gap_score per region
get_nearest_boundary_info <- function(regions_gr, boundaries_gr) {
  dist_hits <- distanceToNearest(regions_gr, boundaries_gr, ignore.strand = TRUE)

  n <- length(regions_gr)
  result <- tibble(
    distance = rep(NA_real_, n),
    nearest_type = rep(NA_character_, n),
    nearest_enriched_in = rep(NA_character_, n),
    nearest_gap_score = rep(NA_real_, n)
  )

  if (length(dist_hits) > 0) {
    qi <- queryHits(dist_hits)
    si <- subjectHits(dist_hits)
    result$distance[qi] <- mcols(dist_hits)$distance
    result$nearest_type[qi] <- as.character(mcols(boundaries_gr)$type[si])
    result$nearest_enriched_in[qi] <- as.character(mcols(boundaries_gr)$enriched_in[si])
    result$nearest_gap_score[qi] <- mcols(boundaries_gr)$gap_score[si]
  }

  return(result)
}

#' Q1: Compare shared vs non-shared anchor proximity to differential boundaries
#' Wilcoxon rank-sum, Fisher's exact at 4 thresholds, permutation test
#'
#' @param shared_gr GRanges of shared anchors
#' @param control_gr GRanges of non-shared control anchors
#' @param boundaries_gr GRanges of differential boundaries
#' @param boundaries_df tibble of differential boundaries (for permutation)
#' @return list with distance tables and test results
q1_proximity_comparison <- function(shared_gr, control_gr, boundaries_gr, boundaries_df) {
  cat("\n========================================\n")
  cat("[Q1] Shared vs non-shared anchor proximity to differential TAD boundaries\n")
  cat("========================================\n")

  # Compute distances
  shared_dist <- compute_distances_to_boundaries(shared_gr, boundaries_gr)
  control_dist <- compute_distances_to_boundaries(control_gr, boundaries_gr)

  # Build per-anchor distance tables with nearest boundary info
  shared_info <- get_nearest_boundary_info(shared_gr, boundaries_gr)
  shared_dist_df <- tibble(
    anchor_id = as.character(mcols(shared_gr)$anchor_id),
    chr = as.character(seqnames(shared_gr)),
    start = start(shared_gr),
    end = end(shared_gr),
    n_lost_loops = mcols(shared_gr)$n_lost_loops,
    n_gained_loops = mcols(shared_gr)$n_gained_loops,
    distance_to_boundary = shared_info$distance,
    nearest_boundary_type = shared_info$nearest_type,
    nearest_boundary_enriched_in = shared_info$nearest_enriched_in,
    nearest_boundary_gap_score = shared_info$nearest_gap_score,
    set = "shared"
  )

  control_dist_df <- tibble(
    anchor_id = paste0("ctrl_", seq_along(control_dist)),
    chr = as.character(seqnames(control_gr)),
    start = start(control_gr),
    end = end(control_gr),
    n_lost_loops = NA_integer_,
    n_gained_loops = NA_integer_,
    distance_to_boundary = control_dist,
    nearest_boundary_type = NA_character_,
    nearest_boundary_enriched_in = NA_character_,
    nearest_boundary_gap_score = NA_real_,
    set = "non-shared"
  )

  # Add nearest boundary info for control set too
  control_info <- get_nearest_boundary_info(control_gr, boundaries_gr)
  control_dist_df$nearest_boundary_type <- control_info$nearest_type
  control_dist_df$nearest_boundary_enriched_in <- control_info$nearest_enriched_in
  control_dist_df$nearest_boundary_gap_score <- control_info$nearest_gap_score

  # Wilcoxon rank-sum test
  shared_vals <- shared_dist[!is.na(shared_dist)]
  control_vals <- control_dist[!is.na(control_dist)]

  wilcox_result <- wilcox.test(shared_vals, control_vals, exact = FALSE)
  n_total <- length(shared_vals) + length(control_vals)
  z_approx <- qnorm(wilcox_result$p.value / 2)
  effect_r <- abs(z_approx) / sqrt(n_total)

  cat(sprintf("\n  Wilcoxon rank-sum test:\n"))
  cat(sprintf("    Shared anchors: n=%d, median=%s bp\n",
              length(shared_vals), format(median(shared_vals), big.mark = ",")))
  cat(sprintf("    Non-shared anchors: n=%d, median=%s bp\n",
              length(control_vals), format(median(control_vals), big.mark = ",")))
  cat(sprintf("    W=%.0f, p=%.4g, effect size r=%.3f\n",
              wilcox_result$statistic, wilcox_result$p.value, effect_r))

  # Fisher's exact at 4 thresholds
  fisher_results <- list()
  for (thresh in DISTANCE_THRESHOLDS) {
    shared_near <- sum(shared_vals <= thresh, na.rm = TRUE)
    shared_far <- length(shared_vals) - shared_near
    control_near <- sum(control_vals <= thresh, na.rm = TRUE)
    control_far <- length(control_vals) - control_near

    contingency <- matrix(c(shared_near, shared_far, control_near, control_far),
                          nrow = 2, byrow = TRUE,
                          dimnames = list(c("shared", "non-shared"), c("near", "far")))
    ft <- fisher.test(contingency)

    fisher_results[[as.character(thresh)]] <- tibble(
      threshold_bp = thresh,
      threshold_label = paste0(thresh / 1000, "kb"),
      shared_near = shared_near,
      shared_total = length(shared_vals),
      shared_prop = shared_near / length(shared_vals),
      control_near = control_near,
      control_total = length(control_vals),
      control_prop = control_near / length(control_vals),
      odds_ratio = as.numeric(ft$estimate),
      or_ci_low = ft$conf.int[1],
      or_ci_high = ft$conf.int[2],
      p_value = ft$p.value
    )

    cat(sprintf("  Fisher at %dkb: shared=%.1f%% vs control=%.1f%%, OR=%.2f (%.2f-%.2f), p=%.4g\n",
                thresh / 1000,
                shared_near / length(shared_vals) * 100,
                control_near / length(control_vals) * 100,
                as.numeric(ft$estimate), ft$conf.int[1], ft$conf.int[2], ft$p.value))
  }
  fisher_df <- bind_rows(fisher_results)

  # Permutation test at 50kb threshold
  cat("\n  Running permutation test (n=1000) at 50kb threshold...\n")
  observed_prop <- sum(shared_vals <= 50000) / length(shared_vals)
  perm_results <- run_anchor_permutation(shared_gr, boundaries_df, boundaries_gr,
                                         threshold = 50000, n_perm = N_PERMUTATIONS)

  # Combine all stats
  proximity_stats <- tibble(
    test = c("Wilcoxon", rep("Fisher", length(DISTANCE_THRESHOLDS)), "Permutation_50kb"),
    detail = c("shared_vs_nonshared",
               paste0("threshold_", DISTANCE_THRESHOLDS / 1000, "kb"),
               "50kb"),
    statistic = c(as.numeric(wilcox_result$statistic),
                  fisher_df$odds_ratio,
                  perm_results$fold_enrichment),
    p_value = c(wilcox_result$p.value,
                fisher_df$p_value,
                perm_results$p_value),
    shared_median = c(median(shared_vals), rep(NA, length(DISTANCE_THRESHOLDS) + 1)),
    control_median = c(median(control_vals), rep(NA, length(DISTANCE_THRESHOLDS) + 1)),
    effect_size = c(effect_r, fisher_df$odds_ratio, perm_results$fold_enrichment)
  )

  return(list(
    shared_dist_df = shared_dist_df,
    control_dist_df = control_dist_df,
    proximity_stats = proximity_stats,
    fisher_df = fisher_df,
    permutation = perm_results,
    wilcoxon = list(
      p_value = wilcox_result$p.value,
      statistic = as.numeric(wilcox_result$statistic),
      shared_median = median(shared_vals),
      control_median = median(control_vals),
      effect_r = effect_r
    )
  ))
}

#' Permutation test: randomize anchor positions matching chr/size distribution
#' @param shared_gr GRanges of shared anchors
#' @param boundaries_df tibble of differential boundaries
#' @param boundaries_gr GRanges of differential boundaries
#' @param threshold Distance threshold (bp)
#' @param n_perm Number of permutations
#' @return list with observed, null_distribution, p_value, fold_enrichment
run_anchor_permutation <- function(shared_gr, boundaries_df, boundaries_gr,
                                   threshold = 50000, n_perm = 1000) {
  # Observed proportion of shared anchors near boundaries
  shared_dist <- compute_distances_to_boundaries(shared_gr, boundaries_gr)
  observed_prop <- mean(shared_dist <= threshold, na.rm = TRUE)

  # Get chromosome + size distribution of shared anchors for matched sampling
  shared_df <- as.data.frame(shared_gr) %>%
    mutate(size = width)
  chr_counts <- table(as.character(seqnames(shared_gr)))

  # Genome-wide chromosome sizes (mm10 autosomes) for sampling positions
  chr_sizes <- c(
    chr1 = 195471971, chr2 = 182113224, chr3 = 160039680, chr4 = 156508116,
    chr5 = 151834684, chr6 = 149736546, chr7 = 145441459, chr8 = 129401213,
    chr9 = 124595110, chr10 = 130694993, chr11 = 122082543, chr12 = 120129022,
    chr13 = 120421639, chr14 = 124902244, chr15 = 104043685, chr16 = 98207768,
    chr17 = 94987271, chr18 = 90702639, chr19 = 61431566
  )

  null_props <- numeric(n_perm)
  median_anchor_size <- median(width(shared_gr))

  for (i in seq_len(n_perm)) {
    # For each chromosome, sample random positions matching count
    random_regions_list <- list()
    for (chrom in names(chr_counts)) {
      n_this_chr <- as.integer(chr_counts[chrom])
      max_start <- chr_sizes[chrom] - median_anchor_size
      if (is.na(max_start) || max_start < 1) next
      random_starts <- sample.int(max_start, size = n_this_chr, replace = FALSE)
      random_regions_list[[chrom]] <- GRanges(
        seqnames = chrom,
        ranges = IRanges(start = random_starts, width = median_anchor_size)
      )
    }

    if (length(random_regions_list) == 0) next
    random_gr <- do.call(c, unname(random_regions_list))

    random_dist <- compute_distances_to_boundaries(random_gr, boundaries_gr)
    null_props[i] <- mean(random_dist <= threshold, na.rm = TRUE)
  }

  # One-tailed p-value (enrichment)
  p_value <- (sum(null_props >= observed_prop) + 1) / (n_perm + 1)
  fold_enrichment <- observed_prop / mean(null_props)

  cat(sprintf("  Permutation test at %dkb:\n", threshold / 1000))
  cat(sprintf("    Observed proportion: %.3f\n", observed_prop))
  cat(sprintf("    Null mean: %.3f (sd: %.3f)\n", mean(null_props), sd(null_props)))
  cat(sprintf("    Fold enrichment: %.2f\n", fold_enrichment))
  cat(sprintf("    p-value (one-tailed): %.4g\n", p_value))

  return(list(
    observed = observed_prop,
    null_mean = mean(null_props),
    null_sd = sd(null_props),
    null_distribution = null_props,
    fold_enrichment = fold_enrichment,
    p_value = p_value,
    threshold = threshold
  ))
}

# ==============================================================================
# 6. Q2: MERGE TYPE ENRICHMENT AT SHARED ANCHORS NEAR BOUNDARIES
# ==============================================================================

#' Q2: Test whether Merge-type boundaries are enriched near shared anchors
#' Chi-square GOF + per-type Fisher's exact (BH-corrected)
#'
#' @param shared_dist_df tibble with per-shared-anchor boundary info (from Q1)
#' @param boundaries_df tibble of all differential boundaries
#' @param threshold Distance threshold (bp) for "near boundary"
#' @return list with chi-square GOF and per-type Fisher results
q2_merge_enrichment <- function(shared_dist_df, boundaries_df, threshold = 50000) {
  cat("\n========================================\n")
  cat("[Q2] Boundary type enrichment near shared anchors\n")
  cat("========================================\n")

  # Filter shared anchors that are near a boundary
  near_shared <- shared_dist_df %>%
    filter(!is.na(distance_to_boundary) & distance_to_boundary <= threshold & !is.na(nearest_boundary_type))

  cat(sprintf("  Shared anchors within %dkb of a differential boundary: %d / %d\n",
              threshold / 1000, nrow(near_shared), nrow(shared_dist_df)))

  if (nrow(near_shared) < 5) {
    cat("  WARNING: Too few shared anchors near boundaries for type analysis.\n")
    return(list(
      gof_test = NULL,
      per_type_fisher = tibble(),
      near_shared = near_shared,
      status = "insufficient_data"
    ))
  }

  # Observed boundary type distribution near shared anchors
  observed_counts <- table(near_shared$nearest_boundary_type)
  cat("  Observed boundary types near shared anchors:\n")
  for (bt in names(observed_counts)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", bt, observed_counts[bt],
                observed_counts[bt] / sum(observed_counts) * 100))
  }

  # Expected distribution: overall differential boundary type proportions
  overall_counts <- table(boundaries_df$Type)
  # Restrict to types present in observed
  common_types <- intersect(names(observed_counts), names(overall_counts))
  if (length(common_types) < 2) {
    cat("  WARNING: Fewer than 2 boundary types shared — cannot run GOF.\n")
    return(list(gof_test = NULL, per_type_fisher = tibble(), near_shared = near_shared, status = "too_few_types"))
  }

  # Chi-square goodness-of-fit
  obs_vec <- observed_counts[common_types]
  expected_probs <- overall_counts[common_types] / sum(overall_counts[common_types])

  gof <- chisq.test(obs_vec, p = as.numeric(expected_probs), rescale.p = TRUE)
  cat(sprintf("\n  Chi-square GOF (observed vs expected type distribution):\n"))
  cat(sprintf("    X² = %.2f, df = %d, p = %.4g\n", gof$statistic, gof$parameter, gof$p.value))

  # Print expected vs observed
  expected_vec <- sum(obs_vec) * expected_probs
  comparison <- tibble(
    type = common_types,
    observed = as.integer(obs_vec),
    expected = round(as.numeric(expected_vec), 1),
    obs_prop = as.numeric(obs_vec) / sum(obs_vec),
    exp_prop = as.numeric(expected_probs),
    ratio = as.numeric(obs_vec) / as.numeric(expected_vec)
  )
  cat("\n  Type comparison (observed / expected):\n")
  for (i in seq_len(nrow(comparison))) {
    r <- comparison[i, ]
    cat(sprintf("    %s: %d / %.1f (ratio=%.2f)\n", r$type, r$observed, r$expected, r$ratio))
  }

  # Per-type Fisher's exact: is type X enriched among near-shared vs all boundaries?
  fisher_results <- list()
  for (btype in common_types) {
    # 2x2: (shared-near boundaries of type X vs other) vs (all boundaries of type X vs other)
    n_shared_this <- sum(near_shared$nearest_boundary_type == btype)
    n_shared_other <- nrow(near_shared) - n_shared_this
    n_all_this <- as.integer(overall_counts[btype])
    n_all_other <- sum(overall_counts) - n_all_this

    contingency <- matrix(
      c(n_shared_this, n_shared_other, n_all_this, n_all_other),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("near_shared", "all_boundaries"), c(btype, "other"))
    )

    ft <- fisher.test(contingency)
    fisher_results[[btype]] <- tibble(
      boundary_type = btype,
      n_near_shared = n_shared_this,
      n_near_shared_total = nrow(near_shared),
      prop_near_shared = n_shared_this / nrow(near_shared),
      n_all = n_all_this,
      n_all_total = sum(overall_counts),
      prop_all = n_all_this / sum(overall_counts),
      odds_ratio = as.numeric(ft$estimate),
      or_ci_low = ft$conf.int[1],
      or_ci_high = ft$conf.int[2],
      p_value = ft$p.value
    )
  }

  fisher_df <- bind_rows(fisher_results) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)

  cat("\n  Per-type Fisher's exact (BH-corrected):\n")
  for (i in seq_len(nrow(fisher_df))) {
    r <- fisher_df[i, ]
    sig <- if (r$p_adj < 0.05) " *" else ""
    cat(sprintf("    %s: OR=%.2f (%.2f-%.2f), p=%.4g, padj=%.4g%s\n",
                r$boundary_type, r$odds_ratio, r$or_ci_low, r$or_ci_high,
                r$p_value, r$p_adj, sig))
  }

  return(list(
    gof_test = gof,
    comparison = comparison,
    per_type_fisher = fisher_df,
    near_shared = near_shared,
    status = "success"
  ))
}

# ==============================================================================
# 7. Q3: LOST vs GAINED PARTNER ANCHOR PROXIMITY
# ==============================================================================

#' Q3: For each loop at a shared anchor, identify the non-shared partner anchor
#' and compare distances to nearest boundary for lost vs gained loop partners.
#'
#' @param shared_loops tibble of shared anchor loops (604 loops)
#' @param shared_anchors_gr GRanges of shared anchors
#' @param boundaries_gr GRanges of differential boundaries
#' @param tolerance bp tolerance for identifying which end is the shared anchor
#' @return list with partner distance table and test results
q3_partner_proximity <- function(shared_loops, shared_anchors_gr, boundaries_gr, tolerance = CONTROL_EXCLUSION_TOLERANCE) {
  cat("\n========================================\n")
  cat("[Q3] Lost vs gained loop partner anchor proximity to boundaries\n")
  cat("========================================\n")

  # For each loop, determine which anchor is the shared one and which is the partner
  partner_data <- list()
  n_both_shared <- 0
  n_assigned <- 0

  for (i in seq_len(nrow(shared_loops))) {
    loop <- shared_loops[i, ]

    a1_gr <- GRanges(seqnames = loop$anchor1_chr,
                     ranges = IRanges(start = loop$anchor1_start, end = loop$anchor1_end))
    a2_gr <- GRanges(seqnames = loop$anchor2_chr,
                     ranges = IRanges(start = loop$anchor2_start, end = loop$anchor2_end))

    a1_shared <- any(countOverlaps(a1_gr, shared_anchors_gr, maxgap = tolerance) > 0)
    a2_shared <- any(countOverlaps(a2_gr, shared_anchors_gr, maxgap = tolerance) > 0)

    if (a1_shared && a2_shared) {
      n_both_shared <- n_both_shared + 1
      next
    }

    # Identify partner anchor (the non-shared one)
    if (a1_shared && !a2_shared) {
      partner_gr <- a2_gr
      partner_chr <- loop$anchor2_chr
      partner_start <- loop$anchor2_start
      partner_end <- loop$anchor2_end
    } else if (a2_shared && !a1_shared) {
      partner_gr <- a1_gr
      partner_chr <- loop$anchor1_chr
      partner_start <- loop$anchor1_start
      partner_end <- loop$anchor1_end
    } else {
      # Neither anchor overlaps a shared anchor — shouldn't happen, skip
      next
    }

    # Distance to nearest boundary
    dist_hit <- distanceToNearest(partner_gr, boundaries_gr, ignore.strand = TRUE)
    partner_dist <- if (length(dist_hit) > 0) mcols(dist_hit)$distance else NA_real_
    nearest_idx <- if (length(dist_hit) > 0) subjectHits(dist_hit) else NA_integer_
    nearest_type <- if (!is.na(nearest_idx)) as.character(mcols(boundaries_gr)$type[nearest_idx]) else NA_character_
    nearest_enriched <- if (!is.na(nearest_idx)) as.character(mcols(boundaries_gr)$enriched_in[nearest_idx]) else NA_character_

    n_assigned <- n_assigned + 1
    partner_data[[n_assigned]] <- tibble(
      loop_id = loop$loop_id,
      direction = loop$direction,
      partner_chr = partner_chr,
      partner_start = partner_start,
      partner_end = partner_end,
      loop_distance = loop$loop_distance,
      partner_dist_to_boundary = partner_dist,
      nearest_boundary_type = nearest_type,
      nearest_boundary_enriched_in = nearest_enriched
    )
  }

  partner_df <- bind_rows(partner_data)

  cat(sprintf("  Total shared anchor loops: %d\n", nrow(shared_loops)))
  cat(sprintf("  Both ends shared (excluded): %d\n", n_both_shared))
  cat(sprintf("  Partner anchors analyzed: %d\n", nrow(partner_df)))

  # Wilcoxon: lost-loop partners vs gained-loop partners
  lost_partners <- partner_df %>% filter(direction == "down_in_mutant") %>%
    pull(partner_dist_to_boundary) %>% na.omit()
  gained_partners <- partner_df %>% filter(direction == "up_in_mutant") %>%
    pull(partner_dist_to_boundary) %>% na.omit()

  cat(sprintf("  Lost-loop partners: n=%d, median=%s bp\n",
              length(lost_partners), format(median(lost_partners), big.mark = ",")))
  cat(sprintf("  Gained-loop partners: n=%d, median=%s bp\n",
              length(gained_partners), format(median(gained_partners), big.mark = ",")))

  if (length(lost_partners) >= 2 && length(gained_partners) >= 2) {
    wilcox_result <- wilcox.test(lost_partners, gained_partners, exact = FALSE)
    n_total <- length(lost_partners) + length(gained_partners)
    z_approx <- qnorm(wilcox_result$p.value / 2)
    effect_r <- abs(z_approx) / sqrt(n_total)

    cat(sprintf("  Wilcoxon: W=%.0f, p=%.4g, r=%.3f\n",
                wilcox_result$statistic, wilcox_result$p.value, effect_r))

    test_result <- tibble(
      test = "Wilcoxon_partner_dist",
      lost_n = length(lost_partners),
      lost_median = median(lost_partners),
      gained_n = length(gained_partners),
      gained_median = median(gained_partners),
      statistic = as.numeric(wilcox_result$statistic),
      p_value = wilcox_result$p.value,
      effect_size_r = effect_r
    )
  } else {
    cat("  Insufficient data for Wilcoxon test.\n")
    test_result <- tibble(
      test = "Wilcoxon_partner_dist",
      lost_n = length(lost_partners),
      lost_median = if (length(lost_partners) > 0) median(lost_partners) else NA_real_,
      gained_n = length(gained_partners),
      gained_median = if (length(gained_partners) > 0) median(gained_partners) else NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      effect_size_r = NA_real_
    )
  }

  return(list(
    partner_df = partner_df,
    test_result = test_result,
    n_both_shared = n_both_shared
  ))
}

# ==============================================================================
# 8. Q4: DIRECTION CONCORDANCE AT SHARED ANCHORS NEAR BOUNDARIES
# ==============================================================================

#' Q4: At shared anchors near boundaries, does loop direction match boundary enrichment?
#' Concordant: lost loop + ctrl-enriched boundary, gained loop + mut-enriched boundary
#'
#' @param shared_loops tibble of shared anchor loops
#' @param shared_dist_df tibble with per-shared-anchor boundary info
#' @param shared_anchors_gr GRanges of shared anchors
#' @param boundaries_gr GRanges of differential boundaries
#' @param threshold Distance threshold (bp)
#' @param tolerance bp tolerance for anchor overlap
#' @return list with concordance results
q4_direction_concordance <- function(shared_loops, shared_dist_df, shared_anchors_gr,
                                     boundaries_gr, threshold = 50000, tolerance = CONTROL_EXCLUSION_TOLERANCE) {
  cat("\n========================================\n")
  cat("[Q4] Direction concordance at shared anchors near boundaries\n")
  cat("========================================\n")

  # Identify shared anchors near boundaries
  near_anchors <- shared_dist_df %>%
    filter(!is.na(distance_to_boundary) & distance_to_boundary <= threshold)

  cat(sprintf("  Shared anchors within %dkb of boundary: %d\n", threshold / 1000, nrow(near_anchors)))

  if (nrow(near_anchors) < 5) {
    cat("  WARNING: Too few shared anchors near boundaries for concordance analysis.\n")
    return(list(concordance_df = tibble(), test = NULL, status = "insufficient_data"))
  }

  near_anchors_gr <- shared_anchors_to_granges(near_anchors %>%
    rename(n_lost_loops_x = n_lost_loops, n_gained_loops_x = n_gained_loops) %>%
    rename(n_lost_loops = n_lost_loops_x, n_gained_loops = n_gained_loops_x))

  # For each loop at a near-boundary shared anchor, get loop direction and nearest boundary enrichment
  concordance_data <- list()
  n_assigned <- 0

  for (i in seq_len(nrow(shared_loops))) {
    loop <- shared_loops[i, ]

    a1_gr <- GRanges(seqnames = loop$anchor1_chr,
                     ranges = IRanges(start = loop$anchor1_start, end = loop$anchor1_end))
    a2_gr <- GRanges(seqnames = loop$anchor2_chr,
                     ranges = IRanges(start = loop$anchor2_start, end = loop$anchor2_end))

    # Check if either anchor overlaps a near-boundary shared anchor
    a1_near <- any(countOverlaps(a1_gr, near_anchors_gr, maxgap = tolerance) > 0)
    a2_near <- any(countOverlaps(a2_gr, near_anchors_gr, maxgap = tolerance) > 0)

    if (!a1_near && !a2_near) next

    # Get boundary enrichment for the shared anchor that is near a boundary
    if (a1_near) {
      anchor_for_boundary <- a1_gr
    } else {
      anchor_for_boundary <- a2_gr
    }

    dist_hit <- distanceToNearest(anchor_for_boundary, boundaries_gr, ignore.strand = TRUE)
    if (length(dist_hit) == 0) next
    boundary_idx <- subjectHits(dist_hit)
    boundary_enriched <- as.character(mcols(boundaries_gr)$enriched_in[boundary_idx])
    boundary_dist <- mcols(dist_hit)$distance
    if (is.na(boundary_enriched) || boundary_dist > threshold) next

    # Determine concordance
    loop_dir <- loop$direction
    concordant <- (loop_dir == "down_in_mutant" & boundary_enriched == "Matrix 1") |
                  (loop_dir == "up_in_mutant" & boundary_enriched == "Matrix 2")

    n_assigned <- n_assigned + 1
    concordance_data[[n_assigned]] <- tibble(
      loop_id = loop$loop_id,
      direction = loop_dir,
      boundary_enriched_in = boundary_enriched,
      boundary_distance = boundary_dist,
      concordant = concordant
    )
  }

  concordance_df <- bind_rows(concordance_data)

  cat(sprintf("  Loops at near-boundary shared anchors: %d\n", nrow(concordance_df)))

  if (nrow(concordance_df) < 5) {
    cat("  WARNING: Too few loops for concordance analysis.\n")
    return(list(concordance_df = concordance_df, test = NULL, status = "insufficient_data"))
  }

  # Build contingency table
  contingency <- table(
    loop_direction = concordance_df$direction,
    boundary_enriched = concordance_df$boundary_enriched_in
  )

  cat("  Contingency table:\n")
  print(contingency)

  # Chi-square test (with simulation for small cells)
  if (all(dim(contingency) >= 2)) {
    chi_test <- chisq.test(contingency, simulate.p.value = TRUE, B = 2000)
    cat(sprintf("  Chi-square: X²=%.2f, p=%.4g (simulated)\n", chi_test$statistic, chi_test$p.value))
  } else {
    chi_test <- NULL
  }

  # Concordance rate
  n_concordant <- sum(concordance_df$concordant)
  n_total <- nrow(concordance_df)
  concordance_rate <- n_concordant / n_total

  # Binomial test (is concordance different from 50%?)
  binom_test <- binom.test(n_concordant, n_total, p = 0.5)

  cat(sprintf("  Concordance rate: %.1f%% (%d/%d)\n", concordance_rate * 100, n_concordant, n_total))
  cat(sprintf("  Binomial test vs 50%%: p=%.4g\n", binom_test$p.value))
  cat(sprintf("  Compare to population-level concordance: 69.6%%\n"))

  concordance_summary <- tibble(
    threshold_kb = threshold / 1000,
    n_loops = n_total,
    n_concordant = n_concordant,
    n_discordant = n_total - n_concordant,
    concordance_rate = concordance_rate,
    binom_p_value = binom_test$p.value,
    chi_sq_stat = if (!is.null(chi_test)) chi_test$statistic else NA_real_,
    chi_sq_p = if (!is.null(chi_test)) chi_test$p.value else NA_real_,
    population_concordance = 0.696
  )

  return(list(
    concordance_df = concordance_df,
    contingency = contingency,
    chi_test = chi_test,
    concordance_summary = concordance_summary,
    status = "success"
  ))
}

# ==============================================================================
# 9. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Plot 01: Shared vs non-shared distance violin
create_proximity_violin <- function(shared_dist_df, control_dist_df, wilcoxon_result, label) {
  plot_df <- bind_rows(
    shared_dist_df %>% select(set, distance_to_boundary),
    control_dist_df %>% select(set, distance_to_boundary)
  ) %>%
    filter(!is.na(distance_to_boundary)) %>%
    mutate(
      log_distance = log10(distance_to_boundary + 1),
      set_label = case_when(
        set == "shared" ~ sprintf("Shared anchors\n(n=%d)", sum(set == "shared")),
        set == "non-shared" ~ sprintf("Non-shared anchors\n(n=%d)", sum(set == "non-shared"))
      )
    )

  # Recalculate labels with correct counts
  n_shared <- sum(plot_df$set == "shared")
  n_nonshared <- sum(plot_df$set == "non-shared")
  plot_df$set_label <- ifelse(
    plot_df$set == "shared",
    sprintf("Shared anchors\n(n=%d)", n_shared),
    sprintf("Non-shared anchors\n(n=%d)", n_nonshared)
  )

  p <- ggplot(plot_df, aes(x = set_label, y = log_distance, fill = set)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = c("shared" = "#d73027", "non-shared" = "#91bfdb"), guide = "none") +
    scale_y_continuous(
      name = "Distance to nearest differential TAD boundary",
      breaks = c(3, 4, 5, 6, 7),
      labels = c("1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb")
    ) +
    labs(
      title = "Shared Anchor Proximity to Differential TAD Boundaries",
      subtitle = sprintf("%s | Wilcoxon p = %.3g | Shared median = %s, Non-shared median = %s",
                        label,
                        wilcoxon_result$p_value,
                        format(wilcoxon_result$shared_median, big.mark = ","),
                        format(wilcoxon_result$control_median, big.mark = ",")),
      x = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 11)
    )

  return(p)
}

#' Plot 02: Fisher's exact OR forest plot across 4 thresholds
create_fisher_forest_plot <- function(fisher_df, label) {
  plot_df <- fisher_df %>%
    mutate(threshold_label = factor(threshold_label,
                                    levels = paste0(DISTANCE_THRESHOLDS / 1000, "kb")))

  p <- ggplot(plot_df, aes(x = odds_ratio, y = threshold_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = or_ci_low, xmax = or_ci_high), height = 0.2, linewidth = 0.8) +
    geom_point(aes(color = p_value < 0.05), size = 4) +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey60"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "NS"),
                       name = "Significance") +
    geom_text(aes(label = sprintf("OR=%.2f\np=%.3g", odds_ratio, p_value)),
              hjust = -0.15, size = 3) +
    scale_x_log10() +
    labs(
      title = "Shared Anchor Boundary Enrichment (Fisher's Exact)",
      subtitle = sprintf("%s | Shared vs non-shared anchors", label),
      x = "Odds Ratio (log scale)",
      y = "Distance threshold"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14))

  return(p)
}

#' Plot 03: Permutation test null distribution
create_permutation_plot <- function(perm_result, label) {
  null_df <- tibble(value = perm_result$null_distribution)

  p <- ggplot(null_df, aes(x = value)) +
    geom_histogram(bins = 50, fill = "grey70", color = "grey40", alpha = 0.7) +
    geom_vline(xintercept = perm_result$observed, color = "#d73027", linewidth = 1.5, linetype = "solid") +
    annotate("text",
             x = perm_result$observed,
             y = Inf, vjust = 2, hjust = 1.1,
             label = sprintf("Observed: %.3f\nFold: %.2f\np = %.4g",
                            perm_result$observed, perm_result$fold_enrichment, perm_result$p_value),
             color = "#d73027", fontface = "bold", size = 3.5) +
    labs(
      title = "Permutation Test: Shared Anchor Boundary Proximity",
      subtitle = sprintf("%s | %d permutations at %dkb threshold",
                        label, length(perm_result$null_distribution), perm_result$threshold / 1000),
      x = sprintf("Proportion of anchors within %dkb of boundary", perm_result$threshold / 1000),
      y = "Count (permutations)"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14))

  return(p)
}

#' Plot 04: Boundary type comparison (expected vs observed)
create_type_comparison_plot <- function(q2_result, label) {
  if (is.null(q2_result$comparison) || nrow(q2_result$comparison) == 0) {
    return(NULL)
  }

  plot_df <- q2_result$comparison %>%
    pivot_longer(cols = c(obs_prop, exp_prop), names_to = "source", values_to = "proportion") %>%
    mutate(source_label = ifelse(source == "obs_prop", "Near shared anchors", "All differential boundaries"))

  p <- ggplot(plot_df, aes(x = type, y = proportion, fill = source_label)) +
    geom_col(position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Near shared anchors" = "#d73027", "All differential boundaries" = "grey70"),
                      name = "Source") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "Boundary Type Distribution: Shared Anchors vs Overall",
      subtitle = sprintf("%s | Chi-sq GOF p = %.4g",
                        label,
                        if (!is.null(q2_result$gof_test)) q2_result$gof_test$p.value else NA),
      x = "Boundary Type",
      y = "Proportion"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  return(p)
}

#' Plot 05: Lost vs gained partner distance violin
create_partner_violin <- function(q3_result, label) {
  plot_df <- q3_result$partner_df %>%
    filter(!is.na(partner_dist_to_boundary)) %>%
    mutate(
      log_distance = log10(partner_dist_to_boundary + 1),
      direction_label = ifelse(direction == "down_in_mutant", "Lost-loop partners", "Gained-loop partners")
    )

  test <- q3_result$test_result

  p <- ggplot(plot_df, aes(x = direction_label, y = log_distance, fill = direction)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = DIRECTION_COLORS, guide = "none") +
    scale_y_continuous(
      name = "Partner anchor distance to nearest boundary",
      breaks = c(3, 4, 5, 6, 7),
      labels = c("1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb")
    ) +
    labs(
      title = "Partner Anchor Proximity to Differential TAD Boundaries",
      subtitle = sprintf("%s | Wilcoxon p = %.3g | Lost median = %s, Gained median = %s",
                        label,
                        test$p_value,
                        format(test$lost_median, big.mark = ","),
                        format(test$gained_median, big.mark = ",")),
      x = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 11)
    )

  return(p)
}

#' Plot 06: Direction concordance heatmap (2x2 tile)
create_concordance_heatmap <- function(q4_result, label) {
  if (is.null(q4_result$contingency) || q4_result$status != "success") {
    return(NULL)
  }

  cont_df <- as.data.frame(q4_result$contingency) %>%
    rename(count = Freq) %>%
    mutate(
      direction_label = ifelse(loop_direction == "down_in_mutant", "Lost loops", "Gained loops"),
      boundary_label = ifelse(boundary_enriched == "Matrix 1", "Ctrl-enriched\nboundary", "Mut-enriched\nboundary"),
      concordant = (loop_direction == "down_in_mutant" & boundary_enriched == "Matrix 1") |
                   (loop_direction == "up_in_mutant" & boundary_enriched == "Matrix 2")
    )

  conc_rate <- q4_result$concordance_summary$concordance_rate

  p <- ggplot(cont_df, aes(x = boundary_label, y = direction_label, fill = count)) +
    geom_tile(color = "white", linewidth = 2) +
    geom_text(aes(label = count), size = 8, fontface = "bold") +
    # Mark concordant cells
    geom_tile(data = cont_df %>% filter(concordant),
              aes(x = boundary_label, y = direction_label),
              fill = NA, color = "#2ca02c", linewidth = 2) +
    scale_fill_gradient(low = "white", high = "#d73027", name = "Count") +
    labs(
      title = "Direction Concordance: Loop Direction x Boundary Enrichment",
      subtitle = sprintf("%s | Concordance = %.1f%% | Binom p = %.4g | Population = 69.6%%",
                        label, conc_rate * 100,
                        q4_result$concordance_summary$binom_p_value),
      x = "Nearest Boundary Enrichment",
      y = "Loop Direction"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      panel.grid = element_blank()
    )

  return(p)
}

#' Plot 07: Combined summary panel
create_combined_panel <- function(plots, label) {
  # Filter out NULLs
  valid <- plots[!sapply(plots, is.null)]
  if (length(valid) < 2) return(NULL)

  # Arrange in 2-column layout
  combined <- wrap_plots(valid, ncol = 2) +
    plot_annotation(
      title = "Shared Anchor x TAD Boundary Cross-Reference",
      subtitle = label,
      theme = theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13)
      )
    )

  return(combined)
}

# ==============================================================================
# 10. REPORT GENERATION
# ==============================================================================

#' Write summary report as plain text
generate_report <- function(q1, q2, q3, q4, config, output_dir) {
  report_path <- file.path(config$output_dir_tsvs, "shared_boundary_summary_report.txt")  # Original: file.path(output_dir, "summary_report.txt")
  label <- config$label

  lines <- c(
    "==========================================================================",
    sprintf("SHARED ANCHOR x TAD BOUNDARY CROSS-REFERENCE ANALYSIS - %s", toupper(label)),
    "==========================================================================",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "--------------------------------------------------------------------------",
    "Q1: ARE SHARED ANCHORS CLOSER TO DIFFERENTIAL TAD BOUNDARIES THAN NON-SHARED?",
    "--------------------------------------------------------------------------",
    sprintf("Shared anchors: n=%d, median distance = %s bp",
            nrow(q1$shared_dist_df), format(q1$wilcoxon$shared_median, big.mark = ",")),
    sprintf("Non-shared control anchors: n=%d, median distance = %s bp",
            nrow(q1$control_dist_df), format(q1$wilcoxon$control_median, big.mark = ",")),
    sprintf("Wilcoxon rank-sum: W=%.0f, p=%.4g, effect size r=%.3f",
            q1$wilcoxon$statistic, q1$wilcoxon$p_value, q1$wilcoxon$effect_r),
    "",
    "Fisher's exact test at distance thresholds:",
    paste(capture.output(print(as.data.frame(q1$fisher_df %>%
      select(threshold_label, shared_prop, control_prop, odds_ratio, p_value)))),
      collapse = "\n"),
    "",
    sprintf("Permutation test (50kb, n=%d): observed=%.3f, null mean=%.3f, fold=%.2f, p=%.4g",
            length(q1$permutation$null_distribution),
            q1$permutation$observed, q1$permutation$null_mean,
            q1$permutation$fold_enrichment, q1$permutation$p_value),
    ""
  )

  # Q2
  lines <- c(lines,
    "--------------------------------------------------------------------------",
    "Q2: IS MERGE-TYPE BOUNDARY ENRICHED NEAR SHARED ANCHORS?",
    "--------------------------------------------------------------------------"
  )
  if (q2$status == "success") {
    lines <- c(lines,
      sprintf("Shared anchors within 50kb of boundary: %d", nrow(q2$near_shared)),
      sprintf("Chi-square GOF: X²=%.2f, df=%d, p=%.4g",
              q2$gof_test$statistic, q2$gof_test$parameter, q2$gof_test$p.value),
      "",
      "Type comparison (observed / expected):",
      paste(capture.output(print(as.data.frame(q2$comparison %>% select(type, observed, expected, ratio)))),
            collapse = "\n"),
      "",
      "Per-type Fisher's exact (BH-corrected):",
      paste(capture.output(print(as.data.frame(q2$per_type_fisher %>%
        select(boundary_type, odds_ratio, p_value, p_adj)))),
        collapse = "\n")
    )
  } else {
    lines <- c(lines, sprintf("  Status: %s (insufficient data)", q2$status))
  }

  # Q3
  lines <- c(lines, "",
    "--------------------------------------------------------------------------",
    "Q3: ARE LOST-LOOP PARTNER ANCHORS CLOSER TO BOUNDARIES THAN GAINED-LOOP PARTNERS?",
    "--------------------------------------------------------------------------",
    sprintf("Partner anchors analyzed: %d (excluded %d both-shared)", nrow(q3$partner_df), q3$n_both_shared),
    sprintf("Lost-loop partners: n=%d, median=%s bp",
            q3$test_result$lost_n, format(q3$test_result$lost_median, big.mark = ",")),
    sprintf("Gained-loop partners: n=%d, median=%s bp",
            q3$test_result$gained_n, format(q3$test_result$gained_median, big.mark = ",")),
    sprintf("Wilcoxon: W=%.0f, p=%.4g, r=%.3f",
            q3$test_result$statistic, q3$test_result$p_value, q3$test_result$effect_size_r)
  )

  # Q4
  lines <- c(lines, "",
    "--------------------------------------------------------------------------",
    "Q4: DOES BOUNDARY ENRICHMENT DIRECTION MATCH LOOP DIRECTION AT SHARED ANCHORS?",
    "--------------------------------------------------------------------------"
  )
  if (q4$status == "success") {
    cs <- q4$concordance_summary
    lines <- c(lines,
      sprintf("Loops at shared anchors within %dkb of boundary: %d", cs$threshold_kb, cs$n_loops),
      "",
      "Contingency table:",
      paste(capture.output(print(q4$contingency)), collapse = "\n"),
      "",
      sprintf("Concordance rate: %.1f%% (%d concordant, %d discordant)",
              cs$concordance_rate * 100, cs$n_concordant, cs$n_discordant),
      sprintf("Binomial test vs 50%%: p=%.4g", cs$binom_p_value),
      sprintf("Chi-square: X²=%.2f, p=%.4g", cs$chi_sq_stat, cs$chi_sq_p),
      sprintf("Population-level concordance (all 2,910 loops): %.1f%%", cs$population_concordance * 100)
    )
  } else {
    lines <- c(lines, sprintf("  Status: %s", q4$status))
  }

  lines <- c(lines, "",
    "==========================================================================",
    "END OF REPORT",
    "=========================================================================="
  )

  writeLines(lines, report_path)
  cat(sprintf("\nReport written to: %s\n", report_path))
}

# ==============================================================================
# 11. MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Run the complete shared anchor x boundary analysis for one timepoint
#' @param config Named list from TIMEPOINT_CONFIG
run_analysis <- function(config) {
  cat(sprintf("\n*** Running analysis for %s ***\n\n", config$label))

  # Create output directories
  tables_dir <- config$output_dir_tsvs  # Original: file.path(config$output_dir, "tables")
  plots_dir <- config$output_dir_plots  # Original: file.path(config$output_dir, "plots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Load data ---
  cat("[Loading data]\n")
  shared_anchors <- load_shared_anchors(config$shared_anchors_file)
  shared_loops <- load_shared_anchor_loops(config$shared_loops_file)
  boundaries_df <- load_boundaries(config$boundaries_file)
  all_loops <- load_all_loops(config$all_loops_file)

  # Convert to GRanges
  shared_gr <- shared_anchors_to_granges(shared_anchors)
  boundaries_gr <- boundaries_to_granges(boundaries_df)

  # Build control set
  control_gr <- build_nonshared_control_anchors(all_loops, shared_gr)

  # --- Q1: Proximity comparison ---
  q1 <- q1_proximity_comparison(shared_gr, control_gr, boundaries_gr, boundaries_df)

  # Save Q1 tables
  write_tsv(q1$shared_dist_df, file.path(tables_dir, "shared_boundary_shared_anchor_boundary_distances.tsv"))  # Original: tables/shared_anchor_boundary_distances.tsv
  write_tsv(q1$control_dist_df, file.path(tables_dir, "shared_boundary_nonshared_anchor_boundary_distances.tsv"))  # Original: tables/nonshared_anchor_boundary_distances.tsv
  write_tsv(q1$proximity_stats, file.path(tables_dir, "shared_boundary_proximity_comparison_stats.tsv"))  # Original: tables/proximity_comparison_stats.tsv
  write_tsv(q1$fisher_df, file.path(tables_dir, "shared_boundary_fisher_threshold_results.tsv"))  # Original: tables/fisher_threshold_results.tsv

  # Save permutation null distribution
  perm_df <- tibble(
    iteration = seq_along(q1$permutation$null_distribution),
    null_proportion = q1$permutation$null_distribution
  )
  perm_summary <- tibble(
    observed = q1$permutation$observed,
    null_mean = q1$permutation$null_mean,
    null_sd = q1$permutation$null_sd,
    fold_enrichment = q1$permutation$fold_enrichment,
    p_value = q1$permutation$p_value,
    threshold_bp = q1$permutation$threshold
  )
  write_tsv(bind_rows(perm_summary, tibble(note = "Null distribution below")),
            file.path(tables_dir, "shared_boundary_permutation_results.tsv"))  # Original: tables/permutation_results.tsv

  cat("  Q1 tables saved.\n")

  # --- Q2: Merge enrichment ---
  q2 <- q2_merge_enrichment(q1$shared_dist_df, boundaries_df, threshold = 50000)

  # Save Q2 tables
  if (q2$status == "success") {
    write_tsv(q2$per_type_fisher, file.path(tables_dir, "shared_boundary_merge_enrichment_results.tsv"))  # Original: tables/merge_enrichment_results.tsv
    write_tsv(q2$comparison, file.path(tables_dir, "shared_boundary_boundary_type_comparison.tsv"))  # Original: tables/boundary_type_comparison.tsv
  }
  cat("  Q2 tables saved.\n")

  # --- Q3: Partner proximity ---
  q3 <- q3_partner_proximity(shared_loops, shared_gr, boundaries_gr)

  # Save Q3 tables
  write_tsv(q3$partner_df, file.path(tables_dir, "shared_boundary_partner_distances.tsv"))  # Original: tables/partner_distances.tsv
  write_tsv(q3$test_result, file.path(tables_dir, "shared_boundary_partner_distance_test.tsv"))  # Original: tables/partner_distance_test.tsv
  cat("  Q3 tables saved.\n")

  # --- Q4: Direction concordance ---
  q4 <- q4_direction_concordance(shared_loops, q1$shared_dist_df, shared_gr,
                                  boundaries_gr, threshold = 50000)

  # Save Q4 tables
  if (q4$status == "success") {
    write_tsv(q4$concordance_df, file.path(tables_dir, "shared_boundary_direction_concordance.tsv"))  # Original: tables/direction_concordance.tsv
    write_tsv(q4$concordance_summary, file.path(tables_dir, "shared_boundary_concordance_summary.tsv"))  # Original: tables/concordance_summary.tsv
  }
  cat("  Q4 tables saved.\n")

  # --- Generate plots ---
  cat("\n[Generating plots]\n")

  p1 <- create_proximity_violin(q1$shared_dist_df, q1$control_dist_df, q1$wilcoxon, config$label)
  save_multiformat_ggplot(p1, file.path(plots_dir, "01_proximity_violin"), width = 8, height = 6)

  p2 <- create_fisher_forest_plot(q1$fisher_df, config$label)
  save_multiformat_ggplot(p2, file.path(plots_dir, "02_fisher_forest_plot"), width = 9, height = 5)

  p3 <- create_permutation_plot(q1$permutation, config$label)
  save_multiformat_ggplot(p3, file.path(plots_dir, "03_permutation_test"), width = 8, height = 5)

  p4 <- create_type_comparison_plot(q2, config$label)
  if (!is.null(p4)) {
    save_multiformat_ggplot(p4, file.path(plots_dir, "04_boundary_type_comparison"), width = 9, height = 6)
  }

  p5 <- create_partner_violin(q3, config$label)
  save_multiformat_ggplot(p5, file.path(plots_dir, "05_partner_distance_violin"), width = 8, height = 6)

  p6 <- create_concordance_heatmap(q4, config$label)
  if (!is.null(p6)) {
    save_multiformat_ggplot(p6, file.path(plots_dir, "06_concordance_heatmap"), width = 7, height = 5)
  }

  # Combined panel
  all_plots <- list(p1, p2, p3, p4, p5, p6)
  p7 <- create_combined_panel(all_plots, config$label)
  if (!is.null(p7)) {
    save_multiformat_ggplot(p7, file.path(plots_dir, "07_combined_summary"), width = 18, height = 16)
  }

  # --- Generate report ---
  generate_report(q1, q2, q3, q4, config, config$output_dir_tsvs)  # Original: config$output_dir

  cat(sprintf("\n*** %s analysis complete ***\n", config$label))
  cat(sprintf("Output directory: %s\n", config$output_dir))
}

# ==============================================================================
# 12. CLI PARSING + EXECUTION
# ==============================================================================

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Parse --timepoint argument
  timepoint <- "late"  # Default
  if (length(args) >= 2) {
    tp_idx <- which(args == "--timepoint")
    if (length(tp_idx) > 0 && tp_idx < length(args)) {
      timepoint <- args[tp_idx + 1]
    }
  }

  stopifnot(timepoint %in% c("late", "early", "both"))

  cat(sprintf("Timepoint: %s\n", timepoint))
  cat(sprintf("Working directory: %s\n\n", BASE_DIR))

  if (timepoint == "both") {
    for (tp in c("late", "early")) {
      run_analysis(TIMEPOINT_CONFIG[[tp]])
      cat("\n")
    }
  } else {
    run_analysis(TIMEPOINT_CONFIG[[timepoint]])
  }

  cat("\nDone. Total time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
}

main()
