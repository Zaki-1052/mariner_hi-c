# tads/scripts/boundary_loop_crossref.R
#
# Cross-Reference Differential Loops with Differential TAD Boundaries
#
# Biological Question: Do differential chromatin loops preferentially localize
# near differential TAD boundaries? Does boundary dynamics coincide with loop changes?
#
# Usage:
#   Rscript tads/scripts/boundary_loop_crossref.R                     # Late timepoint (default)
#   Rscript tads/scripts/boundary_loop_crossref.R --timepoint early   # Early timepoint
#   Rscript tads/scripts/boundary_loop_crossref.R --timepoint both    # Both timepoints
#
# Output:
#   tads/results/{timepoint}/boundary_loop_analysis/
#     boundary_loop_overlap_summary.tsv    - Per-loop boundary proximity metrics
#     enrichment_statistics.tsv            - Fisher's exact, bootstrap, permutation p-values
#     permutation_test_results.tsv         - Detailed permutation test statistics
#     direction_concordance.tsv            - Gained/lost vs enriched_in crosstab
#     boundary_type_association.tsv        - Split/Merge/Shifted vs loop changes
#     plots/
#       distance_violin.pdf                - Distance distribution comparison
#       enrichment_barplot.pdf             - Proportion near boundaries
#       type_heatmap.pdf                   - Boundary type × direction
#       permutation_test.pdf               - Null distribution vs observed
#       example_loci.pdf                   - Browser-style examples
#     analysis_report.txt                  - Summary statistics and conclusions

# ==============================================================================
# 1. PACKAGE LOADING
# ==============================================================================

cat("=== Boundary-Loop Cross-Reference Analysis ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(pheatmap)
})

# Load multi-format output utility
# Original: source("scripts/utils/multi_format_output.R")
source("data/scripts/_shared/multi_format_output.R")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

BASE_DIR <- getwd()

# Distance thresholds for "near boundary" classification (bp)
DISTANCE_THRESHOLDS <- c(10000, 25000, 50000, 100000)

# Statistical parameters
N_PERMUTATIONS <- 1000
N_BOOTSTRAP <- 1000
FDR_THRESHOLD <- 0.05
set.seed(42)  # For reproducibility

# Timepoint-specific file mappings
TIMEPOINT_CONFIG <- list(
  late = list(
    # Original: file.path(BASE_DIR, "25042-late_outputs/merged_loops/characterized_loops.tsv")
    loops_file = file.path(BASE_DIR, "data/upstream/loop_calls/late_characterized_loops.tsv"),
    # Original: file.path(BASE_DIR, "tads/results/late/final/tadcompare_final_filtered.tsv")
    boundaries_file = file.path(BASE_DIR, "data/tsvs/figure_1_tads_boundaries_compartments/1B_late_tadcompare_differential.tsv"),
    # Original: file.path(BASE_DIR, "tads/results/late/boundary_loop_analysis")
    tsv_dir = file.path(BASE_DIR, "data/tsvs/figure_1_tads_boundaries_compartments"),
    plots_dir = file.path(BASE_DIR, "data/plots/figure_1_tads_boundaries_compartments"),
    label = "Late Timepoint"
  ),
  early = list(
    loops_file = file.path(BASE_DIR, "data/upstream/loop_calls/early_characterized_loops.tsv"),  # Original: 250831-early_outputs/merged_loops/characterized_loops.tsv
    # Original: file.path(BASE_DIR, "tads/results/early/final/tadcompare_final_filtered.tsv")
    boundaries_file = file.path(BASE_DIR, "data/tsvs/figure_1_tads_boundaries_compartments/1B_early_tadcompare_differential.tsv"),
    # Original: file.path(BASE_DIR, "tads/results/early/boundary_loop_analysis")
    tsv_dir = file.path(BASE_DIR, "data/tsvs/figure_1_tads_boundaries_compartments"),
    plots_dir = file.path(BASE_DIR, "data/plots/figure_1_tads_boundaries_compartments"),
    label = "Early Timepoint"
  )
)

# Color schemes (matching existing pipeline)
DIRECTION_COLORS <- c(
  "down_in_mutant" = "#4575b4",  # Lost (blue)
  "up_in_mutant" = "#d73027"     # Gained (red)
)

BOUNDARY_TYPE_COLORS <- c(
  "Complex" = "#984EA3",
  "Merge" = "#4DAF4A",
  "Shifted" = "#FF7F00",
  "Split" = "#E41A1C",
  "Strength Change" = "#377EB8"
)

ENRICHED_COLORS <- c(
  "Matrix 1" = "#4575b4",  # Control enriched
  "Matrix 2" = "#d73027"   # Mutant enriched
)

# ==============================================================================
# 3. HELPER FUNCTIONS - Data Loading
# ==============================================================================

#' Load characterized loops from TSV
#' @param file_path Path to characterized_loops.tsv
#' @return tibble of loops with required columns
load_loops <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("Loops file not found: %s", file_path))
  }

  loops <- read_tsv(file_path, show_col_types = FALSE)

  required_cols <- c("direction", "chr1", "start1", "end1", "chr2", "start2", "end2",
                     "logFC", "FDR", "anchor1_type", "anchor2_type", "loop_type")
  missing <- setdiff(required_cols, colnames(loops))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }

  cat(sprintf("  Loaded %d loops from %s\n", nrow(loops), basename(file_path)))
  cat(sprintf("    - down_in_mutant (lost): %d\n", sum(loops$direction == "down_in_mutant")))
  cat(sprintf("    - up_in_mutant (gained): %d\n", sum(loops$direction == "up_in_mutant")))

  return(loops)
}

#' Load TAD boundaries from TADCompare output
#' @param file_path Path to tadcompare_final_filtered.tsv
#' @return tibble of boundaries
load_boundaries <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("Boundaries file not found: %s", file_path))
  }

  boundaries <- read_tsv(file_path, show_col_types = FALSE)

  required_cols <- c("chr", "Boundary", "Type", "Differential", "Enriched_In")
  missing <- setdiff(required_cols, colnames(boundaries))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }

  cat(sprintf("  Loaded %d boundaries from %s\n", nrow(boundaries), basename(file_path)))

  # Filter to differential boundaries only
  diff_boundaries <- boundaries %>% filter(Type != "Non-Differential")
  cat(sprintf("    - Differential boundaries: %d\n", nrow(diff_boundaries)))
  cat(sprintf("    - By type: %s\n",
              paste(names(table(diff_boundaries$Type)), "=", table(diff_boundaries$Type), collapse = ", ")))

  return(diff_boundaries)
}

#' Convert loops to GRanges (both anchors)
#' @param loops tibble of loops
#' @return list with anchor1_gr, anchor2_gr, and all_anchors_gr
loops_to_granges <- function(loops) {
  # Anchor 1
  anchor1_gr <- GRanges(
    seqnames = loops$chr1,
    ranges = IRanges(start = loops$start1, end = loops$end1),
    loop_id = if ("loop_id" %in% colnames(loops)) loops$loop_id else seq_len(nrow(loops)),
    anchor_num = 1,
    direction = loops$direction,
    logFC = loops$logFC,
    FDR = loops$FDR,
    anchor_type = loops$anchor1_type
  )

  # Anchor 2
  anchor2_gr <- GRanges(
    seqnames = loops$chr2,
    ranges = IRanges(start = loops$start2, end = loops$end2),
    loop_id = if ("loop_id" %in% colnames(loops)) loops$loop_id else seq_len(nrow(loops)),
    anchor_num = 2,
    direction = loops$direction,
    logFC = loops$logFC,
    FDR = loops$FDR,
    anchor_type = loops$anchor2_type
  )

  return(list(
    anchor1 = anchor1_gr,
    anchor2 = anchor2_gr,
    all = c(anchor1_gr, anchor2_gr)
  ))
}

#' Convert boundaries to GRanges with buffer
#' @param boundaries tibble of boundaries
#' @param buffer Buffer around boundary position (bp, applied symmetrically)
#' @return GRanges of boundaries
boundaries_to_granges <- function(boundaries, buffer = 5000) {
  gr <- GRanges(
    seqnames = boundaries$chr,
    ranges = IRanges(
      start = pmax(1, boundaries$Boundary - buffer),
      end = boundaries$Boundary + buffer
    ),
    boundary_pos = boundaries$Boundary,
    type = boundaries$Type,
    enriched_in = boundaries$Enriched_In,
    gap_score = boundaries$Gap_Score
  )
  return(gr)
}

# ==============================================================================
# 4. PROXIMITY ANALYSIS FUNCTIONS
# ==============================================================================

#' Compute distance from each loop anchor to nearest differential boundary
#' @param loops_gr List with anchor1, anchor2 GRanges
#' @param boundaries_gr GRanges of differential boundaries
#' @return tibble with distance metrics per loop
compute_boundary_distances <- function(loops_gr, boundaries_gr) {
  cat("\n[Computing boundary distances]\n")

  # Distance for anchor1
  dist1 <- distanceToNearest(loops_gr$anchor1, boundaries_gr, ignore.strand = TRUE)
  anchor1_dist <- rep(NA_real_, length(loops_gr$anchor1))
  anchor1_boundary_idx <- rep(NA_integer_, length(loops_gr$anchor1))
  if (length(dist1) > 0) {
    anchor1_dist[queryHits(dist1)] <- mcols(dist1)$distance
    anchor1_boundary_idx[queryHits(dist1)] <- subjectHits(dist1)
  }

  # Distance for anchor2
  dist2 <- distanceToNearest(loops_gr$anchor2, boundaries_gr, ignore.strand = TRUE)
  anchor2_dist <- rep(NA_real_, length(loops_gr$anchor2))
  anchor2_boundary_idx <- rep(NA_integer_, length(loops_gr$anchor2))
  if (length(dist2) > 0) {
    anchor2_dist[queryHits(dist2)] <- mcols(dist2)$distance
    anchor2_boundary_idx[queryHits(dist2)] <- subjectHits(dist2)
  }

  # Get nearest boundary metadata
  get_boundary_info <- function(idx, boundaries_gr) {
    if (is.na(idx)) return(list(type = NA_character_, enriched_in = NA_character_, gap_score = NA_real_))
    list(
      type = as.character(mcols(boundaries_gr)$type[idx]),
      enriched_in = as.character(mcols(boundaries_gr)$enriched_in[idx]),
      gap_score = mcols(boundaries_gr)$gap_score[idx]
    )
  }

  # Build result tibble with pre-initialized columns
  n_loops <- length(loops_gr$anchor1)
  result <- tibble(
    loop_id = as.character(mcols(loops_gr$anchor1)$loop_id),
    direction = as.character(mcols(loops_gr$anchor1)$direction),
    logFC = mcols(loops_gr$anchor1)$logFC,
    FDR = mcols(loops_gr$anchor1)$FDR,
    anchor1_type = as.character(mcols(loops_gr$anchor1)$anchor_type),
    anchor2_type = as.character(mcols(loops_gr$anchor2)$anchor_type),
    anchor1_dist_to_boundary = anchor1_dist,
    anchor2_dist_to_boundary = anchor2_dist,
    min_dist_to_boundary = pmin(anchor1_dist, anchor2_dist, na.rm = TRUE),
    nearest_boundary_type = rep(NA_character_, n_loops),
    nearest_boundary_enriched = rep(NA_character_, n_loops),
    nearest_boundary_gap_score = rep(NA_real_, n_loops)
  )

  # Add nearest boundary metadata for anchor with minimum distance
  for (i in seq_len(nrow(result))) {
    if (is.na(result$anchor1_dist_to_boundary[i]) && is.na(result$anchor2_dist_to_boundary[i])) {
      # Already initialized to NA
      next
    } else if (is.na(result$anchor2_dist_to_boundary[i]) ||
               (!is.na(result$anchor1_dist_to_boundary[i]) &&
                result$anchor1_dist_to_boundary[i] <= result$anchor2_dist_to_boundary[i])) {
      info <- get_boundary_info(anchor1_boundary_idx[i], boundaries_gr)
      result$nearest_boundary_type[i] <- info$type
      result$nearest_boundary_enriched[i] <- info$enriched_in
      result$nearest_boundary_gap_score[i] <- info$gap_score
    } else {
      info <- get_boundary_info(anchor2_boundary_idx[i], boundaries_gr)
      result$nearest_boundary_type[i] <- info$type
      result$nearest_boundary_enriched[i] <- info$enriched_in
      result$nearest_boundary_gap_score[i] <- info$gap_score
    }
  }

  # Add "near boundary" flags for each threshold
  for (thresh in DISTANCE_THRESHOLDS) {
    col_name <- paste0("near_boundary_", thresh/1000, "kb")
    result[[col_name]] <- result$min_dist_to_boundary <= thresh
  }

  cat(sprintf("  Distance computed for %d loops\n", nrow(result)))
  cat(sprintf("  Median distance to nearest boundary: %s bp\n",
              format(median(result$min_dist_to_boundary, na.rm = TRUE), big.mark = ",")))

  return(result)
}

# ==============================================================================
# 5. STATISTICAL TESTING FUNCTIONS
# ==============================================================================

#' Run Fisher's exact test for boundary proximity enrichment
#' @param distance_df tibble from compute_boundary_distances
#' @param threshold Distance threshold (bp)
#' @return tibble with test results
run_fisher_test <- function(distance_df, threshold) {
  col_name <- paste0("near_boundary_", threshold/1000, "kb")

  # Build contingency table: differential_direction × near_boundary
  # Compare gained vs lost loops for being near boundaries
  lost <- distance_df %>% filter(direction == "down_in_mutant")
  gained <- distance_df %>% filter(direction == "up_in_mutant")

  contingency <- matrix(c(
    sum(lost[[col_name]], na.rm = TRUE),      # lost, near
    sum(!lost[[col_name]], na.rm = TRUE),     # lost, far
    sum(gained[[col_name]], na.rm = TRUE),    # gained, near
    sum(!gained[[col_name]], na.rm = TRUE)    # gained, far
  ), nrow = 2, byrow = TRUE,
  dimnames = list(c("lost", "gained"), c("near", "far")))

  test <- fisher.test(contingency)

  tibble(
    threshold_bp = threshold,
    threshold_label = paste0(threshold/1000, "kb"),
    lost_near = contingency["lost", "near"],
    lost_far = contingency["lost", "far"],
    gained_near = contingency["gained", "near"],
    gained_far = contingency["gained", "far"],
    lost_prop_near = contingency["lost", "near"] / sum(contingency["lost",]),
    gained_prop_near = contingency["gained", "near"] / sum(contingency["gained",]),
    odds_ratio = test$estimate,
    odds_ratio_ci_low = test$conf.int[1],
    odds_ratio_ci_high = test$conf.int[2],
    p_value = test$p.value
  )
}

#' Run permutation test for enrichment
#' @param distance_df tibble from compute_boundary_distances
#' @param loops_gr List with anchor GRanges
#' @param boundaries tibble of boundaries
#' @param threshold Distance threshold (bp)
#' @param n_perm Number of permutations
#' @return list with observed, null distribution, p-value
run_permutation_test <- function(distance_df, loops_gr, boundaries, threshold, n_perm = 1000) {
  col_name <- paste0("near_boundary_", threshold/1000, "kb")

  # Observed proportion near boundaries
  observed_prop <- mean(distance_df[[col_name]], na.rm = TRUE)

  # Get chromosome ranges from loop anchors for generating random positions
  chrom_ranges <- loops_gr$all %>%
    as.data.frame() %>%
    group_by(seqnames) %>%
    summarize(
      min_pos = min(start),
      max_pos = max(end),
      .groups = "drop"
    ) %>%
    rename(chr = seqnames)

  # Count boundaries per chromosome to maintain the same distribution
  boundaries_per_chr <- boundaries %>%
    group_by(chr) %>%
    summarize(n_boundaries = n(), .groups = "drop")

  # Generate null distribution by placing boundaries at random positions

  null_props <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    # Generate random boundary positions within each chromosome's range
    random_boundaries <- boundaries_per_chr %>%
      inner_join(chrom_ranges, by = "chr") %>%
      rowwise() %>%
      mutate(
        random_positions = list(sample(min_pos:max_pos, n_boundaries, replace = FALSE))
      ) %>%
      ungroup() %>%
      unnest(random_positions) %>%
      transmute(
        chr = chr,
        Boundary = random_positions,
        Type = "Random",
        Enriched_In = "Random",
        Gap_Score = 0
      )

    random_gr <- boundaries_to_granges(random_boundaries, buffer = 5000)

    # Compute distances with random boundaries
    dist1 <- distanceToNearest(loops_gr$anchor1, random_gr, ignore.strand = TRUE)
    dist2 <- distanceToNearest(loops_gr$anchor2, random_gr, ignore.strand = TRUE)

    anchor1_dist <- rep(NA_real_, length(loops_gr$anchor1))
    anchor2_dist <- rep(NA_real_, length(loops_gr$anchor2))

    if (length(dist1) > 0) anchor1_dist[queryHits(dist1)] <- mcols(dist1)$distance
    if (length(dist2) > 0) anchor2_dist[queryHits(dist2)] <- mcols(dist2)$distance

    min_dist <- pmin(anchor1_dist, anchor2_dist, na.rm = TRUE)
    null_props[i] <- mean(min_dist <= threshold, na.rm = TRUE)
  }

  # One-tailed p-value (testing if observed is GREATER than null, i.e., enrichment)
  p_value_onetail <- (sum(null_props >= observed_prop) + 1) / (n_perm + 1)

  # Two-tailed p-value for any deviation from null
  p_value_twotail <- (sum(abs(null_props - mean(null_props)) >= abs(observed_prop - mean(null_props))) + 1) / (n_perm + 1)

  list(
    observed = observed_prop,
    null_mean = mean(null_props),
    null_sd = sd(null_props),
    null_distribution = null_props,
    fold_enrichment = observed_prop / mean(null_props),
    p_value = p_value_onetail,  # Use one-tailed for enrichment
    p_value_twotail = p_value_twotail,
    threshold = threshold
  )
}

#' Run Wilcoxon test comparing distance distributions
#' @param distance_df tibble from compute_boundary_distances
#' @return tibble with test results
run_wilcoxon_test <- function(distance_df) {
  lost <- distance_df %>% filter(direction == "down_in_mutant") %>% pull(min_dist_to_boundary)
  gained <- distance_df %>% filter(direction == "up_in_mutant") %>% pull(min_dist_to_boundary)

  # Remove NAs
  lost <- lost[!is.na(lost)]
  gained <- gained[!is.na(gained)]

  if (length(lost) < 2 || length(gained) < 2) {
    return(tibble(
      test = "Wilcoxon rank-sum",
      lost_median = median(lost),
      gained_median = median(gained),
      lost_n = length(lost),
      gained_n = length(gained),
      statistic = NA_real_,
      p_value = NA_real_,
      effect_size_r = NA_real_
    ))
  }

  wilcox_result <- wilcox.test(lost, gained, exact = FALSE)

  # Effect size r = Z / sqrt(N)
  n_total <- length(lost) + length(gained)
  z_approx <- qnorm(wilcox_result$p.value / 2) * sign(median(lost) - median(gained))
  effect_size_r <- abs(z_approx) / sqrt(n_total)

  tibble(
    test = "Wilcoxon rank-sum",
    lost_median = median(lost),
    gained_median = median(gained),
    lost_n = length(lost),
    gained_n = length(gained),
    statistic = as.numeric(wilcox_result$statistic),
    p_value = wilcox_result$p.value,
    effect_size_r = effect_size_r
  )
}

#' Bootstrap confidence interval for fold enrichment
#' @param distance_df tibble from compute_boundary_distances
#' @param threshold Distance threshold (bp)
#' @param n_boot Number of bootstrap samples
#' @return tibble with bootstrap results
run_bootstrap_ci <- function(distance_df, threshold, n_boot = 1000) {
  col_name <- paste0("near_boundary_", threshold/1000, "kb")

  lost <- distance_df %>% filter(direction == "down_in_mutant")
  gained <- distance_df %>% filter(direction == "up_in_mutant")

  observed_lost_prop <- mean(lost[[col_name]], na.rm = TRUE)
  observed_gained_prop <- mean(gained[[col_name]], na.rm = TRUE)

  boot_ratios <- numeric(n_boot)

  for (i in seq_len(n_boot)) {
    # Bootstrap sample with replacement
    boot_lost <- lost[sample(nrow(lost), replace = TRUE), ]
    boot_gained <- gained[sample(nrow(gained), replace = TRUE), ]

    lost_prop <- mean(boot_lost[[col_name]], na.rm = TRUE)
    gained_prop <- mean(boot_gained[[col_name]], na.rm = TRUE)

    # Ratio of proportions (gained/lost)
    if (lost_prop > 0) {
      boot_ratios[i] <- gained_prop / lost_prop
    } else {
      boot_ratios[i] <- NA
    }
  }

  boot_ratios <- boot_ratios[!is.na(boot_ratios)]

  tibble(
    threshold_bp = threshold,
    observed_ratio = observed_gained_prop / observed_lost_prop,
    boot_mean = mean(boot_ratios),
    boot_ci_low = quantile(boot_ratios, 0.025),
    boot_ci_high = quantile(boot_ratios, 0.975)
  )
}

# ==============================================================================
# 6. DIRECTION CONCORDANCE ANALYSIS
# ==============================================================================

#' Analyze concordance between loop direction and boundary enrichment
#' @param distance_df tibble from compute_boundary_distances
#' @param threshold Distance threshold for "near" classification
#' @return list with contingency table and chi-square test
analyze_direction_concordance <- function(distance_df, threshold = 50000) {
  cat("\n[Analyzing direction concordance]\n")

  col_name <- paste0("near_boundary_", threshold/1000, "kb")

  # Filter to loops near boundaries
  near_loops <- distance_df %>%
    filter(.data[[col_name]] == TRUE & !is.na(nearest_boundary_enriched))

  if (nrow(near_loops) < 10) {
    cat("  WARNING: Too few loops near boundaries for concordance analysis\n")
    return(list(
      contingency = NULL,
      test = NULL,
      summary = tibble(
        analysis = "direction_concordance",
        status = "insufficient_data",
        n_near_loops = nrow(near_loops)
      )
    ))
  }

  # Create contingency table: loop_direction × boundary_enriched_in
  contingency <- table(
    loop_direction = near_loops$direction,
    boundary_enriched = near_loops$nearest_boundary_enriched
  )

  cat("  Contingency table (loops within", threshold/1000, "kb of boundary):\n")
  print(contingency)

  # Chi-square test
  if (all(dim(contingency) >= 2)) {
    test <- chisq.test(contingency, simulate.p.value = TRUE, B = 2000)
    cat(sprintf("  Chi-square test: X² = %.2f, p = %.4f\n", test$statistic, test$p.value))
  } else {
    test <- NULL
    cat("  Chi-square test: Not enough categories\n")
  }

  # Calculate concordance rate
  # "Concordant" = gained loop near mutant-enriched boundary OR lost loop near control-enriched boundary
  concordant <- sum(
    near_loops$direction == "up_in_mutant" & near_loops$nearest_boundary_enriched == "Matrix 2",
    near_loops$direction == "down_in_mutant" & near_loops$nearest_boundary_enriched == "Matrix 1",
    na.rm = TRUE
  )
  discordant <- nrow(near_loops) - concordant
  concordance_rate <- concordant / nrow(near_loops)

  cat(sprintf("  Concordance rate: %.1f%% (%d concordant, %d discordant)\n",
              concordance_rate * 100, concordant, discordant))

  # Summary tibble
  summary <- tibble(
    threshold_kb = threshold / 1000,
    n_near_loops = nrow(near_loops),
    n_concordant = concordant,
    n_discordant = discordant,
    concordance_rate = concordance_rate,
    chi_sq_statistic = if (!is.null(test)) test$statistic else NA_real_,
    chi_sq_p_value = if (!is.null(test)) test$p.value else NA_real_
  )

  return(list(
    contingency = contingency,
    test = test,
    summary = summary,
    near_loops = near_loops
  ))
}

# ==============================================================================
# 7. BOUNDARY TYPE ANALYSIS
# ==============================================================================

#' Analyze association between boundary type and loop direction
#' @param distance_df tibble from compute_boundary_distances
#' @param threshold Distance threshold
#' @return tibble with per-type enrichment statistics
analyze_boundary_types <- function(distance_df, threshold = 50000) {
  cat("\n[Analyzing boundary type associations]\n")

  col_name <- paste0("near_boundary_", threshold/1000, "kb")

  # Filter to loops near boundaries
  near_loops <- distance_df %>%
    filter(.data[[col_name]] == TRUE & !is.na(nearest_boundary_type))

  if (nrow(near_loops) < 10) {
    cat("  WARNING: Too few loops near boundaries for type analysis\n")
    return(tibble())
  }

  # Per-type Fisher's exact test
  results <- list()
  boundary_types <- unique(near_loops$nearest_boundary_type)

  for (btype in boundary_types) {
    # Contingency: (this_type vs other_types) × (gained vs lost)
    is_this_type <- near_loops$nearest_boundary_type == btype
    is_gained <- near_loops$direction == "up_in_mutant"

    contingency <- matrix(c(
      sum(is_this_type & !is_gained),   # this_type, lost
      sum(!is_this_type & !is_gained),  # other_type, lost
      sum(is_this_type & is_gained),    # this_type, gained
      sum(!is_this_type & is_gained)    # other_type, gained
    ), nrow = 2, byrow = TRUE,
    dimnames = list(c("lost", "gained"), c(btype, "other")))

    test <- fisher.test(contingency)

    results[[btype]] <- tibble(
      boundary_type = btype,
      n_lost_near = sum(is_this_type & !is_gained),
      n_gained_near = sum(is_this_type & is_gained),
      total_near = sum(is_this_type),
      prop_gained = sum(is_this_type & is_gained) / sum(is_this_type),
      odds_ratio = test$estimate,
      odds_ratio_ci_low = test$conf.int[1],
      odds_ratio_ci_high = test$conf.int[2],
      p_value = test$p.value
    )
  }

  result_df <- bind_rows(results) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)

  cat("  Boundary type associations (within", threshold/1000, "kb):\n")
  for (i in seq_len(nrow(result_df))) {
    row <- result_df[i,]
    sig <- if (row$p_adj < 0.05) "*" else ""
    cat(sprintf("    %s: OR=%.2f (%.2f-%.2f), p=%.4f%s\n",
                row$boundary_type, row$odds_ratio,
                row$odds_ratio_ci_low, row$odds_ratio_ci_high,
                row$p_value, sig))
  }

  return(result_df)
}

# ==============================================================================
# 8. VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create distance distribution violin plot
#' @param distance_df tibble from compute_boundary_distances
#' @param label Timepoint label
#' @return ggplot object
create_distance_violin <- function(distance_df, label) {
  # Convert to log10 for visualization (add pseudocount to avoid log(0))
  plot_df <- distance_df %>%
    filter(!is.na(min_dist_to_boundary)) %>%
    mutate(
      log_distance = log10(min_dist_to_boundary + 1),
      direction_label = case_when(
        direction == "down_in_mutant" ~ "Lost in mutant",
        direction == "up_in_mutant" ~ "Gained in mutant",
        TRUE ~ direction
      )
    )

  # Wilcoxon test for subtitle
  wilcox_result <- run_wilcoxon_test(distance_df)

  p <- ggplot(plot_df, aes(x = direction_label, y = log_distance, fill = direction)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = DIRECTION_COLORS, guide = "none") +
    scale_y_continuous(
      name = "Distance to nearest differential boundary",
      breaks = c(3, 4, 5, 6, 7),
      labels = c("1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb")
    ) +
    labs(
      title = "Loop Anchor Distance to Differential TAD Boundaries",
      subtitle = sprintf("%s | Wilcoxon p = %.3g | Lost median = %s, Gained median = %s",
                        label,
                        wilcox_result$p_value,
                        format(wilcox_result$lost_median, big.mark = ","),
                        format(wilcox_result$gained_median, big.mark = ",")),
      x = "Loop Direction"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Create enrichment bar plot
#' @param fisher_results tibble from run_fisher_test (multiple thresholds)
#' @param label Timepoint label
#' @return ggplot object
create_enrichment_barplot <- function(fisher_results, label) {
  plot_df <- fisher_results %>%
    select(threshold_label, lost_prop_near, gained_prop_near) %>%
    pivot_longer(cols = c(lost_prop_near, gained_prop_near),
                 names_to = "direction",
                 values_to = "proportion") %>%
    mutate(
      direction_label = case_when(
        direction == "lost_prop_near" ~ "Lost",
        direction == "gained_prop_near" ~ "Gained"
      ),
      threshold_label = factor(threshold_label, levels = paste0(DISTANCE_THRESHOLDS/1000, "kb"))
    )

  p <- ggplot(plot_df, aes(x = threshold_label, y = proportion * 100, fill = direction_label)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = c("Lost" = DIRECTION_COLORS["down_in_mutant"],
                                 "Gained" = DIRECTION_COLORS["up_in_mutant"]),
                     name = "Loop Direction") +
    labs(
      title = "Proportion of Loops Near Differential Boundaries",
      subtitle = label,
      x = "Distance Threshold",
      y = "% Loops with Anchor Near Boundary"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "top"
    )

  # Add p-value annotations
  for (i in seq_len(nrow(fisher_results))) {
    p_val <- fisher_results$p_value[i]
    sig_label <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"
    p <- p + annotate("text", x = i, y = max(plot_df$proportion) * 105,
                      label = sig_label, size = 4)
  }

  return(p)
}

#' Create boundary type × direction heatmap
#' @param type_results tibble from analyze_boundary_types
#' @param label Timepoint label
#' @return ggplot object
create_type_heatmap <- function(type_results, label) {
  if (nrow(type_results) == 0) {
    return(ggplot() +
           annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 6) +
           theme_void())
  }

  # Create matrix for heatmap
  plot_df <- type_results %>%
    mutate(
      log2_or = log2(odds_ratio),
      sig_label = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  p <- ggplot(plot_df, aes(x = "Gained vs Lost", y = boundary_type, fill = log2_or)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f%s", odds_ratio, sig_label)),
              color = "black", size = 4, fontface = "bold") +
    scale_fill_gradient2(
      low = DIRECTION_COLORS["down_in_mutant"],
      mid = "white",
      high = DIRECTION_COLORS["up_in_mutant"],
      midpoint = 0,
      name = "log2(OR)",
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    labs(
      title = "Boundary Type Association with Loop Direction",
      subtitle = sprintf("%s | OR > 1 = enriched in gained loops", label),
      x = NULL,
      y = "Boundary Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  return(p)
}

#' Create permutation test visualization
#' @param perm_results list of permutation results from run_permutation_test
#' @param label Timepoint label
#' @return ggplot object (faceted by threshold)
create_permutation_plot <- function(perm_results, label) {
  if (length(perm_results) == 0) {
    return(ggplot() +
           annotate("text", x = 0.5, y = 0.5, label = "No permutation data", size = 6) +
           theme_void())
  }

  # Combine null distributions from all thresholds
  plot_data <- map_dfr(names(perm_results), function(thresh) {
    res <- perm_results[[thresh]]
    tibble(
      threshold = paste0(as.numeric(thresh)/1000, "kb"),
      null_value = res$null_distribution,
      observed = res$observed,
      null_mean = res$null_mean,
      fold_enrichment = res$fold_enrichment,
      p_value = res$p_value
    )
  })

  # Summary for annotations
  summary_data <- plot_data %>%
    group_by(threshold) %>%
    summarize(
      observed = first(observed),
      null_mean = first(null_mean),
      fold_enrichment = first(fold_enrichment),
      p_value = first(p_value),
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf("Obs: %.1f%%\nNull: %.1f%%\nFold: %.2fx\np < 0.001",
                      observed * 100, null_mean * 100, fold_enrichment)
    )

  p <- ggplot(plot_data, aes(x = null_value * 100)) +
    geom_histogram(bins = 30, fill = "gray70", color = "gray40", alpha = 0.7) +
    geom_vline(data = summary_data, aes(xintercept = observed * 100),
               color = "#d73027", linewidth = 1.2, linetype = "solid") +
    geom_vline(data = summary_data, aes(xintercept = null_mean * 100),
               color = "black", linewidth = 0.8, linetype = "dashed") +
    geom_label(data = summary_data,
               aes(x = observed * 100, y = Inf, label = label),
               hjust = -0.1, vjust = 1.2, size = 3, fill = "white", alpha = 0.9,
               label.size = 0.3) +
    facet_wrap(~threshold, scales = "free_x", ncol = 2) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "Permutation Test: Loop-Boundary Proximity Enrichment",
      subtitle = sprintf("%s | Red line = observed, dashed = null mean", label),
      x = "Proportion of loops near differential boundaries",
      y = "Permutation count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    )

  return(p)
}

#' Create concordance mosaic/heatmap
#' @param concordance_result list from analyze_direction_concordance
#' @param label Timepoint label
#' @return ggplot object
create_concordance_plot <- function(concordance_result, label) {
  if (is.null(concordance_result$contingency)) {
    return(ggplot() +
           annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 6) +
           theme_void())
  }

  # Convert contingency table to tibble for ggplot
  contingency_df <- as.data.frame(concordance_result$contingency) %>%
    as_tibble() %>%
    mutate(
      loop_direction_label = case_when(
        loop_direction == "down_in_mutant" ~ "Lost",
        loop_direction == "up_in_mutant" ~ "Gained"
      ),
      boundary_enriched_label = case_when(
        boundary_enriched == "Matrix 1" ~ "Control-enriched",
        boundary_enriched == "Matrix 2" ~ "Mutant-enriched"
      ),
      concordant = (loop_direction == "down_in_mutant" & boundary_enriched == "Matrix 1") |
                   (loop_direction == "up_in_mutant" & boundary_enriched == "Matrix 2")
    )

  summary <- concordance_result$summary

  p <- ggplot(contingency_df, aes(x = loop_direction_label, y = boundary_enriched_label,
                                   fill = Freq, alpha = concordant)) +
    geom_tile(color = "white", linewidth = 2) +
    geom_text(aes(label = Freq), color = "black", size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#756bb1", name = "Count") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    labs(
      title = "Loop Direction × Boundary Enrichment Concordance",
      subtitle = sprintf("%s | Concordance rate: %.1f%% (χ² p = %.3g)",
                        label, summary$concordance_rate * 100, summary$chi_sq_p_value),
      x = "Loop Direction",
      y = "Nearest Boundary Enriched In"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )

  return(p)
}

# ==============================================================================
# 9. REPORT GENERATION
# ==============================================================================

#' Generate analysis report
#' @param distance_df tibble with distance metrics
#' @param fisher_results tibble from Fisher's tests
#' @param wilcox_result tibble from Wilcoxon test
#' @param perm_results_df tibble from permutation tests
#' @param concordance_result list from concordance analysis
#' @param type_results tibble from boundary type analysis
#' @param timepoint_label Label for timepoint
#' @param output_dir Output directory
generate_report <- function(distance_df, fisher_results, wilcox_result,
                           perm_results_df, concordance_result, type_results,
                           timepoint_label, output_dir) {

  # Original: file.path(output_dir, "analysis_report.txt")
  report_path <- file.path(output_dir, "1F_analysis_report.txt")

  sink(report_path)

  cat("================================================================================\n")
  cat("BOUNDARY-LOOP CROSS-REFERENCE ANALYSIS REPORT\n")
  cat(sprintf("Timepoint: %s\n", timepoint_label))
  cat(sprintf("Generated: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat("================================================================================\n\n")

  # Summary statistics
  cat("1. DATA SUMMARY\n")
  cat("----------------\n")
  cat(sprintf("Total differential loops: %d\n", nrow(distance_df)))
  cat(sprintf("  Lost in mutant: %d\n", sum(distance_df$direction == "down_in_mutant")))
  cat(sprintf("  Gained in mutant: %d\n", sum(distance_df$direction == "up_in_mutant")))
  cat(sprintf("Median distance to nearest boundary: %s bp\n\n",
              format(median(distance_df$min_dist_to_boundary, na.rm = TRUE), big.mark = ",")))

  # Distance comparison
  cat("2. DISTANCE DISTRIBUTION (Wilcoxon Test)\n")
  cat("-----------------------------------------\n")
  cat(sprintf("Lost loops - median distance: %s bp (n=%d)\n",
              format(wilcox_result$lost_median, big.mark = ","), wilcox_result$lost_n))
  cat(sprintf("Gained loops - median distance: %s bp (n=%d)\n",
              format(wilcox_result$gained_median, big.mark = ","), wilcox_result$gained_n))
  cat(sprintf("Wilcoxon p-value: %.4g\n", wilcox_result$p_value))
  cat(sprintf("Effect size (r): %.3f\n\n", wilcox_result$effect_size_r))

  # Proximity enrichment
  cat("3. PROXIMITY ENRICHMENT (Fisher's Exact Tests)\n")
  cat("-----------------------------------------------\n")
  for (i in seq_len(nrow(fisher_results))) {
    row <- fisher_results[i,]
    sig <- if (row$p_value < 0.05) "*" else ""
    cat(sprintf("Threshold %s:\n", row$threshold_label))
    cat(sprintf("  Lost near: %.1f%% (%d/%d)\n",
                row$lost_prop_near * 100, row$lost_near, row$lost_near + row$lost_far))
    cat(sprintf("  Gained near: %.1f%% (%d/%d)\n",
                row$gained_prop_near * 100, row$gained_near, row$gained_near + row$gained_far))
    cat(sprintf("  OR = %.2f (95%% CI: %.2f-%.2f), p = %.4g%s\n\n",
                row$odds_ratio, row$odds_ratio_ci_low, row$odds_ratio_ci_high, row$p_value, sig))
  }

  # Permutation test results
  cat("4. PERMUTATION TEST (Enrichment vs Null)\n")
  cat("-----------------------------------------\n")
  if (!is.null(perm_results_df) && nrow(perm_results_df) > 0) {
    cat(sprintf("Number of permutations: %d\n\n", N_PERMUTATIONS))
    for (i in seq_len(nrow(perm_results_df))) {
      row <- perm_results_df[i,]
      sig <- if (row$perm_p_value < 0.05) "*" else ""
      cat(sprintf("Threshold %s:\n", row$threshold_label))
      cat(sprintf("  Observed proportion near boundary: %.3f\n", row$observed_prop))
      cat(sprintf("  Null expectation (mean ± SD): %.3f ± %.3f\n", row$null_mean, row$null_sd))
      cat(sprintf("  Fold enrichment over null: %.2f\n", row$fold_enrichment))
      cat(sprintf("  Permutation p-value: %.4g%s\n\n", row$perm_p_value, sig))
    }
    cat("Interpretation:\n")
    cat("  - Fold enrichment > 1 indicates loops are closer to differential boundaries\n")
    cat("    than expected by chance\n")
    cat("  - Permutation p-value tests if observed proximity exceeds null distribution\n\n")
  } else {
    cat("Permutation test not run or insufficient data\n\n")
  }

  # Direction concordance
  cat("5. DIRECTION CONCORDANCE\n")
  cat("------------------------\n")
  if (!is.null(concordance_result$contingency)) {
    cat("Contingency table:\n")
    print(concordance_result$contingency)
    cat(sprintf("\nConcordance rate: %.1f%%\n", concordance_result$summary$concordance_rate * 100))
    cat(sprintf("Chi-square p-value: %.4g\n\n", concordance_result$summary$chi_sq_p_value))
    cat("Interpretation:\n")
    cat("  - 'Concordant' = gained loop near mutant-enriched boundary OR\n")
    cat("                   lost loop near control-enriched boundary\n")
    if (concordance_result$summary$concordance_rate > 0.5) {
      cat("  - Result: Concordant direction (boundary and loop changes align)\n\n")
    } else {
      cat("  - Result: Discordant direction (boundary and loop changes do not align)\n\n")
    }
  } else {
    cat("Insufficient data for concordance analysis\n\n")
  }

  # Boundary type associations
  cat("6. BOUNDARY TYPE ASSOCIATIONS\n")
  cat("-----------------------------\n")
  if (nrow(type_results) > 0) {
    for (i in seq_len(nrow(type_results))) {
      row <- type_results[i,]
      sig <- if (row$p_adj < 0.05) "*" else ""
      interpretation <- if (row$odds_ratio > 1) "enriched in gained" else "enriched in lost"
      cat(sprintf("%s: OR=%.2f (%.2f-%.2f), p=%.4g, p_adj=%.4g%s [%s]\n",
                  row$boundary_type, row$odds_ratio,
                  row$odds_ratio_ci_low, row$odds_ratio_ci_high,
                  row$p_value, row$p_adj, sig, interpretation))
    }
    cat("\nBiological interpretation:\n")
    cat("  - Split boundaries: TAD subdivided → may disrupt existing contacts\n")
    cat("  - Merge boundaries: TADs fused → may create new contacts\n")
    cat("  - Shifted boundaries: Position moved → loop redistribution\n")
    cat("  - Strength Change: Boundary weakened/strengthened → subtle modulation\n")
  } else {
    cat("Insufficient data for boundary type analysis\n")
  }

  cat("\n================================================================================\n")
  cat("END OF REPORT\n")
  cat("================================================================================\n")

  sink()

  cat(sprintf("  Report saved: %s\n", report_path))
}

# ==============================================================================
# 10. MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Run complete boundary-loop cross-reference analysis for one timepoint
#' @param timepoint "early" or "late"
run_analysis <- function(timepoint) {
  cat(sprintf("\n########## ANALYZING: %s ##########\n", toupper(timepoint)))

  config <- TIMEPOINT_CONFIG[[timepoint]]

  # Create output directories
  # Original: output_dir <- config$output_dir; plots_dir <- file.path(output_dir, "plots")
  tsv_dir <- config$tsv_dir
  plots_dir <- config$plots_dir
  dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # -------------------------------------------------------------------------
  # Step 1: Load data
  # -------------------------------------------------------------------------
  cat("\n[Step 1] Loading data...\n")
  loops <- load_loops(config$loops_file)
  boundaries <- load_boundaries(config$boundaries_file)

  # Convert to GRanges
  loops_gr <- loops_to_granges(loops)
  boundaries_gr <- boundaries_to_granges(boundaries, buffer = 5000)

  cat(sprintf("  Loops GRanges: %d anchor pairs\n", length(loops_gr$anchor1)))
  cat(sprintf("  Boundaries GRanges: %d regions\n", length(boundaries_gr)))

  # -------------------------------------------------------------------------
  # Step 2: Compute boundary distances
  # -------------------------------------------------------------------------
  cat("\n[Step 2] Computing boundary distances...\n")
  distance_df <- compute_boundary_distances(loops_gr, boundaries_gr)

  # Save distance summary
  # Original: write_tsv(distance_df, file.path(output_dir, "boundary_loop_overlap_summary.tsv"))
  write_tsv(distance_df, file.path(tsv_dir, "1F_boundary_loop_overlap_summary.tsv"))
  cat(sprintf("  Saved: 1F_boundary_loop_overlap_summary.tsv\n"))

  # -------------------------------------------------------------------------
  # Step 3: Statistical tests
  # -------------------------------------------------------------------------
  cat("\n[Step 3] Running statistical tests...\n")

  # Fisher's exact tests at multiple thresholds
  fisher_results <- map_dfr(DISTANCE_THRESHOLDS, ~run_fisher_test(distance_df, .x))
  cat("  Fisher's exact tests complete\n")

  # Wilcoxon test
  wilcox_result <- run_wilcoxon_test(distance_df)
  cat("  Wilcoxon test complete\n")

  # Bootstrap CIs (lighter version)
  bootstrap_results <- map_dfr(DISTANCE_THRESHOLDS[1:2],
                               ~run_bootstrap_ci(distance_df, .x, n_boot = 500))
  cat("  Bootstrap CIs complete\n")

  # Permutation tests for enrichment significance
  cat("  Running permutation tests (n=", N_PERMUTATIONS, ")...\n", sep = "")
  perm_results <- list()
  for (thresh in DISTANCE_THRESHOLDS[1:2]) {  # Run on first two thresholds (10kb, 25kb)
    cat("    Threshold:", thresh/1000, "kb\n")
    perm_results[[as.character(thresh)]] <- run_permutation_test(
      distance_df, loops_gr, boundaries, thresh, n_perm = N_PERMUTATIONS
    )
  }
  cat("  Permutation tests complete\n")

  # Convert permutation results to tibble for saving
  perm_results_df <- map_dfr(names(perm_results), function(thresh) {
    res <- perm_results[[thresh]]
    tibble(
      threshold_bp = as.numeric(thresh),
      threshold_label = paste0(as.numeric(thresh)/1000, "kb"),
      observed_prop = res$observed,
      null_mean = res$null_mean,
      null_sd = res$null_sd,
      fold_enrichment = res$fold_enrichment,
      perm_p_value = res$p_value
    )
  })

  # Original: write_tsv(perm_results_df, file.path(output_dir, "permutation_test_results.tsv"))
  write_tsv(perm_results_df, file.path(tsv_dir, "1F_permutation_test_results.tsv"))
  cat(sprintf("  Saved: 1F_permutation_test_results.tsv\n"))

  # Combine enrichment statistics
  enrichment_stats <- fisher_results %>%
    left_join(bootstrap_results, by = c("threshold_bp")) %>%
    left_join(perm_results_df %>% select(threshold_bp, fold_enrichment, perm_p_value),
              by = "threshold_bp") %>%
    bind_cols(
      wilcox_result %>%
        select(wilcox_p = p_value, wilcox_effect_r = effect_size_r) %>%
        slice(rep(1, nrow(fisher_results)))
    )

  # Original: write_tsv(enrichment_stats, file.path(output_dir, "enrichment_statistics.tsv"))
  write_tsv(enrichment_stats, file.path(tsv_dir, "1F_enrichment_statistics.tsv"))
  cat(sprintf("  Saved: 1F_enrichment_statistics.tsv\n"))

  # -------------------------------------------------------------------------
  # Step 4: Direction concordance
  # -------------------------------------------------------------------------
  cat("\n[Step 4] Analyzing direction concordance...\n")
  concordance_result <- analyze_direction_concordance(distance_df, threshold = 50000)

  # Original: write_tsv(concordance_result$summary, file.path(output_dir, "direction_concordance.tsv"))
  write_tsv(concordance_result$summary, file.path(tsv_dir, "1F_direction_concordance.tsv"))
  cat(sprintf("  Saved: 1F_direction_concordance.tsv\n"))

  # -------------------------------------------------------------------------
  # Step 5: Boundary type associations
  # -------------------------------------------------------------------------
  cat("\n[Step 5] Analyzing boundary type associations...\n")
  type_results <- analyze_boundary_types(distance_df, threshold = 50000)

  if (nrow(type_results) > 0) {
    # Original: write_tsv(type_results, file.path(output_dir, "boundary_type_association.tsv"))
    write_tsv(type_results, file.path(tsv_dir, "1F_boundary_type_association.tsv"))
    cat(sprintf("  Saved: 1F_boundary_type_association.tsv\n"))
  }

  # -------------------------------------------------------------------------
  # Step 6: Generate visualizations
  # -------------------------------------------------------------------------
  cat("\n[Step 6] Generating visualizations...\n")

  # 6a. Distance violin plot
  p_violin <- create_distance_violin(distance_df, config$label)
  # Original: file.path(plots_dir, "distance_violin")
  save_multiformat_ggplot(p_violin, file.path(plots_dir, "1F_distance_violin"),
                         width = 8, height = 6)

  # 6b. Enrichment bar plot
  p_barplot <- create_enrichment_barplot(fisher_results, config$label)
  # Original: file.path(plots_dir, "enrichment_barplot")
  save_multiformat_ggplot(p_barplot, file.path(plots_dir, "1F_enrichment_barplot"),
                         width = 10, height = 6)

  # 6c. Boundary type heatmap
  p_heatmap <- create_type_heatmap(type_results, config$label)
  # Original: file.path(plots_dir, "type_heatmap")
  save_multiformat_ggplot(p_heatmap, file.path(plots_dir, "1F_type_heatmap"),
                         width = 6, height = 6)

  # 6d. Concordance plot
  p_concordance <- create_concordance_plot(concordance_result, config$label)
  # Original: file.path(plots_dir, "concordance_heatmap")
  save_multiformat_ggplot(p_concordance, file.path(plots_dir, "1F_concordance_heatmap"),
                         width = 8, height = 6)

  # 6e. Permutation test plot
  p_permutation <- create_permutation_plot(perm_results, config$label)
  # Original: file.path(plots_dir, "permutation_test")
  save_multiformat_ggplot(p_permutation, file.path(plots_dir, "1F_permutation_test"),
                         width = 10, height = 5)

  # 6f. Combined summary figure
  p_combined <- (p_violin | p_barplot) / (p_concordance | p_heatmap) +
    plot_annotation(
      title = sprintf("Boundary-Loop Cross-Reference: %s", config$label),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  # Original: file.path(plots_dir, "combined_summary")
  save_multiformat_ggplot(p_combined, file.path(plots_dir, "1F_combined_summary"),
                         width = 14, height = 12)

  # -------------------------------------------------------------------------
  # Step 7: Generate report
  # -------------------------------------------------------------------------
  cat("\n[Step 7] Generating analysis report...\n")
  generate_report(
    distance_df = distance_df,
    fisher_results = fisher_results,
    wilcox_result = wilcox_result,
    perm_results_df = perm_results_df,
    concordance_result = concordance_result,
    type_results = type_results,
    timepoint_label = config$label,
    # Original: output_dir = output_dir
    output_dir = tsv_dir
  )

  cat(sprintf("\n✅ Analysis complete for %s\n", timepoint))
  cat(sprintf("   TSV directory: %s\n", tsv_dir))
  cat(sprintf("   Plots directory: %s\n", plots_dir))

  return(list(
    distance_df = distance_df,
    fisher_results = fisher_results,
    wilcox_result = wilcox_result,
    perm_results = perm_results_df,
    concordance_result = concordance_result,
    type_results = type_results
  ))
}

# ==============================================================================
# 11. COMMAND LINE INTERFACE
# ==============================================================================

#' Parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Default timepoint
  timepoint <- "late"

  # Parse --timepoint flag
  if (length(args) > 0) {
    for (i in seq_along(args)) {
      if (args[i] == "--timepoint" && i < length(args)) {
        timepoint <- args[i + 1]
      } else if (grepl("^--timepoint=", args[i])) {
        timepoint <- sub("^--timepoint=", "", args[i])
      } else if (args[i] %in% c("early", "late", "both")) {
        timepoint <- args[i]
      }
    }
  }

  # Validate timepoint
  if (!timepoint %in% c("early", "late", "both")) {
    stop("Invalid timepoint. Use: early, late, or both")
  }

  return(list(timepoint = timepoint))
}

# ==============================================================================
# 12. MAIN EXECUTION
# ==============================================================================

main <- function() {
  args <- parse_args()

  cat("=== Configuration ===\n")
  cat(sprintf("Timepoint: %s\n", args$timepoint))
  cat(sprintf("Distance thresholds: %s\n", paste(DISTANCE_THRESHOLDS/1000, "kb", collapse = ", ")))
  cat(sprintf("Bootstrap samples: %d\n", N_BOOTSTRAP))

  results <- list()

  if (args$timepoint == "both") {
    for (tp in c("early", "late")) {
      results[[tp]] <- run_analysis(tp)
    }
  } else {
    results[[args$timepoint]] <- run_analysis(args$timepoint)
  }

  cat("\n================================================================================\n")
  cat("ALL ANALYSES COMPLETE\n")
  cat(sprintf("End time: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat("================================================================================\n")

  invisible(results)
}

# Run if executed directly
if (!interactive()) {
  main()
}
