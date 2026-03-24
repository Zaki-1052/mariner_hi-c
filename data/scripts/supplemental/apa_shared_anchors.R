# scripts/apa_shared_anchors.R
#
# Aggregate Peak Analysis (APA) for Shared Anchor Loop Subsets
#
# Performs APA on loops from shared anchor analysis:
#   1. Lost long-range loops at shared anchors (141 loops, >1.2 Mb)
#   2. Gained short-range loops at shared anchors (157 loops, <370 kb)
#
# Expected results:
#   - Lost loops: stronger enrichment in CONTROL (loops are lost in mutant)
#   - Gained loops: stronger enrichment in MUTANT (loops are gained)
#
# @usage Rscript scripts/apa_shared_anchors.R [--timepoint late|early]
# @param timepoint Optional: "late" (default) or "early"
#
# @author Zakir Alibhai
# @date 2026-01-28

# ============================================================================
# Setup and Configuration
# ============================================================================

suppressPackageStartupMessages({
  library(mariner)
  library(InteractionSet)
  library(strawr)
  library(HDF5Array)
  library(SummarizedExperiment)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(yaml)
  library(scales)
})

# Load multi-format output utility for PDF + SVG + JPEG output
source("scripts/utils/multi_format_output.R")

# Load paths configuration
cat("\n========================================\n")
cat("APA Analysis: Shared Anchor Loop Subsets\n")
cat("========================================\n\n")

config <- yaml::read_yaml("config/paths_config.yaml")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
TIMEPOINT <- "late"  # Default

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      TIMEPOINT <- args[i + 1]
    }
  }
}

cat(sprintf("Timepoint: %s\n", TIMEPOINT))

# APA parameters
RESOLUTION <- 10000  # 10kb - matches input loop binning
BUFFER_BP <- 50000   # 50kb window around each loop
BUFFER_BINS <- floor(BUFFER_BP / RESOLUTION)  # 5 bins on each side

cat(sprintf("Resolution: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
cat(sprintf("Window size: %d bins (%d kb total)\n",
            2 * BUFFER_BINS + 1, (2 * BUFFER_BINS + 1) * RESOLUTION / 1000))

# .hic file paths from config
HIC_FILES <- c(
  ctrl_M1 = config$hic_files$ctrl_M1,
  ctrl_M2 = config$hic_files$ctrl_M2,
  ctrl_M3 = config$hic_files$ctrl_M3,
  mut_M1 = config$hic_files$mut_M1,
  mut_M2 = config$hic_files$mut_M2,
  mut_M3 = config$hic_files$mut_M3
)

# Sample groups for statistical comparison
GROUPS <- factor(c("ctrl", "ctrl", "ctrl", "mut", "mut", "mut"))
CTRL_INDICES <- which(GROUPS == "ctrl")
MUT_INDICES <- which(GROUPS == "mut")

# Input files
INPUT_DIR <- file.path("output/shared_anchor_analysis", TIMEPOINT, "apa_subsets")
INPUT_FILES <- list(
  lost_longrange = file.path(INPUT_DIR, "shared_lost_longrange.rds"),
  gained_shortrange = file.path(INPUT_DIR, "shared_gained_shortrange.rds")
)

# Output directory
OUTPUT_BASE <- "outputs/apa_shared_anchors"
OUTPUT_DIR <- file.path(OUTPUT_BASE, TIMEPOINT)
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

# ============================================================================
# Core Functions (adapted from apa_analysis.R)
# ============================================================================

#' Check .hic file availability and normalization methods
check_hic_files <- function(hic_files) {
  files_exist <- all(sapply(hic_files, file.exists))

  if (!files_exist) {
    missing <- names(hic_files)[!sapply(hic_files, file.exists)]
    cat("  Warning: Missing .hic files:\n")
    for (f in missing) {
      cat(sprintf("    - %s: %s\n", f, hic_files[f]))
    }
    return(list(files_exist = FALSE, available_norm = NULL))
  }

  # Check normalization methods
  available_norms <- tryCatch({
    strawr::readHicNormTypes(hic_files[1])
  }, error = function(e) {
    c("NONE")
  })

  # Prefer KR, fallback to VC, then NONE
  norm_choice <- if ("KR" %in% available_norms) {
    "KR"
  } else if ("VC" %in% available_norms) {
    "VC"
  } else {
    "NONE"
  }

  cat(sprintf("  All .hic files found\n"))
  cat(sprintf("  Using normalization: %s\n", norm_choice))

  return(list(files_exist = TRUE, available_norm = norm_choice))
}

#' Extract Hi-C matrices around loop positions for APA
extract_apa_matrices <- function(loops, hic_files, resolution, buffer, norm = "KR",
                                  output_dir, subset_name) {

  cat(sprintf("  Extracting matrices: %d loops x %d samples\n", length(loops), length(hic_files)))
  cat(sprintf("  Buffer: %d bins (%dkb window)\n", buffer, (buffer * 2 + 1) * resolution / 1000))

  # Create temporary HDF5 directory
  hdf5_dir <- file.path(output_dir, sprintf("temp_hdf5_%s", subset_name))
  if (dir.exists(hdf5_dir)) {
    unlink(hdf5_dir, recursive = TRUE)
  }
  dir.create(hdf5_dir, recursive = TRUE)

  hdf5_file <- file.path(hdf5_dir, "apa_matrices.h5")

  # Extract matrices using mariner
  tryCatch({
    pixels <- pullHicMatrices(
      x = loops,
      files = hic_files,
      binSize = resolution,
      h5File = hdf5_file,
      half = "both",
      norm = norm,
      matrix = "observed",
      blockSize = 1e6,
      onDisk = TRUE,
      compressionLevel = 1
    )

    cat(sprintf("  Extraction complete: %s\n", hdf5_file))
    return(pixels)

  }, error = function(e) {
    cat(sprintf("  Error during extraction: %s\n", e$message))
    return(NULL)
  })
}

#' Aggregate matrices across loops to create mean contact map
aggregate_apa_matrices <- function(pixels, method = "mean") {

  if (is.null(pixels)) {
    return(NULL)
  }

  tryCatch({
    count_array <- counts(pixels)
    dims <- dim(count_array)
    n_bins <- dims[1]
    n_loops <- dims[3]
    n_samples <- dims[4]

    cat(sprintf("  Aggregating: %d loops x %d samples -> %dx%d matrix per sample\n",
                n_loops, n_samples, n_bins, n_bins))

    # Aggregate per sample
    agg_matrices <- array(NA, dim = c(n_bins, n_bins, n_samples))

    for (j in 1:n_samples) {
      sample_matrices <- as.array(count_array[, , , j, drop = FALSE])

      if (method == "mean") {
        agg_matrices[, , j] <- apply(sample_matrices, c(1, 2), mean, na.rm = TRUE)
      } else if (method == "median") {
        agg_matrices[, , j] <- apply(sample_matrices, c(1, 2), median, na.rm = TRUE)
      }
    }

    cat(sprintf("  Aggregation complete using %s\n", method))
    return(agg_matrices)

  }, error = function(e) {
    cat(sprintf("  Error during aggregation: %s\n", e$message))
    return(NULL)
  })
}

#' Calculate P2LL enrichment scores (center / local background)
calculate_enrichment_scores <- function(pixels, sample_names, groups) {

  if (is.null(pixels)) {
    return(NULL)
  }

  tryCatch({
    count_array <- counts(pixels)
    dims <- dim(count_array)
    n_bins <- dims[1]
    n_loops <- dims[3]
    n_samples <- dims[4]

    cat(sprintf("  Calculating P2LL enrichment for %d loops\n", n_loops))

    # Center pixel position
    center_idx <- ceiling(n_bins / 2)

    # Corner pixels for background (4 corners of matrix)
    corner_indices <- list(
      c(1, 1), c(1, n_bins),
      c(n_bins, 1), c(n_bins, n_bins)
    )

    # Initialize results
    enrichment_df <- data.frame(
      loop_id = rep(1:n_loops, each = n_samples),
      sample = rep(sample_names, times = n_loops),
      group = rep(groups, times = n_loops),
      center_value = NA_real_,
      background_mean = NA_real_,
      enrichment = NA_real_
    )

    row_idx <- 1
    for (i in 1:n_loops) {
      for (j in 1:n_samples) {
        mat <- as.matrix(count_array[, , i, j])

        # Extract center pixel
        center_val <- mat[center_idx, center_idx]

        # Extract corner pixels for background
        corner_vals <- sapply(corner_indices, function(idx) mat[idx[1], idx[2]])
        bg_mean <- mean(corner_vals, na.rm = TRUE)

        # Calculate enrichment (P2LL)
        enrichment <- if (!is.na(center_val) && !is.na(bg_mean) && bg_mean > 0) {
          center_val / bg_mean
        } else {
          NA_real_
        }

        enrichment_df$center_value[row_idx] <- center_val
        enrichment_df$background_mean[row_idx] <- bg_mean
        enrichment_df$enrichment[row_idx] <- enrichment

        row_idx <- row_idx + 1
      }
    }

    valid_enrichment <- enrichment_df$enrichment[!is.na(enrichment_df$enrichment) &
                                                   is.finite(enrichment_df$enrichment)]

    cat(sprintf("  Enrichment calculation complete\n"))
    cat(sprintf("    Mean enrichment: %.2f (range: %.2f - %.2f)\n",
                mean(valid_enrichment, na.rm = TRUE),
                min(valid_enrichment, na.rm = TRUE),
                max(valid_enrichment, na.rm = TRUE)))

    return(enrichment_df)

  }, error = function(e) {
    cat(sprintf("  Error calculating enrichment: %s\n", e$message))
    return(NULL)
  })
}

#' Generate APA heatmap for a single condition
plot_apa_heatmap <- function(agg_matrix, condition_name, sample_indices, resolution,
                              center_label = "Loop") {

  if (is.null(agg_matrix)) {
    return(NULL)
  }

  # Average across specified samples
  mean_matrix <- apply(agg_matrix[, , sample_indices, drop = FALSE], c(1, 2), mean, na.rm = TRUE)

  # Prepare data for ggplot
  n_bins <- nrow(mean_matrix)
  buffer_kb <- (n_bins - 1) / 2 * resolution / 1000

  df <- expand.grid(
    x = 1:n_bins,
    y = 1:n_bins
  )
  df$value <- as.vector(mean_matrix)

  # Create axis labels
  center_idx <- ceiling(n_bins / 2)
  axis_breaks <- c(1, center_idx, n_bins)
  axis_labels <- c(sprintf("-%dkb", buffer_kb), center_label, sprintf("+%dkb", buffer_kb))

  # Plot
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis(name = "Contacts", option = "magma") +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = sprintf("APA: %s", condition_name),
      x = "Anchor 2",
      y = "Anchor 1"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )

  return(p)
}

#' Generate difference heatmap (comparison - reference)
plot_difference_heatmap <- function(agg_matrix, reference_indices, comparison_indices,
                                     resolution, title = "Difference: Mutant - Control") {

  if (is.null(agg_matrix)) {
    return(NULL)
  }

  # Calculate mean matrices
  ref_mean <- apply(agg_matrix[, , reference_indices, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  comp_mean <- apply(agg_matrix[, , comparison_indices, drop = FALSE], c(1, 2), mean, na.rm = TRUE)

  # Difference (comparison - reference)
  diff_matrix <- comp_mean - ref_mean

  # Prepare data
  n_bins <- nrow(diff_matrix)
  buffer_kb <- (n_bins - 1) / 2 * resolution / 1000

  df <- expand.grid(
    x = 1:n_bins,
    y = 1:n_bins
  )
  df$value <- as.vector(diff_matrix)

  # Axis labels
  center_idx <- ceiling(n_bins / 2)
  axis_breaks <- c(1, center_idx, n_bins)
  axis_labels <- c(sprintf("-%dkb", buffer_kb), "Loop", sprintf("+%dkb", buffer_kb))

  # Determine symmetric limits
  max_abs <- max(abs(df$value), na.rm = TRUE)
  limit <- min(max_abs, 5)  # Cap at +/- 5

  # Plot with RdBu diverging scale
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      name = "Difference",
      low = muted("blue"),
      mid = "white",
      high = muted("red"),
      midpoint = 0,
      limits = c(-limit, limit),
      oob = squish
    ) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = title,
      x = "Anchor 2",
      y = "Anchor 1"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )

  return(p)
}

#' Generate box plot comparing enrichment between conditions with statistics
plot_enrichment_comparison <- function(enrichment_df, loop_type, expected_higher = "ctrl") {

  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
    return(NULL)
  }

  # Filter valid enrichment scores
  df_valid <- enrichment_df %>%
    filter(!is.na(enrichment), is.finite(enrichment))

  if (nrow(df_valid) == 0) {
    cat("  Warning: No valid enrichment scores for plotting\n")
    return(NULL)
  }

  # Statistical test
  ctrl_vals <- df_valid$enrichment[df_valid$group == "ctrl"]
  mut_vals <- df_valid$enrichment[df_valid$group == "mut"]

  test_result <- wilcox.test(ctrl_vals, mut_vals, alternative = "two.sided")

  # Effect size (Cliff's Delta)
  cliff_delta <- (sum(outer(mut_vals, ctrl_vals, `-`) > 0) -
                   sum(outer(mut_vals, ctrl_vals, `-`) < 0)) /
                  (length(mut_vals) * length(ctrl_vals))

  # Determine expected outcome
  expected_text <- ifelse(expected_higher == "ctrl",
                          "Expected: Control > Mutant (loops lost in mutant)",
                          "Expected: Mutant > Control (loops gained in mutant)")

  # Actual outcome
  mean_ctrl <- mean(ctrl_vals, na.rm = TRUE)
  mean_mut <- mean(mut_vals, na.rm = TRUE)
  actual_higher <- ifelse(mean_ctrl > mean_mut, "Control", "Mutant")

  # Create plot
  p <- ggplot(df_valid, aes(x = group, y = enrichment, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("ctrl" = "steelblue", "mut" = "firebrick")) +
    scale_x_discrete(labels = c("ctrl" = "Control", "mut" = "BAP1-KO")) +
    labs(
      title = sprintf("P2LL Enrichment: %s", loop_type),
      subtitle = sprintf("Wilcoxon p = %.2e | Cliff's delta = %.3f\n%s higher (mean: ctrl=%.2f, mut=%.2f)",
                        test_result$p.value, cliff_delta, actual_higher, mean_ctrl, mean_mut),
      x = "",
      y = "Pixel Enrichment (P2LL)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )

  return(list(
    plot = p,
    test = test_result,
    cliff_delta = cliff_delta,
    n_ctrl = length(ctrl_vals),
    n_mut = length(mut_vals),
    mean_ctrl = mean_ctrl,
    mean_mut = mean_mut
  ))
}

#' Validate loop resolutions and report statistics
#' Detects mixed-resolution loops that require normalization via assignToBins()
#' @param loops GInteractions object
#' @param name Character label for logging
#' @return TRUE if valid (>0 loops), FALSE otherwise
validate_loop_resolutions <- function(loops, name) {
  if (length(loops) == 0) {
    cat(sprintf("  Warning: No loops in %s\n", name))
    return(FALSE)
  }

  # Check anchor widths to detect mixed resolutions
  widths1 <- width(anchors(loops, type = "first"))
  widths2 <- width(anchors(loops, type = "second"))
  unique_widths <- unique(c(widths1, widths2))

  if (length(unique_widths) > 1) {
    cat(sprintf("  Note: Mixed resolutions detected in %s\n", name))
    cat(sprintf("    Unique anchor widths: %s bp\n",
                paste(sort(unique_widths), collapse = ", ")))
    cat(sprintf("    Will normalize to target resolution (%d bp) via assignToBins()\n", RESOLUTION))
  } else {
    cat(sprintf("  Anchor width: %d bp (single resolution)\n", unique_widths[1]))
  }

  return(TRUE)
}

#' Run complete APA analysis for a loop subset
run_apa_analysis <- function(loops, subset_name, output_dir, hic_files,
                              norm_method, expected_enriched = "ctrl") {

  cat(sprintf("\n--- Processing: %s ---\n", subset_name))

  # Create subset output directory
  subset_output_dir <- file.path(output_dir, subset_name)
  if (!dir.exists(subset_output_dir)) {
    dir.create(subset_output_dir, recursive = TRUE)
  }

  n_loops <- length(loops)
  cat(sprintf("  Loop count: %d\n", n_loops))

  if (n_loops < 10) {
    cat(sprintf("  Skipping: Only %d loops (minimum 10 required)\n", n_loops))
    return(FALSE)
  }

  # Step 1: Normalize loops to target resolution (fixes mixed-resolution input)
  # This is CRITICAL: input loops from shared_anchor_analysis.R may have mixed
  # resolutions (e.g., 10kb and 25kb) because they originate from multi-resolution
  # edgeR analysis. pixelsToMatrices() requires uniform anchor widths.
  cat("\n  Step 1: Normalizing loops to target resolution...\n")
  loops_binned <- assignToBins(
    x = loops,
    binSize = RESOLUTION,
    pos1 = "center",
    pos2 = "center"
  )
  cat(sprintf("    Binned %d loops to %d bp resolution\n", length(loops_binned), RESOLUTION))

  # Validate post-binning loop count
  if (length(loops_binned) < 10) {
    cat(sprintf("  Skipping: Only %d loops after binning (minimum 10 required)\n",
                length(loops_binned)))
    return(FALSE)
  }

  # Step 2: Prepare loops with buffer for matrix extraction
  cat("\n  Step 2: Preparing loops with buffer for matrix extraction...\n")
  loops_buffered <- pixelsToMatrices(
    x = loops_binned,
    buffer = BUFFER_BINS
  )

  # Step 3: Extract matrices
  cat("\n  Step 3: Extracting Hi-C matrices...\n")
  pixels <- extract_apa_matrices(
    loops_buffered,
    hic_files,
    RESOLUTION,
    BUFFER_BINS,
    norm_method,
    subset_output_dir,
    subset_name
  )

  if (is.null(pixels)) {
    cat("  Matrix extraction failed\n")
    return(FALSE)
  }

  # Step 4: Aggregate matrices
  cat("\n  Step 4: Aggregating matrices...\n")
  agg_matrices <- aggregate_apa_matrices(pixels, method = "mean")

  if (is.null(agg_matrices)) {
    cat("  Aggregation failed\n")
    return(FALSE)
  }

  # Step 5: Calculate enrichment scores
  cat("\n  Step 5: Calculating enrichment scores...\n")
  enrichment_df <- calculate_enrichment_scores(pixels, names(hic_files), GROUPS)

  if (!is.null(enrichment_df)) {
    enrichment_file <- file.path(subset_output_dir, "enrichment_scores.tsv")
    write.table(enrichment_df, enrichment_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", enrichment_file))
  }

  # Step 6: Generate visualizations
  cat("\n  Step 6: Generating visualizations...\n")

  # 6a. Control aggregate heatmap
  p_ctrl <- plot_apa_heatmap(agg_matrices, "Control", CTRL_INDICES, RESOLUTION)
  if (!is.null(p_ctrl)) {
    save_multiformat_ggplot(p_ctrl, file.path(subset_output_dir, "aggregate_heatmap_ctrl"),
                            width = 6, height = 5)
  }

  # 6b. Mutant aggregate heatmap
  p_mut <- plot_apa_heatmap(agg_matrices, "BAP1-KO (Mutant)", MUT_INDICES, RESOLUTION)
  if (!is.null(p_mut)) {
    save_multiformat_ggplot(p_mut, file.path(subset_output_dir, "aggregate_heatmap_mut"),
                            width = 6, height = 5)
  }

  # 6c. Difference heatmap
  p_diff <- plot_difference_heatmap(agg_matrices, CTRL_INDICES, MUT_INDICES, RESOLUTION,
                                     title = "Difference: BAP1-KO - Control")
  if (!is.null(p_diff)) {
    save_multiformat_ggplot(p_diff, file.path(subset_output_dir, "difference_heatmap"),
                            width = 6, height = 5)
  }

  # 6d. Enrichment comparison box plot
  loop_type_label <- gsub("_", " ", subset_name)
  loop_type_label <- tools::toTitleCase(loop_type_label)

  enrichment_result <- plot_enrichment_comparison(enrichment_df, loop_type_label, expected_enriched)

  if (!is.null(enrichment_result)) {
    save_multiformat_ggplot(enrichment_result$plot,
                            file.path(subset_output_dir, "enrichment_boxplot"),
                            width = 6, height = 5)

    # Save statistical test results
    test_file <- file.path(subset_output_dir, "statistical_tests.tsv")
    test_df <- data.frame(
      subset = subset_name,
      test = "Wilcoxon rank-sum",
      statistic = enrichment_result$test$statistic,
      p_value = enrichment_result$test$p.value,
      cliff_delta = enrichment_result$cliff_delta,
      n_ctrl = enrichment_result$n_ctrl,
      n_mut = enrichment_result$n_mut,
      mean_ctrl = enrichment_result$mean_ctrl,
      mean_mut = enrichment_result$mean_mut,
      expected_enriched = expected_enriched,
      actual_enriched = ifelse(enrichment_result$mean_ctrl > enrichment_result$mean_mut,
                               "ctrl", "mut")
    )
    write.table(test_df, test_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Statistical test: p = %.2e\n", enrichment_result$test$p.value))
  }

  # Cleanup HDF5 temporary files
  hdf5_dir <- file.path(subset_output_dir, sprintf("temp_hdf5_%s", subset_name))
  if (dir.exists(hdf5_dir)) {
    unlink(hdf5_dir, recursive = TRUE)
  }

  cat(sprintf("  All outputs saved to: %s\n", subset_output_dir))

  return(TRUE)
}

# ============================================================================
# Main Execution
# ============================================================================

main <- function() {

  # Check .hic file availability
  cat("\nChecking .hic files...\n")
  hic_check <- check_hic_files(HIC_FILES)

  if (!hic_check$files_exist) {
    cat("\nCannot proceed without .hic files\n")
    cat("Expected locations:\n")
    for (f in names(HIC_FILES)) {
      cat(sprintf("  %s: %s\n", f, HIC_FILES[f]))
    }
    quit(status = 1)
  }

  norm_method <- hic_check$available_norm

  # Check input files
  cat("\nChecking input files...\n")
  for (name in names(INPUT_FILES)) {
    if (!file.exists(INPUT_FILES[[name]])) {
      cat(sprintf("  Missing: %s\n", INPUT_FILES[[name]]))
      quit(status = 1)
    }
    cat(sprintf("  Found: %s\n", INPUT_FILES[[name]]))
  }

  # Summary counters
  success_count <- 0
  total_analyses <- 0
  results_summary <- list()

  # -------------------------
  # Process: Lost Long-Range Loops
  # -------------------------
  cat("\n========================================\n")
  cat("Processing: Lost Long-Range Loops\n")
  cat("(Expected enrichment: CONTROL > Mutant)\n")
  cat("========================================\n")

  loops_lost <- readRDS(INPUT_FILES$lost_longrange)
  cat(sprintf("Loaded %d loops from %s\n", length(loops_lost), INPUT_FILES$lost_longrange))
  validate_loop_resolutions(loops_lost, "lost_longrange")

  total_analyses <- total_analyses + 1
  if (run_apa_analysis(loops_lost, "lost_longrange", OUTPUT_DIR, HIC_FILES,
                       norm_method, expected_enriched = "ctrl")) {
    success_count <- success_count + 1
    results_summary$lost_longrange <- "SUCCESS"
  } else {
    results_summary$lost_longrange <- "FAILED"
  }

  # -------------------------
  # Process: Gained Short-Range Loops
  # -------------------------
  cat("\n========================================\n")
  cat("Processing: Gained Short-Range Loops\n")
  cat("(Expected enrichment: MUTANT > Control)\n")
  cat("========================================\n")

  loops_gained <- readRDS(INPUT_FILES$gained_shortrange)
  cat(sprintf("Loaded %d loops from %s\n", length(loops_gained), INPUT_FILES$gained_shortrange))
  validate_loop_resolutions(loops_gained, "gained_shortrange")

  total_analyses <- total_analyses + 1
  if (run_apa_analysis(loops_gained, "gained_shortrange", OUTPUT_DIR, HIC_FILES,
                       norm_method, expected_enriched = "mut")) {
    success_count <- success_count + 1
    results_summary$gained_shortrange <- "SUCCESS"
  } else {
    results_summary$gained_shortrange <- "FAILED"
  }

  # ============================================================================
  # Final Summary
  # ============================================================================

  cat("\n========================================\n")
  cat("APA Analysis Complete\n")
  cat("========================================\n")
  cat(sprintf("Successful analyses: %d / %d\n", success_count, total_analyses))
  cat(sprintf("Results directory: %s\n", OUTPUT_DIR))

  # Write summary report
  summary_file <- file.path(OUTPUT_DIR, "summary_report.txt")
  sink(summary_file)
  cat("APA Analysis Summary Report: Shared Anchor Loop Subsets\n")
  cat(sprintf("Generated: %s\n", Sys.time()))
  cat(sprintf("Timepoint: %s\n\n", TIMEPOINT))

  cat("Analysis Parameters:\n")
  cat(sprintf("  Resolution: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
  cat(sprintf("  Window: %d bins (%d kb)\n", 2 * BUFFER_BINS + 1, (2 * BUFFER_BINS + 1) * RESOLUTION / 1000))
  cat(sprintf("  Normalization: %s\n\n", norm_method))

  cat("Results:\n")
  cat(sprintf("  Lost Long-Range: %s (%d loops)\n",
              results_summary$lost_longrange, length(loops_lost)))
  cat(sprintf("  Gained Short-Range: %s (%d loops)\n",
              results_summary$gained_shortrange, length(loops_gained)))

  cat("\nExpected Biological Interpretation:\n")
  cat("  - Lost loops: Should show STRONGER signal in Control (loops lost in mutant)\n")
  cat("  - Gained loops: Should show STRONGER signal in Mutant (loops gained)\n")

  cat("\nOutput Structure:\n")
  cat(sprintf("  %s/\n", OUTPUT_DIR))
  cat("    lost_longrange/\n")
  cat("      aggregate_heatmap_ctrl.{pdf,svg,jpg}\n")
  cat("      aggregate_heatmap_mut.{pdf,svg,jpg}\n")
  cat("      difference_heatmap.{pdf,svg,jpg}\n")
  cat("      enrichment_boxplot.{pdf,svg,jpg}\n")
  cat("      enrichment_scores.tsv\n")
  cat("      statistical_tests.tsv\n")
  cat("    gained_shortrange/\n")
  cat("      (same structure)\n")
  cat("    summary_report.txt\n")
  sink()

  cat(sprintf("\nSummary saved: %s\n\n", summary_file))

  # Print statistical results
  cat("Statistical Results:\n")
  for (subset in c("lost_longrange", "gained_shortrange")) {
    test_file <- file.path(OUTPUT_DIR, subset, "statistical_tests.tsv")
    if (file.exists(test_file)) {
      test_df <- read.table(test_file, sep = "\t", header = TRUE)
      cat(sprintf("\n  %s:\n", subset))
      cat(sprintf("    p-value: %.2e\n", test_df$p_value))
      cat(sprintf("    Cliff's delta: %.3f\n", test_df$cliff_delta))
      cat(sprintf("    Mean enrichment - Control: %.2f, Mutant: %.2f\n",
                  test_df$mean_ctrl, test_df$mean_mut))
      cat(sprintf("    Expected enriched: %s, Actual enriched: %s\n",
                  test_df$expected_enriched, test_df$actual_enriched))
    }
  }

  cat("\n")
  return(invisible(NULL))
}

# Run main function
main()
