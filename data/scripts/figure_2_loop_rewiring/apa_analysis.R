# scripts/apa_analysis.R

#' Aggregate Peak Analysis (APA) for Differential Chromatin Loops
#'
#' Performs comprehensive APA analysis across multiple resolutions comparing
#' Hi-C contact enrichment between BAP1-KO and control conditions.
#'
#' Pipeline: Extract Hi-C matrices → Aggregate → Calculate enrichment → Visualize
#'
#' @usage Rscript scripts/apa_analysis.R [--resolution RESOLUTION] [--loops LOOPS]
#' @param resolution Optional: specific resolution (5000, 10000, or 25000). Default: all three
#' @param loops Optional: "merged", "resolution_specific", or "both". Default: both
#'
#' @examples
#' # Run for all resolutions and loop sets (default)
#' Rscript scripts/apa_analysis.R
#'
#' # Run only 5kb resolution with merged loops
#' Rscript scripts/apa_analysis.R --resolution 5000 --loops merged

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
})

# Load multi-format output utility for PDF + SVG + JPEG output
source("data/scripts/_shared/multi_format_output.R")  # Original: source("scripts/utils/multi_format_output.R")

# Load paths configuration
cat("\nLoading paths configuration...\n")
config <- yaml::read_yaml("config/paths_config.yaml")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
TARGET_RESOLUTION <- NULL
TARGET_LOOPS <- "both"  # "merged", "resolution_specific", or "both"

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--resolution" && i < length(args)) {
      TARGET_RESOLUTION <- as.numeric(args[i + 1])
    }
    if (args[i] == "--loops" && i < length(args)) {
      TARGET_LOOPS <- args[i + 1]
    }
  }
}

# Define resolutions to process
RESOLUTIONS <- if (!is.null(TARGET_RESOLUTION)) {
  TARGET_RESOLUTION
} else {
  c(5000, 10000, 25000)
}

# .hic file paths from config (match extract_counts.R pattern)
# TODO: not in data/ — .hic files are HPC-only, paths read from config
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

# Output directories
# Original: OUTPUT_BASE <- "outputs/apa_results"
APA_PLOT_DIR <- "data/plots/figure_2_loop_rewiring/apa"
APA_TSV_DIR  <- "data/tsvs/figure_2_loop_rewiring"
if (!dir.exists(APA_PLOT_DIR)) {
  dir.create(APA_PLOT_DIR, recursive = TRUE)
}
if (!dir.exists(APA_TSV_DIR)) {
  dir.create(APA_TSV_DIR, recursive = TRUE)
}

cat("\n========================================\n")
cat("Aggregate Peak Analysis (APA) Pipeline\n")
cat("========================================\n\n")
cat(sprintf("Resolutions to process: %s\n", paste(RESOLUTIONS/1000, "kb", collapse = ", ")))
cat(sprintf("Loop sets: %s\n", TARGET_LOOPS))
cat(sprintf("Plot directory: %s\n", APA_PLOT_DIR))
cat(sprintf("TSV directory: %s\n\n", APA_TSV_DIR))

# ============================================================================
# Core Functions
# ============================================================================

#' Calculate buffer size for 50kb window at given resolution
#'
#' @param resolution Bin resolution in base pairs
#' @return Number of bins needed for 50kb window
get_buffer_size <- function(resolution) {
  buffer_bp <- 50000  # 50kb window
  buffer_bins <- floor(buffer_bp / resolution)
  return(buffer_bins)
}

#' Check .hic file availability and normalization methods
#'
#' @param hic_files Named vector of .hic file paths
#' @return List with files_exist (logical) and available_norm (character vector)
check_hic_files <- function(hic_files) {
  files_exist <- all(sapply(hic_files, file.exists))

  if (!files_exist) {
    missing <- names(hic_files)[!sapply(hic_files, file.exists)]
    cat("  ⚠ Warning: Missing .hic files:\n")
    for (f in missing) {
      cat(sprintf("    - %s: %s\n", f, hic_files[f]))
    }
    return(list(files_exist = FALSE, available_norm = NULL))
  }

  # Check normalization methods (use first file as representative)
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

  cat(sprintf("  ✓ All .hic files found\n"))
  cat(sprintf("  ✓ Using normalization: %s\n\n", norm_choice))

  return(list(files_exist = TRUE, available_norm = norm_choice))
}

#' Extract Hi-C matrices around loop positions for APA
#'
#' @param loops GInteractions object with loop coordinates
#' @param hic_files Named vector of .hic file paths
#' @param resolution Bin resolution
#' @param buffer Number of bins to extend around loop
#' @param norm Normalization method ("KR", "VC", or "NONE")
#' @return InteractionArray with extracted matrices
extract_apa_matrices <- function(loops, hic_files, resolution, buffer, norm = "KR") {

  cat(sprintf("  Extracting matrices: %d loops × %d samples\n", length(loops), length(hic_files)))
  cat(sprintf("  Buffer: %d bins (%dkb window)\n", buffer, (buffer * 2 + 1) * resolution / 1000))

  # Create temporary HDF5 directory for on-disk storage
  hdf5_dir <- file.path(APA_PLOT_DIR, sprintf("temp_hdf5_res_%dkb", resolution/1000))  # Original: file.path(OUTPUT_BASE, ...)
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

    cat(sprintf("  ✓ Extraction complete: %s\n", hdf5_file))
    return(pixels)

  }, error = function(e) {
    cat(sprintf("  ✗ Error during extraction: %s\n", e$message))
    return(NULL)
  })
}

#' Aggregate matrices across loops to create mean contact map
#'
#' @param pixels InteractionArray from extract_apa_matrices
#' @param method Aggregation method ("mean" or "median")
#' @return Matrix of aggregated contacts (bins × bins)
aggregate_apa_matrices <- function(pixels, method = "mean") {

  if (is.null(pixels)) {
    return(NULL)
  }

  # Use mariner's aggHicMatrices if available, otherwise manual aggregation
  tryCatch({
    # Get dimensions
    count_array <- counts(pixels)
    dims <- dim(count_array)
    n_bins <- dims[1]  # Should equal dims[2] (square matrix)
    n_loops <- dims[3]
    n_samples <- dims[4]

    cat(sprintf("  Aggregating: %d loops × %d samples → %d×%d matrix per sample\n",
                n_loops, n_samples, n_bins, n_bins))

    # Aggregate per sample
    agg_matrices <- array(NA, dim = c(n_bins, n_bins, n_samples))

    for (j in 1:n_samples) {
      # Convert DelayedArray slice to regular array for apply() with vector MARGIN
      sample_matrices <- as.array(count_array[, , , j, drop = FALSE])

      if (method == "mean") {
        agg_matrices[, , j] <- apply(sample_matrices, c(1, 2), mean, na.rm = TRUE)
      } else if (method == "median") {
        agg_matrices[, , j] <- apply(sample_matrices, c(1, 2), median, na.rm = TRUE)
      }
    }

    cat(sprintf("  ✓ Aggregation complete using %s\n", method))
    return(agg_matrices)

  }, error = function(e) {
    cat(sprintf("  ✗ Error during aggregation: %s\n", e$message))
    return(NULL)
  })
}

#' Calculate P2LL enrichment scores (center / local background)
#'
#' @param pixels InteractionArray from extract_apa_matrices
#' @return Data frame with per-loop enrichment scores per sample
calculate_enrichment_scores <- function(pixels) {

  if (is.null(pixels)) {
    return(NULL)
  }

  tryCatch({
    count_array <- counts(pixels)
    dims <- dim(count_array)
    n_bins <- dims[1]
    n_loops <- dims[3]
    n_samples <- dims[4]

    cat(sprintf("  Calculating P2LL enrichment scores for %d loops\n", n_loops))

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
      sample = rep(names(HIC_FILES)[1:n_samples], times = n_loops),
      group = rep(GROUPS, times = n_loops),
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

    # Filter outliers (>3 SD from mean)
    valid_enrichment <- enrichment_df$enrichment[!is.na(enrichment_df$enrichment) &
                                                   is.finite(enrichment_df$enrichment)]
    if (length(valid_enrichment) > 0) {
      mean_enrich <- mean(valid_enrichment)
      sd_enrich <- sd(valid_enrichment)
      outlier_thresh <- 3

      n_outliers <- sum(abs(enrichment_df$enrichment - mean_enrich) > outlier_thresh * sd_enrich,
                        na.rm = TRUE)
      if (n_outliers > 0) {
        cat(sprintf("  ⚠ Flagged %d outlier enrichment scores (>3 SD)\n", n_outliers))
      }
    }

    cat(sprintf("  ✓ Enrichment calculation complete\n"))
    cat(sprintf("    Mean enrichment: %.2f (range: %.2f - %.2f)\n",
                mean(valid_enrichment, na.rm = TRUE),
                min(valid_enrichment, na.rm = TRUE),
                max(valid_enrichment, na.rm = TRUE)))

    return(enrichment_df)

  }, error = function(e) {
    cat(sprintf("  ✗ Error calculating enrichment: %s\n", e$message))
    return(NULL)
  })
}

#' Generate APA heatmap for a single condition
#'
#' @param agg_matrix Aggregated matrix (bins × bins × samples)
#' @param condition_name Name for plot title ("Control" or "BAP1-KO")
#' @param sample_indices Which samples to average (numeric vector)
#' @param resolution Bin resolution for axis labels
#' @return ggplot object
plot_apa_heatmap <- function(agg_matrix, condition_name, sample_indices, resolution) {

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
  axis_labels <- c(sprintf("-%dkb", buffer_kb), "3'", sprintf("+%dkb", buffer_kb))

  # Plot
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis(name = "µ Contacts", option = "magma") +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = sprintf("APA: %s", condition_name),
      x = "",
      y = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )

  return(p)
}

#' Generate difference heatmap (Mutant - Control)
#'
#' @param agg_matrix Aggregated matrix (bins × bins × samples)
#' @param ctrl_indices Control sample indices
#' @param mut_indices Mutant sample indices
#' @param resolution Bin resolution for axis labels
#' @return ggplot object
plot_difference_heatmap <- function(agg_matrix, ctrl_indices, mut_indices, resolution) {

  if (is.null(agg_matrix)) {
    return(NULL)
  }

  # Calculate mean matrices
  ctrl_mean <- apply(agg_matrix[, , ctrl_indices, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  mut_mean <- apply(agg_matrix[, , mut_indices, drop = FALSE], c(1, 2), mean, na.rm = TRUE)

  # Difference
  diff_matrix <- mut_mean - ctrl_mean

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
  axis_labels <- c(sprintf("-%dkb", buffer_kb), "3'", sprintf("+%dkb", buffer_kb))

  # Plot with RdBu diverging scale
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      name = "Difference",
      low = scales::muted("blue"),
      mid = "white",
      high = scales::muted("red"),
      midpoint = 0,
      limits = c(-5, 5),
      oob = scales::squish
    ) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = "Difference: BAP1-KO - Control",
      x = "",
      y = ""
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
#'
#' @param enrichment_df Data frame from calculate_enrichment_scores
#' @param direction Loop direction ("Up-regulated" or "Down-regulated")
#' @return List with ggplot object and test results
plot_enrichment_comparison <- function(enrichment_df, direction = "Up-regulated") {

  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
    return(NULL)
  }

  # Filter valid enrichment scores
  df_valid <- enrichment_df %>%
    filter(!is.na(enrichment), is.finite(enrichment))

  if (nrow(df_valid) == 0) {
    cat("  ⚠ No valid enrichment scores for plotting\n")
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

  # Create plot
  p <- ggplot(df_valid, aes(x = group, y = enrichment, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
    scale_fill_manual(values = c("ctrl" = "grey60", "mut" = "firebrick")) +
    scale_x_discrete(labels = c("ctrl" = "Control", "mut" = "BAP1-KO")) +
    labs(
      title = sprintf("Pixel Enrichment: %s Loops", direction),
      subtitle = sprintf("Wilcoxon p = %.2e, Cliff's Δ = %.3f",
                        test_result$p.value, cliff_delta),
      x = "",
      y = "Pixel Enrichment (P2LL)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  return(list(
    plot = p,
    test = test_result,
    cliff_delta = cliff_delta,
    n_ctrl = length(ctrl_vals),
    n_mut = length(mut_vals)
  ))
}

#' Generate per-replicate heatmaps (6 panels)
#'
#' @param agg_matrix Aggregated matrix (bins × bins × samples)
#' @param resolution Bin resolution for axis labels
#' @param sample_names Names of samples
#' @return Combined ggplot object
plot_replicate_heatmaps <- function(agg_matrix, resolution, sample_names) {

  if (is.null(agg_matrix)) {
    return(NULL)
  }

  n_samples <- dim(agg_matrix)[3]
  n_bins <- dim(agg_matrix)[1]
  buffer_kb <- (n_bins - 1) / 2 * resolution / 1000

  # Prepare data for all samples
  plot_list <- list()

  for (i in 1:n_samples) {
    mat <- agg_matrix[, , i]

    df <- expand.grid(
      x = 1:n_bins,
      y = 1:n_bins
    )
    df$value <- as.vector(mat)

    center_idx <- ceiling(n_bins / 2)
    axis_breaks <- c(1, center_idx, n_bins)
    axis_labels <- c(sprintf("-%dkb", buffer_kb), "3'", sprintf("+%dkb", buffer_kb))

    p <- ggplot(df, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_viridis(name = "Contacts", option = "magma") +
      scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
      scale_y_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
      coord_fixed() +
      labs(title = sample_names[i]) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        panel.grid = element_blank(),
        legend.key.height = unit(0.5, "cm")
      )

    plot_list[[i]] <- p
  }

  # Combine plots using patchwork if available, otherwise gridExtra
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plot_list, ncol = 3)
  } else {
    combined <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)
  }

  return(combined)
}

#' Main APA analysis orchestrator
#'
#' @param resolution Bin resolution
#' @param loops_gi GInteractions object with loops
#' @param loop_set_name Name for output directory ("merged_loops" or "resolution_specific")
#' @param direction Loop direction ("up" or "down")
#' @param hic_files Named vector of .hic file paths
#' @param norm Normalization method
#' @return Logical indicating success
run_apa_for_loop_set <- function(resolution, loops_gi, loop_set_name, direction,
                                  hic_files, norm) {

  cat(sprintf("\n--- Processing: %s %s-regulated loops at %dkb ---\n",
              loop_set_name, direction, resolution/1000))

  # Create output directories
  # Original: output_dir <- file.path(OUTPUT_BASE, sprintf("res_%dkb", resolution/1000), loop_set_name, sprintf("%s_loops", direction))
  plot_subdir <- file.path(APA_PLOT_DIR, sprintf("res_%dkb", resolution/1000), loop_set_name, sprintf("%s_loops", direction))
  if (!dir.exists(plot_subdir)) {
    dir.create(plot_subdir, recursive = TRUE)
  }
  if (!dir.exists(APA_TSV_DIR)) {
    dir.create(APA_TSV_DIR, recursive = TRUE)
  }
  # Prefix for flat TSV filenames
  tsv_prefix <- sprintf("2C_%dkb_%s_%s", resolution/1000, loop_set_name, direction)

  # Check minimum loop count
  n_loops <- length(loops_gi)
  if (n_loops < 10) {
    cat(sprintf("  ⚠ Skipping: Only %d loops (minimum 10 required)\n", n_loops))
    return(FALSE)
  }

  cat(sprintf("  Loop count: %d\n", n_loops))

  # Step 1: Extract matrices
  buffer <- get_buffer_size(resolution)
  pixels <- extract_apa_matrices(loops_gi, hic_files, resolution, buffer, norm)

  if (is.null(pixels)) {
    cat("  ✗ Matrix extraction failed\n")
    return(FALSE)
  }

  # Step 2: Aggregate matrices
  agg_matrices <- aggregate_apa_matrices(pixels, method = "mean")

  if (is.null(agg_matrices)) {
    cat("  ✗ Aggregation failed\n")
    return(FALSE)
  }

  # Step 3: Calculate enrichment scores
  enrichment_df <- calculate_enrichment_scores(pixels)

  if (!is.null(enrichment_df)) {
    enrichment_file <- file.path(APA_TSV_DIR, sprintf("%s_enrichment_scores.tsv", tsv_prefix))  # Original: file.path(output_dir, "enrichment_scores.tsv")
    write.table(enrichment_df, enrichment_file, sep = "\t",
                quote = FALSE, row.names = FALSE)
    cat(sprintf("  ✓ Saved: %s\n", enrichment_file))
  }

  # Step 4: Generate plots
  cat("\n  Generating visualizations...\n")

  # 4a. Control aggregate heatmap
  ctrl_indices <- which(GROUPS == "ctrl")
  mut_indices <- which(GROUPS == "mut")

  p_ctrl <- plot_apa_heatmap(agg_matrices, "Control", ctrl_indices, resolution)
  if (!is.null(p_ctrl)) {
    save_multiformat_ggplot(p_ctrl, file.path(plot_subdir, "2C_aggregate_heatmap_ctrl"), width = 6, height = 5)  # Original: file.path(output_dir, "aggregate_heatmap_ctrl")
  }

  # 4b. Mutant aggregate heatmap
  p_mut <- plot_apa_heatmap(agg_matrices, "BAP1-KO", mut_indices, resolution)
  if (!is.null(p_mut)) {
    save_multiformat_ggplot(p_mut, file.path(plot_subdir, "2C_aggregate_heatmap_mut"), width = 6, height = 5)  # Original: file.path(output_dir, "aggregate_heatmap_mut")
  }

  # 4c. Difference heatmap
  p_diff <- plot_difference_heatmap(agg_matrices, ctrl_indices, mut_indices, resolution)
  if (!is.null(p_diff)) {
    save_multiformat_ggplot(p_diff, file.path(plot_subdir, "2C_difference_heatmap"), width = 6, height = 5)  # Original: file.path(output_dir, "difference_heatmap")
  }

  # 4d. Enrichment box plot
  direction_label <- ifelse(direction == "up", "Up-regulated", "Down-regulated")
  enrichment_result <- plot_enrichment_comparison(enrichment_df, direction_label)

  if (!is.null(enrichment_result)) {
    save_multiformat_ggplot(enrichment_result$plot, file.path(plot_subdir, "2C_enrichment_boxplot"), width = 6, height = 5)  # Original: file.path(output_dir, "enrichment_boxplot")

    # Save statistical test results
    test_file <- file.path(APA_TSV_DIR, sprintf("%s_statistical_tests.tsv", tsv_prefix))  # Original: file.path(output_dir, "statistical_tests.tsv")
    test_df <- data.frame(
      test = "Wilcoxon rank-sum",
      statistic = enrichment_result$test$statistic,
      p_value = enrichment_result$test$p.value,
      cliff_delta = enrichment_result$cliff_delta,
      n_ctrl = enrichment_result$n_ctrl,
      n_mut = enrichment_result$n_mut
    )
    write.table(test_df, test_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  ✓ Statistical test: p = %.2e\n", enrichment_result$test$p.value))
  }

  # 4e. Per-replicate heatmaps
  p_replicates <- plot_replicate_heatmaps(agg_matrices, resolution, names(hic_files))
  if (!is.null(p_replicates)) {
    save_multiformat_ggplot(p_replicates, file.path(plot_subdir, "2C_replicate_heatmaps"), width = 12, height = 8)  # Original: file.path(output_dir, "replicate_heatmaps")
  }

  cat(sprintf("  ✓ All outputs saved to: %s (plots), %s (TSVs)\n", plot_subdir, APA_TSV_DIR))

  # Cleanup HDF5 temporary files
  hdf5_dir <- file.path(APA_PLOT_DIR, sprintf("temp_hdf5_res_%dkb", resolution/1000))  # Original: file.path(OUTPUT_BASE, ...)
  if (dir.exists(hdf5_dir)) {
    unlink(hdf5_dir, recursive = TRUE)
  }

  return(TRUE)
}

# ============================================================================
# Main Execution
# ============================================================================

main <- function() {

  # Check .hic file availability
  hic_check <- check_hic_files(HIC_FILES)
  if (!hic_check$files_exist) {
    cat("\n✗ Cannot proceed without .hic files\n")
    cat("  Expected locations:\n")
    for (f in names(HIC_FILES)) {
      cat(sprintf("    %s: %s\n", f, HIC_FILES[f]))
    }
    quit(status = 1)
  }

  norm_method <- hic_check$available_norm

  # Summary counters
  success_count <- 0
  total_analyses <- 0

  # Process each resolution
  for (resolution in RESOLUTIONS) {

    cat(sprintf("\n========================================\n"))
    cat(sprintf("Resolution: %dkb\n", resolution/1000))
    cat(sprintf("========================================\n"))

    res_kb <- resolution / 1000

    # -------------------------
    # MERGED LOOPS
    # -------------------------
    if (TARGET_LOOPS %in% c("both", "merged")) {

      cat("\n## Processing MERGED loops ##\n")

      merged_file <- "outputs/merged_loops/non_redundant_loops.rds"  # TODO: not in data/

      if (!file.exists(merged_file)) {
        cat(sprintf("  ⚠ Merged loops file not found: %s\n", merged_file))
        cat("    Run downstream_analysis.R first\n")
      } else {
        loops_gi <- readRDS(merged_file)

        # Bin merged loops to current resolution and buffer for matrix extraction
        buffer <- get_buffer_size(resolution)
        loops_gi <- assignToBins(
          x = loops_gi,
          binSize = resolution,
          pos1 = "center",
          pos2 = "center"
        )
        loops_gi <- pixelsToMatrices(
          x = loops_gi,
          buffer = buffer
        )

        # Separate by direction
        up_loops <- loops_gi[mcols(loops_gi)$logFC > 0]
        down_loops <- loops_gi[mcols(loops_gi)$logFC < 0]

        cat(sprintf("  Total merged loops: %d (up: %d, down: %d)\n",
                    length(loops_gi), length(up_loops), length(down_loops)))

        # Process up-regulated
        total_analyses <- total_analyses + 1
        if (run_apa_for_loop_set(resolution, up_loops, "merged_loops", "up",
                                   HIC_FILES, norm_method)) {
          success_count <- success_count + 1
        }

        # Process down-regulated
        total_analyses <- total_analyses + 1
        if (run_apa_for_loop_set(resolution, down_loops, "merged_loops", "down",
                                   HIC_FILES, norm_method)) {
          success_count <- success_count + 1
        }
      }
    }

    # -------------------------
    # RESOLUTION-SPECIFIC LOOPS
    # -------------------------
    if (TARGET_LOOPS %in% c("both", "resolution_specific")) {

      cat("\n## Processing RESOLUTION-SPECIFIC loops ##\n")

      binned_file <- file.path("outputs", sprintf("res_%dkb", res_kb), "03_binned.rds")  # TODO: not in data/
      results_file <- file.path("outputs", sprintf("edgeR_results_res_%dkb", res_kb),
                                "primary_analysis", "all_results_primary.tsv")  # TODO: not in data/

      if (!file.exists(binned_file) || !file.exists(results_file)) {
        cat(sprintf("  ⚠ Required files not found:\n"))
        cat(sprintf("    Binned: %s\n", binned_file))
        cat(sprintf("    Results: %s\n", results_file))
      } else {
        # Load binned (unbuffered) GInteractions - coordinates match edgeR results
        loops_gi <- readRDS(binned_file)

        # Load results to get significance and logFC
        results_df <- read.table(results_file, sep = "\t", header = TRUE)

        # Match loops with results using coordinate-based matching
        # This handles cases where edgeR filtering removes some loops

        # Create coordinate IDs for results
        results_df$coord_id <- paste(
          results_df$chr1, results_df$start1, results_df$end1,
          results_df$chr2, results_df$start2, results_df$end2,
          sep = "_"
        )

        # Create coordinate IDs for GInteractions
        gi_coords <- paste(
          as.character(seqnames(anchors(loops_gi, "first"))),
          start(anchors(loops_gi, "first")),
          end(anchors(loops_gi, "first")),
          as.character(seqnames(anchors(loops_gi, "second"))),
          start(anchors(loops_gi, "second")),
          end(anchors(loops_gi, "second")),
          sep = "_"
        )

        # Match results to GI coordinates
        match_idx <- match(results_df$coord_id, gi_coords)

        if (any(is.na(match_idx))) {
          cat(sprintf("  ⚠ Warning: %d results could not be matched to GI\n",
                      sum(is.na(match_idx))))
          # Remove unmatched results
          results_df <- results_df[!is.na(match_idx), ]
          match_idx <- match_idx[!is.na(match_idx)]
        }

        # Subset GI to matched loops and add metadata
        loops_gi_matched <- loops_gi[match_idx]
        mcols(loops_gi_matched)$logFC <- results_df$logFC
        mcols(loops_gi_matched)$FDR <- results_df$FDR
        mcols(loops_gi_matched)$significant <- results_df$FDR < 0.05

        cat(sprintf("  Matched %d loops from GI (%d) to results (%d)\n",
                    length(loops_gi_matched), length(loops_gi), nrow(results_df)))

        # Filter to significant only
        sig_loops <- loops_gi_matched[mcols(loops_gi_matched)$significant]

        # Separate by direction
        up_loops <- sig_loops[mcols(sig_loops)$logFC > 0]
        down_loops <- sig_loops[mcols(sig_loops)$logFC < 0]

        cat(sprintf("  Significant loops: %d (up: %d, down: %d)\n",
                    length(sig_loops), length(up_loops), length(down_loops)))

        # Process up-regulated
        total_analyses <- total_analyses + 1
        if (run_apa_for_loop_set(resolution, up_loops, "resolution_specific", "up",
                                   HIC_FILES, norm_method)) {
          success_count <- success_count + 1
        }

        # Process down-regulated
        total_analyses <- total_analyses + 1
        if (run_apa_for_loop_set(resolution, down_loops, "resolution_specific", "down",
                                   HIC_FILES, norm_method)) {
          success_count <- success_count + 1
        }
      }
    }
  }

  # ============================================================================
  # Final Summary
  # ============================================================================

  cat("\n========================================\n")
  cat("APA Analysis Complete\n")
  cat("========================================\n")
  cat(sprintf("Successful analyses: %d / %d\n", success_count, total_analyses))
  cat(sprintf("Plot directory: %s\n", APA_PLOT_DIR))
  cat(sprintf("TSV directory: %s\n", APA_TSV_DIR))

  # Write summary report
  summary_file <- file.path(APA_TSV_DIR, "2C_apa_summary_report.txt")  # Original: file.path(OUTPUT_BASE, "summary_report.txt")
  sink(summary_file)
  cat("APA Analysis Summary Report\n")
  cat(sprintf("Generated: %s\n\n", Sys.time()))
  cat(sprintf("Resolutions processed: %s\n",
              paste(RESOLUTIONS/1000, "kb", collapse = ", ")))
  cat(sprintf("Loop sets: %s\n", TARGET_LOOPS))
  cat(sprintf("Normalization: %s\n", norm_method))
  cat(sprintf("Successful analyses: %d / %d\n\n", success_count, total_analyses))

  cat("Output structure:\n")
  cat(sprintf("  Plots: %s/\n", APA_PLOT_DIR))
  cat(sprintf("  TSVs:  %s/\n", APA_TSV_DIR))
  for (res in RESOLUTIONS) {
    cat(sprintf("    res_%dkb/\n", res/1000))
    if (TARGET_LOOPS %in% c("both", "merged")) {
      cat("      merged_loops/{up,down}_loops/\n")
    }
    if (TARGET_LOOPS %in% c("both", "resolution_specific")) {
      cat("      resolution_specific/{up,down}_loops/\n")
    }
  }

  cat("\nFiles generated per analysis:\n")
  cat("  - aggregate_heatmap_ctrl.pdf\n")
  cat("  - aggregate_heatmap_mut.pdf\n")
  cat("  - difference_heatmap.pdf\n")
  cat("  - enrichment_boxplot.pdf\n")
  cat("  - replicate_heatmaps.pdf\n")
  cat("  - enrichment_scores.tsv\n")
  cat("  - statistical_tests.tsv\n")
  sink()

  cat(sprintf("\n✓ Summary saved: %s\n\n", summary_file))

  return(invisible(NULL))
}

# Run main function
main()
