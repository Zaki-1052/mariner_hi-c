# scripts/apa_se_analysis.R
#
# APA analysis for SE-DEG contact frequency
#
# Measures Hi-C contact enrichment between DEG promoters and superenhancers
# vs. DEG promoters and regular enhancers vs. invariant gene promoters and SEs.
# Compares control vs BAP1-KO conditions across biological replicates.
#
# Requires HPC execution (.hic files on Expanse).
#
# Usage:
#   Rscript scripts/apa_se_analysis.R --timepoint late
#   Rscript scripts/apa_se_analysis.R --timepoint late --resolution 10000
#   Rscript scripts/apa_se_analysis.R --timepoint both
#
# Input:
#   - BEDPE pairs from se_deg_proximity.R (apa_inputs/)
#   - Gene table with gene_class from se_deg_proximity.R
#   - SE pair metadata (se_pairs_with_k27ac.tsv) for K27ac sub-stratification
#   - .hic files (6 samples, paths from config/paths_config.yaml)
#
# Output:
#   loops/output/superenhancer_analysis/1a_se_deg_proximity/apa_results/{timepoint}/

suppressPackageStartupMessages({
  library(mariner)
  library(InteractionSet)
  library(strawr)
  library(HDF5Array)
  library(SummarizedExperiment)
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(yaml)
})

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

BASE_DIR <- getwd()

source(file.path(BASE_DIR, "loops/scripts/utils/multi_format_output.R"))

HIC_ROOT <- "/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic"

TIMEPOINT_HIC_MAP <- list(
  late  = "250402",
  early = "250831"
)

SAMPLE_NAMES <- c("ctrl_M1", "ctrl_M2", "ctrl_M3", "mut_M1", "mut_M2", "mut_M3")

build_hic_paths <- function(timepoint) {
  tp_dir <- TIMEPOINT_HIC_MAP[[timepoint]]
  if (is.null(tp_dir)) stop(sprintf("Unknown timepoint: %s", timepoint))
  paths <- file.path(HIC_ROOT, tp_dir, paste0(SAMPLE_NAMES, ".hic"))
  names(paths) <- SAMPLE_NAMES
  paths
}

GROUPS <- factor(c("ctrl", "ctrl", "ctrl", "mut", "mut", "mut"))
CTRL_INDICES <- which(GROUPS == "ctrl")
MUT_INDICES <- which(GROUPS == "mut")

APA_INPUT_DIR <- file.path(BASE_DIR,
  "loops/output/superenhancer_analysis/1a_se_deg_proximity")

PAIR_CLASSES <- list(
  deg_se = list(
    file_pattern = "%s_deg_se_pairs.bedpe",
    label = "DEG-SE",
    description = "DEG promoters x Superenhancers"
  ),
  invariant_se = list(
    file_pattern = "%s_invariant_se_pairs.bedpe",
    label = "Invariant-SE",
    description = "Invariant gene promoters x Superenhancers"
  )
)

BUFFER_BP <- 50000
MIN_PAIRS <- 10

# =============================================================================
# 2. ARGUMENT PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
TIMEPOINT <- "late"
TARGET_RESOLUTION <- 10000

i <- 1
while (i <= length(args)) {
  if (args[i] == "--timepoint" && i < length(args)) {
    TIMEPOINT <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--resolution" && i < length(args)) {
    TARGET_RESOLUTION <- as.numeric(args[i + 1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

# =============================================================================
# 3. CORE FUNCTIONS
# =============================================================================

check_hic_files <- function(hic_files) {
  files_exist <- sapply(hic_files, file.exists)
  if (!all(files_exist)) {
    missing <- names(hic_files)[!files_exist]
    cat("Missing .hic files:\n")
    for (f in missing) cat(sprintf("  - %s: %s\n", f, hic_files[f]))
    stop("Cannot proceed without .hic files")
  }

  norms <- tryCatch(strawr::readHicNormTypes(hic_files[1]),
                     error = function(e) "NONE")
  norm <- if ("KR" %in% norms) "KR" else if ("VC" %in% norms) "VC" else "NONE"
  cat(sprintf("  All .hic files found. Normalization: %s\n", norm))
  return(norm)
}

load_bedpe_as_gi <- function(bedpe_file, resolution) {
  if (!file.exists(bedpe_file)) {
    stop(sprintf("BEDPE not found: %s", bedpe_file))
  }

  bedpe <- read.table(bedpe_file, header = FALSE, sep = "\t",
                      col.names = c("chr1", "x1", "x2", "chr2", "y1", "y2"))

  if (nrow(bedpe) == 0) stop(sprintf("Empty BEDPE: %s", bedpe_file))

  gi <- GInteractions(
    anchor1 = GRanges(seqnames = bedpe$chr1,
                      ranges = IRanges(start = bedpe$x1, end = bedpe$x2)),
    anchor2 = GRanges(seqnames = bedpe$chr2,
                      ranges = IRanges(start = bedpe$y1, end = bedpe$y2))
  )

  gi <- gi[seqnames(anchors(gi, "first")) == seqnames(anchors(gi, "second"))]

  gi <- assignToBins(gi, binSize = resolution, pos1 = "center", pos2 = "center")

  buffer_bins <- floor(BUFFER_BP / resolution)
  gi <- pixelsToMatrices(gi, buffer = buffer_bins)

  cat(sprintf("  Loaded %d pairs from %s (binned to %dkb)\n",
              length(gi), basename(bedpe_file), resolution / 1000))
  return(gi)
}

extract_and_aggregate <- function(gi, hic_files, resolution, norm,
                                   output_dir, pair_label) {
  n_pairs <- length(gi)
  if (n_pairs < MIN_PAIRS) {
    cat(sprintf("  Skipping %s: only %d pairs (min %d)\n",
                pair_label, n_pairs, MIN_PAIRS))
    return(NULL)
  }

  buffer_bins <- floor(BUFFER_BP / resolution)
  matrix_size <- 2 * buffer_bins + 1

  cat(sprintf("  Extracting: %d pairs x %d samples (%dx%d matrices)\n",
              n_pairs, length(hic_files), matrix_size, matrix_size))

  hdf5_dir <- file.path(output_dir, "temp_hdf5")
  dir.create(hdf5_dir, recursive = TRUE, showWarnings = FALSE)
  hdf5_file <- file.path(hdf5_dir, sprintf("%s_matrices.h5", pair_label))
  on.exit(unlink(hdf5_dir, recursive = TRUE), add = TRUE)

  pixels <- tryCatch({
    pullHicMatrices(
      x = gi,
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
  }, error = function(e) {
    cat(sprintf("  Extraction failed: %s\n", e$message))
    return(NULL)
  })

  if (is.null(pixels)) return(NULL)

  count_array <- counts(pixels)
  dims <- dim(count_array)
  n_bins <- dims[1]
  n_loops <- dims[3]
  n_samples <- dims[4]

  cat(sprintf("  Aggregating %d pairs across %d samples...\n", n_loops, n_samples))

  agg_matrices <- array(NA, dim = c(n_bins, n_bins, n_samples))
  for (j in seq_len(n_samples)) {
    sample_data <- as.array(count_array[, , , j, drop = FALSE])
    agg_matrices[, , j] <- apply(sample_data, c(1, 2), mean, na.rm = TRUE)
  }

  center_idx <- ceiling(n_bins / 2)
  corner_indices <- list(c(1, 1), c(1, n_bins), c(n_bins, 1), c(n_bins, n_bins))

  enrichment_df <- data.frame(
    pair_class = pair_label,
    sample = names(hic_files),
    group = as.character(GROUPS),
    center_value = NA_real_,
    background_mean = NA_real_,
    enrichment = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(n_samples)) {
    mat <- agg_matrices[, , j]
    center_val <- mat[center_idx, center_idx]
    corner_vals <- sapply(corner_indices, function(idx) mat[idx[1], idx[2]])
    bg <- mean(corner_vals, na.rm = TRUE)
    enrichment_df$center_value[j] <- center_val
    enrichment_df$background_mean[j] <- bg
    enrichment_df$enrichment[j] <- if (!is.na(bg) && bg > 0) center_val / bg else NA_real_
  }

  per_loop_enrichment <- data.frame(
    pair_class = rep(pair_label, n_loops * n_samples),
    loop_id = rep(seq_len(n_loops), each = n_samples),
    sample = rep(names(hic_files), times = n_loops),
    group = rep(as.character(GROUPS), times = n_loops),
    center_value = NA_real_,
    background_mean = NA_real_,
    enrichment = NA_real_,
    stringsAsFactors = FALSE
  )

  row_idx <- 1
  for (i in seq_len(n_loops)) {
    for (j in seq_len(n_samples)) {
      mat <- as.matrix(count_array[, , i, j])
      cv <- mat[center_idx, center_idx]
      bg <- mean(sapply(corner_indices, function(idx) mat[idx[1], idx[2]]), na.rm = TRUE)
      per_loop_enrichment$center_value[row_idx] <- cv
      per_loop_enrichment$background_mean[row_idx] <- bg
      per_loop_enrichment$enrichment[row_idx] <- if (!is.na(bg) && bg > 0) cv / bg else NA_real_
      row_idx <- row_idx + 1
    }
  }

  cat(sprintf("  Mean aggregate P2LL — ctrl: %.2f, mut: %.2f\n",
              mean(enrichment_df$enrichment[enrichment_df$group == "ctrl"], na.rm = TRUE),
              mean(enrichment_df$enrichment[enrichment_df$group == "mut"], na.rm = TRUE)))

  return(list(
    agg_matrices = agg_matrices,
    enrichment_agg = enrichment_df,
    enrichment_per_loop = per_loop_enrichment,
    n_pairs = n_pairs,
    resolution = resolution
  ))
}

# =============================================================================
# 4. PLOTTING FUNCTIONS
# =============================================================================

plot_apa_heatmap <- function(agg_matrix, condition_name, sample_indices, resolution) {
  mean_matrix <- apply(agg_matrix[, , sample_indices, drop = FALSE],
                       c(1, 2), mean, na.rm = TRUE)

  n_bins <- nrow(mean_matrix)
  buffer_kb <- (n_bins - 1) / 2 * resolution / 1000
  center_idx <- ceiling(n_bins / 2)

  df <- expand.grid(x = 1:n_bins, y = 1:n_bins)
  df$value <- as.vector(mean_matrix)

  axis_breaks <- c(1, center_idx, n_bins)
  axis_labels <- c(sprintf("-%dkb", buffer_kb), "center", sprintf("+%dkb", buffer_kb))

  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis(name = "Contacts", option = "magma") +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    coord_fixed() +
    labs(title = condition_name, x = "SE anchor", y = "Gene TSS") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid = element_blank())
}

plot_diff_heatmap <- function(agg_matrix, resolution) {
  ctrl_mean <- apply(agg_matrix[, , CTRL_INDICES, drop = FALSE],
                     c(1, 2), mean, na.rm = TRUE)
  mut_mean <- apply(agg_matrix[, , MUT_INDICES, drop = FALSE],
                    c(1, 2), mean, na.rm = TRUE)
  diff_matrix <- mut_mean - ctrl_mean

  n_bins <- nrow(diff_matrix)
  buffer_kb <- (n_bins - 1) / 2 * resolution / 1000
  center_idx <- ceiling(n_bins / 2)

  df <- expand.grid(x = 1:n_bins, y = 1:n_bins)
  df$value <- as.vector(diff_matrix)

  axis_breaks <- c(1, center_idx, n_bins)
  axis_labels <- c(sprintf("-%dkb", buffer_kb), "center", sprintf("+%dkb", buffer_kb))

  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(name = "Diff", low = scales::muted("blue"),
                         mid = "white", high = scales::muted("red"),
                         midpoint = 0, oob = scales::squish) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels, expand = c(0, 0)) +
    coord_fixed() +
    labs(title = "BAP1-KO - Control", x = "SE anchor", y = "Gene TSS") +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid = element_blank())
}

plot_enrichment_comparison <- function(per_loop_df, pair_classes, resolution) {
  df <- per_loop_df %>%
    dplyr::filter(!is.na(enrichment), is.finite(enrichment))

  if (nrow(df) == 0) return(NULL)

  p <- ggplot(df, aes(x = group, y = enrichment, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.1, size = 0.3) +
    facet_wrap(~pair_class, scales = "free_y") +
    scale_fill_manual(values = c("ctrl" = "grey60", "mut" = "firebrick")) +
    scale_x_discrete(labels = c("ctrl" = "Control", "mut" = "BAP1-KO")) +
    labs(
      title = sprintf("SE-Gene Contact Enrichment (%dkb)", resolution / 1000),
      x = "", y = "Pixel Enrichment (P2LL)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(face = "bold"))

  return(p)
}

plot_deg_vs_invariant_enrichment <- function(per_loop_all) {
  sample_means <- per_loop_all %>%
    dplyr::filter(!is.na(enrichment), is.finite(enrichment)) %>%
    dplyr::group_by(pair_class, group) %>%
    dplyr::summarise(
      mean_enrichment = mean(enrichment, na.rm = TRUE),
      se = sd(enrichment, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  if (nrow(sample_means) == 0) return(NULL)

  p <- ggplot(sample_means, aes(x = pair_class, y = mean_enrichment, fill = group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean_enrichment - se, ymax = mean_enrichment + se),
                  position = position_dodge(width = 0.7), width = 0.2) +
    scale_fill_manual(values = c("ctrl" = "grey60", "mut" = "firebrick"),
                      labels = c("ctrl" = "Control", "mut" = "BAP1-KO")) +
    labs(
      title = "SE Contact Enrichment: DEGs vs Invariant Genes",
      x = "", y = "Mean P2LL Enrichment", fill = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top")

  return(p)
}

# =============================================================================
# 5. MAIN ANALYSIS
# =============================================================================

run_timepoint <- function(timepoint, resolution) {
  cat(sprintf("\n========== APA SE Analysis: %s at %dkb ==========\n\n",
              toupper(timepoint), resolution / 1000))

  output_dir <- file.path(APA_INPUT_DIR, "apa_results", timepoint,
                           sprintf("res_%dkb", resolution / 1000))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  hic_files <- build_hic_paths(timepoint)
  norm <- check_hic_files(hic_files)

  all_per_loop <- list()
  all_agg <- list()
  results_list <- list()

  for (pc_name in names(PAIR_CLASSES)) {
    pc <- PAIR_CLASSES[[pc_name]]
    bedpe_file <- file.path(APA_INPUT_DIR, "apa_inputs",
                             sprintf(pc$file_pattern, timepoint))

    if (!file.exists(bedpe_file)) {
      cat(sprintf("  BEDPE not found: %s — skipping\n", basename(bedpe_file)))
      next
    }

    cat(sprintf("\n--- %s: %s ---\n", pc$label, pc$description))

    gi <- load_bedpe_as_gi(bedpe_file, resolution)
    result <- extract_and_aggregate(gi, hic_files, resolution, norm,
                                     output_dir, pc$label)

    if (is.null(result)) next

    results_list[[pc_name]] <- result
    all_per_loop[[pc_name]] <- result$enrichment_per_loop
    all_agg[[pc_name]] <- result$enrichment_agg

    pair_dir <- file.path(output_dir, pc_name)
    dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)

    p_ctrl <- plot_apa_heatmap(result$agg_matrices, sprintf("%s: Control", pc$label),
                                CTRL_INDICES, resolution)
    save_multiformat_ggplot(p_ctrl, file.path(pair_dir, "apa_heatmap_ctrl"),
                            width = 6, height = 5)

    p_mut <- plot_apa_heatmap(result$agg_matrices, sprintf("%s: BAP1-KO", pc$label),
                               MUT_INDICES, resolution)
    save_multiformat_ggplot(p_mut, file.path(pair_dir, "apa_heatmap_mut"),
                            width = 6, height = 5)

    p_diff <- plot_diff_heatmap(result$agg_matrices, resolution)
    save_multiformat_ggplot(p_diff, file.path(pair_dir, "apa_heatmap_diff"),
                            width = 6, height = 5)

    write.table(result$enrichment_agg,
                file.path(pair_dir, "aggregate_enrichment.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(result$enrichment_per_loop,
                file.path(pair_dir, "per_loop_enrichment.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    cat(sprintf("  Saved outputs to %s/\n", pair_dir))
  }

  if (length(all_per_loop) == 0) {
    cat("\n  No pair classes produced results — skipping comparison plots\n")
    return(invisible(NULL))
  }

  # --- Combined comparison plots ---
  cat("\n--- Generating comparison plots ---\n")

  per_loop_combined <- bind_rows(all_per_loop)

  p_comparison <- plot_enrichment_comparison(per_loop_combined, names(results_list), resolution)
  if (!is.null(p_comparison)) {
    save_multiformat_ggplot(p_comparison, file.path(output_dir, "enrichment_comparison"),
                            width = 10, height = 6)
  }

  p_deg_vs_inv <- plot_deg_vs_invariant_enrichment(per_loop_combined)
  if (!is.null(p_deg_vs_inv)) {
    save_multiformat_ggplot(p_deg_vs_inv, file.path(output_dir, "deg_vs_invariant_enrichment"),
                            width = 7, height = 6)
  }

  # --- Statistical tests ---
  cat("\n--- Statistical tests ---\n")
  stats_lines <- c(
    sprintf("APA SE-DEG Contact Frequency Analysis"),
    sprintf("Timepoint: %s", timepoint),
    sprintf("Resolution: %dkb", resolution / 1000),
    sprintf("Date: %s\n", Sys.time())
  )

  for (pc_name in names(results_list)) {
    r <- results_list[[pc_name]]
    pc <- PAIR_CLASSES[[pc_name]]

    ctrl_enrich <- r$enrichment_per_loop %>%
      dplyr::filter(group == "ctrl", !is.na(enrichment), is.finite(enrichment)) %>%
      dplyr::pull(enrichment)
    mut_enrich <- r$enrichment_per_loop %>%
      dplyr::filter(group == "mut", !is.na(enrichment), is.finite(enrichment)) %>%
      dplyr::pull(enrichment)

    if (length(ctrl_enrich) >= 3 && length(mut_enrich) >= 3) {
      wt <- wilcox.test(ctrl_enrich, mut_enrich)
      cliff_d <- (sum(outer(mut_enrich, ctrl_enrich, `-`) > 0) -
                   sum(outer(mut_enrich, ctrl_enrich, `-`) < 0)) /
                  (length(mut_enrich) * length(ctrl_enrich))

      stats_lines <- c(stats_lines,
        sprintf("--- %s (%d pairs) ---", pc$label, r$n_pairs),
        sprintf("  Ctrl median P2LL: %.3f", median(ctrl_enrich)),
        sprintf("  Mut  median P2LL: %.3f", median(mut_enrich)),
        sprintf("  Wilcoxon p = %.2e", wt$p.value),
        sprintf("  Cliff's delta = %.3f\n", cliff_d)
      )
      cat(sprintf("  %s: Wilcoxon p = %.2e, Cliff's d = %.3f\n",
                  pc$label, wt$p.value, cliff_d))
    }
  }

  if (length(results_list) >= 2) {
    deg_enrich <- results_list[["deg_se"]]$enrichment_per_loop %>%
      dplyr::filter(!is.na(enrichment), is.finite(enrichment))
    inv_enrich <- results_list[["invariant_se"]]$enrichment_per_loop %>%
      dplyr::filter(!is.na(enrichment), is.finite(enrichment))

    if (nrow(deg_enrich) >= 3 && nrow(inv_enrich) >= 3) {
      wt_class <- wilcox.test(deg_enrich$enrichment, inv_enrich$enrichment)
      stats_lines <- c(stats_lines,
        "--- DEG-SE vs Invariant-SE (all samples pooled) ---",
        sprintf("  DEG-SE median P2LL: %.3f", median(deg_enrich$enrichment)),
        sprintf("  Inv-SE median P2LL: %.3f", median(inv_enrich$enrichment)),
        sprintf("  Wilcoxon p = %.2e\n", wt_class$p.value)
      )
    }
  }

  writeLines(stats_lines, file.path(output_dir, "apa_statistics.txt"))
  cat(sprintf("\n  Wrote statistics to %s\n", file.path(output_dir, "apa_statistics.txt")))

  cat(sprintf("\n========== Done: %s ==========\n", toupper(timepoint)))
}

# =============================================================================
# 6. DISPATCH
# =============================================================================

if (TIMEPOINT == "both") {
  for (tp in c("late", "early")) {
    run_timepoint(tp, TARGET_RESOLUTION)
  }
} else {
  run_timepoint(TIMEPOINT, TARGET_RESOLUTION)
}
