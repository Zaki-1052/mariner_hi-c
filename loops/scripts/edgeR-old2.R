# scripts/07_edgeR_analysis.R
# Differential chromatin loop analysis using edgeR with fixed dispersion
# Pipeline: mariner Hi-C extraction → edgeR statistical testing
# Constraint: No biological replicates (n=1 per condition)

# =============================================================================
# Setup & Package Loading
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(InteractionSet)
  library(GenomicRanges)
  library(ggplot2)
  library(parallel)
})

cat("\n=================================================================\n")
cat("edgeR Differential Loop Analysis - No Replicate Design\n")
cat("=================================================================\n\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
#cat("edgeR version:", packageVersion("edgeR"), "\n")
cat("R version:", R.version.string, "\n\n")

# =============================================================================
# Configuration Loading
# =============================================================================

load_config <- function(config_path = "config/edgeR_params.yaml") {
  cat("Loading configuration from:", config_path, "\n")
  
  if (!file.exists(config_path)) {
    stop("Configuration file not found: ", config_path)
  }
  
  config <- yaml::read_yaml(config_path)
  cat("✓ Configuration loaded successfully\n\n")
  return(config)
}

# =============================================================================
# Helper Functions - Validation & Utilities
# =============================================================================

validate_input_files <- function(config) {
  cat("Validating input files...\n")
  
  count_matrix_path <- config$paths$count_matrix
  binned_coords_path <- config$paths$binned_coords
  
  if (!file.exists(count_matrix_path)) {
    stop("Count matrix not found: ", count_matrix_path)
  }
  
  if (!file.exists(binned_coords_path)) {
    stop("Binned coordinates not found: ", binned_coords_path)
  }
  
  cat("✓ Count matrix found:", count_matrix_path, "\n")
  cat("✓ Coordinates found:", binned_coords_path, "\n")
  cat("✓ Input files validated\n\n")
}

create_output_directories <- function(config) {
  cat("Creating output directory structure...\n")
  
  base_dir <- config$paths$output_base
  subdirs <- c(
    config$paths$tables_dir,
    config$paths$plots_dir,
    config$paths$reports_dir,
    config$paths$checkpoint_dir
  )
  
  for (subdir in subdirs) {
    dir_path <- file.path(base_dir, subdir)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      cat("  Created:", dir_path, "\n")
    }
  }
  
  cat("✓ Output directories ready\n\n")
}

save_checkpoint <- function(object, name, config) {
  checkpoint_dir <- file.path(config$paths$output_base, config$paths$checkpoint_dir)
  filepath <- file.path(checkpoint_dir, paste0(name, ".rds"))
  saveRDS(object, filepath)
  cat("  Checkpoint saved:", filepath, "\n")
}

# =============================================================================
# Data Loading & Validation
# =============================================================================

load_count_matrix <- function(config) {
  cat("Loading count matrix...\n")
  
  count_matrix <- readRDS(config$paths$count_matrix)
  
  # Validate structure
  if (!is.matrix(count_matrix)) {
    stop("Count matrix is not a matrix object")
  }
  
  if (ncol(count_matrix) != 2) {
    stop("Expected 2 columns (ctrl, mut), found: ", ncol(count_matrix))
  }
  
  expected_colnames <- config$experimental_design$sample_names
  if (!all(colnames(count_matrix) == expected_colnames)) {
    stop("Column names mismatch. Expected: ", paste(expected_colnames, collapse=", "))
  }
  
  # Data integrity checks
  if (any(is.na(count_matrix))) {
    stop("NA values detected in count matrix")
  }
  
  if (any(count_matrix < 0)) {
    stop("Negative values detected in count matrix")
  }
  
  cat("✓ Count matrix loaded successfully\n")
  cat("  Dimensions:", nrow(count_matrix), "loops ×", ncol(count_matrix), "samples\n")
  cat("  Column names:", paste(colnames(count_matrix), collapse=", "), "\n")
  cat("  Value range: [", min(count_matrix), ",", max(count_matrix), "]\n")
  cat("  Total counts:", paste(colSums(count_matrix), collapse=", "), "\n\n")
  
  return(count_matrix)
}

extract_genomic_coordinates <- function(config) {
  cat("Extracting genomic coordinates from binned GInteractions...\n")
  
  gi <- readRDS(config$paths$binned_coords)
  
  if (!inherits(gi, "GInteractions")) {
    stop("Expected GInteractions object, got: ", class(gi))
  }
  
  # Extract anchor coordinates
  anchor1 <- anchors(gi, type="first")
  anchor2 <- anchors(gi, type="second")
  
  # Build annotation dataframe
  coords <- data.frame(
    loop_id = paste0("loop_", seq_len(length(gi))),
    chr1 = as.character(seqnames(anchor1)),
    start1 = start(anchor1),
    end1 = end(anchor1),
    chr2 = as.character(seqnames(anchor2)),
    start2 = start(anchor2),
    end2 = end(anchor2),
    stringsAsFactors = FALSE
  )
  
  # Create coordinate string for easier reading
  coords$coordinate <- paste0(
    coords$chr1, ":", coords$start1, "-", coords$end1,
    "_",
    coords$chr2, ":", coords$start2, "-", coords$end2
  )
  
  cat("✓ Coordinates extracted successfully\n")
  cat("  Number of loops:", nrow(coords), "\n")
  cat("  Chromosomes (anchor1):", paste(unique(coords$chr1), collapse=", "), "\n\n")
  
  return(coords)
}

# =============================================================================
# DGEList Construction
# =============================================================================

create_dgelist <- function(counts, coords, config) {
  cat("Constructing DGEList object...\n")
  
  # Verify counts and coords alignment
  if (nrow(counts) != nrow(coords)) {
    stop("Row count mismatch: counts (", nrow(counts), ") vs coords (", nrow(coords), ")")
  }
  
  # Set rownames to loop_ids
  rownames(counts) <- coords$loop_id
  
  # Create group factor
  group <- factor(config$experimental_design$group_levels, 
                  levels = config$experimental_design$group_levels)
  
  # Create DGEList with genomic annotations
  y <- DGEList(
    counts = counts,
    group = group,
    genes = coords
  )
  
  cat("✓ DGEList created successfully\n")
  cat("  Samples:", paste(colnames(y), collapse=", "), "\n")
  cat("  Groups:", paste(y$samples$group, collapse=", "), "\n")
  cat("  Library sizes:", paste(y$samples$lib.size, collapse=", "), "\n")
  cat("  Library size ratio (ctrl/mut):", 
      round(y$samples$lib.size[1] / y$samples$lib.size[2], 3), "\n\n")
  
  save_checkpoint(y, "01_raw_dgelist", config)
  
  return(y)
}

# =============================================================================
# Filtering
# =============================================================================

filter_low_count_loops <- function(y, config) {
  cat("Filtering low-count loops...\n")
  
  n_before <- nrow(y)
  
  # Use edgeR's filterByExpr
  keep <- filterByExpr(
    y, 
    group = y$samples$group,
    min.count = config$filtering$min_count,
    min.total.count = config$filtering$min_total_count
  )
  
  n_keep <- sum(keep)
  filter_rate <- 1 - (n_keep / n_before)
  
  cat("  Loops before filtering:", n_before, "\n")
  cat("  Loops passing filter:", n_keep, "\n")
  cat("  Filter rate:", round(filter_rate * 100, 2), "%\n")
  
  # Check if filter rate is too high
  if (filter_rate > config$filtering$max_filter_rate) {
    warning("Filter rate exceeds threshold (", 
            round(filter_rate * 100, 1), "% > ",
            round(config$filtering$max_filter_rate * 100, 1), "%)")
    
    if (config$filtering$manual_fallback$enabled) {
      cat("  Applying fallback manual filtering...\n")
      cpm_vals <- cpm(y)
      keep_manual <- rowSums(cpm_vals >= config$filtering$manual_fallback$min_cpm) >= 
                     config$filtering$manual_fallback$min_samples
      keep <- keep_manual
      n_keep <- sum(keep)
      cat("  Loops after fallback filter:", n_keep, "\n")
    }
  }
  
  # Check minimum loops retained
  if (n_keep < config$filtering$min_loops_retained) {
    stop("Too few loops retained (", n_keep, " < ", 
         config$filtering$min_loops_retained, ")")
  }
  
  # Apply filter and recalculate library sizes
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  cat("✓ Filtering complete\n")
  cat("  Final loop count:", nrow(y), "\n\n")
  
  save_checkpoint(y, "02_filtered_dgelist", config)
  
  return(y)
}

# =============================================================================
# Normalization
# =============================================================================

normalize_libraries <- function(y, config) {
  cat("Calculating normalization factors (TMM)...\n")
  
  y <- calcNormFactors(
    y, 
    method = config$normalization$method,
    logratioTrim = config$normalization$logratioTrim,
    sumTrim = config$normalization$sumTrim,
    doWeighting = config$normalization$doWeighting
  )
  
  norm_factors <- y$samples$norm.factors
  eff_lib_sizes <- y$samples$lib.size * norm_factors
  
  cat("✓ Normalization complete\n")
  cat("  Normalization factors:", paste(round(norm_factors, 4), collapse=", "), "\n")
  cat("  Effective library sizes:", paste(round(eff_lib_sizes), collapse=", "), "\n")
  
  # Validate norm factors are in reasonable range
  expected_range <- config$normalization$norm_factor_range
  if (any(norm_factors < expected_range[1] | norm_factors > expected_range[2])) {
    warning("Normalization factors outside expected range [", 
            expected_range[1], ", ", expected_range[2], "]: ",
            paste(round(norm_factors, 4), collapse=", "))
  } else {
    cat("  Norm factors within expected range ✓\n")
  }
  
  cat("\n")
  
  save_checkpoint(y, "03_normalized_dgelist", config)
  
  return(y)
}

# =============================================================================
# Primary Differential Analysis (Fixed BCV)
# =============================================================================

run_primary_analysis <- function(y, config) {
  cat("=================================================================\n")
  cat("PRIMARY ANALYSIS: BCV = ", config$dispersion$primary_bcv, "\n")
  cat("=================================================================\n\n")
  
  bcv <- config$dispersion$primary_bcv
  dispersion <- bcv^2
  
  cat("Running differential analysis with fixed dispersion...\n")
  cat("  BCV:", bcv, "\n")
  cat("  Dispersion:", dispersion, "\n")
  cat("  Justification:", config$dispersion$justification, "\n\n")
  
  # Create design matrix
  design <- model.matrix(~ y$samples$group)
  colnames(design) <- c("Intercept", "mut_vs_ctrl")
  
  cat("Design matrix:\n")
  print(design)
  cat("\n")
  
  # Fit GLM with fixed dispersion
  cat("Fitting GLM with fixed dispersion...\n")
  fit <- glmFit(
    y, 
    design, 
    dispersion = dispersion,
    prior.count = config$testing$prior_count
  )
  cat("✓ GLM fitting complete\n\n")
  
  # Perform likelihood ratio test
  cat("Performing likelihood ratio test (coefficient 2: mut_vs_ctrl)...\n")
  lrt <- glmLRT(fit, coef = config$testing$test_coefficient)
  cat("✓ LRT complete\n\n")
  
  # Extract results
  results <- topTags(
    lrt, 
    n = Inf,
    adjust.method = config$testing$adjust_method,
    sort.by = "PValue"
  )$table
  
  # Add derived columns
  results$fold_change <- 2^results$logFC
  results$abs_logFC <- abs(results$logFC)
  
  # Categorize significance
  results$significance_category <- ifelse(
    results$FDR < config$thresholds$fdr_stringent,
    "Stringent",
    ifelse(
      results$FDR < config$thresholds$fdr_primary,
      "Significant",
      "NS"
    )
  )
  
  # Determine direction
  results$direction <- ifelse(
    results$FDR >= config$thresholds$fdr_primary,
    "unchanged",
    ifelse(
      results$logFC > 0,
      "up_in_mut",
      "down_in_mut"
    )
  )
  
  # Summary statistics
  cat("Results Summary:\n")
  cat("  Total loops tested:", nrow(results), "\n")
  cat("  Significant (FDR < ", config$thresholds$fdr_primary, "):", 
      sum(results$FDR < config$thresholds$fdr_primary), "\n")
  cat("  Stringent (FDR < ", config$thresholds$fdr_stringent, "):", 
      sum(results$FDR < config$thresholds$fdr_stringent), "\n")
  cat("  Up in mutant:", sum(results$direction == "up_in_mut"), "\n")
  cat("  Down in mutant:", sum(results$direction == "down_in_mut"), "\n")
  cat("  |logFC| > 1:", sum(results$abs_logFC > 1 & results$FDR < config$thresholds$fdr_primary), "\n")
  cat("  |logFC| > 2:", sum(results$abs_logFC > 2 & results$FDR < config$thresholds$fdr_primary), "\n\n")
  
  # Save results
  save_checkpoint(results, paste0("04_primary_results_BCV", bcv), config)
  
  result_obj <- list(
    fit = fit,
    lrt = lrt,
    results = results,
    bcv = bcv
  )
  
  return(result_obj)
}

# =============================================================================
# Sensitivity Analysis
# =============================================================================

run_sensitivity_analysis <- function(y, config) {
  cat("=================================================================\n")
  cat("SENSITIVITY ANALYSIS: Testing Multiple BCV Values\n")
  cat("=================================================================\n\n")
  
  bcv_values <- config$dispersion$sensitivity_bcvs
  labels <- config$dispersion$sensitivity_labels
  
  cat("Testing BCV values:", paste(bcv_values, collapse=", "), "\n\n")
  
  # Create design matrix (same for all)
  design <- model.matrix(~ y$samples$group)
  colnames(design) <- c("Intercept", "mut_vs_ctrl")
  
  # Function to run analysis for one BCV
  analyze_bcv <- function(bcv, label) {
    cat("  BCV =", bcv, "(", label, ")\n")
    
    dispersion <- bcv^2
    
    # Fit and test
    fit <- glmFit(y, design, dispersion = dispersion, 
                  prior.count = config$testing$prior_count)
    lrt <- glmLRT(fit, coef = config$testing$test_coefficient)
    
    results <- topTags(lrt, n = Inf, 
                       adjust.method = config$testing$adjust_method,
                       sort.by = "none")$table  # Keep original order
    
    return(list(
      bcv = bcv,
      label = label,
      n_sig_005 = sum(results$FDR < 0.05),
      n_sig_001 = sum(results$FDR < 0.01),
      sig_loops = rownames(results)[results$FDR < 0.05]
    ))
  }
  
  # Run analyses (parallel if enabled)
  if (config$computation$use_parallel && config$computation$n_cores > 1) {
    cat("Running in parallel with", config$computation$n_cores, "cores...\n")
    sensitivity_results <- mclapply(
      seq_along(bcv_values),
      function(i) analyze_bcv(bcv_values[i], labels[i]),
      mc.cores = config$computation$n_cores
    )
  } else {
    sensitivity_results <- lapply(
      seq_along(bcv_values),
      function(i) analyze_bcv(bcv_values[i], labels[i])
    )
  }
  
  # Summarize results
  cat("\n")
  cat("Sensitivity Analysis Summary:\n")
  cat(sprintf("%-12s %-8s %-15s %-15s\n", "BCV", "Label", "Sig (FDR<0.05)", "Sig (FDR<0.01)"))
  cat(strrep("-", 55), "\n")
  
  for (res in sensitivity_results) {
    cat(sprintf("%-12.2f %-8s %-15d %-15d\n", 
                res$bcv, res$label, res$n_sig_005, res$n_sig_001))
  }
  
  # Calculate overlap
  cat("\nOverlap analysis (FDR < 0.05):\n")
  all_sig_loops <- lapply(sensitivity_results, function(x) x$sig_loops)
  core_loops <- Reduce(intersect, all_sig_loops)
  cat("  Core significant loops (all BCVs agree):", length(core_loops), "\n")
  
  # Any significant in at least one
  any_sig_loops <- Reduce(union, all_sig_loops)
  cat("  Any significant (at least one BCV):", length(any_sig_loops), "\n\n")
  
  save_checkpoint(sensitivity_results, "05_sensitivity_analysis", config)
  
  return(sensitivity_results)
}

# =============================================================================
# Output Generation
# =============================================================================

write_results_tables <- function(primary_results, config) {
  cat("Writing results tables...\n")
  
  results <- primary_results$results
  bcv <- primary_results$bcv
  
  # Full results table
  output_prefix <- config$output_naming$prefix
  tables_dir <- file.path(config$paths$output_base, config$paths$tables_dir)
  
  full_table_path <- file.path(tables_dir, 
                                sprintf("%s_results_full_BCV%.1f.tsv", output_prefix, bcv))
  write.table(results, full_table_path, 
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Full results:", full_table_path, "\n")
  
  # Significant subset (FDR < 0.05)
  sig_results <- results[results$FDR < config$thresholds$fdr_primary, ]
  if (nrow(sig_results) > 0) {
    sig_table_path <- file.path(tables_dir, 
                                 sprintf("%s_significant_FDR0.05_BCV%.1f.tsv", 
                                         output_prefix, bcv))
    write.table(sig_results, sig_table_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Significant (FDR<0.05):", sig_table_path, "\n")
  }
  
  # Stringent subset (FDR < 0.01)
  stringent_results <- results[results$FDR < config$thresholds$fdr_stringent, ]
  if (nrow(stringent_results) > 0) {
    stringent_table_path <- file.path(tables_dir, 
                                       sprintf("%s_stringent_FDR0.01_BCV%.1f.tsv", 
                                               output_prefix, bcv))
    write.table(stringent_results, stringent_table_path, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Stringent (FDR<0.01):", stringent_table_path, "\n")
  }
  
  cat("✓ Results tables written\n\n")
}

# =============================================================================
# Visualization
# =============================================================================

plot_ma <- function(primary_results, config) {
  cat("Generating MA plot...\n")
  
  results <- primary_results$results
  bcv <- primary_results$bcv
  
  # Prepare data
  plot_data <- data.frame(
    baseMean = results$logCPM,
    logFC = results$logFC,
    significant = results$FDR < config$thresholds$fdr_primary
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = baseMean, y = logFC)) +
    geom_point(aes(color = significant), 
               alpha = config$visualization$ma_plot$alpha,
               size = config$visualization$ma_plot$point_size) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(
      values = c("FALSE" = config$visualization$ma_plot$nonsig_color,
                 "TRUE" = config$visualization$ma_plot$sig_color),
      labels = c("NS", "FDR < 0.05")
    ) +
    labs(
      title = sprintf("MA Plot (BCV = %.1f)", bcv),
      x = "Average log-CPM",
      y = "log2 Fold Change (mut vs ctrl)",
      color = "Significance"
    ) +
    theme_bw(base_size = config$visualization$font_size) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Save plot
  plots_dir <- file.path(config$paths$output_base, config$paths$plots_dir)
  plot_path <- file.path(plots_dir, 
                          sprintf("%s_MA_plot_BCV%.1f.pdf", 
                                  config$output_naming$prefix, bcv))
  
  ggsave(plot_path, p, 
         width = config$visualization$figure_width,
         height = config$visualization$figure_height,
         dpi = config$visualization$dpi)
  
  cat("  Saved:", plot_path, "\n\n")
  
  return(p)
}

plot_volcano <- function(primary_results, config) {
  cat("Generating volcano plot...\n")
  
  results <- primary_results$results
  bcv <- primary_results$bcv
  
  # Prepare data with categories
  plot_data <- data.frame(
    logFC = results$logFC,
    negLog10FDR = -log10(results$FDR),
    category = ifelse(
      results$FDR >= config$thresholds$fdr_primary,
      "NS",
      ifelse(results$logFC > 0, "Up", "Down")
    )
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(color = category), 
               alpha = config$visualization$volcano_plot$alpha,
               size = config$visualization$volcano_plot$point_size) +
    geom_vline(xintercept = c(-config$visualization$volcano_plot$logfc_threshold,
                               config$visualization$volcano_plot$logfc_threshold),
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(config$visualization$volcano_plot$fdr_threshold),
               linetype = "dashed", color = "gray50") +
    scale_color_manual(
      values = c(
        "Up" = config$visualization$volcano_plot$colors$up,
        "Down" = config$visualization$volcano_plot$colors$down,
        "NS" = config$visualization$volcano_plot$colors$ns
      ),
      labels = c("Down in mut", "NS", "Up in mut")
    ) +
    labs(
      title = sprintf("Volcano Plot (BCV = %.1f)", bcv),
      x = "log2 Fold Change (mut vs ctrl)",
      y = "-log10(FDR)",
      color = "Direction"
    ) +
    theme_bw(base_size = config$visualization$font_size) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Save plot
  plots_dir <- file.path(config$paths$output_base, config$paths$plots_dir)
  plot_path <- file.path(plots_dir, 
                          sprintf("%s_volcano_plot_BCV%.1f.pdf", 
                                  config$output_naming$prefix, bcv))
  
  ggsave(plot_path, p, 
         width = config$visualization$figure_width,
         height = config$visualization$figure_height,
         dpi = config$visualization$dpi)
  
  cat("  Saved:", plot_path, "\n\n")
  
  return(p)
}

plot_sensitivity <- function(sensitivity_results, config) {
  cat("Generating sensitivity analysis plot...\n")
  
  # Extract data
  plot_data <- data.frame(
    bcv = sapply(sensitivity_results, function(x) x$bcv),
    label = sapply(sensitivity_results, function(x) x$label),
    n_sig_005 = sapply(sensitivity_results, function(x) x$n_sig_005),
    n_sig_001 = sapply(sensitivity_results, function(x) x$n_sig_001)
  )
  
  # Reshape for plotting
  plot_data_long <- data.frame(
    bcv = rep(plot_data$bcv, 2),
    label = rep(plot_data$label, 2),
    threshold = rep(c("FDR < 0.05", "FDR < 0.01"), each = nrow(plot_data)),
    n_sig = c(plot_data$n_sig_005, plot_data$n_sig_001)
  )
  
  # Create plot
  p <- ggplot(plot_data_long, aes(x = bcv, y = n_sig, color = threshold, group = threshold)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = plot_data$bcv, 
                       labels = paste0(plot_data$bcv, "\n(", plot_data$label, ")")) +
    scale_color_manual(values = c("FDR < 0.05" = "#E74C3C", "FDR < 0.01" = "#3498DB")) +
    labs(
      title = "Sensitivity Analysis: Effect of BCV on Differential Calls",
      x = "Biological Coefficient of Variation (BCV)",
      y = "Number of Significant Loops",
      color = "Threshold"
    ) +
    theme_bw(base_size = config$visualization$font_size) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Save plot
  plots_dir <- file.path(config$paths$output_base, config$paths$plots_dir)
  plot_path <- file.path(plots_dir, 
                          sprintf("%s_sensitivity_analysis.pdf", 
                                  config$output_naming$prefix))
  
  ggsave(plot_path, p, 
         width = config$visualization$figure_width,
         height = config$visualization$figure_height,
         dpi = config$visualization$dpi)
  
  cat("  Saved:", plot_path, "\n\n")
  
  return(p)
}

# =============================================================================
# Summary Report
# =============================================================================

generate_summary_report <- function(primary_results, sensitivity_results, config) {
  cat("=================================================================\n")
  cat("ANALYSIS SUMMARY\n")
  cat("=================================================================\n\n")
  
  results <- primary_results$results
  bcv <- primary_results$bcv
  
  cat("Analysis Parameters:\n")
  cat("  Primary BCV:", bcv, "\n")
  cat("  Dispersion:", bcv^2, "\n")
  cat("  Design:", config$testing$design_formula, "\n")
  cat("  Test type:", config$testing$test_type, "\n")
  cat("  Multiple testing:", config$testing$adjust_method, "\n\n")
  
  cat("Dataset:\n")
  cat("  Total loops analyzed:", nrow(results), "\n")
  cat("  Samples:", paste(config$experimental_design$sample_names, collapse=" vs "), "\n\n")
  
  cat("Primary Results (BCV = ", bcv, "):\n", sep="")
  cat("  Significant (FDR < 0.05):", sum(results$FDR < 0.05), "\n")
  cat("  Stringent (FDR < 0.01):", sum(results$FDR < 0.01), "\n")
  cat("  Up-regulated in mutant:", sum(results$direction == "up_in_mut"), "\n")
  cat("  Down-regulated in mutant:", sum(results$direction == "down_in_mut"), "\n")
  cat("  Median |logFC| (significant):", 
      round(median(results$abs_logFC[results$FDR < 0.05]), 3), "\n")
  cat("  Max |logFC|:", round(max(results$abs_logFC), 3), "\n\n")
  
  cat("Effect Size Distribution (FDR < 0.05):\n")
  sig_results <- results[results$FDR < 0.05, ]
  if (nrow(sig_results) > 0) {
    cat("  |logFC| > 0.5:", sum(sig_results$abs_logFC > 0.5), "\n")
    cat("  |logFC| > 1.0:", sum(sig_results$abs_logFC > 1.0), "\n")
    cat("  |logFC| > 1.5:", sum(sig_results$abs_logFC > 1.5), "\n")
    cat("  |logFC| > 2.0:", sum(sig_results$abs_logFC > 2.0), "\n")
  }
  cat("\n")
  
  cat("Sensitivity Analysis:\n")
  for (res in sensitivity_results) {
    cat(sprintf("  BCV %.2f (%s): %d significant (FDR<0.05)\n",
                res$bcv, res$label, res$n_sig_005))
  }
  
  # Calculate core set
  all_sig_loops <- lapply(sensitivity_results, function(x) x$sig_loops)
  core_loops <- Reduce(intersect, all_sig_loops)
  cat("\n  Core differential loops (all BCVs agree):", length(core_loops), "\n")
  
  cat("\n")
  cat("Output Files:\n")
  cat("  Tables:", file.path(config$paths$output_base, config$paths$tables_dir), "\n")
  cat("  Plots:", file.path(config$paths$output_base, config$paths$plots_dir), "\n")
  cat("  Checkpoints:", file.path(config$paths$output_base, config$paths$checkpoint_dir), "\n")
  
  cat("\n=================================================================\n")
  cat("Analysis Complete: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("=================================================================\n\n")
}

# =============================================================================
# Main Execution Pipeline
# =============================================================================

main <- function() {
  set.seed(42)  # For reproducibility
  
  # Load configuration
  config <- load_config("config/edgeR_params.yaml")
  
  # Setup
  validate_input_files(config)
  create_output_directories(config)
  
  # Load data
  count_matrix <- load_count_matrix(config)
  coords <- extract_genomic_coordinates(config)
  
  # Create DGEList
  y <- create_dgelist(count_matrix, coords, config)
  
  # Filtering
  y <- filter_low_count_loops(y, config)
  
  # Normalization
  y <- normalize_libraries(y, config)
  
  # Primary analysis
  primary_results <- run_primary_analysis(y, config)
  
  # Sensitivity analysis
  sensitivity_results <- run_sensitivity_analysis(y, config)
  
  # Generate outputs
  write_results_tables(primary_results, config)
  
  # Visualizations
  plot_ma(primary_results, config)
  plot_volcano(primary_results, config)
  plot_sensitivity(sensitivity_results, config)
  
  # Summary report
  generate_summary_report(primary_results, sensitivity_results, config)
  
  cat("All analyses completed successfully!\n")
}

# Execute main pipeline
main()
