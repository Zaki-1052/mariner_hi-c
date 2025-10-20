# scripts/07_edgeR_analysis.R
# Differential chromatin loop analysis using edgeR
# Critical constraint: No biological replicates (n=1 per condition)
# Approach: Fixed dispersion (BCV=0.4) + sensitivity analysis

suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(GenomicRanges)
  library(InteractionSet)
  library(ggplot2)
  library(data.table)
  library(RColorBrewer)
})

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Load and validate YAML configuration
load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  config <- yaml::read_yaml(config_path)
  cat("✓ Configuration loaded\n")
  return(config)
}

#' Expand environment variables in paths
expand_path <- function(path) {
  # Simple ${var} expansion
  path <- gsub("\\$\\{base_dir\\}", 
               "/expanse/lustre/projects/csd940/zalibhai/mariner", path)
  path <- gsub("\\$\\{input_dir\\}", 
               "/expanse/lustre/projects/csd940/zalibhai/mariner/outputs/full", path)
  path <- gsub("\\$\\{output_base\\}", 
               "/expanse/lustre/projects/csd940/zalibhai/mariner/outputs/differential", path)
  path <- gsub("\\$\\{qc_dir\\}", 
               "/expanse/lustre/projects/csd940/zalibhai/mariner/outputs/qc_report", path)
  return(path)
}

#' Create output directories
create_output_dirs <- function(config) {
  dirs <- c(
    expand_path(config$paths$output_tables),
    expand_path(config$paths$output_plots),
    expand_path(config$paths$output_reports),
    expand_path(config$paths$logs_dir)
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

#' Extract genomic coordinates from GInteractions
extract_coordinates <- function(gi) {
  coords <- data.frame(
    chr1 = as.character(seqnames(anchors(gi, type = "first"))),
    start1 = start(anchors(gi, type = "first")),
    end1 = end(anchors(gi, type = "first")),
    chr2 = as.character(seqnames(anchors(gi, type = "second"))),
    start2 = start(anchors(gi, type = "second")),
    end2 = end(anchors(gi, type = "second")),
    stringsAsFactors = FALSE
  )
  
  coords$coordinate <- paste0(
    coords$chr1, ":", coords$start1, "-", coords$end1,
    "_",
    coords$chr2, ":", coords$start2, "-", coords$end2
  )
  
  return(coords)
}

#' Run differential analysis for single BCV value
run_differential_analysis <- function(y, design, bcv, coef = 2) {
  dispersion <- bcv^2
  
  # Fit GLM with fixed dispersion
  fit <- glmFit(y, design, dispersion = dispersion, prior.count = 0.125)
  
  # Likelihood ratio test
  lrt <- glmLRT(fit, coef = coef)
  
  # Extract results
  results <- topTags(lrt, n = Inf, sort.by = "none")$table
  
  return(list(fit = fit, lrt = lrt, results = results))
}

#' Classify loops by significance
classify_significance <- function(results, fdr_05 = 0.05, fdr_01 = 0.01, fc_threshold = 1) {
  results$significance_category <- "NS"
  results$significance_category[results$FDR < fdr_05] <- "FDR_05"
  results$significance_category[results$FDR < fdr_01] <- "FDR_01"
  
  results$direction <- "unchanged"
  sig_idx <- results$FDR < fdr_05
  results$direction[sig_idx & results$logFC > 0] <- "up_in_mut"
  results$direction[sig_idx & results$logFC < 0] <- "down_in_mut"
  
  results$high_fc <- abs(results$logFC) > fc_threshold
  
  return(results)
}

################################################################################
# MAIN ANALYSIS PIPELINE
################################################################################

main <- function() {
  
  cat("================================================================================\n")
  cat("edgeR Differential Loop Analysis\n")
  cat("Pipeline: Mariner extraction → edgeR differential analysis\n")
  cat("Critical constraint: No biological replicates (n=1 per condition)\n")
  cat("================================================================================\n\n")
  
  # Set reproducibility
  set.seed(42)
  
  #-----------------------------------------------------------------------------
  # 1. CONFIGURATION & SETUP
  #-----------------------------------------------------------------------------
  
  cat("STEP 1: Loading configuration\n")
  cat("--------------------------------------------------------------------------------\n")
  
  config_path <- "config/analysis_params.yaml"
  config <- load_config(config_path)
  create_output_dirs(config)
  
  # Extract key parameters
  bcv_primary <- config$edger$dispersion$primary_analysis$bcv
  bcv_sensitivity <- config$edger$dispersion$sensitivity_analysis$bcv_values
  fdr_threshold <- config$edger$testing$fdr_threshold
  fdr_strict <- config$edger$testing$fdr_strict
  
  cat("\nPrimary analysis parameters:\n")
  cat("  BCV:", bcv_primary, "\n")
  cat("  Dispersion:", bcv_primary^2, "\n")
  cat("  FDR threshold:", fdr_threshold, "\n")
  cat("  Sensitivity BCV range:", paste(bcv_sensitivity, collapse = ", "), "\n\n")
  
  #-----------------------------------------------------------------------------
  # 2. INPUT VALIDATION
  #-----------------------------------------------------------------------------
  
  cat("STEP 2: Loading and validating input data\n")
  cat("--------------------------------------------------------------------------------\n")
  
  counts_path <- expand_path(config$paths$counts_matrix)
  coords_path <- expand_path(config$paths$binned_interactions)
  
  cat("Loading count matrix from:", counts_path, "\n")
  counts <- readRDS(counts_path)
  
  cat("Loading coordinates from:", coords_path, "\n")
  binned_gi <- readRDS(coords_path)
  
  # Validate dimensions
  expected_loops <- config$qc$pre_checks$expected_loop_count
  if (nrow(counts) != expected_loops) {
    warning("Unexpected loop count: ", nrow(counts), " (expected ", expected_loops, ")")
  }
  
  # Validate no NAs or negatives
  if (any(is.na(counts))) {
    stop("Count matrix contains NA values")
  }
  if (any(counts < 0)) {
    stop("Count matrix contains negative values")
  }
  
  # Validate column names
  expected_cols <- config$qc$pre_checks$verify_column_names
  if (!all(colnames(counts) == expected_cols)) {
    stop("Unexpected column names. Expected: ", paste(expected_cols, collapse = ", "))
  }
  
  cat("\n✓ Input validation passed\n")
  cat("  Dimensions:", nrow(counts), "loops ×", ncol(counts), "samples\n")
  cat("  Sample names:", paste(colnames(counts), collapse = ", "), "\n")
  cat("  Count range: [", min(counts), ", ", max(counts), "]\n\n")
  
  #-----------------------------------------------------------------------------
  # 3. COORDINATE RECOVERY
  #-----------------------------------------------------------------------------
  
  cat("STEP 3: Extracting genomic coordinates\n")
  cat("--------------------------------------------------------------------------------\n")
  
  loop_coords <- extract_coordinates(binned_gi)
  loop_coords$loop_id <- paste0("loop_", seq_len(nrow(loop_coords)))
  rownames(counts) <- loop_coords$loop_id
  
  cat("✓ Coordinates extracted for", nrow(loop_coords), "loops\n")
  cat("  Example coordinate:", loop_coords$coordinate[1], "\n\n")
  
  #-----------------------------------------------------------------------------
  # 4. DGELIST CONSTRUCTION
  #-----------------------------------------------------------------------------
  
  cat("STEP 4: Creating DGEList object\n")
  cat("--------------------------------------------------------------------------------\n")
  
  # Create group factor
  group <- factor(colnames(counts), levels = c("ctrl", "mut"))
  
  # Build DGEList with gene annotation
  y <- DGEList(
    counts = counts,
    group = group,
    genes = loop_coords[, c("loop_id", "chr1", "start1", "end1", 
                             "chr2", "start2", "end2", "coordinate")]
  )
  
  cat("✓ DGEList created\n")
  cat("  Library sizes:\n")
  print(y$samples)
  cat("\n")
  
  # Save checkpoint
  checkpoint_file <- file.path(expand_path(config$paths$output_tables), 
                               "07_checkpoint_initial_dge.rds")
  saveRDS(y, checkpoint_file)
  cat("✓ Checkpoint saved:", checkpoint_file, "\n\n")
  
  #-----------------------------------------------------------------------------
  # 5. FILTERING PIPELINE
  #-----------------------------------------------------------------------------
  
  cat("STEP 5: Filtering low-count loops\n")
  cat("--------------------------------------------------------------------------------\n")
  
  # Calculate pre-filter stats
  n_loops_initial <- nrow(y)
  
  # Apply filterByExpr
  keep <- filterByExpr(y, group = group)
  n_kept <- sum(keep)
  filter_rate <- 1 - (n_kept / n_loops_initial)
  
  cat("Filtering results:\n")
  cat("  Initial loops:", n_loops_initial, "\n")
  cat("  Loops retained:", n_kept, "\n")
  cat("  Filter rate:", sprintf("%.1f%%", filter_rate * 100), "\n")
  
  # Check if filter rate is too high
  max_filter_rate <- config$qc$post_filter_checks$max_filter_rate
  if (filter_rate > max_filter_rate) {
    warning("High filter rate (", sprintf("%.1f%%", filter_rate * 100), 
            ") exceeds threshold (", max_filter_rate * 100, "%)")
  }
  
  # Check minimum loops retained
  min_loops <- config$qc$post_filter_checks$min_loops_retained
  if (n_kept < min_loops) {
    warning("Few loops retained (", n_kept, ") below recommended minimum (", min_loops, ")")
  }
  
  # Apply filter
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  cat("\n✓ Filtering complete\n")
  cat("  Final loop count:", nrow(y), "\n")
  cat("  Updated library sizes:", paste(y$samples$lib.size, collapse = ", "), "\n\n")
  
  # Save filtered checkpoint
  checkpoint_file <- file.path(expand_path(config$paths$output_tables), 
                               "07_checkpoint_filtered_dge.rds")
  saveRDS(y, checkpoint_file)
  cat("✓ Checkpoint saved:", checkpoint_file, "\n\n")
  
  #-----------------------------------------------------------------------------
  # 6. NORMALIZATION
  #-----------------------------------------------------------------------------
  
  cat("STEP 6: TMM normalization\n")
  cat("--------------------------------------------------------------------------------\n")
  
  y <- calcNormFactors(y, method = "TMM")
  
  cat("Normalization factors:\n")
  print(y$samples)
  
  # Validate norm.factors range
  valid_range <- config$edger$normalization$validate_range
  if (any(y$samples$norm.factors < valid_range[1]) || 
      any(y$samples$norm.factors > valid_range[2])) {
    warning("Normalization factors outside expected range [", 
            paste(valid_range, collapse = ", "), "]")
  }
  
  cat("\n✓ Normalization complete\n")
  cat("  Effective library sizes:", 
      paste(round(y$samples$lib.size * y$samples$norm.factors), collapse = ", "), "\n\n")
  
  # Save normalized checkpoint
  checkpoint_file <- file.path(expand_path(config$paths$output_tables), 
                               "07_checkpoint_normalized_dge.rds")
  saveRDS(y, checkpoint_file)
  cat("✓ Checkpoint saved:", checkpoint_file, "\n\n")
  
  #-----------------------------------------------------------------------------
  # 7. PRIMARY DIFFERENTIAL ANALYSIS
  #-----------------------------------------------------------------------------
  
  cat("STEP 7: Primary differential analysis (BCV =", bcv_primary, ")\n")
  cat("--------------------------------------------------------------------------------\n")
  
  # Create design matrix
  design <- model.matrix(~group)
  colnames(design) <- c("Intercept", "mut_vs_ctrl")
  
  cat("Design matrix:\n")
  print(design)
  cat("\n")
  
  # Run analysis
  cat("Running GLM + LRT with fixed dispersion =", bcv_primary^2, "\n")
  analysis_primary <- run_differential_analysis(
    y = y, 
    design = design, 
    bcv = bcv_primary, 
    coef = 2
  )
  
  # Extract and annotate results
  results_primary <- analysis_primary$results
  results_primary <- cbind(y$genes, results_primary)
  results_primary <- classify_significance(
    results_primary,
    fdr_05 = fdr_threshold,
    fdr_01 = fdr_strict,
    fc_threshold = config$edger$testing$logfc_threshold
  )
  
  # Calculate summary statistics
  n_sig_05 <- sum(results_primary$FDR < fdr_threshold)
  n_sig_01 <- sum(results_primary$FDR < fdr_strict)
  n_up <- sum(results_primary$direction == "up_in_mut")
  n_down <- sum(results_primary$direction == "down_in_mut")
  n_high_fc <- sum(results_primary$high_fc & results_primary$FDR < fdr_threshold)
  
  cat("\n✓ Primary analysis complete\n")
  cat("\nSummary statistics:\n")
  cat("  Total loops analyzed:", nrow(results_primary), "\n")
  cat("  Significant (FDR < 0.05):", n_sig_05, 
      sprintf("(%.1f%%)", 100 * n_sig_05 / nrow(results_primary)), "\n")
  cat("  Significant (FDR < 0.01):", n_sig_01, 
      sprintf("(%.1f%%)", 100 * n_sig_01 / nrow(results_primary)), "\n")
  cat("  Up in mutant:", n_up, "\n")
  cat("  Down in mutant:", n_down, "\n")
  cat("  High fold-change (|logFC| > 1):", n_high_fc, "\n")
  
  if (n_sig_05 > 0) {
    sig_results <- results_primary[results_primary$FDR < fdr_threshold, ]
    cat("  Median |logFC| (significant):", 
        sprintf("%.3f", median(abs(sig_results$logFC))), "\n")
    cat("  Range logFC (significant): [", 
        sprintf("%.3f", min(sig_results$logFC)), ", ",
        sprintf("%.3f", max(sig_results$logFC)), "]\n")
  }
  cat("\n")
  
  # Save primary results
  results_file <- file.path(expand_path(config$paths$output_tables), 
                           "07_differential_loops_primary_BCV0.4.tsv")
  write.table(results_primary, results_file, 
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("✓ Primary results saved:", results_file, "\n")
  
  results_rds <- file.path(expand_path(config$paths$output_tables), 
                          "07_differential_loops_primary_BCV0.4.rds")
  saveRDS(results_primary, results_rds)
  cat("✓ Primary results (RDS) saved:", results_rds, "\n\n")
  
  #-----------------------------------------------------------------------------
  # 8. SENSITIVITY ANALYSIS
  #-----------------------------------------------------------------------------
  
  cat("STEP 8: Sensitivity analysis across BCV range\n")
  cat("--------------------------------------------------------------------------------\n")
  
  sensitivity_results <- list()
  sensitivity_summary <- data.frame(
    bcv = numeric(),
    dispersion = numeric(),
    n_significant_fdr05 = integer(),
    n_significant_fdr01 = integer(),
    median_logfc = numeric(),
    n_up = integer(),
    n_down = integer(),
    stringsAsFactors = FALSE
  )
  
  for (bcv in bcv_sensitivity) {
    cat("\nAnalyzing BCV =", bcv, "(dispersion =", bcv^2, ")\n")
    
    analysis <- run_differential_analysis(y, design, bcv, coef = 2)
    results <- analysis$results
    results <- cbind(y$genes, results)
    results <- classify_significance(results, fdr_05 = fdr_threshold, fdr_01 = fdr_strict)
    
    # Store results
    sensitivity_results[[paste0("BCV_", bcv)]] <- results
    
    # Summary stats
    n_sig_05 <- sum(results$FDR < fdr_threshold)
    n_sig_01 <- sum(results$FDR < fdr_strict)
    sig_loops <- results[results$FDR < fdr_threshold, ]
    median_lfc <- if (nrow(sig_loops) > 0) median(abs(sig_loops$logFC)) else NA
    
    sensitivity_summary <- rbind(sensitivity_summary, data.frame(
      bcv = bcv,
      dispersion = bcv^2,
      n_significant_fdr05 = n_sig_05,
      n_significant_fdr01 = n_sig_01,
      median_logfc = median_lfc,
      n_up = sum(results$direction == "up_in_mut"),
      n_down = sum(results$direction == "down_in_mut"),
      stringsAsFactors = FALSE
    ))
    
    cat("  Significant loops (FDR < 0.05):", n_sig_05, "\n")
  }
  
  cat("\n✓ Sensitivity analysis complete\n\n")
  cat("Sensitivity summary:\n")
  print(sensitivity_summary)
  cat("\n")
  
  # Identify core differential set (significant across all BCVs)
  sig_lists <- lapply(sensitivity_results, function(res) {
    res$loop_id[res$FDR < fdr_threshold]
  })
  core_differential <- Reduce(intersect, sig_lists)
  
  cat("Core differential loops (significant across all BCV values):", 
      length(core_differential), "\n\n")
  
  # Save sensitivity results
  sensitivity_file <- file.path(expand_path(config$paths$output_tables), 
                                "07_sensitivity_analysis.rds")
  saveRDS(list(
    results = sensitivity_results,
    summary = sensitivity_summary,
    core_differential = core_differential
  ), sensitivity_file)
  cat("✓ Sensitivity analysis saved:", sensitivity_file, "\n\n")
  
  #-----------------------------------------------------------------------------
  # 9. QUALITY CONTROL PLOTS
  #-----------------------------------------------------------------------------
  
  cat("STEP 9: Generating diagnostic plots\n")
  cat("--------------------------------------------------------------------------------\n")
  
  plot_dir <- expand_path(config$paths$output_plots)
  colors <- config$visualization$colors
  
  # MA Plot
  cat("Creating MA plot...\n")
  pdf(file.path(plot_dir, "07_MA_plot.pdf"), width = 8, height = 6)
  
  results_primary$color <- colors$significance$not_sig
  results_primary$color[results_primary$FDR < fdr_threshold] <- 
    colors$significance$fdr_05
  results_primary$color[results_primary$FDR < fdr_strict] <- 
    colors$significance$fdr_01
  
  plot(results_primary$logCPM, results_primary$logFC,
       pch = 16, cex = 0.5, col = results_primary$color,
       xlab = "Average log2 CPM",
       ylab = "log2 Fold Change (mut vs ctrl)",
       main = paste0("MA Plot (BCV = ", bcv_primary, ")"))
  abline(h = c(-1, 0, 1), lty = 2, col = "gray50")
  legend("topright", 
         legend = c("Not significant", 
                   paste0("FDR < ", fdr_threshold), 
                   paste0("FDR < ", fdr_strict)),
         col = c(colors$significance$not_sig, 
                colors$significance$fdr_05, 
                colors$significance$fdr_01),
         pch = 16, cex = 0.8)
  
  dev.off()
  
  # Volcano Plot
  cat("Creating volcano plot...\n")
  pdf(file.path(plot_dir, "07_volcano_plot.pdf"), width = 8, height = 6)
  
  plot(results_primary$logFC, -log10(results_primary$PValue),
       pch = 16, cex = 0.5, col = results_primary$color,
       xlab = "log2 Fold Change (mut vs ctrl)",
       ylab = "-log10 P-value",
       main = paste0("Volcano Plot (BCV = ", bcv_primary, ")"))
  abline(v = c(-1, 1), lty = 2, col = "gray50")
  abline(h = -log10(fdr_threshold), lty = 2, col = "red")
  legend("topright", 
         legend = c("Not significant", 
                   paste0("FDR < ", fdr_threshold), 
                   paste0("FDR < ", fdr_strict)),
         col = c(colors$significance$not_sig, 
                colors$significance$fdr_05, 
                colors$significance$fdr_01),
         pch = 16, cex = 0.8)
  
  dev.off()
  
  # Dispersion sensitivity plot
  cat("Creating dispersion sensitivity plot...\n")
  pdf(file.path(plot_dir, "07_dispersion_sensitivity.pdf"), width = 8, height = 6)
  
  plot(sensitivity_summary$bcv, sensitivity_summary$n_significant_fdr05,
       type = "b", pch = 16, col = "blue", lwd = 2,
       xlab = "BCV",
       ylab = "Number of Significant Loops",
       main = "Sensitivity to Dispersion Parameter",
       ylim = c(0, max(sensitivity_summary$n_significant_fdr05) * 1.1))
  lines(sensitivity_summary$bcv, sensitivity_summary$n_significant_fdr01,
        type = "b", pch = 16, col = "red", lwd = 2)
  abline(v = bcv_primary, lty = 2, col = "gray50")
  legend("topright", 
         legend = c(paste0("FDR < ", fdr_threshold), 
                   paste0("FDR < ", fdr_strict),
                   paste0("Primary BCV (", bcv_primary, ")")),
         col = c("blue", "red", "gray50"),
         lty = c(1, 1, 2), lwd = c(2, 2, 1), pch = c(16, 16, NA))
  
  dev.off()
  
  # P-value distribution (QC check)
  cat("Creating p-value distribution plot...\n")
  pdf(file.path(plot_dir, "07_pvalue_distribution.pdf"), width = 8, height = 6)
  
  hist(results_primary$PValue, breaks = 50, 
       col = "lightblue", border = "white",
       xlab = "P-value", main = "P-value Distribution",
       freq = FALSE)
  abline(h = 1, lty = 2, col = "red", lwd = 2)
  
  dev.off()
  
  cat("✓ All diagnostic plots generated in:", plot_dir, "\n\n")
  
  #-----------------------------------------------------------------------------
  # 10. FINAL SUMMARY REPORT
  #-----------------------------------------------------------------------------
  
  cat("STEP 10: Generating summary report\n")
  cat("--------------------------------------------------------------------------------\n")
  
  summary_stats <- list(
    analysis_date = Sys.Date(),
    primary_bcv = bcv_primary,
    total_loops_input = expected_loops,
    loops_after_filtering = nrow(y),
    filter_rate = filter_rate,
    normalization_method = "TMM",
    norm_factors = y$samples$norm.factors,
    effective_lib_sizes = y$samples$lib.size * y$samples$norm.factors,
    n_significant_fdr05 = n_sig_05,
    n_significant_fdr01 = n_sig_01,
    n_up_in_mut = n_up,
    n_down_in_mut = n_down,
    n_high_fc = n_high_fc,
    sensitivity_range = range(sensitivity_summary$n_significant_fdr05),
    core_differential_loops = length(core_differential)
  )
  
  # Save summary
  summary_file <- file.path(expand_path(config$paths$output_reports), 
                           "07_analysis_summary.rds")
  saveRDS(summary_stats, summary_file)
  
  # Print final summary
  cat("\n")
  cat("================================================================================\n")
  cat("ANALYSIS COMPLETE - SUMMARY\n")
  cat("================================================================================\n\n")
  cat("Input:\n")
  cat("  Total loops (before filtering):", summary_stats$total_loops_input, "\n")
  cat("  Loops analyzed (after filtering):", summary_stats$loops_after_filtering, "\n")
  cat("  Filter rate:", sprintf("%.1f%%", summary_stats$filter_rate * 100), "\n\n")
  
  cat("Primary Analysis (BCV =", summary_stats$primary_bcv, "):\n")
  cat("  Significant loops (FDR < 0.05):", summary_stats$n_significant_fdr05, "\n")
  cat("  Significant loops (FDR < 0.01):", summary_stats$n_significant_fdr01, "\n")
  cat("  Up-regulated in mutant:", summary_stats$n_up_in_mut, "\n")
  cat("  Down-regulated in mutant:", summary_stats$n_down_in_mut, "\n")
  cat("  High fold-change (|logFC| > 1):", summary_stats$n_high_fc, "\n\n")
  
  cat("Sensitivity Analysis:\n")
  cat("  BCV range tested:", paste(bcv_sensitivity, collapse = ", "), "\n")
  cat("  Significant loops range:", 
      paste(summary_stats$sensitivity_range, collapse = " - "), "\n")
  cat("  Core differential set:", summary_stats$core_differential_loops, "\n\n")
  
  cat("Output Files:\n")
  cat("  Primary results (TSV):", results_file, "\n")
  cat("  Primary results (RDS):", results_rds, "\n")
  cat("  Sensitivity analysis:", sensitivity_file, "\n")
  cat("  Summary statistics:", summary_file, "\n")
  cat("  Diagnostic plots:", plot_dir, "\n\n")
  
  cat("Session Info:\n")
  print(sessionInfo())
  cat("\n")
  
  cat("================================================================================\n")
  cat("edgeR differential analysis completed successfully!\n")
  cat("Next step: Run 08_compare_methods.R to compare with Hiccups results\n")
  cat("================================================================================\n")
  
  return(invisible(summary_stats))
}

################################################################################
# EXECUTE MAIN PIPELINE
################################################################################

if (!interactive()) {
  tryCatch({
    main()
  }, error = function(e) {
    cat("\n")
    cat("================================================================================\n")
    cat("ERROR: Analysis failed\n")
    cat("================================================================================\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Traceback:\n")
    print(traceback())
    quit(status = 1)
  })
}
