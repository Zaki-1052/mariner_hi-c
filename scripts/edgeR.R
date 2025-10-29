# scripts/01_run_edgeR_analysis.R
# Replicate-aware edgeR differential loop analysis for Mariner Hi-C data
# Modified for biological replicates (n=3 per condition)
# Multi-resolution support: accepts resolution as command-line argument
# Author: Zakir
# Date: 2025-10-20

# =============================================================================
# SETUP AND CONFIGURATION
# =============================================================================

# Parse command-line arguments for resolution
args <- commandArgs(trailingOnly = TRUE)
RESOLUTION <- if (length(args) > 0) as.numeric(args[1]) else 5000

cat("\n========================================\n")
cat("edgeR Differential Loop Analysis\n")
cat(sprintf("RESOLUTION: %d bp (%d kb)\n", RESOLUTION, RESOLUTION/1000))
cat("========================================\n\n")

# Load required libraries
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(GenomicRanges)
  library(InteractionSet)
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# Load configuration
cat("Loading configuration...\n")
config <- yaml::read_yaml("config/edgeR_config.yaml")

# Override config paths to use resolution-specific directories
input_dir <- sprintf("outputs/res_%dkb", RESOLUTION/1000)
config$paths$input$counts_matrix <- file.path(input_dir, "06_counts_matrix.rds")
config$paths$input$coordinates <- file.path(input_dir, "03_binned.rds")
config$paths$input$merged_loops <- file.path(input_dir, "02_merged.rds")
config$paths$input$qc_summary <- file.path(input_dir, "qc_report/qc_report_summary.rds")

# Update output paths to be resolution-specific
config$paths$output$base <- file.path("outputs", sprintf("edgeR_results_res_%dkb", RESOLUTION/1000))
config$paths$output$primary <- file.path(config$paths$output$base, "primary_analysis")
config$paths$output$comparison <- file.path(config$paths$output$base, "hiccups_comparison")
config$paths$output$plots <- file.path(config$paths$output$base, "plots")
config$paths$output$logs <- file.path(config$paths$output$base, "logs")

# Set working directory to base
setwd(config$paths$base_dir)

# Set seed for reproducibility
set.seed(config$runtime$seed)

# Create output directories
cat("Creating output directories...\n")
output_dirs <- c(
  config$paths$output$base,
  config$paths$output$primary,
  config$paths$output$plots,
  config$paths$output$logs
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Start logging
log_file <- file.path(config$paths$output$logs, "analysis_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "message", append = TRUE)

cat("\nSetup complete\n")
cat("   - Config loaded from: config/edgeR_config.yaml\n")
cat("   - Output directory: ", config$paths$output$base, "\n")
cat("   - Random seed: ", config$runtime$seed, "\n\n")

# =============================================================================
# DATA LOADING AND VALIDATION
# =============================================================================

cat("Loading input data...\n")

# Load count matrix
cat("   - Loading count matrix... ")
counts_matrix <- readRDS(config$paths$input$counts_matrix)
cat("(", nrow(counts_matrix), " loops × ", ncol(counts_matrix), " samples)\n", sep = "")

# Validate expected number of samples
expected_samples <- length(config$samples$names)
if (ncol(counts_matrix) != expected_samples) {
  stop(sprintf("ERROR: Expected %d samples but count matrix has %d columns",
               expected_samples, ncol(counts_matrix)))
}

# Load coordinates (GInteractions)
cat("   - Loading loop coordinates... ")
binned_gi <- readRDS(config$paths$input$coordinates)
cat("(", length(binned_gi), " genomic positions)\n", sep = "")

# Load QC summary for shift status (optional - may not exist yet)
shift_status <- NULL
if (file.exists(config$paths$input$qc_summary)) {
  cat("   - Loading QC summary... ")
  qc_summary <- readRDS(config$paths$input$qc_summary)
  if ("shift_analysis" %in% names(qc_summary) &&
      "loop_shift_status" %in% names(qc_summary$shift_analysis)) {
    shift_status <- qc_summary$shift_analysis$loop_shift_status
    cat("(shift data for ", length(shift_status), " loops)\n", sep = "")
  } else {
    cat("(no shift data available)\n")
  }
} else {
  cat("   - QC summary not found (skipping shift analysis)\n")
}

# Validate data integrity
cat("\nValidating data integrity...\n")

# Check dimensions match
stopifnot(
  "Count matrix and coordinates dimension mismatch" =
    nrow(counts_matrix) == length(binned_gi)
)

# Check for NAs
na_count <- sum(is.na(counts_matrix))
stopifnot("Count matrix contains NA values" = na_count == 0)

# Check for negative values
neg_count <- sum(counts_matrix < 0)
stopifnot("Count matrix contains negative values" = neg_count == 0)

# Check for all-zero rows (will be filtered later, but log now)
zero_rows <- sum(rowSums(counts_matrix) == 0)
if (zero_rows > 0) {
  cat("   Warning:", zero_rows, "loops have zero total counts (will be filtered)\n")
}

cat("   ✓ Data integrity checks passed\n")

# =============================================================================
# CREATE GENE ANNOTATION
# =============================================================================

cat("\nCreating genomic annotations...\n")

# Extract coordinates from GInteractions
# anchors(binned_gi) returns a list with first and second anchors
anchor1 <- anchors(binned_gi, type = "first")
anchor2 <- anchors(binned_gi, type = "second")

# Create annotation dataframe
genes_df <- data.frame(
  loop_id = paste0("loop_", seq_len(length(binned_gi))),
  chr1 = as.character(seqnames(anchor1)),
  start1 = start(anchor1),
  end1 = end(anchor1),
  chr2 = as.character(seqnames(anchor2)),
  start2 = start(anchor2),
  end2 = end(anchor2),
  stringsAsFactors = FALSE
)

# Create compact coordinate string for display
genes_df$coord_string <- paste0(
  genes_df$chr1, ":", genes_df$start1, "-", genes_df$end1, "_",
  genes_df$chr2, ":", genes_df$start2, "-", genes_df$end2
)

# Add shift status if available
if (!is.null(shift_status)) {
  genes_df$shift_status <- shift_status
}

# Set rownames to match count matrix
rownames(genes_df) <- rownames(counts_matrix)

cat("   ✓ Created annotation for", nrow(genes_df), "loops\n")
cat("   - Columns:", paste(colnames(genes_df), collapse = ", "), "\n")
if (!is.null(shift_status)) {
  cat("   - Shifted loops:", sum(genes_df$shift_status),
      "(", round(100 * mean(genes_df$shift_status), 1), "%)\n", sep = "")
}

# =============================================================================
# CREATE DGEList OBJECT
# =============================================================================

cat("\nCreating DGEList object...\n")

# Create group factor from config
# config$samples$groups should be c("ctrl", "ctrl", "ctrl", "mut", "mut", "mut")
group <- factor(config$samples$groups, levels = c("ctrl", "mut"))

# Verify sample names match count matrix columns
if (!all(colnames(counts_matrix) == config$samples$names)) {
  cat("   WARNING: Renaming count matrix columns to match config\n")
  colnames(counts_matrix) <- config$samples$names
}

# Create DGEList
y <- DGEList(
  counts = counts_matrix,
  group = group,
  genes = genes_df
)

# Add sample names to samples dataframe
y$samples$sample_name <- config$samples$names

cat("   ✓ DGEList created\n")
cat("   - Total samples:", ncol(y), "\n")
cat("   - Groups:", paste(unique(group), collapse = ", "), "\n")
cat("   - Replicates per group: n =", table(group)[1], "\n")
cat("\n   Per-sample library sizes:\n")
for (i in 1:ncol(y)) {
  cat(sprintf("     %s (%s): %.0f\n",
              y$samples$sample_name[i],
              y$samples$group[i],
              y$samples$lib.size[i]))
}
cat(sprintf("\n   Median library size by group:\n"))
cat(sprintf("     ctrl: %.0f\n", median(y$samples$lib.size[group == "ctrl"])))
cat(sprintf("     mut: %.0f\n", median(y$samples$lib.size[group == "mut"])))

# =============================================================================
# FILTERING LOW-COUNT LOOPS
# =============================================================================

cat("\nFiltering low-count loops...\n")

# Apply filterByExpr with custom parameters
keep <- filterByExpr(
  y,
  group = y$samples$group,
  min.count = config$statistics$filtering$min_count,
  min.total.count = config$statistics$filtering$min_total_count,
  min.prop = config$statistics$filtering$min_prop
)

cat("   - Filtering with min.count =", config$statistics$filtering$min_count, "\n")
cat("   - Before filtering:", nrow(y), "loops\n")
cat("   - Loops passing filter:", sum(keep), "\n")
cat("   - Loops filtered out:", sum(!keep), 
    "(", round(100 * mean(!keep), 1), "%)\n", sep = "")

# Apply filter and recalculate library sizes
y <- y[keep, , keep.lib.sizes = FALSE]

cat("   After filtering:", nrow(y), "loops retained\n")
cat("   - New library sizes: ctrl =", y$samples$lib.size[1], 
    ", mut =", y$samples$lib.size[2], "\n")

# =============================================================================
# NORMALIZATION
# =============================================================================

cat("\nPerforming TMM normalization...\n")

# Calculate normalization factors
y <- normLibSizes(y, method = config$statistics$normalization_method)

cat("   ✓ Normalization complete\n")
cat("   - Method:", config$statistics$normalization_method, "\n")
cat("   - Normalization factors: ctrl =", round(y$samples$norm.factors[1], 4), 
    ", mut =", round(y$samples$norm.factors[2], 4), "\n")
cat("   - Effective library sizes: ctrl =", 
    round(y$samples$lib.size[1] * y$samples$norm.factors[1]), 
    ", mut =", round(y$samples$lib.size[2] * y$samples$norm.factors[2]), "\n")

# =============================================================================
# SAMPLE QUALITY CONTROL - MDS PLOT
# =============================================================================

cat("\nGenerating MDS plot for sample QC...\n")

pdf(file.path(config$paths$output$plots, "mds_plot.pdf"),
    width = 8, height = 6)

# MDS plot colored by group
plotMDS(
  y,
  col = c(rep("blue", 3), rep("red", 3)),
  pch = 16,
  cex = 2,
  main = "MDS Plot - Sample Relationships",
  labels = y$samples$sample_name
)

legend(
  "topright",
  legend = c("Control", "Mutant"),
  col = c("blue", "red"),
  pch = 16,
  pt.cex = 2,
  bty = "n"
)

dev.off()

cat("   ✓ MDS plot saved\n")
cat("   - Check plot to verify:\n")
cat("     • Replicates cluster by condition\n")
cat("     • No obvious outliers\n")
cat("     • Clear separation between groups\n\n")

# =============================================================================
# DISPERSION ESTIMATION (DATA-DRIVEN)
# =============================================================================

cat("Estimating dispersions from biological replicates...\n")

# Design matrix - treatment effect parameterization
design <- model.matrix(~group, data = y$samples)
colnames(design) <- c("Intercept", "MutantEffect")

cat("   Design matrix:\n")
cat("     - Intercept: Baseline (control) level\n")
cat("     - MutantEffect: Difference (mutant - control)\n\n")

# Estimate dispersions with robust method
cat("   Estimating common, trended, and tagwise dispersions...\n")
y <- estimateDisp(y, design, robust = TRUE)

cat("   ✓ Dispersion estimation complete\n")
cat(sprintf("     - Common dispersion: %.4f (BCV = %.3f)\n",
            y$common.dispersion, sqrt(y$common.dispersion)))
cat(sprintf("     - Tagwise dispersion range: %.4f - %.4f\n",
            min(y$tagwise.dispersion), max(y$tagwise.dispersion)))
cat(sprintf("     - Median tagwise BCV: %.3f\n\n",
            median(sqrt(y$tagwise.dispersion))))

# Plot BCV
cat("   Generating BCV plot...\n")
pdf(file.path(config$paths$output$plots, "bcv_plot.pdf"),
    width = 8, height = 6)
plotBCV(y, main = "Biological Coefficient of Variation")
dev.off()
cat("   ✓ BCV plot saved\n\n")

# =============================================================================
# QUASI-LIKELIHOOD GLM FIT (RECOMMENDED METHOD)
# =============================================================================

cat("Fitting quasi-likelihood GLM...\n")
cat("   Method: glmQLFit with robust=TRUE\n")
cat("   Benefits:\n")
cat("     • Accounts for gene-specific variability\n")
cat("     • More rigorous error rate control\n")
cat("     • Robust to outliers\n\n")

# QL fit
fit <- glmQLFit(y, design, robust = TRUE)

cat("   ✓ QL fit complete\n")
#cat(sprintf("     - Prior df: %.1f\n", fit$df.prior))
cat(sprintf("     - Residual df: %d\n", min(fit$df.residual)))

# Plot QL dispersions
cat("\n   Generating QL dispersion plot...\n")
pdf(file.path(config$paths$output$plots, "ql_dispersion_plot.pdf"),
    width = 8, height = 6)
plotQLDisp(fit, main = "Quasi-Likelihood Dispersions")
dev.off()
cat("   ✓ QL dispersion plot saved\n\n")

# =============================================================================
# DIFFERENTIAL TESTING
# =============================================================================

cat("Testing for differential loops...\n")
cat("   Test: Mutant effect (coefficient 2)\n")
cat("   Method: Quasi-likelihood F-test\n\n")

# QL F-test for mutant effect
qlf <- glmQLFTest(fit, coef = 2)

# Extract all results
results <- topTags(qlf, n = Inf, sort.by = "none")$table

cat("   ✓ Testing complete\n\n")

# Add significance classification
results$significant <- results$FDR < config$statistics$fdr_primary
results$exploratory <- results$FDR < config$statistics$fdr_exploratory

# Categorize by effect size
results$category <- "non_significant"
results$category[results$significant & abs(results$logFC) > config$statistics$fold_change_thresholds$strong] <- "strong_differential"
results$category[results$significant & abs(results$logFC) > config$statistics$fold_change_thresholds$moderate &
                 abs(results$logFC) <= config$statistics$fold_change_thresholds$strong] <- "moderate_differential"
results$category[results$significant & abs(results$logFC) <= config$statistics$fold_change_thresholds$moderate] <- "weak_differential"
results$category[!results$significant & results$exploratory] <- "trending"

# Add direction
results$direction <- "unchanged"
results$direction[results$significant & results$logFC > 0] <- "up_in_mutant"
results$direction[results$significant & results$logFC < 0] <- "down_in_mutant"

# Print summary
cat("=== Differential Loop Results ===\n")
n_sig <- sum(results$significant)
n_up <- sum(results$direction == "up_in_mutant")
n_down <- sum(results$direction == "down_in_mutant")

cat(sprintf("   Significant loops (FDR < %.2f): %d (%.1f%%)\n",
            config$statistics$fdr_primary, n_sig,
            100 * n_sig / nrow(results)))
cat(sprintf("     - Up in mutant: %d\n", n_up))
cat(sprintf("     - Down in mutant: %d\n", n_down))
cat(sprintf("\n   By effect size:\n"))
cat(sprintf("     - Strong (|logFC| > 1): %d\n",
            sum(results$category == "strong_differential")))
cat(sprintf("     - Moderate (|logFC| > 0.5): %d\n",
            sum(results$category == "moderate_differential")))
cat(sprintf("     - Weak: %d\n",
            sum(results$category == "weak_differential")))
cat(sprintf("\n   Trending (FDR < 0.10): %d\n\n",
            sum(results$category == "trending")))

# =============================================================================
# SHIFTED LOOP ENRICHMENT ANALYSIS (if shift data available)
# =============================================================================

shifted_tested <- NA
shifted_sig <- NA
nonshifted_tested <- NA
nonshifted_sig <- NA
fisher_result <- NULL

if (!is.null(shift_status)) {
  cat("Analyzing shifted loop enrichment...\n")

  # Calculate enrichment of shifted loops in significant set
  shifted_tested <- sum(results$shift_status)
  shifted_sig <- sum(results$shift_status & results$significant)
  nonshifted_tested <- sum(!results$shift_status)
  nonshifted_sig <- sum(!results$shift_status & results$significant)

  # Fisher's exact test
  contingency <- matrix(
    c(shifted_sig, shifted_tested - shifted_sig,
      nonshifted_sig, nonshifted_tested - nonshifted_sig),
    nrow = 2,
    dimnames = list(
      c("Shifted", "Non-shifted"),
      c("Significant", "Non-significant")
    )
  )

  fisher_result <- fisher.test(contingency)

  cat("   - Shifted loops tested:", shifted_tested, "\n")
  cat("   - Shifted loops significant:", shifted_sig,
      "(", round(100 * shifted_sig / shifted_tested, 1), "%)\n", sep = "")
  cat("   - Non-shifted loops tested:", nonshifted_tested, "\n")
  cat("   - Non-shifted loops significant:", nonshifted_sig,
      "(", round(100 * nonshifted_sig / nonshifted_tested, 1), "%)\n", sep = "")
  cat("   - Fisher's exact test p-value:", format.pval(fisher_result$p.value, digits = 3), "\n")
  cat("   - Odds ratio:", round(fisher_result$estimate, 3), "\n\n")
} else {
  cat("Skipping shifted loop enrichment (no shift data available)\n\n")
}

# =============================================================================
# GENERATE SUMMARY STATISTICS
# =============================================================================

cat("Computing summary statistics...\n")

summary_stats <- list(
  total_loops_input = length(binned_gi),
  total_loops_filtered = sum(!keep),
  total_loops_tested = nrow(results),
  
  # Significance counts
  significant_fdr05 = sum(results$significant),
  significant_fdr10 = sum(results$exploratory),
  
  # Direction
  up_in_mutant = sum(results$direction == "up_in_mutant"),
  down_in_mutant = sum(results$direction == "down_in_mutant"),
  
  # Effect size categories
  strong_differential = sum(results$category == "strong_differential"),
  moderate_differential = sum(results$category == "moderate_differential"),
  weak_differential = sum(results$category == "weak_differential"),
  trending = sum(results$category == "trending"),
  
  # Fold change statistics
  median_abs_logfc_all = median(abs(results$logFC)),
  median_abs_logfc_sig = median(abs(results$logFC[results$significant])),
  max_abs_logfc = max(abs(results$logFC)),
  
  # Shifted loop statistics
  shifted_tested = shifted_tested,
  shifted_significant = shifted_sig,
  shifted_enrichment_pval = fisher_result$p.value,
  shifted_odds_ratio = as.numeric(fisher_result$estimate),
  
  # Library information
  ctrl_lib_size = y$samples$lib.size[1],
  mut_lib_size = y$samples$lib.size[2],
  ctrl_norm_factor = y$samples$norm.factors[1],
  mut_norm_factor = y$samples$norm.factors[2],
  
  # Analysis parameters
  primary_bcv = config$statistics$primary_bcv,
  fdr_threshold = config$statistics$fdr_primary,
  normalization_method = config$statistics$normalization_method
)

cat("   ✓ Summary statistics computed\n\n")

# =============================================================================
# VISUALIZATION
# =============================================================================

cat("Generating visualizations...\n")

# Prepare color mapping
results$plot_color <- config$visualization$colors$non_significant
results$plot_color[results$significant & results$logFC > 0] <- config$visualization$colors$significant_up
results$plot_color[results$significant & results$logFC < 0] <- config$visualization$colors$significant_down

# 1. MA Plot
cat("   - Creating MA plot... ")
pdf(file.path(config$paths$output$plots, "ma_plot_primary.pdf"),
    width = config$visualization$width, 
    height = config$visualization$height)

par(mar = c(5, 5, 4, 2))
plot(
  results$logCPM,
  results$logFC,
  pch = 16,
  cex = config$visualization$ma_plot$point_size,
  col = adjustcolor(results$plot_color, alpha.f = config$visualization$ma_plot$alpha),
  xlab = "Average log2 CPM",
  ylab = "log2 Fold Change (Mutant / Control)",
  main = paste0("MA Plot (BCV = ", config$statistics$primary_bcv, ")")
)

# Add reference lines
abline(h = 0, col = "black", lty = 2, lwd = 1.5)
abline(h = config$visualization$ma_plot$fc_lines, col = "gray40", lty = 2, lwd = 1)

# Add legend
legend(
  "topright",
  legend = c(
    paste0("Up in mutant (n=", sum(results$direction == "up_in_mutant"), ")"),
    paste0("Down in mutant (n=", sum(results$direction == "down_in_mutant"), ")"),
    "Not significant"
  ),
  col = c(
    config$visualization$colors$significant_up,
    config$visualization$colors$significant_down,
    config$visualization$colors$non_significant
  ),
  pch = 16,
  pt.cex = 1.5,
  bty = "n"
)

dev.off()
cat("✓\n")

# 2. Volcano Plot
cat("   - Creating volcano plot... ")
pdf(file.path(config$paths$output$plots, "volcano_plot_primary.pdf"),
    width = config$visualization$width,
    height = config$visualization$height)

par(mar = c(5, 5, 4, 2))
plot(
  results$logFC,
  -log10(results$PValue),
  pch = 16,
  cex = config$visualization$volcano_plot$point_size,
  col = adjustcolor(results$plot_color, alpha.f = config$visualization$volcano_plot$alpha),
  xlab = "log2 Fold Change (Mutant / Control)",
  ylab = "-log10(P-value)",
  main = paste0("Volcano Plot (BCV = ", config$statistics$primary_bcv, ")")
)

# Add reference lines
abline(v = 0, col = "black", lty = 2, lwd = 1.5)
abline(v = c(-1, 1), col = "gray40", lty = 2, lwd = 1)
abline(h = -log10(config$statistics$fdr_primary), col = "red", lty = 2, lwd = 1.5)

# Label top differential loops
top_indices <- order(results$PValue)[1:min(config$visualization$volcano_plot$label_top_n, sum(results$significant))]
if (length(top_indices) > 0) {
  text(
    results$logFC[top_indices],
    -log10(results$PValue[top_indices]),
    labels = results$loop_id[top_indices],
    cex = 0.6,
    pos = 3
  )
}

# Add legend
legend(
  "topright",
  legend = c(
    paste0("Up in mutant (n=", sum(results$direction == "up_in_mutant"), ")"),
    paste0("Down in mutant (n=", sum(results$direction == "down_in_mutant"), ")"),
    paste0("FDR < ", config$statistics$fdr_primary)
  ),
  col = c(
    config$visualization$colors$significant_up,
    config$visualization$colors$significant_down,
    "red"
  ),
  lty = c(NA, NA, 2),
  pch = c(16, 16, NA),
  pt.cex = 1.5,
  lwd = 1.5,
  bty = "n"
)

dev.off()
cat("✓\n")

# 3. decideTests Summary Plot
cat("   - Creating results summary plot... ")
pdf(file.path(config$paths$output$plots, "results_summary.pdf"),
    width = 8, height = 6)

par(mar = c(5, 5, 4, 2))

# Create summary barplot
summary_counts <- c(
  n_up,
  n_down,
  sum(results$category == "strong_differential"),
  sum(results$category == "moderate_differential"),
  sum(results$category == "weak_differential")
)

names_labels <- c(
  "Up in\nmutant",
  "Down in\nmutant",
  "Strong\n(|logFC|>1)",
  "Moderate\n(|logFC|>0.5)",
  "Weak"
)

barplot(
  summary_counts,
  names.arg = names_labels,
  col = c(config$visualization$colors$significant_up,
          config$visualization$colors$significant_down,
          "darkgreen", "forestgreen", "lightgreen"),
  ylab = "Number of Loops",
  main = sprintf("Differential Loop Summary (FDR < %.2f)", config$statistics$fdr_primary),
  las = 1,
  ylim = c(0, max(summary_counts) * 1.2)
)

# Add counts on top of bars
text(
  x = seq_along(summary_counts) * 1.2 - 0.5,
  y = summary_counts + max(summary_counts) * 0.05,
  labels = summary_counts,
  cex = 1.2
)

dev.off()
cat("✓\n")

# 4. Shifted Loop Enrichment Plot (if shift data available)
if (!is.null(shift_status) && !is.null(fisher_result)) {
  cat("   - Creating shifted loop enrichment plot... ")
  pdf(file.path(config$paths$output$plots, "shifted_loop_enrichment.pdf"),
      width = 8, height = 6)

  par(mar = c(5, 5, 4, 2))
  prop_data <- c(
    100 * shifted_sig / shifted_tested,
    100 * nonshifted_sig / nonshifted_tested
  )

  barplot(
    prop_data,
    names.arg = c("Shifted Loops", "Non-shifted Loops"),
    col = c(config$visualization$colors$shifted_loops, "gray70"),
    ylab = "% Significant (FDR < 0.05)",
    main = "Differential Loop Rate by Shift Status",
    ylim = c(0, max(prop_data) * 1.2),
    las = 1
  )

  # Add text labels
  text(
    x = c(0.7, 1.9),
    y = prop_data + max(prop_data) * 0.05,
    labels = paste0(
      round(prop_data, 1), "%\n(",
      c(shifted_sig, nonshifted_sig), "/",
      c(shifted_tested, nonshifted_tested), ")"
    ),
    cex = 0.9
  )

  # Add Fisher's test result
  text(
    x = 1.3,
    y = max(prop_data) * 1.1,
    labels = paste0(
      "Fisher's exact test\n",
      "OR = ", round(fisher_result$estimate, 2),
      ", p = ", format.pval(fisher_result$p.value, digits = 2)
    ),
    cex = 0.8
  )

  dev.off()
  cat("✓\n")
} else {
  cat("   - Skipping shifted loop enrichment plot (no shift data)\n")
}

cat("\n   All plots generated\n\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("Saving results...\n")

# Save DGEList object
saveRDS(y, file.path(config$paths$output$primary, "edgeR_dge_object.rds"))
cat("   ✓ DGEList object saved\n")

# Save primary results (all formats)
write.table(
  results,
  file = file.path(config$paths$output$primary, "all_results_primary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("   ✓ Primary results saved (TSV)\n")

saveRDS(
  results,
  file.path(config$paths$output$primary, "all_results_primary.rds")
)
cat("   ✓ Primary results saved (RDS)\n")

# Save significant loops only
sig_results <- results[results$significant, ]
write.table(
  sig_results,
  file = file.path(config$paths$output$primary, "significant_loops_fdr05.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("   ✓ Significant loops saved (", nrow(sig_results), " loops)\n", sep = "")

# Save top 100 by effect size
top100 <- results[order(abs(results$logFC), decreasing = TRUE)[1:min(100, nrow(results))], ]
write.table(
  top100,
  file = file.path(config$paths$output$primary, "top100_differential.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("   ✓ Top 100 by effect size saved\n")

# Save final results with stringent thresholds (|logFC| > 0.3, FDR < 0.03)
final_results <- results[abs(results$logFC) > 0.3 & results$FDR < 0.03, ]
write.table(
  final_results,
  file = file.path(config$paths$output$primary, "final_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("   ✓ Final results saved (|logFC| > 0.3, FDR < 0.03): ", nrow(final_results), " loops\n", sep = "")

# Save decideTests summary
dt_summary <- summary(decideTests(qlf))
write.table(
  dt_summary,
  file = file.path(config$paths$output$primary, "decideTests_summary.txt"),
  sep = "\t",
  quote = FALSE
)
cat("   ✓ decideTests summary saved\n")

# Save summary statistics
summary_text <- capture.output({
  cat("========================================\n")
  cat("edgeR Differential Loop Analysis Summary\n")
  cat("========================================\n\n")
  cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("INPUT DATA\n")
  cat("----------\n")
  cat("Total loops (input):", summary_stats$total_loops_input, "\n")
  cat("Loops filtered out:", summary_stats$total_loops_filtered, "\n")
  cat("Loops tested:", summary_stats$total_loops_tested, "\n\n")
  
  cat("LIBRARY INFORMATION\n")
  cat("-------------------\n")
  cat("Control library size:", summary_stats$ctrl_lib_size, "\n")
  cat("Mutant library size:", summary_stats$mut_lib_size, "\n")
  cat("Control normalization factor:", round(summary_stats$ctrl_norm_factor, 4), "\n")
  cat("Mutant normalization factor:", round(summary_stats$mut_norm_factor, 4), "\n")
  cat("Normalization method:", summary_stats$normalization_method, "\n\n")
  
  cat("ANALYSIS PARAMETERS\n")
  cat("-------------------\n")
  cat("Experimental design: n=3 replicates per condition\n")
  cat("Statistical method: Quasi-likelihood GLM with robust estimation\n")
  cat("Estimated common BCV:", round(sqrt(y$common.dispersion), 3), "\n")
  cat("Median tagwise BCV:", round(median(sqrt(y$tagwise.dispersion)), 3), "\n")
  cat("Residual degrees of freedom:", min(fit$df.residual), "\n")
  cat("Prior df for QL:", round(fit$df.prior, 1), "\n")
  cat("FDR threshold:", summary_stats$fdr_threshold, "\n")
  cat("Filtering min count:", config$statistics$filtering$min_count, "\n\n")
  
  cat("DIFFERENTIAL LOOPS (FDR < 0.05)\n")
  cat("--------------------------------\n")
  cat("Total significant:", summary_stats$significant_fdr05, 
      "(", round(100 * summary_stats$significant_fdr05 / summary_stats$total_loops_tested, 2), "%)\n", sep = "")
  cat("  - Up in mutant:", summary_stats$up_in_mutant, "\n")
  cat("  - Down in mutant:", summary_stats$down_in_mutant, "\n\n")
  
  cat("EFFECT SIZE CATEGORIES\n")
  cat("----------------------\n")
  cat("Strong (|logFC| > 1):", summary_stats$strong_differential, "\n")
  cat("Moderate (|logFC| > 0.5):", summary_stats$moderate_differential, "\n")
  cat("Weak (|logFC| <= 0.5):", summary_stats$weak_differential, "\n")
  cat("Trending (FDR < 0.10):", summary_stats$trending, "\n\n")
  
  cat("FOLD CHANGE STATISTICS\n")
  cat("----------------------\n")
  cat("Median |logFC| (all):", round(summary_stats$median_abs_logfc_all, 3), "\n")
  cat("Median |logFC| (significant):", round(summary_stats$median_abs_logfc_sig, 3), "\n")
  cat("Max |logFC|:", round(summary_stats$max_abs_logfc, 3), "\n\n")
  
  cat("SHIFTED LOOP ANALYSIS\n")
  cat("---------------------\n")
  cat("Shifted loops tested:", summary_stats$shifted_tested, "\n")
  cat("Shifted loops significant:", summary_stats$shifted_significant,
      "(", round(100 * summary_stats$shifted_significant / summary_stats$shifted_tested, 1), "%)\n", sep = "")
  cat("Enrichment p-value:", format.pval(summary_stats$shifted_enrichment_pval, digits = 3), "\n")
  cat("Odds ratio:", round(summary_stats$shifted_odds_ratio, 3), "\n\n")
})

writeLines(summary_text, file.path(config$paths$output$primary, "summary_statistics.txt"))
cat("   ✓ Summary statistics saved\n")

# Save session info
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(config$paths$output$logs, "session_info.txt"))
cat("   ✓ Session info saved\n\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("========================================\n")
cat("REPLICATE-AWARE ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Experimental Design:\n")
cat("   - Samples: 6 (3 ctrl + 3 mut)\n")
cat("   - Method: Quasi-likelihood GLM with robust estimation\n")
cat("   - Estimated BCV:", round(sqrt(y$common.dispersion), 3), "\n\n")

cat("Key Results:\n")
cat("   - Loops tested:", summary_stats$total_loops_tested, "\n")
cat("   - Significant (FDR < 0.05):", summary_stats$significant_fdr05,
    "(", round(100 * summary_stats$significant_fdr05 / summary_stats$total_loops_tested, 1), "%)\n", sep = "")
cat("   - Up in mutant:", summary_stats$up_in_mutant, "\n")
cat("   - Down in mutant:", summary_stats$down_in_mutant, "\n")
cat("   - Strong differential (|logFC| > 1):", summary_stats$strong_differential, "\n\n")

cat("Improvement over merged analysis:\n")
cat("   - Previous (merged, n=1): 2 significant loops\n")
cat("   - Current (replicates, n=3):", summary_stats$significant_fdr05, "significant loops\n")
cat("   - Fold improvement:", round(summary_stats$significant_fdr05 / 2, 1), "×\n\n")

cat("Output Files:\n")
cat("   Primary results:", config$paths$output$primary, "\n")
cat("   Plots:", config$paths$output$plots, "\n")
cat("   Logs:", config$paths$output$logs, "\n\n")

cat("Next Steps:\n")
cat("   1. Review MDS plot to confirm replicate quality\n")
cat("   2. Check BCV and QL dispersion plots\n")
cat("   3. Examine MA and volcano plots for biological patterns\n")
cat("   4. Compare with Hiccups differential results (38,593 loops)\n\n")

# Close log
sink(type = "message")
close(log_con)

cat("Analysis log saved to:", log_file, "\n")

