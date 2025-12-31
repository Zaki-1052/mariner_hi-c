#!/usr/bin/env Rscript
# stripes/scripts/phase3_edgeR.R
# Phase 3: edgeR Differential Stripe Analysis
#
# Performs differential analysis using quasi-likelihood GLM on stripe counts.
# Key difference from loop pipeline: SKIP FILTERING (retain all stripes).
#
# Usage: Rscript phase3_edgeR.R <timepoint>
#   timepoint: "early" (250831) or "late" (250402)

suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(ggplot2)
  library(dplyr)
})

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript phase3_edgeR.R <timepoint>\n",
       "  timepoint: 'early' or 'late'")
}
TIMEPOINT <- args[1]

if (!TIMEPOINT %in% c("early", "late")) {
  stop("Timepoint must be 'early' or 'late'")
}

cat("\n========================================\n")
cat("Phase 3: edgeR Differential Analysis\n")
cat(sprintf("Timepoint: %s\n", TIMEPOINT))
cat("========================================\n\n")

# ==============================================================================
# LOAD CONFIGURATION
# ==============================================================================
cat("Loading configuration...\n")

config_file <- "config/stripe_config.yaml"
if (!file.exists(config_file)) {
  stop(sprintf("Config file not found: %s\n  Run from project base directory",
               config_file))
}

config <- yaml::read_yaml(config_file)
cat(sprintf("  Config: %s\n", config_file))

# Set seed for reproducibility
set.seed(config$runtime$seed)

# ==============================================================================
# PATH SETUP
# ==============================================================================
resolution_kb <- config$stripe_data$resolutions$primary / 1000
input_dir <- file.path(config$outputs$base_dir, TIMEPOINT,
                       paste0("res_", resolution_kb, "kb"))
output_dir <- file.path(input_dir, "edgeR_results")
plots_dir <- file.path(output_dir, "plots")

# Create output directories
for (dir in c(output_dir, plots_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

cat(sprintf("  Input: %s\n", input_dir))
cat(sprintf("  Output: %s\n", output_dir))

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("\n=== Loading Data ===\n")

# Load count matrix from Phase 2
counts_file <- file.path(input_dir, "02_stripe_counts.rds")
if (!file.exists(counts_file)) {
  stop(sprintf("Phase 2 output not found: %s\n  Run phase2_quantification.R first",
               counts_file))
}
counts_matrix <- readRDS(counts_file)
cat(sprintf("  Count matrix: %d stripes x %d samples\n",
            nrow(counts_matrix), ncol(counts_matrix)))

# Load stripe metadata from Phase 1
stripes_file <- file.path(input_dir, "01_unified_stripes.tsv")
stripes <- read.table(stripes_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
cat(sprintf("  Stripe metadata: %d stripes\n", nrow(stripes)))

# Validate dimensions match
if (nrow(counts_matrix) != nrow(stripes)) {
  stop(sprintf("Dimension mismatch: counts (%d) vs stripes (%d)",
               nrow(counts_matrix), nrow(stripes)))
}

# ==============================================================================
# CREATE GENE ANNOTATION
# ==============================================================================
cat("\n=== Creating Annotations ===\n")

genes_df <- data.frame(
  stripe_id = stripes$stripe_id,
  chr = stripes$chr,
  anchor_x1 = stripes$anchor_x1,
  anchor_x2 = stripes$anchor_x2,
  span_y1 = stripes$span_y1,
  span_y2 = stripes$span_y2,
  direction_type = stripes$direction_type,
  source = stripes$source,
  pval_ctrl = stripes$pval_ctrl,
  pval_mut = stripes$pval_mut,
  n_ctrl_reps = stripes$n_ctrl_reps,
  n_mut_reps = stripes$n_mut_reps,
  in_10kb = stripes$in_10kb,
  confidence = stripes$confidence,
  stringsAsFactors = FALSE
)

# Create coordinate string
genes_df$coord_string <- paste0(
  genes_df$chr, ":", genes_df$anchor_x1, "-", genes_df$anchor_x2,
  "->", genes_df$span_y1, "-", genes_df$span_y2
)

rownames(genes_df) <- genes_df$stripe_id
rownames(counts_matrix) <- genes_df$stripe_id

cat(sprintf("  Annotation created for %d stripes\n", nrow(genes_df)))

# ==============================================================================
# CREATE DGEList
# ==============================================================================
cat("\n=== Creating DGEList ===\n")

# Create group factor
group <- factor(config$samples$groups, levels = c("ctrl", "mut"))

# Verify sample names
if (!all(colnames(counts_matrix) == config$samples$names)) {
  cat("  Warning: Renaming columns to match config\n")
  colnames(counts_matrix) <- config$samples$names
}

# Create DGEList
dge <- DGEList(
  counts = counts_matrix,
  group = group,
  genes = genes_df
)

dge$samples$sample_name <- config$samples$names

cat(sprintf("  DGEList created: %d stripes x %d samples\n",
            nrow(dge), ncol(dge)))
cat(sprintf("  Groups: %s\n", paste(levels(group), collapse = ", ")))
cat(sprintf("  Replicates per group: n = %d\n", sum(group == "ctrl")))

cat("\n  Library sizes:\n")
for (i in 1:ncol(dge)) {
  cat(sprintf("    %s (%s): %.0f\n",
              dge$samples$sample_name[i],
              dge$samples$group[i],
              dge$samples$lib.size[i]))
}

# ==============================================================================
# SKIP FILTERING (KEY DIFFERENCE FROM LOOP PIPELINE)
# ==============================================================================
cat("\n=== Filtering ===\n")

if (config$edger$skip_filtering) {
  cat("  SKIPPING low-count filtering (per design decision)\n")
  cat(sprintf("  Retaining ALL %d stripes for analysis\n", nrow(dge)))
  cat("  Rationale: Small feature set (~150-300 stripes)\n")
} else {
  # Optional filtering if config changes
  keep <- filterByExpr(dge, group = group, min.count = config$edger$min_count)
  cat(sprintf("  Loops passing filter: %d / %d\n", sum(keep), length(keep)))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
}

# ==============================================================================
# NORMALIZATION
# ==============================================================================
cat("\n=== Normalization ===\n")

dge <- calcNormFactors(dge, method = config$edger$normalization_method)

cat(sprintf("  Method: %s\n", config$edger$normalization_method))
cat("  Normalization factors:\n")
for (i in 1:ncol(dge)) {
  cat(sprintf("    %s: %.4f\n",
              dge$samples$sample_name[i],
              dge$samples$norm.factors[i]))
}

# ==============================================================================
# MDS PLOT (Sample QC)
# ==============================================================================
cat("\n=== Generating MDS Plot ===\n")

pdf(file.path(plots_dir, "mds_plot.pdf"), width = 8, height = 6)

plotMDS(
  dge,
  col = c(rep("blue", 3), rep("red", 3)),
  pch = 16,
  cex = 2,
  main = sprintf("MDS Plot - Stripe Quantification (%s)", TIMEPOINT),
  labels = dge$samples$sample_name
)

legend(
  "topright",
  legend = c("Control", "BAP1-KO"),
  col = c("blue", "red"),
  pch = 16,
  pt.cex = 2,
  bty = "n"
)

dev.off()
cat("  Saved: mds_plot.pdf\n")

# ==============================================================================
# DISPERSION ESTIMATION
# ==============================================================================
cat("\n=== Dispersion Estimation ===\n")

# Design matrix
design <- model.matrix(~group, data = dge$samples)
colnames(design) <- c("Intercept", "MutantEffect")

cat("  Design matrix:\n")
cat("    - Intercept: Baseline (control)\n")
cat("    - MutantEffect: Difference (mutant - control)\n")

# Estimate dispersions with robust method
dge <- estimateDisp(dge, design, robust = config$edger$robust_dispersion)

cat(sprintf("\n  Common dispersion: %.4f\n", dge$common.dispersion))
cat(sprintf("  Common BCV: %.3f\n", sqrt(dge$common.dispersion)))
cat(sprintf("  Tagwise BCV range: %.3f - %.3f\n",
            min(sqrt(dge$tagwise.dispersion)),
            max(sqrt(dge$tagwise.dispersion))))
cat(sprintf("  Median tagwise BCV: %.3f\n",
            median(sqrt(dge$tagwise.dispersion))))

# BCV plot
pdf(file.path(plots_dir, "bcv_plot.pdf"), width = 8, height = 6)
plotBCV(dge, main = sprintf("BCV Plot - Stripes (%s)", TIMEPOINT))
dev.off()
cat("  Saved: bcv_plot.pdf\n")

# ==============================================================================
# QUASI-LIKELIHOOD GLM FIT
# ==============================================================================
cat("\n=== Quasi-Likelihood GLM ===\n")

fit <- glmQLFit(dge, design, robust = config$edger$robust_dispersion)

cat(sprintf("  Residual df: %d\n", min(fit$df.residual)))

# QL dispersion plot
pdf(file.path(plots_dir, "ql_dispersion_plot.pdf"), width = 8, height = 6)
plotQLDisp(fit, main = sprintf("QL Dispersion - Stripes (%s)", TIMEPOINT))
dev.off()
cat("  Saved: ql_dispersion_plot.pdf\n")

# ==============================================================================
# DIFFERENTIAL TESTING
# ==============================================================================
cat("\n=== Differential Testing ===\n")

qlf <- glmQLFTest(fit, coef = 2)  # Test mutant effect
results <- topTags(qlf, n = Inf, sort.by = "none")$table

cat("  Test: Mutant effect (positive logFC = higher in mutant)\n")

# Add significance flags
results$significant_FDR05 <- results$FDR < config$edger$fdr_primary
results$significant_FDR10 <- results$FDR < config$edger$fdr_exploratory

# Add direction
results$direction_edgeR <- "unchanged"
results$direction_edgeR[results$significant_FDR05 & results$logFC > 0] <- "up_in_mutant"
results$direction_edgeR[results$significant_FDR05 & results$logFC < 0] <- "down_in_mutant"

# Summary statistics
n_sig <- sum(results$significant_FDR05)
n_up <- sum(results$direction_edgeR == "up_in_mutant")
n_down <- sum(results$direction_edgeR == "down_in_mutant")

cat(sprintf("\n  Results (FDR < %.2f):\n", config$edger$fdr_primary))
cat(sprintf("    Significant: %d / %d (%.1f%%)\n",
            n_sig, nrow(results), 100 * n_sig / nrow(results)))
cat(sprintf("    Up in mutant: %d\n", n_up))
cat(sprintf("    Down in mutant: %d\n", n_down))
cat(sprintf("\n  Exploratory (FDR < %.2f): %d\n",
            config$edger$fdr_exploratory, sum(results$significant_FDR10)))

# ==============================================================================
# VOLCANO PLOT
# ==============================================================================
cat("\n=== Generating Volcano Plot ===\n")

# Color mapping
results$plot_color <- "gray60"
results$plot_color[results$significant_FDR05 & results$logFC > 0] <- "firebrick3"
results$plot_color[results$significant_FDR05 & results$logFC < 0] <- "steelblue3"

pdf(file.path(plots_dir, "volcano_plot.pdf"), width = 10, height = 8)

par(mar = c(5, 5, 4, 2))
plot(
  results$logFC,
  -log10(results$PValue),
  pch = 16,
  cex = 0.8,
  col = adjustcolor(results$plot_color, alpha.f = 0.7),
  xlab = "log2 Fold Change (Mutant / Control)",
  ylab = "-log10(P-value)",
  main = sprintf("Volcano Plot - Differential Stripes (%s)\n%d significant (FDR < %.2f)",
                 TIMEPOINT, n_sig, config$edger$fdr_primary)
)

abline(h = 0, col = "black", lty = 2)
abline(v = 0, col = "black", lty = 2)
abline(v = c(-0.3, 0.3), col = "gray40", lty = 3)

# Highlight by source
points(
  results$logFC[results$source == "control_only"],
  -log10(results$PValue[results$source == "control_only"]),
  pch = 1, cex = 1.2, col = "blue"
)
points(
  results$logFC[results$source == "mutant_only"],
  -log10(results$PValue[results$source == "mutant_only"]),
  pch = 1, cex = 1.2, col = "red"
)

legend(
  "topright",
  legend = c(
    sprintf("Up in mutant (%d)", n_up),
    sprintf("Down in mutant (%d)", n_down),
    "Not significant",
    "Control-only source",
    "Mutant-only source"
  ),
  col = c("firebrick3", "steelblue3", "gray60", "blue", "red"),
  pch = c(16, 16, 16, 1, 1),
  pt.cex = 1.5,
  bty = "n"
)

dev.off()
cat("  Saved: volcano_plot.pdf\n")

# ==============================================================================
# MA PLOT
# ==============================================================================
pdf(file.path(plots_dir, "ma_plot.pdf"), width = 10, height = 8)

par(mar = c(5, 5, 4, 2))
plot(
  results$logCPM,
  results$logFC,
  pch = 16,
  cex = 0.8,
  col = adjustcolor(results$plot_color, alpha.f = 0.7),
  xlab = "Average log2 CPM",
  ylab = "log2 Fold Change (Mutant / Control)",
  main = sprintf("MA Plot - Differential Stripes (%s)", TIMEPOINT)
)

abline(h = 0, col = "black", lty = 2)
abline(h = c(-0.3, 0.3), col = "gray40", lty = 3)

legend(
  "topright",
  legend = c(
    sprintf("Up in mutant (%d)", n_up),
    sprintf("Down in mutant (%d)", n_down),
    "Not significant"
  ),
  col = c("firebrick3", "steelblue3", "gray60"),
  pch = 16,
  pt.cex = 1.5,
  bty = "n"
)

dev.off()
cat("  Saved: ma_plot.pdf\n")

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

# Save DGEList object
saveRDS(dge, file.path(output_dir, "03_dge_object.rds"))
cat("  Saved: 03_dge_object.rds\n")

# Save fit object
saveRDS(fit, file.path(output_dir, "03_fit_object.rds"))
cat("  Saved: 03_fit_object.rds\n")

# Save all results
write.table(
  results,
  file.path(output_dir, "03_all_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
saveRDS(results, file.path(output_dir, "03_all_results.rds"))
cat("  Saved: 03_all_results.tsv/rds\n")

# Save significant results
sig_results <- results[results$significant_FDR05, ]
write.table(
  sig_results,
  file.path(output_dir, "03_significant_FDR05.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat(sprintf("  Saved: 03_significant_FDR05.tsv (%d stripes)\n", nrow(sig_results)))

# Save summary statistics
summary_stats <- list(
  timepoint = TIMEPOINT,
  n_stripes_total = nrow(results),
  n_significant_FDR05 = n_sig,
  n_significant_FDR10 = sum(results$significant_FDR10),
  n_up_in_mutant = n_up,
  n_down_in_mutant = n_down,
  common_dispersion = dge$common.dispersion,
  common_bcv = sqrt(dge$common.dispersion),
  median_tagwise_bcv = median(sqrt(dge$tagwise.dispersion)),
  normalization_method = config$edger$normalization_method,
  skip_filtering = config$edger$skip_filtering,
  fdr_threshold = config$edger$fdr_primary,
  date = Sys.Date()
)

saveRDS(summary_stats, file.path(output_dir, "03_summary_stats.rds"))

# Write summary text
summary_text <- sprintf("
========================================
Phase 3: edgeR Differential Analysis
Timepoint: %s
========================================

INPUT
-----
Total stripes: %d
Filtering: %s

NORMALIZATION
-------------
Method: %s
Normalization factors: %.4f - %.4f

DISPERSION
----------
Common BCV: %.3f
Median tagwise BCV: %.3f

RESULTS (FDR < %.2f)
--------------------
Significant stripes: %d (%.1f%%)
  Up in mutant: %d
  Down in mutant: %d

Exploratory (FDR < %.2f): %d

OUTPUT
------
%s

Next: Rscript scripts/phase4_integration.R %s
",
  TIMEPOINT,
  nrow(results),
  ifelse(config$edger$skip_filtering, "SKIPPED (all stripes retained)", "Applied"),
  config$edger$normalization_method,
  min(dge$samples$norm.factors), max(dge$samples$norm.factors),
  sqrt(dge$common.dispersion),
  median(sqrt(dge$tagwise.dispersion)),
  config$edger$fdr_primary,
  n_sig, 100 * n_sig / nrow(results),
  n_up, n_down,
  config$edger$fdr_exploratory, sum(results$significant_FDR10),
  output_dir,
  TIMEPOINT
)

writeLines(summary_text, file.path(output_dir, "03_summary.txt"))
cat("  Saved: 03_summary.txt\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n=================================\n")
cat("Phase 3 Complete\n")
cat("=================================\n\n")

cat(sprintf("Timepoint: %s\n", TIMEPOINT))
cat(sprintf("Stripes analyzed: %d\n", nrow(results)))
cat(sprintf("Significant (FDR < %.2f): %d (%.1f%%)\n",
            config$edger$fdr_primary, n_sig, 100 * n_sig / nrow(results)))
cat(sprintf("  Up in mutant: %d\n", n_up))
cat(sprintf("  Down in mutant: %d\n", n_down))
cat(sprintf("Common BCV: %.3f\n", sqrt(dge$common.dispersion)))

cat(sprintf("\nOutput: %s\n", output_dir))
cat(sprintf("\nNext: Rscript scripts/phase4_integration.R %s\n\n", TIMEPOINT))
