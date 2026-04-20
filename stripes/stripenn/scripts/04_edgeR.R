#!/usr/bin/env Rscript
# stripes/stripenn/scripts/04_edgeR.R
# Stage 4: edgeR differential stripe analysis using O_Sum_added from
# stripenn score output. TMM + robust QL-GLM with skip-filtering.
# Ported from stripes/scripts/phase3_edgeR.R.
#
# Usage: Rscript 04_edgeR.R <timepoint> <resolution_bp>
#   <timepoint>     : 250831 | 250402
#   <resolution_bp> : 5000 | 10000

suppressPackageStartupMessages({
  library(edgeR)
  library(yaml)
  library(ggplot2)
  library(dplyr)
  library(svglite)
})

save_multiformat <- function(plot_code, base_path, width, height, dpi = 300) {
  pdf(paste0(base_path, ".pdf"), width = width, height = height)
  tryCatch(eval(plot_code), finally = dev.off())
  svglite(paste0(base_path, ".svg"), width = width, height = height)
  tryCatch(eval(plot_code), finally = dev.off())
  jpeg(paste0(base_path, ".jpg"), width = width * dpi, height = height * dpi,
       res = dpi, quality = 95)
  tryCatch(eval(plot_code), finally = dev.off())
  cat(sprintf("  Saved: %s.{pdf,svg,jpg}\n", basename(base_path)))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_edgeR.R <timepoint> <resolution_bp>")
}
TIMEPOINT <- args[1]
RES <- as.numeric(args[2])
RES_KB <- RES / 1000

CODE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"

cat("\n========================================\n")
cat("Stage 4: edgeR Differential Analysis\n")
cat(sprintf("Timepoint:  %s\n", TIMEPOINT))
cat(sprintf("Resolution: %d bp (%g kb)\n", RES, RES_KB))
cat("========================================\n\n")

config_file <- file.path(CODE_DIR, "config", "stripenn_config.yaml")
config <- yaml::read_yaml(config_file)
cat(sprintf("Config loaded: %s\n", config_file))

set.seed(config$stripenn_compute$seed)

base_dir <- file.path(DATA_DIR, "outputs", TIMEPOINT, paste0("res_", RES_KB, "kb"))
scores_dir <- file.path(base_dir, "scores")
output_dir <- file.path(base_dir, "04_edgeR")
plots_dir <- file.path(output_dir, "plots")

for (dir in c(output_dir, plots_dir)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# ==============================================================================
# LOAD DATA — Build count matrix from score files
# ==============================================================================
cat("\n=== Loading Score Data ===\n")

sample_names <- config$samples$replicates
groups <- config$samples$groups

union_file <- file.path(base_dir, "union_stripes.tsv")
if (!file.exists(union_file)) {
  stop(sprintf("Union stripes not found: %s\n  Run 02_build_union.R first", union_file))
}
union_df <- read.delim(union_file, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Union stripes: %d\n", nrow(union_df)))

count_matrix <- matrix(NA_real_, nrow = nrow(union_df), ncol = length(sample_names),
                       dimnames = list(union_df$stripe_id, sample_names))

for (sample in sample_names) {
  score_file <- file.path(scores_dir, paste0(sample, ".scores.tsv"))
  if (!file.exists(score_file)) {
    stop(sprintf("Score file missing: %s", score_file))
  }
  scores <- read.delim(score_file, header = TRUE, stringsAsFactors = FALSE)

  if (nrow(scores) != nrow(union_df)) {
    stop(sprintf("Row count mismatch for %s: %d (scores) vs %d (union)",
                 sample, nrow(scores), nrow(union_df)))
  }

  if (!"O_Sum_added" %in% colnames(scores)) {
    stop(sprintf("O_Sum_added column missing in %s", score_file))
  }

  count_matrix[, sample] <- round(as.numeric(scores$O_Sum_added))
  cat(sprintf("  %s: loaded (median O_Sum=%.0f)\n",
              sample, median(scores$O_Sum_added, na.rm = TRUE)))
}

na_count <- sum(is.na(count_matrix))
if (na_count > 0) {
  cat(sprintf("  WARNING: %d NA values in count matrix, setting to 0\n", na_count))
  count_matrix[is.na(count_matrix)] <- 0
}

neg_count <- sum(count_matrix < 0)
if (neg_count > 0) {
  cat(sprintf("  WARNING: %d negative values in count matrix, setting to 0\n", neg_count))
  count_matrix[count_matrix < 0] <- 0
}

cat(sprintf("\n  Count matrix: %d stripes x %d samples\n",
            nrow(count_matrix), ncol(count_matrix)))

# ==============================================================================
# CHROMOSOME FILTERING
# ==============================================================================
exclude_chroms <- config$filtering$exclude_chromosomes
if (length(exclude_chroms) > 0) {
  keep_rows <- !union_df$chr %in% exclude_chroms
  n_removed <- sum(!keep_rows)
  if (n_removed > 0) {
    cat(sprintf("\n  Removing %d stripes on excluded chromosomes (%s)\n",
                n_removed, paste(exclude_chroms, collapse = ", ")))
    count_matrix <- count_matrix[keep_rows, , drop = FALSE]
    union_df <- union_df[keep_rows, ]
  }
}

cat(sprintf("  After filtering: %d stripes\n", nrow(count_matrix)))

# ==============================================================================
# GENE ANNOTATION
# ==============================================================================
genes_df <- data.frame(
  stripe_id = union_df$stripe_id,
  chr = union_df$chr,
  pos1 = union_df$pos1,
  pos2 = union_df$pos2,
  pos3 = union_df$pos3,
  pos4 = union_df$pos4,
  direction_type = union_df$direction_type,
  source = union_df$source,
  pval_ctrl = union_df$pval_ctrl,
  pval_mut = union_df$pval_mut,
  stripe_length = union_df$stripe_length,
  stringsAsFactors = FALSE
)
genes_df$coord_string <- paste0(
  genes_df$chr, ":", genes_df$pos1, "-", genes_df$pos2,
  "->", genes_df$pos3, "-", genes_df$pos4
)
rownames(genes_df) <- genes_df$stripe_id

# ==============================================================================
# CREATE DGEList
# ==============================================================================
cat("\n=== Creating DGEList ===\n")

group <- factor(groups, levels = c("ctrl", "mut"))

dge <- DGEList(
  counts = count_matrix,
  group = group,
  genes = genes_df
)
dge$samples$sample_name <- sample_names

cat(sprintf("  DGEList: %d stripes x %d samples\n", nrow(dge), ncol(dge)))
cat(sprintf("  Groups: %s\n", paste(levels(group), collapse = ", ")))

cat("\n  Library sizes:\n")
for (i in seq_len(ncol(dge))) {
  cat(sprintf("    %s (%s): %.0f\n",
              dge$samples$sample_name[i],
              as.character(dge$samples$group[i]),
              dge$samples$lib.size[i]))
}

# ==============================================================================
# FILTERING
# ==============================================================================
cat("\n=== Filtering ===\n")

if (config$edger$skip_filtering) {
  cat("  SKIPPING low-count filtering (per design decision)\n")
  cat(sprintf("  Retaining ALL %d stripes for analysis\n", nrow(dge)))
} else {
  keep <- filterByExpr(dge, group = group)
  cat(sprintf("  Stripes passing filter: %d / %d\n", sum(keep), length(keep)))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
}

# ==============================================================================
# NORMALIZATION
# ==============================================================================
cat("\n=== Normalization ===\n")

dge <- calcNormFactors(dge, method = config$edger$normalization_method)

cat(sprintf("  Method: %s\n", config$edger$normalization_method))
cat("  Normalization factors:\n")
for (i in seq_len(ncol(dge))) {
  cat(sprintf("    %s: %.4f\n",
              dge$samples$sample_name[i],
              dge$samples$norm.factors[i]))
}

# ==============================================================================
# MDS PLOT
# ==============================================================================
cat("\n=== Generating MDS Plot ===\n")

save_multiformat(quote({
  plotMDS(
    dge,
    col = c(rep("blue", 3), rep("red", 3)),
    pch = 16, cex = 2,
    main = sprintf("MDS Plot - Stripenn Stripes (%s, %gkb)", TIMEPOINT, RES_KB),
    labels = dge$samples$sample_name
  )
  legend("topright",
         legend = c("Control", "BAP1-KO"),
         col = c("blue", "red"),
         pch = 16, pt.cex = 2, bty = "n")
}), file.path(plots_dir, "mds_plot"), width = 8, height = 6)

# ==============================================================================
# DISPERSION ESTIMATION
# ==============================================================================
cat("\n=== Dispersion Estimation ===\n")

design <- model.matrix(~group, data = dge$samples)
colnames(design) <- c("Intercept", "MutantEffect")

dge <- estimateDisp(dge, design, robust = config$edger$robust_dispersion)

cat(sprintf("  Common dispersion: %.4f\n", dge$common.dispersion))
cat(sprintf("  Common BCV: %.3f\n", sqrt(dge$common.dispersion)))
cat(sprintf("  Tagwise BCV range: %.3f - %.3f\n",
            min(sqrt(dge$tagwise.dispersion)),
            max(sqrt(dge$tagwise.dispersion))))
cat(sprintf("  Median tagwise BCV: %.3f\n",
            median(sqrt(dge$tagwise.dispersion))))

save_multiformat(quote({
  plotBCV(dge, main = sprintf("BCV Plot - Stripenn Stripes (%s, %gkb)", TIMEPOINT, RES_KB))
}), file.path(plots_dir, "bcv_plot"), width = 8, height = 6)

# ==============================================================================
# QUASI-LIKELIHOOD GLM
# ==============================================================================
cat("\n=== Quasi-Likelihood GLM ===\n")

fit <- glmQLFit(dge, design, robust = config$edger$robust_dispersion)

save_multiformat(quote({
  plotQLDisp(fit, main = sprintf("QL Dispersion - Stripenn Stripes (%s, %gkb)", TIMEPOINT, RES_KB))
}), file.path(plots_dir, "ql_dispersion_plot"), width = 8, height = 6)

# ==============================================================================
# DIFFERENTIAL TESTING
# ==============================================================================
cat("\n=== Differential Testing ===\n")

qlf <- glmQLFTest(fit, coef = 2)
results <- topTags(qlf, n = Inf, sort.by = "none")$table

results$significant_FDR05 <- results$FDR < config$edger$fdr_primary
results$significant_FDR10 <- results$FDR < config$edger$fdr_exploratory

results$direction_edgeR <- "unchanged"
results$direction_edgeR[results$significant_FDR05 & results$logFC > 0] <- "up_in_mutant"
results$direction_edgeR[results$significant_FDR05 & results$logFC < 0] <- "down_in_mutant"

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

results$plot_color <- "gray60"
results$plot_color[results$significant_FDR05 & results$logFC > 0] <- "firebrick3"
results$plot_color[results$significant_FDR05 & results$logFC < 0] <- "steelblue3"

save_multiformat(quote({
  par(mar = c(5, 5, 4, 2))
  plot(
    results$logFC, -log10(results$PValue),
    pch = 16, cex = 0.8,
    col = adjustcolor(results$plot_color, alpha.f = 0.7),
    xlab = "log2 Fold Change (Mutant / Control)",
    ylab = "-log10(P-value)",
    main = sprintf("Volcano - Stripenn Stripes (%s, %gkb)\n%d significant (FDR < %.2f)",
                   TIMEPOINT, RES_KB, n_sig, config$edger$fdr_primary)
  )
  abline(h = 0, col = "black", lty = 2)
  abline(v = 0, col = "black", lty = 2)
  abline(v = c(-0.3, 0.3), col = "gray40", lty = 3)

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

  legend("topright",
         legend = c(sprintf("Up in mutant (%d)", n_up),
                    sprintf("Down in mutant (%d)", n_down),
                    "Not significant",
                    "Control-only source",
                    "Mutant-only source"),
         col = c("firebrick3", "steelblue3", "gray60", "blue", "red"),
         pch = c(16, 16, 16, 1, 1), pt.cex = 1.5, bty = "n")
}), file.path(plots_dir, "volcano_plot"), width = 10, height = 8)

# MA plot
save_multiformat(quote({
  par(mar = c(5, 5, 4, 2))
  plot(
    results$logCPM, results$logFC,
    pch = 16, cex = 0.8,
    col = adjustcolor(results$plot_color, alpha.f = 0.7),
    xlab = "Average log2 CPM",
    ylab = "log2 Fold Change (Mutant / Control)",
    main = sprintf("MA Plot - Stripenn Stripes (%s, %gkb)", TIMEPOINT, RES_KB)
  )
  abline(h = 0, col = "black", lty = 2)
  abline(h = c(-0.3, 0.3), col = "gray40", lty = 3)
  legend("topright",
         legend = c(sprintf("Up in mutant (%d)", n_up),
                    sprintf("Down in mutant (%d)", n_down),
                    "Not significant"),
         col = c("firebrick3", "steelblue3", "gray60"),
         pch = 16, pt.cex = 1.5, bty = "n")
}), file.path(plots_dir, "ma_plot"), width = 10, height = 8)

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

saveRDS(dge, file.path(output_dir, "dge_object.rds"))
cat("  Saved: dge_object.rds\n")

saveRDS(fit, file.path(output_dir, "fit_object.rds"))
cat("  Saved: fit_object.rds\n")

write.table(results, file.path(output_dir, "all_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(results, file.path(output_dir, "all_results.rds"))
cat("  Saved: all_results.tsv/rds\n")

sig_results <- results[results$significant_FDR05, ]
write.table(sig_results, file.path(output_dir, "significant_FDR05.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: significant_FDR05.tsv (%d stripes)\n", nrow(sig_results)))

# Save count matrix for reproducibility
write.table(count_matrix, file.path(output_dir, "count_matrix.tsv"),
            sep = "\t", quote = FALSE)
cat("  Saved: count_matrix.tsv\n")

summary_text <- sprintf("
========================================
Stage 4: edgeR Differential Analysis (Stripenn)
Timepoint:  %s
Resolution: %gkb
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
",
  TIMEPOINT, RES_KB,
  nrow(results),
  ifelse(config$edger$skip_filtering, "SKIPPED (all stripes retained)", "Applied"),
  config$edger$normalization_method,
  min(dge$samples$norm.factors), max(dge$samples$norm.factors),
  sqrt(dge$common.dispersion),
  median(sqrt(dge$tagwise.dispersion)),
  config$edger$fdr_primary,
  n_sig, 100 * n_sig / nrow(results),
  n_up, n_down,
  config$edger$fdr_exploratory, sum(results$significant_FDR10)
)

writeLines(summary_text, file.path(output_dir, "summary.txt"))
cat("  Saved: summary.txt\n")

cat("\n========================================\n")
cat("Stage 4 complete.\n")
cat(sprintf("End: %s\n", Sys.time()))
cat("========================================\n")
