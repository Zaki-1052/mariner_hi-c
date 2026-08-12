#!/usr/bin/env Rscript
# stripes/stripenn/scripts/05_integration.R
# Stage 5: Merge union annotations with edgeR statistics, classify direction.
# Ported from stripes/scripts/phase4_integration.R.
#
# Usage: Rscript 05_integration.R <timepoint> <resolution_bp>
#   <timepoint>     : 250831 | 250402
#   <resolution_bp> : 5000 | 10000

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05_integration.R <timepoint> <resolution_bp>")
}
TIMEPOINT <- args[1]
RES <- as.numeric(args[2])
RES_KB <- RES / 1000

CODE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"

cat("\n========================================\n")
cat("Stage 5: Integration & Classification\n")
cat(sprintf("Timepoint:  %s\n", TIMEPOINT))
cat(sprintf("Resolution: %d bp (%g kb)\n", RES, RES_KB))
cat("========================================\n\n")

config_file <- file.path(CODE_DIR, "config", "stripenn_config.yaml")
config <- yaml::read_yaml(config_file)

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("=== Loading Data ===\n")

base_dir <- file.path(DATA_DIR, "outputs", TIMEPOINT, paste0("res_", RES_KB, "kb"))

union_file <- file.path(base_dir, "union_stripes.tsv")
if (!file.exists(union_file)) {
  stop(sprintf("Union stripes not found: %s", union_file))
}
unified <- read.delim(union_file, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Union stripes: %d\n", nrow(unified)))

edger_file <- file.path(base_dir, "04_edgeR", "all_results.tsv")
if (!file.exists(edger_file)) {
  stop(sprintf("edgeR results not found: %s\n  Run 04_edgeR.R first", edger_file))
}
edger_results <- read.delim(edger_file, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  edgeR results: %d stripes\n", nrow(edger_results)))

# ==============================================================================
# MERGE
# ==============================================================================
cat("\n=== Merging Data ===\n")

edger_cols <- c("stripe_id", "logFC", "logCPM", "F", "PValue", "FDR",
                "significant_FDR05", "significant_FDR10", "direction_edgeR")

edger_available <- intersect(edger_cols, colnames(edger_results))

final <- merge(
  unified,
  edger_results[, edger_available],
  by = "stripe_id",
  all.x = TRUE
)

cat(sprintf("  Merged: %d stripes\n", nrow(final)))

unmatched <- sum(is.na(final$logFC))
if (unmatched > 0) {
  cat(sprintf("  WARNING: %d stripes have no edgeR statistics\n", unmatched))
}

# ==============================================================================
# DIRECTION CLASSIFICATION
# ==============================================================================
cat("\n=== Direction Classification ===\n")

fdr_thresh <- config$classification$fdr_threshold
logfc_thresh <- config$classification$logFC_threshold
high_fdr_thresh <- config$classification$tiered$high_fdr
medium_logfc_thresh <- config$classification$tiered$medium_logfc

cat(sprintf("  FDR threshold:    %.2f\n", fdr_thresh))
cat(sprintf("  logFC threshold:  %.2f\n", logfc_thresh))
cat(sprintf("  Tiered: high_fdr=%.2f, medium_logfc=%.2f\n",
            high_fdr_thresh, medium_logfc_thresh))

final$direction <- case_when(
  final$source == "control_only" ~ "lost",
  final$source == "mutant_only" ~ "gained",
  final$source == "shared" & final$FDR < fdr_thresh & final$logFC > logfc_thresh ~ "strengthened",
  final$source == "shared" & final$FDR < fdr_thresh & final$logFC < -logfc_thresh ~ "weakened",
  final$source == "shared" ~ "unchanged"
)

final$direction_confidence <- case_when(
  final$source == "control_only" & final$FDR < high_fdr_thresh & final$logFC < 0 ~ "high",
  final$source == "control_only" & final$logFC < -medium_logfc_thresh ~ "medium",
  final$source == "control_only" ~ "low",

  final$source == "mutant_only" & final$FDR < high_fdr_thresh & final$logFC > 0 ~ "high",
  final$source == "mutant_only" & final$logFC > medium_logfc_thresh ~ "medium",
  final$source == "mutant_only" ~ "low",

  final$source == "shared" & final$FDR < 0.05 ~ "high",
  final$source == "shared" & final$FDR < 0.10 ~ "medium",
  TRUE ~ "low"
)

final$direction_consistent <- case_when(
  final$source == "control_only" & final$logFC < 0 ~ TRUE,
  final$source == "control_only" & final$logFC >= 0 ~ FALSE,
  final$source == "mutant_only" & final$logFC > 0 ~ TRUE,
  final$source == "mutant_only" & final$logFC <= 0 ~ FALSE,
  final$source == "shared" ~ NA
)

cat("\nDirection classification:\n")
print(table(final$direction))

cat("\nSource x Direction:\n")
print(table(final$source, final$direction))

# ==============================================================================
# DERIVED METRICS
# ==============================================================================
cat("\n=== Adding Derived Metrics ===\n")

final$effect_category <- case_when(
  abs(final$logFC) > 1.0 ~ "strong",
  abs(final$logFC) > 0.5 ~ "moderate",
  abs(final$logFC) > 0.3 ~ "weak",
  TRUE ~ "minimal"
)

# ==============================================================================
# BEDPE OUTPUT (Juicebox)
# ==============================================================================
cat("\n=== Creating BEDPE Output ===\n")

bedpe_df <- data.frame(
  chr1 = final$chr,
  x1 = final$pos1,
  x2 = final$pos2,
  chr2 = final$chr,
  y1 = final$pos3,
  y2 = final$pos4,
  name = final$stripe_id,
  score = -log10(final$FDR + 1e-300),
  strand1 = ".",
  strand2 = ".",
  color = case_when(
    final$direction == "lost" & final$direction_confidence == "high" ~ "0,0,255",
    final$direction == "lost" & final$direction_confidence == "medium" ~ "100,100,255",
    final$direction == "lost" ~ "180,180,255",
    final$direction == "gained" & final$direction_confidence == "high" ~ "255,0,0",
    final$direction == "gained" & final$direction_confidence == "medium" ~ "255,100,100",
    final$direction == "gained" ~ "255,180,180",
    final$direction == "strengthened" ~ "255,165,0",
    final$direction == "weakened" ~ "100,149,237",
    final$direction == "unchanged" ~ "128,128,128",
    TRUE ~ "200,200,200"
  ),
  logFC = final$logFC,
  FDR = final$FDR,
  direction = final$direction,
  source = final$source,
  stringsAsFactors = FALSE
)

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

write.table(final, file.path(base_dir, "05_final_differential.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(final, file.path(base_dir, "05_final_differential.rds"))
cat(sprintf("  Saved: 05_final_differential.tsv/rds (%d stripes)\n", nrow(final)))

write.table(bedpe_df, file.path(base_dir, "05_final_differential.bedpe"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 05_final_differential.bedpe\n")

for (dir_cat in c("lost", "gained", "strengthened", "weakened", "unchanged")) {
  subset <- final[final$direction == dir_cat, ]
  if (nrow(subset) > 0) {
    outfile <- file.path(base_dir, paste0("05_stripes_", dir_cat, ".tsv"))
    write.table(subset, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("  Saved: 05_stripes_%s.tsv (%d stripes)\n", dir_cat, nrow(subset)))
  }
}

# ==============================================================================
# SUMMARY
# ==============================================================================
n_lost <- sum(final$direction == "lost")
n_gained <- sum(final$direction == "gained")
n_str <- sum(final$direction == "strengthened")
n_wk <- sum(final$direction == "weakened")
n_unch <- sum(final$direction == "unchanged")
n_consistent <- sum(final$direction_consistent, na.rm = TRUE)
n_inconsistent <- sum(!final$direction_consistent, na.rm = TRUE)
pct_consistent <- if ((n_consistent + n_inconsistent) > 0) {
  round(100 * n_consistent / (n_consistent + n_inconsistent), 1)
} else { NA }

summary_text <- sprintf("
========================================
Stage 5: Integration & Classification (Stripenn)
Timepoint:  %s
Resolution: %gkb
========================================

TOTAL: %d stripes

BY DIRECTION:
  Lost:          %d
  Gained:        %d
  Strengthened:  %d
  Weakened:      %d
  Unchanged:     %d

CONFIDENCE TIERS (Lost):
  High:   %d   Medium: %d   Low: %d

CONFIDENCE TIERS (Gained):
  High:   %d   Medium: %d   Low: %d

DIRECTIONAL CONSISTENCY:
  Consistent:   %d
  Inconsistent: %d
  Pct consistent: %s%%

EFFECT SIZE (significant stripes):
  Strong (|logFC|>1):   %d
  Moderate (0.5-1):     %d
  Weak (0.3-0.5):       %d
  Minimal (<0.3):       %d

MEDIAN logFC: %.3f (all), %.3f (significant)
",
  TIMEPOINT, RES_KB,
  nrow(final),
  n_lost, n_gained, n_str, n_wk, n_unch,
  sum(final$direction == "lost" & final$direction_confidence == "high"),
  sum(final$direction == "lost" & final$direction_confidence == "medium"),
  sum(final$direction == "lost" & final$direction_confidence == "low"),
  sum(final$direction == "gained" & final$direction_confidence == "high"),
  sum(final$direction == "gained" & final$direction_confidence == "medium"),
  sum(final$direction == "gained" & final$direction_confidence == "low"),
  n_consistent, n_inconsistent,
  ifelse(is.na(pct_consistent), "N/A", as.character(pct_consistent)),
  sum(final$effect_category == "strong" & final$significant_FDR05, na.rm = TRUE),
  sum(final$effect_category == "moderate" & final$significant_FDR05, na.rm = TRUE),
  sum(final$effect_category == "weak" & final$significant_FDR05, na.rm = TRUE),
  sum(final$effect_category == "minimal" & final$significant_FDR05, na.rm = TRUE),
  median(final$logFC, na.rm = TRUE),
  median(final$logFC[final$significant_FDR05], na.rm = TRUE)
)

writeLines(summary_text, file.path(base_dir, "05_summary.txt"))
cat("  Saved: 05_summary.txt\n")

saveRDS(list(
  timepoint = TIMEPOINT,
  resolution = RES,
  total = nrow(final),
  n_lost = n_lost, n_gained = n_gained,
  n_strengthened = n_str, n_weakened = n_wk, n_unchanged = n_unch,
  pct_consistent = pct_consistent,
  date = Sys.Date()
), file.path(base_dir, "05_summary_stats.rds"))

cat("\n========================================\n")
cat("Stage 5 complete.\n")
cat(sprintf("End: %s\n", Sys.time()))
cat("========================================\n")
