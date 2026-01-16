#!/usr/bin/env Rscript
# stripes/scripts/phase4_integration.R
# Phase 4: Integration & Direction Classification
#
# Merges Phase 1 annotations with Phase 3 statistics and classifies
# each stripe as lost/gained/strengthened/weakened/unchanged/ambiguous.
#
# Usage: Rscript phase4_integration.R <timepoint>
#   timepoint: "early" (250831) or "late" (250402)

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
})

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript phase4_integration.R <timepoint>\n",
       "  timepoint: 'early' or 'late'")
}
TIMEPOINT <- args[1]

if (!TIMEPOINT %in% c("early", "late")) {
  stop("Timepoint must be 'early' or 'late'")
}

cat("\n========================================\n")
cat("Phase 4: Integration & Classification\n")
cat(sprintf("Timepoint: %s\n", TIMEPOINT))
cat("========================================\n\n")

# ==============================================================================
# LOAD CONFIGURATION
# ==============================================================================
cat("Loading configuration...\n")

config_file <- "stripes/config/stripe_config.yaml"
if (!file.exists(config_file)) {
  stop(sprintf("Config file not found: %s", config_file))
}

config <- yaml::read_yaml(config_file)
cat(sprintf("  Config: %s\n", config_file))

# ==============================================================================
# PATH SETUP
# ==============================================================================
resolution_kb <- config$stripe_data$resolutions$primary / 1000
base_dir <- file.path(config$outputs$base_dir, TIMEPOINT)
edger_dir <- file.path(base_dir, "edgeR_results")

cat(sprintf("  Base dir: %s\n", base_dir))
cat(sprintf("  edgeR dir: %s\n", edger_dir))

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("\n=== Loading Data ===\n")

# Phase 1 unified stripes
phase1_file <- file.path(base_dir, "01_unified_stripes.tsv")
if (!file.exists(phase1_file)) {
  stop(sprintf("Phase 1 output not found: %s", phase1_file))
}
unified <- read.table(phase1_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
cat(sprintf("  Phase 1: %d stripes\n", nrow(unified)))

# Phase 3 edgeR results
phase3_file <- file.path(edger_dir, "03_all_results.tsv")
if (!file.exists(phase3_file)) {
  stop(sprintf("Phase 3 output not found: %s\n  Run phase3_edgeR.R first",
               phase3_file))
}
edger_results <- read.table(phase3_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
cat(sprintf("  Phase 3: %d stripes with edgeR statistics\n", nrow(edger_results)))

# ==============================================================================
# MERGE ANNOTATIONS AND STATISTICS
# ==============================================================================
cat("\n=== Merging Data ===\n")

# Select columns from edgeR results
edger_cols <- c("stripe_id", "logFC", "logCPM", "F", "PValue", "FDR",
                "significant_FDR05", "significant_FDR10", "direction_edgeR")

# Merge by stripe_id
final <- merge(
  unified,
  edger_results[, edger_cols],
  by = "stripe_id",
  all.x = TRUE
)

cat(sprintf("  Merged: %d stripes\n", nrow(final)))

# Check for unmatched stripes
unmatched <- sum(is.na(final$logFC))
if (unmatched > 0) {
  cat(sprintf("  WARNING: %d stripes have no edgeR statistics\n", unmatched))
}

# ==============================================================================
# DIRECTION CLASSIFICATION
# ==============================================================================
cat("\n=== Direction Classification ===\n")

# Get thresholds from config
fdr_thresh <- config$classification$fdr_threshold
logfc_thresh <- config$classification$logFC_threshold

cat(sprintf("  FDR threshold: %.2f\n", fdr_thresh))
cat(sprintf("  logFC threshold: %.2f\n", logfc_thresh))

# Tiered classification logic
# Detection (source) is PRIMARY evidence; FDR + direction provide confidence tiers
# This replaces the original FDR-dependent logic that produced 0 differential stripes

# Get tiered thresholds from config (or use defaults)
high_fdr_thresh <- if (!is.null(config$classification$tiered$high_fdr)) {
  config$classification$tiered$high_fdr
} else {
  0.10  # Default: FDR < 0.1 for high confidence
}
medium_logfc_thresh <- if (!is.null(config$classification$tiered$medium_logfc)) {
  config$classification$tiered$medium_logfc
} else {
  0.2  # Default: |logFC| > 0.2 for medium confidence
}

cat(sprintf("  Tiered thresholds: high_fdr=%.2f, medium_logfc=%.2f\n",
            high_fdr_thresh, medium_logfc_thresh))

# Direction classification: all condition-specific stripes get classified
final$direction <- case_when(
  # === LOST (control_only stripes) ===
  # All control_only stripes are "lost" - detection is primary evidence
  final$source == "control_only" ~ "lost",

  # === GAINED (mutant_only stripes) ===
  # All mutant_only stripes are "gained" - detection is primary evidence
  final$source == "mutant_only" ~ "gained",

  # === SHARED stripes ===
  # Strengthened: FDR significant + positive logFC above threshold
  final$source == "shared" & final$FDR < fdr_thresh & final$logFC > logfc_thresh ~ "strengthened",
  # Weakened: FDR significant + negative logFC below threshold
  final$source == "shared" & final$FDR < fdr_thresh & final$logFC < -logfc_thresh ~ "weakened",
  # Unchanged: shared without significant change
  final$source == "shared" ~ "unchanged"
)

# Add confidence tier as separate column
final$direction_confidence <- case_when(
  # Lost confidence tiers
  final$source == "control_only" & final$FDR < high_fdr_thresh & final$logFC < 0 ~ "high",
  final$source == "control_only" & final$logFC < -medium_logfc_thresh ~ "medium",
  final$source == "control_only" ~ "low",

  # Gained confidence tiers
  final$source == "mutant_only" & final$FDR < high_fdr_thresh & final$logFC > 0 ~ "high",
  final$source == "mutant_only" & final$logFC > medium_logfc_thresh ~ "medium",
  final$source == "mutant_only" ~ "low",

  # Shared stripes confidence (based on FDR only)
  final$source == "shared" & final$FDR < 0.05 ~ "high",
  final$source == "shared" & final$FDR < 0.10 ~ "medium",
  TRUE ~ "low"
)

# Flag directional consistency for condition-specific stripes
# TRUE = logFC direction matches expected for source
# FALSE = logFC direction opposite to expected
# NA = shared stripes (not applicable)
final$direction_consistent <- case_when(
  final$source == "control_only" & final$logFC < 0 ~ TRUE,
  final$source == "control_only" & final$logFC >= 0 ~ FALSE,
  final$source == "mutant_only" & final$logFC > 0 ~ TRUE,
  final$source == "mutant_only" & final$logFC <= 0 ~ FALSE,
  final$source == "shared" ~ NA
)

# Print direction summary
cat("\nDirection classification:\n")
print(table(final$direction))

# Cross-tabulate source vs direction
cat("\nSource x Direction:\n")
print(table(final$source, final$direction))

# ==============================================================================
# ADD DERIVED METRICS
# ==============================================================================
cat("\n=== Adding Derived Metrics ===\n")

# Stripe length (span_y2 - span_y1)
final$stripe_length <- final$span_y2 - final$span_y1

# Anchor width
final$anchor_width <- final$anchor_x2 - final$anchor_x1

# Effect size category
final$effect_category <- case_when(
  abs(final$logFC) > 1.0 ~ "strong",
  abs(final$logFC) > 0.5 ~ "moderate",
  abs(final$logFC) > 0.3 ~ "weak",
  TRUE ~ "minimal"
)

cat("  Added: stripe_length, anchor_width, effect_category\n")

# Summary of effect categories for significant stripes
cat("\nEffect categories (significant stripes):\n")
print(table(final$effect_category[final$significant_FDR05]))

# ==============================================================================
# CREATE BEDPE OUTPUT
# ==============================================================================
cat("\n=== Creating BEDPE Output ===\n")

# BEDPE format for Juicebox visualization
bedpe_df <- data.frame(
  chr1 = final$chr,
  x1 = final$anchor_x1,
  x2 = final$anchor_x2,
  chr2 = final$chr,
  y1 = final$span_y1,
  y2 = final$span_y2,
  name = final$stripe_id,
  score = -log10(final$FDR + 1e-300),  # Transform FDR to score
  strand1 = ".",
  strand2 = ".",
  # Color gradient by direction and confidence
  color = case_when(
    # Lost stripes: blue gradient (dark=high conf, light=low conf)
    final$direction == "lost" &
      final$direction_confidence == "high" ~ "0,0,255",
    final$direction == "lost" &
      final$direction_confidence == "medium" ~ "100,100,255",
    final$direction == "lost" &
      final$direction_confidence == "low" ~ "180,180,255",
    # Gained stripes: red gradient (dark=high conf, light=low conf)
    final$direction == "gained" &
      final$direction_confidence == "high" ~ "255,0,0",
    final$direction == "gained" &
      final$direction_confidence == "medium" ~ "255,100,100",
    final$direction == "gained" &
      final$direction_confidence == "low" ~ "255,180,180",
    # Shared stripes
    final$direction == "strengthened" ~ "255,165,0",     # Orange
    final$direction == "weakened" ~ "100,149,237",       # Cornflower blue
    final$direction == "unchanged" ~ "128,128,128",      # Gray
    TRUE ~ "200,200,200"
  ),
  logFC = final$logFC,
  FDR = final$FDR,
  direction = final$direction,
  source = final$source,
  stringsAsFactors = FALSE
)

cat(sprintf("  BEDPE created: %d entries\n", nrow(bedpe_df)))

# ==============================================================================
# GENERATE SUMMARY STATISTICS
# ==============================================================================
cat("\n=== Summary Statistics ===\n")

summary_stats <- list(
  timepoint = TIMEPOINT,
  total_stripes = nrow(final),

  # By source
  n_control_only = sum(final$source == "control_only"),
  n_mutant_only = sum(final$source == "mutant_only"),
  n_shared = sum(final$source == "shared"),

  # By direction (now all classified, no ambiguous)
  n_lost = sum(final$direction == "lost"),
  n_gained = sum(final$direction == "gained"),
  n_strengthened = sum(final$direction == "strengthened"),
  n_weakened = sum(final$direction == "weakened"),
  n_unchanged = sum(final$direction == "unchanged"),

  # Tiered confidence breakdown for lost/gained
  n_lost_high = sum(final$direction == "lost" &
                    final$direction_confidence == "high"),
  n_lost_medium = sum(final$direction == "lost" &
                      final$direction_confidence == "medium"),
  n_lost_low = sum(final$direction == "lost" &
                   final$direction_confidence == "low"),
  n_gained_high = sum(final$direction == "gained" &
                      final$direction_confidence == "high"),
  n_gained_medium = sum(final$direction == "gained" &
                        final$direction_confidence == "medium"),
  n_gained_low = sum(final$direction == "gained" &
                     final$direction_confidence == "low"),

  # Directional consistency metrics
  n_direction_consistent = sum(final$direction_consistent, na.rm = TRUE),
  n_direction_inconsistent = sum(!final$direction_consistent, na.rm = TRUE),
  pct_direction_consistent = round(
    mean(final$direction_consistent, na.rm = TRUE) * 100, 1
  ),

  # By edgeR significance
  n_significant_FDR05 = sum(final$significant_FDR05),
  n_significant_FDR10 = sum(final$significant_FDR10),

  # By Phase 1 confidence (detection quality)
  n_high_detection_conf = sum(final$confidence == "high"),
  n_medium_detection_conf = sum(final$confidence == "medium"),
  n_low_detection_conf = sum(final$confidence == "low"),

  # Fold change stats
  median_logFC = median(final$logFC, na.rm = TRUE),
  median_logFC_sig = median(final$logFC[final$significant_FDR05], na.rm = TRUE),
  max_logFC = max(final$logFC, na.rm = TRUE),
  min_logFC = min(final$logFC, na.rm = TRUE),

  # Stripe geometry
  median_stripe_length = median(final$stripe_length),
  median_anchor_width = median(final$anchor_width),

  # Parameters used
  fdr_threshold = fdr_thresh,
  logfc_threshold = logfc_thresh,
  high_fdr_threshold = high_fdr_thresh,
  medium_logfc_threshold = medium_logfc_thresh,

  date = Sys.Date()
)

# Print key stats
cat(sprintf("\nTotal stripes: %d\n", summary_stats$total_stripes))
cat(sprintf("\nBy source:\n"))
cat(sprintf("  control_only: %d\n", summary_stats$n_control_only))
cat(sprintf("  mutant_only: %d\n", summary_stats$n_mutant_only))
cat(sprintf("  shared: %d\n", summary_stats$n_shared))
cat(sprintf("\nBy direction (tiered confidence):\n"))
cat(sprintf("  LOST: %d total\n", summary_stats$n_lost))
cat(sprintf("    - High:   %d (FDR<0.1 + direction match)\n",
            summary_stats$n_lost_high))
cat(sprintf("    - Medium: %d (direction match, logFC<-0.2)\n",
            summary_stats$n_lost_medium))
cat(sprintf("    - Low:    %d (detection only)\n", summary_stats$n_lost_low))
cat(sprintf("  GAINED: %d total\n", summary_stats$n_gained))
cat(sprintf("    - High:   %d (FDR<0.1 + direction match)\n",
            summary_stats$n_gained_high))
cat(sprintf("    - Medium: %d (direction match, logFC>0.2)\n",
            summary_stats$n_gained_medium))
cat(sprintf("    - Low:    %d (detection only)\n", summary_stats$n_gained_low))
cat(sprintf("  strengthened: %d\n", summary_stats$n_strengthened))
cat(sprintf("  weakened: %d\n", summary_stats$n_weakened))
cat(sprintf("  unchanged: %d\n", summary_stats$n_unchanged))
cat(sprintf("\nDirectional consistency: %.1f%% (%d/%d)\n",
            summary_stats$pct_direction_consistent,
            summary_stats$n_direction_consistent,
            summary_stats$n_direction_consistent +
              summary_stats$n_direction_inconsistent))

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

output_dir <- base_dir

# Final integrated results (TSV and RDS)
final_tsv <- file.path(output_dir, "04_final_differential_stripes.tsv")
final_rds <- file.path(output_dir, "04_final_differential_stripes.rds")

write.table(final, final_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(final, final_rds)
cat(sprintf("  Saved: 04_final_differential_stripes.tsv/rds (%d stripes)\n",
            nrow(final)))

# BEDPE for Juicebox
bedpe_file <- file.path(output_dir, "04_final_differential_stripes.bedpe")
write.table(bedpe_df, bedpe_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
cat(sprintf("  Saved: 04_final_differential_stripes.bedpe\n"))

# Summary statistics
saveRDS(summary_stats, file.path(output_dir, "04_summary_stats.rds"))

# Summary text file
summary_text <- sprintf("
================================================================================
DIFFERENTIAL STRIPE ANALYSIS - FINAL RESULTS (Tiered Classification)
================================================================================
Timepoint: %s (%s)
Date: %s

STRIPE COUNTS
-------------
Total stripes analyzed: %d

By detection source:
  control_only:  %3d (detected in ctrl merged only)
  mutant_only:   %3d (detected in mut merged only)
  shared:        %3d (detected in both conditions)

DIFFERENTIAL DIRECTION (Tiered Confidence)
------------------------------------------
LOST (control_only stripes): %d total
  - High confidence:   %3d (FDR < %.2f + logFC < 0)
  - Medium confidence: %3d (logFC < -%.1f, any FDR)
  - Low confidence:    %3d (detection only)

GAINED (mutant_only stripes): %d total
  - High confidence:   %3d (FDR < %.2f + logFC > 0)
  - Medium confidence: %3d (logFC > %.1f, any FDR)
  - Low confidence:    %3d (detection only)

SHARED stripes:
  strengthened:  %3d (FDR < %.2f, logFC > %.1f)
  weakened:      %3d (FDR < %.2f, logFC < -%.1f)
  unchanged:     %3d

DIRECTIONAL CONSISTENCY
-----------------------
Consistent (logFC matches source): %d (%.1f%%)
Inconsistent:                      %d

Note: Directional consistency measures whether quantitative logFC
agrees with detection-based classification. Low consistency (~50%%)
may indicate noisy detection or weak biological signal.

edgeR SIGNIFICANCE
------------------
FDR < %.2f: %d stripes
FDR < %.2f: %d stripes

DETECTION CONFIDENCE (from Phase 1)
------------------------------------
High:   %d (replicate support + 10kb validation)
Medium: %d
Low:    %d

FOLD CHANGE STATISTICS
----------------------
Median logFC (all): %.3f
Median logFC (sig): %.3f
Range: %.3f to %.3f

STRIPE GEOMETRY
---------------
Median stripe length: %.0f bp
Median anchor width: %.0f bp

OUTPUTS
-------
%s/04_final_differential_stripes.tsv  - Full results
%s/04_final_differential_stripes.rds  - R object
%s/04_final_differential_stripes.bedpe - Juicebox visualization

================================================================================
",
  TIMEPOINT,
  config$stripe_data$timepoints[[TIMEPOINT]],
  Sys.Date(),

  summary_stats$total_stripes,

  summary_stats$n_control_only,
  summary_stats$n_mutant_only,
  summary_stats$n_shared,

  # LOST section
  summary_stats$n_lost,
  summary_stats$n_lost_high, high_fdr_thresh,
  summary_stats$n_lost_medium, medium_logfc_thresh,
  summary_stats$n_lost_low,

  # GAINED section
  summary_stats$n_gained,
  summary_stats$n_gained_high, high_fdr_thresh,
  summary_stats$n_gained_medium, medium_logfc_thresh,
  summary_stats$n_gained_low,

  # SHARED section
  summary_stats$n_strengthened, fdr_thresh, logfc_thresh,
  summary_stats$n_weakened, fdr_thresh, logfc_thresh,
  summary_stats$n_unchanged,

  # Directional consistency
  summary_stats$n_direction_consistent, summary_stats$pct_direction_consistent,
  summary_stats$n_direction_inconsistent,

  # edgeR significance
  fdr_thresh, summary_stats$n_significant_FDR05,
  config$edger$fdr_exploratory, summary_stats$n_significant_FDR10,

  # Detection confidence
  summary_stats$n_high_detection_conf,
  summary_stats$n_medium_detection_conf,
  summary_stats$n_low_detection_conf,

  # Fold change stats
  summary_stats$median_logFC,
  summary_stats$median_logFC_sig,
  summary_stats$min_logFC, summary_stats$max_logFC,

  # Geometry
  summary_stats$median_stripe_length,
  summary_stats$median_anchor_width,

  # Output paths
  output_dir, output_dir, output_dir
)

writeLines(summary_text, file.path(output_dir, "04_summary.txt"))
cat("  Saved: 04_summary.txt\n")

# Save differential stripes by direction
for (dir_cat in c("lost", "gained", "strengthened", "weakened")) {
  subset_df <- final[final$direction == dir_cat, ]
  if (nrow(subset_df) > 0) {
    subset_file <- file.path(output_dir,
                             sprintf("04_stripes_%s.tsv", dir_cat))
    write.table(subset_df, subset_file, sep = "\t",
                quote = FALSE, row.names = FALSE)
    cat(sprintf("  Saved: 04_stripes_%s.tsv (%d stripes)\n",
                dir_cat, nrow(subset_df)))
  }
}

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n=================================\n")
cat("Phase 4 Complete\n")
cat("=================================\n\n")

cat(sprintf("Timepoint: %s\n", TIMEPOINT))
cat(sprintf("Total stripes: %d\n", nrow(final)))
cat(sprintf("\nDifferential stripes (tiered):\n"))
cat(sprintf("  Lost: %d (H:%d M:%d L:%d)\n",
            summary_stats$n_lost, summary_stats$n_lost_high,
            summary_stats$n_lost_medium, summary_stats$n_lost_low))
cat(sprintf("  Gained: %d (H:%d M:%d L:%d)\n",
            summary_stats$n_gained, summary_stats$n_gained_high,
            summary_stats$n_gained_medium, summary_stats$n_gained_low))
cat(sprintf("  Strengthened: %d\n", summary_stats$n_strengthened))
cat(sprintf("  Weakened: %d\n", summary_stats$n_weakened))
cat(sprintf("\nDirectional consistency: %.1f%%\n",
            summary_stats$pct_direction_consistent))

cat(sprintf("\nOutput: %s\n", output_dir))
cat("\nPipeline complete! Review results in:\n")
cat(sprintf("  %s/04_final_differential_stripes.tsv\n", output_dir))
cat(sprintf("  %s/04_final_differential_stripes.bedpe (for Juicebox)\n\n", output_dir))
