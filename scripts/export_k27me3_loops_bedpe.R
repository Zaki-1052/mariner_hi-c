# scripts/export_k27me3_loops_bedpe.R
# Export K27me3-K27me3 loops as BEDPE files
#
# Outputs:
#   1. ALL K27me3-K27me3 loops (both anchors H3K27me3+) for early and late
#   2. Differential K27me3-K27me3 loops (from polycomb_shared_anchor analysis)
#
# Usage:
#   Rscript scripts/export_k27me3_loops_bedpe.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

cat("=== Export K27me3-K27me3 Loops as BEDPE ===\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Configuration
INPUT_FILES <- list(
  late = "25042-late_outputs/merged_loops/merged_all_results.tsv",
  early = "250831-early_outputs/merged_loops/merged_all_results.tsv"
)

PEAK_FILES <- list(
  late = "peaks/beds/H3K27me3CerebellumLate1.bed",
  early = "peaks/beds/H3K27me3CerebellumEarly1.bed"
)

DIFFERENTIAL_FILES <- list(
  late = "output/shared_anchor_analysis/late/polycomb_specific/tables/both_polycomb_shared_loops.tsv"
)

OUTPUT_DIR <- "output/k27me3_loops_bedpe"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: Load H3K27me3 peaks
load_k27me3_peaks <- function(bed_path) {
  if (!file.exists(bed_path)) {
    stop(sprintf("Peak file not found: %s", bed_path))
  }
  df <- read.table(bed_path, sep = "\t", header = FALSE)
  GRanges(seqnames = df$V1, ranges = IRanges(start = df$V2, end = df$V3))
}

# Helper: Convert loops to BEDPE format
write_bedpe <- function(loops_df, output_path) {
  # Standard BEDPE format: chr1, start1, end1, chr2, start2, end2, name, score, strand1, strand2
  bedpe <- loops_df %>%
    mutate(
      name = if ("loop_id" %in% names(.)) loop_id else paste0("loop_", row_number()),
      score = if ("logFC" %in% names(.)) round(abs(logFC) * 100, 0) else 0,
      strand1 = ".",
      strand2 = "."
    ) %>%
    select(chr1, start1, end1, chr2, start2, end2, name, score, strand1, strand2)

  write.table(bedpe, output_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Wrote %d loops to: %s\n", nrow(bedpe), output_path))
}

# Process each timepoint for ALL K27me3-K27me3 loops
for (timepoint in c("late", "early")) {
  cat(sprintf("\n=== Processing %s timepoint ===\n", toupper(timepoint)))

  input_file <- INPUT_FILES[[timepoint]]
  peak_file <- PEAK_FILES[[timepoint]]

  if (!file.exists(input_file)) {
    cat(sprintf("  Skipping: Input file not found: %s\n", input_file))
    next
  }

  # Load all loops
  cat("  Loading all loops...\n")
  loops <- read_tsv(input_file, show_col_types = FALSE)
  cat(sprintf("  Total loops: %d\n", nrow(loops)))

  # Load K27me3 peaks
  cat("  Loading H3K27me3 peaks...\n")
  k27me3_gr <- load_k27me3_peaks(peak_file)
  cat(sprintf("  H3K27me3 peaks: %d\n", length(k27me3_gr)))

  # Create GRanges for anchors
  anchor1_gr <- GRanges(
    seqnames = loops$chr1,
    ranges = IRanges(start = loops$start1, end = loops$end1)
  )
  anchor2_gr <- GRanges(
    seqnames = loops$chr2,
    ranges = IRanges(start = loops$start2, end = loops$end2)
  )

  # Compute overlaps
  loops$anchor1_K27me3 <- countOverlaps(anchor1_gr, k27me3_gr) > 0
  loops$anchor2_K27me3 <- countOverlaps(anchor2_gr, k27me3_gr) > 0

  # Filter to K27me3-K27me3 (both anchors)
  loops_k27me3_both <- loops %>%
    filter(anchor1_K27me3 & anchor2_K27me3)

  cat(sprintf("  K27me3-K27me3 loops (both anchors): %d\n", nrow(loops_k27me3_both)))

  # Export as BEDPE
  output_file <- file.path(OUTPUT_DIR, sprintf("all_k27me3_k27me3_%s.bedpe", timepoint))
  write_bedpe(loops_k27me3_both, output_file)

  # Also export as TSV with full annotations
  tsv_file <- file.path(OUTPUT_DIR, sprintf("all_k27me3_k27me3_%s.tsv", timepoint))
  write_tsv(loops_k27me3_both, tsv_file)
  cat(sprintf("  Wrote TSV: %s\n", tsv_file))
}

# Export differential K27me3-K27me3 loops (late only, from polycomb analysis)
cat("\n=== Exporting Differential K27me3-K27me3 Loops (Late) ===\n")
diff_file <- DIFFERENTIAL_FILES$late
if (file.exists(diff_file)) {
  diff_loops <- read_tsv(diff_file, show_col_types = FALSE)
  cat(sprintf("  Differential K27me3-K27me3 loops: %d\n", nrow(diff_loops)))

  output_file <- file.path(OUTPUT_DIR, "differential_k27me3_k27me3_late.bedpe")
  write_bedpe(diff_loops, output_file)

  # Also save TSV copy
  tsv_file <- file.path(OUTPUT_DIR, "differential_k27me3_k27me3_late.tsv")
  file.copy(diff_file, tsv_file, overwrite = TRUE)
  cat(sprintf("  Copied TSV: %s\n", tsv_file))
} else {
  cat(sprintf("  WARNING: Differential file not found: %s\n", diff_file))
}

cat("\n=== Export Complete ===\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Files:\n")
system(sprintf("ls -la %s", OUTPUT_DIR))
