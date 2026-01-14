# scripts/05_filter_blacklist.R
# Post-filter TADCompare results to remove boundaries in blacklisted regions
# TADCompare does not have native blacklist support - this is post-hoc filtering

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)

# =============================================================================
# Configuration
# =============================================================================

RESOLUTION <- 10000
TIMEPOINTS <- c("early", "late")
BASE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"

# Two blacklist files: lab-specific and ENCODE standard (union both)
BLACKLIST_LAB <- file.path(BASE_DIR, "250123blacklist.bed")
BLACKLIST_ENCODE <- file.path(BASE_DIR, "mm10-blacklist.v2.bed")

cat("=== TADCompare Blacklist Filtering ===\n")
cat("Base directory:", BASE_DIR, "\n")
cat("Lab blacklist:", BLACKLIST_LAB, "\n")
cat("ENCODE blacklist:", BLACKLIST_ENCODE, "\n")
cat("Resolution:", RESOLUTION, "bp\n")
cat("Timepoints:", paste(TIMEPOINTS, collapse = ", "), "\n\n")

# =============================================================================
# Load and union blacklists
# =============================================================================

# Check files exist
if (!file.exists(BLACKLIST_LAB)) {
  stop("Lab blacklist file not found: ", BLACKLIST_LAB)
}
if (!file.exists(BLACKLIST_ENCODE)) {
  stop("ENCODE blacklist file not found: ", BLACKLIST_ENCODE)
}

cat("Loading blacklist regions...\n")

# Load both blacklists
blacklist_lab <- import.bed(BLACKLIST_LAB)
blacklist_encode <- import.bed(BLACKLIST_ENCODE)
cat("  Lab blacklist:", length(blacklist_lab), "regions\n")
cat("  ENCODE blacklist:", length(blacklist_encode), "regions\n")

# Create union (merge overlapping regions from both)
blacklist <- reduce(c(blacklist_lab, blacklist_encode))
cat("  Combined union:", length(blacklist), "regions (after merging overlaps)\n")
cat("  Total blacklisted bases:", sum(width(blacklist)), "bp\n\n")

# =============================================================================
# Filter function
# =============================================================================

filter_boundaries_by_blacklist <- function(results, blacklist, resolution) {
  # Convert boundaries to GRanges (boundary bin spans resolution)
  boundaries_gr <- GRanges(
    seqnames = results$chr,
    ranges = IRanges(
      start = results$Boundary,
      end = results$Boundary + resolution - 1
    )
  )

  # Find boundaries NOT overlapping blacklist (any overlap = remove)
  overlaps <- countOverlaps(boundaries_gr, blacklist)
  keep <- overlaps == 0

  return(keep)
}

# =============================================================================
# Process each timepoint
# =============================================================================

all_summaries <- list()

for (timepoint in TIMEPOINTS) {
  cat("Processing", timepoint, "timepoint...\n")

  # Define paths
  input_file <- file.path(BASE_DIR, "results", timepoint, "final",
                          "tadcompare_final_annotated.tsv")
  output_file <- file.path(BASE_DIR, "results", timepoint, "final",
                           "tadcompare_final_filtered.tsv")
  removed_file <- file.path(BASE_DIR, "results", timepoint, "final",
                            "blacklist_removed_boundaries.tsv")
  summary_file <- file.path(BASE_DIR, "results", timepoint, "final",
                            "blacklist_filter_summary.txt")

  # Check input exists
  if (!file.exists(input_file)) {
    cat("  WARNING: Input file not found:", input_file, "\n")
    cat("  Skipping", timepoint, "timepoint\n\n")
    next
  }

  # Load results
  results <- read_tsv(input_file, show_col_types = FALSE)
  n_before <- nrow(results)
  cat("  Loaded", n_before, "boundaries\n")

  # Apply blacklist filter
  keep <- filter_boundaries_by_blacklist(results, blacklist, RESOLUTION)

  # Split into filtered and removed
  results_filtered <- results[keep, ]
  results_removed <- results[!keep, ]
  n_after <- nrow(results_filtered)
  n_removed <- n_before - n_after

  cat("  Removed", n_removed, "boundaries (",
      round(100 * n_removed / n_before, 1), "%)\n")
  cat("  Remaining", n_after, "boundaries\n")

  # Save filtered results
  write_tsv(results_filtered, output_file)
  cat("  Saved filtered results to:", basename(output_file), "\n")

  # Save removed boundaries for inspection
  if (nrow(results_removed) > 0) {
    write_tsv(results_removed, removed_file)
    cat("  Saved removed boundaries to:", basename(removed_file), "\n")
  }

  # Calculate detailed statistics
  n_diff_before <- sum(results$Differential == "Differential")
  n_diff_removed <- sum(!keep & results$Differential == "Differential")
  n_diff_after <- sum(results_filtered$Differential == "Differential")

  n_robust_before <- sum(results$robustness != "Neither")
  n_robust_removed <- sum(!keep & results$robustness != "Neither")
  n_robust_after <- sum(results_filtered$robustness != "Neither")

  # Chromosome breakdown of removed
  if (n_removed > 0) {
    chr_counts <- sort(table(results_removed$chr), decreasing = TRUE)
    top_chrs <- paste(
      paste0(names(chr_counts)[1:min(5, length(chr_counts))], " (",
             chr_counts[1:min(5, length(chr_counts))], ")"),
      collapse = ", "
    )
  } else {
    top_chrs <- "None"
  }

  # Gap Score distribution of removed
  if (n_removed > 0) {
    removed_gap_scores <- results_removed$Gap_Score
    gap_score_stats <- sprintf(
      "min=%.2f, median=%.2f, max=%.2f",
      min(removed_gap_scores), median(removed_gap_scores), max(removed_gap_scores)
    )
  } else {
    gap_score_stats <- "N/A"
  }

  # Write summary
  summary_text <- sprintf(
"Blacklist Filtering Summary - %s timepoint
==========================================
Date: %s
Lab blacklist: %s (%d regions)
ENCODE blacklist: %s (%d regions)
Combined union: %d regions (after merging overlaps)
Resolution: %d bp
Overlap threshold: Any overlap (conservative)

BOUNDARIES
----------
Before filtering: %d
Removed: %d (%.1f%%)
After filtering: %d

DIFFERENTIAL BOUNDARIES
-----------------------
Before filtering: %d
Removed: %d (%.1f%% of differential)
After filtering: %d

HIGH-CONFIDENCE (robustness != Neither)
---------------------------------------
Before filtering: %d
Removed: %d (%.1f%% of high-confidence)
After filtering: %d

CHROMOSOMES WITH MOST REMOVED
-----------------------------
%s

GAP SCORES OF REMOVED BOUNDARIES
--------------------------------
%s

OUTPUT FILES
------------
Filtered results: %s
Removed boundaries: %s
",
    timepoint,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    BLACKLIST_LAB, length(blacklist_lab),
    BLACKLIST_ENCODE, length(blacklist_encode),
    length(blacklist),
    RESOLUTION,
    n_before, n_removed, 100 * n_removed / n_before, n_after,
    n_diff_before, n_diff_removed,
    ifelse(n_diff_before > 0, 100 * n_diff_removed / n_diff_before, 0),
    n_diff_after,
    n_robust_before, n_robust_removed,
    ifelse(n_robust_before > 0, 100 * n_robust_removed / n_robust_before, 0),
    n_robust_after,
    top_chrs,
    gap_score_stats,
    output_file,
    ifelse(n_removed > 0, removed_file, "N/A (none removed)")
  )

  writeLines(summary_text, summary_file)
  cat("  Saved summary to:", basename(summary_file), "\n\n")

  # Store for final summary
  all_summaries[[timepoint]] <- list(
    n_before = n_before,
    n_removed = n_removed,
    n_after = n_after,
    n_diff_removed = n_diff_removed,
    pct_removed = 100 * n_removed / n_before
  )
}

# =============================================================================
# Final summary
# =============================================================================

cat("=== Final Summary ===\n")
for (tp in names(all_summaries)) {
  s <- all_summaries[[tp]]
  cat(sprintf("%s: %d -> %d boundaries (removed %d, %.1f%%)\n",
              tp, s$n_before, s$n_after, s$n_removed, s$pct_removed))
}
cat("\nDone!\n")
