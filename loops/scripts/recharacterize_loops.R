#!/usr/bin/env Rscript
# scripts/recharacterize_loops.R
# Re-annotate characterized_loops.tsv with proper timepoint-specific ChIP-seq peaks
# Author: Zakir Alibhai
# Date: 2026-01-28
#
# Purpose:
#   Replace stale ChIP-seq annotations (2 marks, 4 anchor categories) with
#   correct timepoint-specific annotations (5 marks, 7 anchor categories).
#   The old annotations used P12_ctrl_H3K27ac_early_peaks.bed for BOTH timepoints
#   and only had H3K27ac + H3K4me1. The corrected annotations use the proper
#   peaks/beds/ files matched to each timepoint.
#
# What changes:
#   9 columns REPLACED (same names, correct data):
#     anchor1/2_H3K27ac_overlap, anchor1/2_H3K4me1_overlap,
#     anchor1/2_type (7 categories), anchor1/2_is_promoter, loop_type
#   6 columns ADDED (new marks):
#     anchor1/2_H3K27me3_overlap, anchor1/2_H3K4me3_overlap,
#     anchor1/2_Bivalent_Promoter_overlap
#
# Usage:
#   Rscript scripts/recharacterize_loops.R                     # Both timepoints
#   Rscript scripts/recharacterize_loops.R --timepoint late    # Late only
#   Rscript scripts/recharacterize_loops.R --timepoint early   # Early only
#   Rscript scripts/recharacterize_loops.R --dry-run           # Preview without writing
#   Rscript scripts/recharacterize_loops.R --no-backup         # Skip .bak creation

# =============================================================================
# SOURCE HELPER FUNCTIONS
# =============================================================================
# Gets: load_chip_peaks(), annotate_chip_overlaps_extended(),
#        classify_anchor_type_extended(), classify_loop_type_extended(),
#        PEAK_FILES config, ANCHOR_TYPE_ORDER, libraries (GenomicRanges, etc.)

source("scripts/annotate_loops_extended.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Map timepoints to input/output paths and peak set keys
TIMEPOINT_CONFIG <- list(
  early = list(
    characterized = "250831-early_outputs/merged_loops/characterized_loops.tsv",
    non_redundant = "250831-early_outputs/merged_loops/non_redundant_loops.tsv",
    peak_key      = "early"
  ),
  late = list(
    characterized = "25042-late_outputs/merged_loops/characterized_loops.tsv",
    non_redundant = "25042-late_outputs/merged_loops/non_redundant_loops.tsv",
    peak_key      = "late"
  )
)

# Stale columns to drop before recomputing
STALE_COLUMNS <- c(
  "anchor1_H3K27ac_overlap", "anchor1_H3K4me1_overlap",
  "anchor1_type", "anchor1_is_promoter",
  "anchor2_H3K27ac_overlap", "anchor2_H3K4me1_overlap",
  "anchor2_type", "anchor2_is_promoter",
  "loop_type"
)

# Columns for non_redundant_loops.tsv (matches downstream_analysis.R select())
NR_COLUMNS <- c(
  "loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
  "logFC", "logCPM", "FDR", "PValue", "significant",
  "resolution_kb", "n_overlaps", "loop_distance", "loop_type",
  "direction", "category"
)

# =============================================================================
# BACKUP UTILITY
# =============================================================================

#' Create a timestamped backup of a file
#'
#' @param file_path Path to file to back up
#' @param skip_backup If TRUE, skip backup creation
#' @return Backup path (or NULL if skipped)
create_backup <- function(file_path, skip_backup = FALSE) {
  if (skip_backup) {
    cat(sprintf("  Skipping backup for: %s\n", file_path))
    return(NULL)
  }

  bak_path <- paste0(file_path, ".bak")

  # Timestamped backup if .bak already exists
  if (file.exists(bak_path)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    bak_path <- paste0(file_path, ".bak.", timestamp)
  }

  file.copy(file_path, bak_path)
  cat(sprintf("  Backup: %s\n", bak_path))
  bak_path
}

# =============================================================================
# CORE RECHARACTERIZATION
# =============================================================================

#' Recharacterize one timepoint's loop files
#'
#' Drops stale ChIP-seq columns, recomputes from proper peak files,
#' and overwrites characterized_loops.tsv and non_redundant_loops.tsv.
#'
#' @param timepoint "early" or "late"
#' @param config TIMEPOINT_CONFIG entry for this timepoint
#' @param no_backup If TRUE, skip .bak file creation
#' @param dry_run If TRUE, show what would change without writing
#' @return Invisible data.frame of corrected characterized loops
recharacterize_one_timepoint <- function(timepoint, config, no_backup = FALSE,
                                         dry_run = FALSE) {
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("RECHARACTERIZING: %s timepoint\n", toupper(timepoint)))
  cat(sprintf("%s\n\n", paste(rep("=", 60), collapse = "")))

  char_file <- config$characterized
  nr_file   <- config$non_redundant
  peak_key  <- config$peak_key

  # --- Validate inputs ---
  if (!file.exists(char_file)) {
    stop(sprintf("characterized_loops.tsv not found: %s", char_file))
  }
  if (!file.exists(nr_file)) {
    stop(sprintf("non_redundant_loops.tsv not found: %s", nr_file))
  }

  cat(sprintf("  Input:  %s\n", char_file))
  cat(sprintf("  Peaks:  %s set\n\n", peak_key))

  # --- Step 1: Read characterized_loops.tsv ---
  cat("Step 1: Reading characterized_loops.tsv...\n")
  df <- read.table(char_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded %d loops x %d columns\n", nrow(df), ncol(df)))

  # Verify stale columns exist
  missing <- setdiff(STALE_COLUMNS, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing expected stale columns: %s\n  File may have already been recharacterized.",
                 paste(missing, collapse = ", ")))
  }

  # Capture old values for comparison
  old_anchor1_type <- df$anchor1_type
  old_anchor2_type <- df$anchor2_type
  old_loop_type    <- df$loop_type

  # --- Step 2: Drop stale columns ---
  cat("Step 2: Dropping 9 stale ChIP-seq columns...\n")
  stale_idx <- which(colnames(df) %in% STALE_COLUMNS)
  cat(sprintf("  Dropping columns at positions: %s\n", paste(stale_idx, collapse = ", ")))
  df_clean <- df[, -stale_idx, drop = FALSE]
  cat(sprintf("  Remaining: %d columns\n\n", ncol(df_clean)))

  # --- Step 3: Build anchor GRanges ---
  cat("Step 3: Building anchor GRanges...\n")

  anchor1_gr <- GRanges(
    seqnames = df_clean$chr1,
    ranges = IRanges(start = df_clean$start1, end = df_clean$end1)
  )
  anchor2_gr <- GRanges(
    seqnames = df_clean$chr2,
    ranges = IRanges(start = df_clean$start2, end = df_clean$end2)
  )

  cat(sprintf("  Anchor1: %d regions\n", length(anchor1_gr)))
  cat(sprintf("  Anchor2: %d regions\n\n", length(anchor2_gr)))

  # --- Step 4: Load ChIP-seq peaks ---
  cat(sprintf("Step 4: Loading ChIP-seq peaks (%s)...\n", peak_key))

  peak_files <- PEAK_FILES[[peak_key]]
  k27ac_gr   <- load_chip_peaks(peak_files$h3k27ac,  "H3K27ac")
  k27me3_gr  <- load_chip_peaks(peak_files$h3k27me3, "H3K27me3")
  k4me1_gr   <- load_chip_peaks(peak_files$h3k4me1,  "H3K4me1")
  k4me3_gr   <- load_chip_peaks(peak_files$h3k4me3,  "H3K4me3")
  bivalent_gr <- load_chip_peaks(peak_files$bivalent, "Bivalent (K4me3+K27me3)")
  cat("\n")

  # --- Step 5: Compute ChIP-seq overlaps ---
  cat("Step 5: Computing ChIP-seq overlaps...\n")

  anchor1_chip <- annotate_chip_overlaps_extended(anchor1_gr, k27ac_gr, k27me3_gr,
                                                   k4me1_gr, k4me3_gr, bivalent_gr)
  anchor2_chip <- annotate_chip_overlaps_extended(anchor2_gr, k27ac_gr, k27me3_gr,
                                                   k4me1_gr, k4me3_gr, bivalent_gr)

  for (mark in c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "Bivalent_Promoter")) {
    col <- paste0(mark, "_overlap")
    cat(sprintf("  Anchor1 %s+: %d (%.1f%%)\n", mark,
                sum(anchor1_chip[[col]]), 100 * mean(anchor1_chip[[col]])))
    cat(sprintf("  Anchor2 %s+: %d (%.1f%%)\n", mark,
                sum(anchor2_chip[[col]]), 100 * mean(anchor2_chip[[col]])))
  }
  cat("\n")

  # --- Step 6: Classify anchor types (reuse existing distance_to_tss) ---
  cat("Step 6: Classifying anchor types (7 categories)...\n")
  cat("  Reusing existing anchor1/2_distance_to_tss (gene annotation, not ChIP-dependent)\n")

  anchor1_type <- classify_anchor_type_extended(
    anchor1_chip$H3K27ac_overlap,
    anchor1_chip$H3K27me3_overlap,
    anchor1_chip$H3K4me1_overlap,
    anchor1_chip$H3K4me3_overlap,
    anchor1_chip$Bivalent_Promoter_overlap,
    df_clean$anchor1_distance_to_tss,
    tss_threshold = 2000
  )

  anchor2_type <- classify_anchor_type_extended(
    anchor2_chip$H3K27ac_overlap,
    anchor2_chip$H3K27me3_overlap,
    anchor2_chip$H3K4me1_overlap,
    anchor2_chip$H3K4me3_overlap,
    anchor2_chip$Bivalent_Promoter_overlap,
    df_clean$anchor2_distance_to_tss,
    tss_threshold = 2000
  )

  cat("\n  Anchor1 type distribution:\n")
  for (type in ANCHOR_TYPE_ORDER) {
    count <- sum(anchor1_type == type)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", type, count, 100 * count / length(anchor1_type)))
  }
  cat("\n  Anchor2 type distribution:\n")
  for (type in ANCHOR_TYPE_ORDER) {
    count <- sum(anchor2_type == type)
    cat(sprintf("    %-20s: %5d (%.1f%%)\n", type, count, 100 * count / length(anchor2_type)))
  }
  cat("\n")

  # --- Step 7: Classify loop types ---
  cat("Step 7: Classifying loop types...\n")

  loop_type <- classify_loop_type_extended(anchor1_type, anchor2_type)

  loop_type_counts <- sort(table(loop_type), decreasing = TRUE)
  cat(sprintf("  %d unique loop types (was ~10, now up to 28)\n\n", length(loop_type_counts)))
  for (i in seq_along(loop_type_counts)) {
    cat(sprintf("    %-35s: %5d (%.1f%%)\n",
                names(loop_type_counts)[i],
                loop_type_counts[i],
                100 * loop_type_counts[i] / length(loop_type)))
  }
  cat("\n")

  # --- Step 8: Attach corrected columns ---
  cat("Step 8: Attaching corrected columns...\n")

  # Replaced columns (same names as originals, at original positions)
  df_clean$anchor1_H3K27ac_overlap  <- anchor1_chip$H3K27ac_overlap
  df_clean$anchor1_H3K4me1_overlap  <- anchor1_chip$H3K4me1_overlap
  df_clean$anchor1_type             <- anchor1_type
  df_clean$anchor1_is_promoter      <- anchor1_type == "Active_Promoter"
  df_clean$anchor2_H3K27ac_overlap  <- anchor2_chip$H3K27ac_overlap
  df_clean$anchor2_H3K4me1_overlap  <- anchor2_chip$H3K4me1_overlap
  df_clean$anchor2_type             <- anchor2_type
  df_clean$anchor2_is_promoter      <- anchor2_type == "Active_Promoter"
  df_clean$loop_type                <- loop_type

  # New columns (appended at end)
  df_clean$anchor1_H3K27me3_overlap          <- anchor1_chip$H3K27me3_overlap
  df_clean$anchor2_H3K27me3_overlap          <- anchor2_chip$H3K27me3_overlap
  df_clean$anchor1_H3K4me3_overlap           <- anchor1_chip$H3K4me3_overlap
  df_clean$anchor2_H3K4me3_overlap           <- anchor2_chip$H3K4me3_overlap
  df_clean$anchor1_Bivalent_Promoter_overlap <- anchor1_chip$Bivalent_Promoter_overlap
  df_clean$anchor2_Bivalent_Promoter_overlap <- anchor2_chip$Bivalent_Promoter_overlap

  cat(sprintf("  Final: %d loops x %d columns (was %d)\n", nrow(df_clean), ncol(df_clean), ncol(df)))

  # --- Step 9: Before/after comparison ---
  cat("\nStep 9: Before/after comparison...\n")

  # Anchor type changes
  old_types <- sort(unique(old_anchor1_type))
  new_types <- sort(unique(anchor1_type))
  cat(sprintf("  Old anchor categories (%d): %s\n", length(old_types), paste(old_types, collapse = ", ")))
  cat(sprintf("  New anchor categories (%d): %s\n", length(new_types), paste(new_types, collapse = ", ")))

  # Loop type changes
  old_lt <- sort(unique(old_loop_type))
  new_lt <- sort(unique(loop_type))
  cat(sprintf("  Old loop types (%d): %s\n", length(old_lt), paste(old_lt, collapse = ", ")))
  cat(sprintf("  New loop types (%d): %s\n", length(new_lt), paste(new_lt, collapse = ", ")))

  # Count how many loops changed type
  changed <- sum(old_loop_type != loop_type)
  cat(sprintf("  Loops that changed type: %d / %d (%.1f%%)\n\n",
              changed, length(loop_type), 100 * changed / length(loop_type)))

  if (dry_run) {
    cat("  [DRY RUN] Would write corrected files but --dry-run specified. Skipping.\n\n")
    return(invisible(df_clean))
  }

  # --- Step 10: Write corrected characterized_loops.tsv ---
  cat("Step 10: Writing corrected files...\n")

  create_backup(char_file, skip_backup = no_backup)
  write.table(df_clean, char_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Wrote: %s (%d rows x %d cols)\n", char_file, nrow(df_clean), ncol(df_clean)))

  # --- Step 11: Rebuild non_redundant_loops.tsv ---
  cat("\nStep 11: Rebuilding non_redundant_loops.tsv...\n")

  # Read existing NR to preserve row count and verify consistency
  nr_df <- read.table(nr_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Old NR: %d rows x %d cols\n", nrow(nr_df), ncol(nr_df)))

  # Sanity check: row counts must match
  if (nrow(nr_df) != nrow(df_clean)) {
    stop(sprintf("Row count mismatch: characterized=%d, non_redundant=%d",
                 nrow(df_clean), nrow(nr_df)))
  }

  # Rebuild NR from corrected characterized data
  nr_cols_present <- intersect(NR_COLUMNS, colnames(df_clean))
  nr_missing <- setdiff(NR_COLUMNS, colnames(df_clean))
  if (length(nr_missing) > 0) {
    warning(sprintf("NR columns missing from characterized: %s", paste(nr_missing, collapse = ", ")))
  }

  nr_new <- df_clean[, nr_cols_present, drop = FALSE]

  create_backup(nr_file, skip_backup = no_backup)
  write.table(nr_new, nr_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Wrote: %s (%d rows x %d cols)\n\n", nr_file, nrow(nr_new), ncol(nr_new)))

  cat(sprintf("Recharacterization complete for %s timepoint.\n", toupper(timepoint)))
  invisible(df_clean)
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  timepoint  <- "both"
  no_backup  <- FALSE
  dry_run    <- FALSE

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]

    if (arg == "--help" || arg == "-h") {
      cat("Usage: Rscript scripts/recharacterize_loops.R [OPTIONS]\n\n")
      cat("Re-annotate characterized_loops.tsv with proper timepoint-specific ChIP-seq peaks.\n\n")
      cat("Options:\n")
      cat("  --timepoint TIMEPOINT  Which timepoint to process: early, late, or both (default: both)\n")
      cat("  --no-backup            Skip creating .bak backup files\n")
      cat("  --dry-run              Preview changes without writing files\n")
      cat("  --help, -h             Show this help message\n\n")
      cat("Stale columns replaced (same names, corrected data):\n")
      cat("  anchor1/2_H3K27ac_overlap, anchor1/2_H3K4me1_overlap,\n")
      cat("  anchor1/2_type (7 categories), anchor1/2_is_promoter, loop_type\n\n")
      cat("New columns added:\n")
      cat("  anchor1/2_H3K27me3_overlap, anchor1/2_H3K4me3_overlap,\n")
      cat("  anchor1/2_Bivalent_Promoter_overlap\n\n")
      quit(status = 0)

    } else if (arg == "--timepoint") {
      i <- i + 1
      if (i > length(args)) stop("--timepoint requires a value: early, late, or both")
      timepoint <- args[i]
      if (!timepoint %in% c("early", "late", "both")) {
        stop("--timepoint must be 'early', 'late', or 'both'")
      }

    } else if (arg == "--no-backup") {
      no_backup <- TRUE

    } else if (arg == "--dry-run") {
      dry_run <- TRUE

    } else {
      stop(sprintf("Unknown argument: %s (use --help for usage)", arg))
    }

    i <- i + 1
  }

  list(timepoint = timepoint, no_backup = no_backup, dry_run = dry_run)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if (!interactive() && sys.nframe() == 0) {
  args <- parse_arguments()

  cat("\n")
  cat("================================================================\n")
  cat("  Recharacterize Loop Annotations with Proper ChIP-seq Peaks\n")
  cat("================================================================\n")
  cat(sprintf("  Timepoint: %s\n", args$timepoint))
  cat(sprintf("  Dry run:   %s\n", args$dry_run))
  cat(sprintf("  Backup:    %s\n", !args$no_backup))
  cat("================================================================\n")

  # Determine which timepoints to process
  if (args$timepoint == "both") {
    timepoints <- c("early", "late")
  } else {
    timepoints <- args$timepoint
  }

  # Process each timepoint
  for (tp in timepoints) {
    config <- TIMEPOINT_CONFIG[[tp]]
    recharacterize_one_timepoint(tp, config, args$no_backup, args$dry_run)
  }

  # Print downstream re-run instructions
  cat("\n")
  cat("================================================================\n")
  cat("  DOWNSTREAM SCRIPTS TO RE-RUN\n")
  cat("================================================================\n\n")
  cat("The following scripts read from the files that were just updated:\n\n")
  cat("  1. loop_distance_analysis.R (Figures 6-7 use anchor_type / loop_type)\n")
  cat("     Rscript scripts/loop_distance_analysis.R --timepoint both\n\n")
  cat("  2. visualizations.R (loop type pie charts) -- only if it reads\n")
  cat("     from 25042-late_outputs/ or 250831-early_outputs/\n\n")
  cat("  Scripts that are SAFE (recompute their own annotations):\n")
  cat("    - loop_distance_k27me3_filtered.R\n")
  cat("    - chip_distance_analysis.R\n\n")
  cat("================================================================\n")
  cat("  VERIFICATION STEPS\n")
  cat("================================================================\n\n")
  cat("  1. Diff .bak vs new file to confirm only ChIP columns changed:\n")
  cat("     diff <(head -1 FILE.bak) <(head -1 FILE)\n\n")
  cat("  2. Verify row counts preserved:\n")
  cat("     wc -l 25042-late_outputs/merged_loops/characterized_loops.tsv\n")
  cat("     wc -l 250831-early_outputs/merged_loops/characterized_loops.tsv\n\n")
  cat("  3. Spot-check anchor types match annotate_loops_extended.R output\n\n")
  cat("================================================================\n")
}
