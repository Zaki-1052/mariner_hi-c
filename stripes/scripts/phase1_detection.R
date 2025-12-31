#!/usr/bin/env Rscript
# stripes/scripts/phase1_detection.R
# Phase 1: Detection & Union Set Creation for Differential Stripe Analysis
#
# Creates unified stripe set from merged BEDPE files with:
#   - Direction-aware anchor matching (50kb tolerance)
#   - Source classification (control_only, mutant_only, shared)
#   - Replicate support counting for confidence scoring
#   - 10kb resolution validation annotation
#
# Usage: Rscript phase1_detection.R <timepoint>
#   timepoint: "early" (250831) or "late" (250402)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(yaml)
  library(dplyr)
})

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript phase1_detection.R <timepoint>\n",
       "  timepoint: 'early' or 'late'")
}
TIMEPOINT <- args[1]

if (!TIMEPOINT %in% c("early", "late")) {
  stop("Timepoint must be 'early' or 'late'")
}

cat("\n========================================\n")
cat("Phase 1: Detection & Union Set Creation\n")
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
cat(sprintf("  Project base: %s\n", config$project$base_dir))

# ==============================================================================
# HELPER: Build Paths from Config
# ==============================================================================
get_stripe_bedpe_path <- function(resolution_kb, timepoint, sample, config) {
  # Build path: {base_path}/{resolution}kb/{timepoint}/{sample}_stripes.bedpe
  timepoint_dir <- config$stripe_data$timepoints[[timepoint]]
  file.path(
    config$stripe_data$base_path,
    paste0(resolution_kb, "kb"),
    timepoint_dir,
    paste0(sample, "_stripes.bedpe")
  )
}

get_output_dir <- function(timepoint, resolution_kb, config) {
  # Build path: {base_dir}/{timepoint}/
  file.path(
    config$outputs$base_dir,
    timepoint
  )
}

# ==============================================================================
# FUNCTION: Parse Stripe BEDPE Format
# ==============================================================================
parse_stripe_bedpe <- function(filepath, source_name) {
  if (!file.exists(filepath)) {
    warning(sprintf("File not found: %s", filepath))
    return(NULL)
  }

  bedpe <- tryCatch({
    read.table(filepath, header = FALSE, skip = 1, stringsAsFactors = FALSE,
               sep = "\t",
               col.names = c("chr1", "x1", "x2", "chr2", "y1", "y2", "pval"))
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", filepath, e$message))
    return(NULL)
  })

  if (is.null(bedpe) || nrow(bedpe) == 0) return(NULL)

  # Filter to intrachromosomal only
  bedpe <- bedpe[bedpe$chr1 == bedpe$chr2, ]
  if (nrow(bedpe) == 0) return(NULL)

  # Infer stripe direction
  # 3' stripe: span downstream (y1 >= x2)
  # 5' stripe: span upstream (y2 <= x1)
  bedpe$direction <- case_when(
    bedpe$y1 >= bedpe$x2 ~ "3prime",
    bedpe$y2 <= bedpe$x1 ~ "5prime",
    TRUE ~ "ambiguous"
  )

  # Handle ambiguous by midpoint comparison
  ambig_idx <- bedpe$direction == "ambiguous"
  if (any(ambig_idx)) {
    anchor_mid <- (bedpe$x1[ambig_idx] + bedpe$x2[ambig_idx]) / 2
    span_mid <- (bedpe$y1[ambig_idx] + bedpe$y2[ambig_idx]) / 2
    bedpe$direction[ambig_idx] <- ifelse(span_mid > anchor_mid, "3prime", "5prime")
  }

  bedpe$anchor_center <- (bedpe$x1 + bedpe$x2) / 2
  bedpe$source_file <- source_name

  bedpe
}

# ==============================================================================
# FUNCTION: Convert to GRanges for Anchor Matching
# ==============================================================================
stripes_to_granges <- function(stripe_df) {
  if (is.null(stripe_df) || nrow(stripe_df) == 0) return(GRanges())

  GRanges(
    seqnames = stripe_df$chr1,
    ranges = IRanges(start = stripe_df$x1, end = stripe_df$x2),
    direction = stripe_df$direction,
    pval = stripe_df$pval,
    span_y1 = stripe_df$y1,
    span_y2 = stripe_df$y2,
    anchor_center = stripe_df$anchor_center
  )
}

# ==============================================================================
# FUNCTION: Match Anchors with Direction Awareness
# ==============================================================================
match_anchors <- function(gr_ctrl, gr_mut, tolerance) {
  if (length(gr_ctrl) == 0 || length(gr_mut) == 0) {
    return(data.frame(ctrl_idx = integer(0), mut_idx = integer(0),
                      same_direction = logical(0)))
  }

  # Extend for tolerance matching
  gr_ctrl_ext <- resize(gr_ctrl, width(gr_ctrl) + 2 * tolerance, fix = "center")
  gr_mut_ext <- resize(gr_mut, width(gr_mut) + 2 * tolerance, fix = "center")

  overlaps <- findOverlaps(gr_ctrl_ext, gr_mut_ext)
  if (length(overlaps) == 0) {
    return(data.frame(ctrl_idx = integer(0), mut_idx = integer(0),
                      same_direction = logical(0)))
  }

  data.frame(
    ctrl_idx = queryHits(overlaps),
    mut_idx = subjectHits(overlaps),
    same_direction = gr_ctrl$direction[queryHits(overlaps)] ==
                     gr_mut$direction[subjectHits(overlaps)]
  )
}

# ==============================================================================
# FUNCTION: Build Union Set
# ==============================================================================
build_union_set <- function(ctrl_merged, mut_merged, tolerance) {
  gr_ctrl <- stripes_to_granges(ctrl_merged)
  gr_mut <- stripes_to_granges(mut_merged)

  cat(sprintf("  Control merged: %d stripes\n", length(gr_ctrl)))
  cat(sprintf("  Mutant merged: %d stripes\n", length(gr_mut)))

  matches <- match_anchors(gr_ctrl, gr_mut, tolerance)

  # Shared = same anchor AND same direction
  shared <- matches[matches$same_direction, ]
  shared <- shared[!duplicated(shared$ctrl_idx), ]
  shared <- shared[!duplicated(shared$mut_idx), ]

  ctrl_shared <- unique(shared$ctrl_idx)
  mut_shared <- unique(shared$mut_idx)
  ctrl_only <- setdiff(seq_len(length(gr_ctrl)), ctrl_shared)
  mut_only <- setdiff(seq_len(length(gr_mut)), mut_shared)

  cat(sprintf("  Shared: %d, Control-only: %d, Mutant-only: %d\n",
              nrow(shared), length(ctrl_only), length(mut_only)))

  # Build union dataframe
  union_list <- list()

  # Shared stripes
  if (nrow(shared) > 0) {
    for (i in seq_len(nrow(shared))) {
      ci <- shared$ctrl_idx[i]
      mi <- shared$mut_idx[i]
      union_list[[length(union_list) + 1]] <- data.frame(
        chr = as.character(seqnames(gr_ctrl)[ci]),
        anchor_x1 = start(gr_ctrl)[ci], anchor_x2 = end(gr_ctrl)[ci],
        span_y1 = gr_ctrl$span_y1[ci], span_y2 = gr_ctrl$span_y2[ci],
        direction_type = gr_ctrl$direction[ci],
        source = "shared",
        pval_ctrl = gr_ctrl$pval[ci], pval_mut = gr_mut$pval[mi],
        anchor_center = gr_ctrl$anchor_center[ci],
        stringsAsFactors = FALSE
      )
    }
  }

  # Control-only stripes
  for (i in ctrl_only) {
    union_list[[length(union_list) + 1]] <- data.frame(
      chr = as.character(seqnames(gr_ctrl)[i]),
      anchor_x1 = start(gr_ctrl)[i], anchor_x2 = end(gr_ctrl)[i],
      span_y1 = gr_ctrl$span_y1[i], span_y2 = gr_ctrl$span_y2[i],
      direction_type = gr_ctrl$direction[i],
      source = "control_only",
      pval_ctrl = gr_ctrl$pval[i], pval_mut = NA_real_,
      anchor_center = gr_ctrl$anchor_center[i],
      stringsAsFactors = FALSE
    )
  }

  # Mutant-only stripes
  for (i in mut_only) {
    union_list[[length(union_list) + 1]] <- data.frame(
      chr = as.character(seqnames(gr_mut)[i]),
      anchor_x1 = start(gr_mut)[i], anchor_x2 = end(gr_mut)[i],
      span_y1 = gr_mut$span_y1[i], span_y2 = gr_mut$span_y2[i],
      direction_type = gr_mut$direction[i],
      source = "mutant_only",
      pval_ctrl = NA_real_, pval_mut = gr_mut$pval[i],
      anchor_center = gr_mut$anchor_center[i],
      stringsAsFactors = FALSE
    )
  }

  if (length(union_list) == 0) stop("No stripes in union set")

  union_df <- do.call(rbind, union_list)
  union_df <- union_df[order(union_df$chr, union_df$anchor_x1), ]
  union_df$stripe_id <- paste0("stripe_", sprintf("%04d", seq_len(nrow(union_df))))

  union_df
}

# ==============================================================================
# FUNCTION: Count Replicate Support
# ==============================================================================
count_replicate_support <- function(union_df, rep_files_ctrl, rep_files_mut,
                                    tolerance) {
  union_df$n_ctrl_reps <- 0
  union_df$n_mut_reps <- 0

  union_gr <- GRanges(
    seqnames = union_df$chr,
    ranges = IRanges(union_df$anchor_x1, union_df$anchor_x2),
    direction = union_df$direction_type
  )
  union_gr_ext <- resize(union_gr, width(union_gr) + 2 * tolerance, fix = "center")

  # Count ctrl replicate support
  for (rep_file in rep_files_ctrl) {
    rep_df <- parse_stripe_bedpe(rep_file, basename(rep_file))
    if (is.null(rep_df)) next
    rep_gr <- stripes_to_granges(rep_df)
    if (length(rep_gr) == 0) next

    ov <- findOverlaps(union_gr_ext, rep_gr)
    for (j in seq_len(length(ov))) {
      qi <- queryHits(ov)[j]
      si <- subjectHits(ov)[j]
      if (union_df$direction_type[qi] == rep_gr$direction[si]) {
        union_df$n_ctrl_reps[qi] <- union_df$n_ctrl_reps[qi] + 1
      }
    }
  }

  # Count mut replicate support
  for (rep_file in rep_files_mut) {
    rep_df <- parse_stripe_bedpe(rep_file, basename(rep_file))
    if (is.null(rep_df)) next
    rep_gr <- stripes_to_granges(rep_df)
    if (length(rep_gr) == 0) next

    ov <- findOverlaps(union_gr_ext, rep_gr)
    for (j in seq_len(length(ov))) {
      qi <- queryHits(ov)[j]
      si <- subjectHits(ov)[j]
      if (union_df$direction_type[qi] == rep_gr$direction[si]) {
        union_df$n_mut_reps[qi] <- union_df$n_mut_reps[qi] + 1
      }
    }
  }

  # Cap at 3
  union_df$n_ctrl_reps <- pmin(union_df$n_ctrl_reps, 3)
  union_df$n_mut_reps <- pmin(union_df$n_mut_reps, 3)

  union_df
}

# ==============================================================================
# FUNCTION: Annotate 10kb Validation
# ==============================================================================
annotate_10kb_validation <- function(union_df, config, timepoint, tolerance) {
  union_df$in_10kb <- FALSE

  # Get 10kb merged files
  ctrl_10kb <- get_stripe_bedpe_path(10, timepoint, "ctrl_merged", config)
  mut_10kb <- get_stripe_bedpe_path(10, timepoint, "mut_merged", config)

  all_10kb <- list()
  if (file.exists(ctrl_10kb)) {
    df <- parse_stripe_bedpe(ctrl_10kb, "ctrl_10kb")
    if (!is.null(df)) all_10kb[[length(all_10kb) + 1]] <- df
  }
  if (file.exists(mut_10kb)) {
    df <- parse_stripe_bedpe(mut_10kb, "mut_10kb")
    if (!is.null(df)) all_10kb[[length(all_10kb) + 1]] <- df
  }

  if (length(all_10kb) == 0) {
    cat("  Warning: No 10kb files found\n")
    return(union_df)
  }

  all_10kb_df <- do.call(rbind, all_10kb)
  gr_10kb <- stripes_to_granges(all_10kb_df)
  cat(sprintf("  10kb stripes for validation: %d\n", length(gr_10kb)))

  union_gr <- GRanges(
    seqnames = union_df$chr,
    ranges = IRanges(union_df$anchor_x1, union_df$anchor_x2),
    direction = union_df$direction_type
  )
  union_gr_ext <- resize(union_gr, width(union_gr) + 2 * tolerance, fix = "center")

  ov <- findOverlaps(union_gr_ext, gr_10kb)
  for (j in seq_len(length(ov))) {
    qi <- queryHits(ov)[j]
    si <- subjectHits(ov)[j]
    if (union_df$direction_type[qi] == gr_10kb$direction[si]) {
      union_df$in_10kb[qi] <- TRUE
    }
  }

  cat(sprintf("  Validated in 10kb: %d (%.1f%%)\n",
              sum(union_df$in_10kb), 100 * mean(union_df$in_10kb)))

  union_df
}

# ==============================================================================
# FUNCTION: Assign Confidence Scores
# ==============================================================================
assign_confidence <- function(union_df, config) {
  conf <- config$detection$confidence

  union_df$n_reps_relevant <- case_when(
    union_df$source == "control_only" ~ union_df$n_ctrl_reps,
    union_df$source == "mutant_only" ~ union_df$n_mut_reps,
    union_df$source == "shared" ~ pmax(union_df$n_ctrl_reps, union_df$n_mut_reps)
  )

  union_df$pval_min <- pmin(union_df$pval_ctrl, union_df$pval_mut, na.rm = TRUE)

  union_df$confidence <- case_when(
    union_df$n_reps_relevant >= conf$high_min_reps & union_df$in_10kb ~ "high",
    union_df$n_reps_relevant >= conf$medium_min_reps |
      union_df$in_10kb |
      union_df$pval_min < conf$high_pval_threshold ~ "medium",
    TRUE ~ "low"
  )

  union_df
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
timepoint_dir <- config$stripe_data$timepoints[[TIMEPOINT]]
tolerance <- config$detection$anchor_tolerance_bp
resolution_kb <- config$stripe_data$resolutions$primary / 1000

cat(sprintf("\nTimepoint directory: %s\n", timepoint_dir))
cat(sprintf("Anchor tolerance: %d kb\n", tolerance / 1000))
cat(sprintf("Primary resolution: %d kb\n", resolution_kb))

# Build file paths from config
cat("\n=== Loading BEDPE Files ===\n")
ctrl_merged_file <- get_stripe_bedpe_path(resolution_kb, TIMEPOINT,
                                          "ctrl_merged", config)
mut_merged_file <- get_stripe_bedpe_path(resolution_kb, TIMEPOINT,
                                         "mut_merged", config)

cat(sprintf("  ctrl_merged: %s\n", ctrl_merged_file))
cat(sprintf("  mut_merged: %s\n", mut_merged_file))

ctrl_merged <- parse_stripe_bedpe(ctrl_merged_file, "ctrl_merged")
mut_merged <- parse_stripe_bedpe(mut_merged_file, "mut_merged")

if (is.null(ctrl_merged) || is.null(mut_merged)) {
  stop("Failed to load merged BEDPE files")
}

cat(sprintf("\n  Loaded ctrl: %d stripes\n", nrow(ctrl_merged)))
cat(sprintf("  Loaded mut: %d stripes\n", nrow(mut_merged)))

# Individual replicate files
rep_files_ctrl <- sapply(c("M1", "M2", "M3"), function(rep) {
  get_stripe_bedpe_path(resolution_kb, TIMEPOINT, paste0("ctrl_", rep), config)
})
rep_files_mut <- sapply(c("M1", "M2", "M3"), function(rep) {
  get_stripe_bedpe_path(resolution_kb, TIMEPOINT, paste0("mut_", rep), config)
})

# Build union set
cat("\n=== Building Union Set ===\n")
union_df <- build_union_set(ctrl_merged, mut_merged, tolerance)
cat(sprintf("Total union: %d stripes\n", nrow(union_df)))

# Count replicate support
cat("\n=== Counting Replicate Support ===\n")
union_df <- count_replicate_support(union_df, rep_files_ctrl, rep_files_mut,
                                    tolerance)

# 10kb validation
cat("\n=== 10kb Validation ===\n")
union_df <- annotate_10kb_validation(union_df, config, TIMEPOINT, tolerance)

# Confidence scoring
cat("\n=== Confidence Scoring ===\n")
union_df <- assign_confidence(union_df, config)

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

output_dir <- get_output_dir(TIMEPOINT, resolution_kb, config)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("  Created: %s\n", output_dir))
}

# Reorder columns
output_cols <- c("stripe_id", "chr", "anchor_x1", "anchor_x2",
                 "span_y1", "span_y2", "direction_type", "source",
                 "pval_ctrl", "pval_mut", "n_ctrl_reps", "n_mut_reps",
                 "in_10kb", "confidence", "anchor_center")
union_df <- union_df[, output_cols]

# Save RDS and TSV
rds_file <- file.path(output_dir, "01_unified_stripes.rds")
tsv_file <- file.path(output_dir, "01_unified_stripes.tsv")

saveRDS(union_df, rds_file)
write.table(union_df, tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("  Saved: %s\n", rds_file))
cat(sprintf("  Saved: %s\n", tsv_file))

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n=================================\n")
cat("Phase 1 Complete\n")
cat("=================================\n\n")

cat(sprintf("Timepoint: %s (%s)\n", TIMEPOINT, timepoint_dir))
cat(sprintf("Total stripes: %d\n\n", nrow(union_df)))

cat("By source:\n")
print(table(union_df$source))

cat("\nBy direction:\n")
print(table(union_df$direction_type))

cat("\nBy confidence:\n")
print(table(union_df$confidence))

cat(sprintf("\nOutput: %s\n", output_dir))
cat(sprintf("\nNext: Rscript scripts/phase2_quantification.R %s\n\n", TIMEPOINT))
