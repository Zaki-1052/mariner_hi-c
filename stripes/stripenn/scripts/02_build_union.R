#!/usr/bin/env Rscript
# stripes/stripenn/scripts/02_build_union.R
# Stage 2: Build union stripe set from ctrl/mut merged stripenn calls.
# Ported from stripes/scripts/phase1_detection.R — adapted for stripenn's
# result_filtered.tsv format and simplified (no replicate counting or 10kb
# validation; those are handled by Stages 3/6 respectively).
#
# Usage: Rscript 02_build_union.R <timepoint> <resolution_bp>
#   <timepoint>     : 250831 | 250402
#   <resolution_bp> : 5000 | 10000

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(yaml)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 02_build_union.R <timepoint> <resolution_bp>")
}
TIMEPOINT <- args[1]
RES <- as.numeric(args[2])
RES_KB <- RES / 1000

CODE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"

cat("\n========================================\n")
cat("Stage 2: Build Union Stripe Set\n")
cat(sprintf("Timepoint:  %s\n", TIMEPOINT))
cat(sprintf("Resolution: %d bp (%g kb)\n", RES, RES_KB))
cat("========================================\n\n")

config_file <- file.path(CODE_DIR, "config", "stripenn_config.yaml")
if (!file.exists(config_file)) {
  stop(sprintf("Config file not found: %s", config_file))
}
config <- yaml::read_yaml(config_file)
cat(sprintf("Config loaded: %s\n", config_file))

tolerance <- config$detection$anchor_tolerance_bp
exclude_chroms <- config$filtering$exclude_chromosomes
cat(sprintf("Anchor tolerance: %d bp (%d kb)\n", tolerance, tolerance / 1000))
cat(sprintf("Exclude chromosomes: %s\n", paste(exclude_chroms, collapse = ", ")))

# ==============================================================================
# FUNCTION: Parse stripenn result_filtered.tsv
# ==============================================================================
parse_stripenn_result <- function(filepath, source_name) {
  if (!file.exists(filepath)) {
    stop(sprintf("File not found: %s", filepath))
  }

  df <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE)
  cat(sprintf("  Read %d stripes from %s\n", nrow(df), basename(filepath)))

  if (nrow(df) == 0) return(NULL)

  expected_cols <- c("chr", "pos1", "pos2", "chr2", "pos3", "pos4",
                     "length", "width", "Mean", "maxpixel", "pvalue", "Stripiness")
  missing <- setdiff(expected_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns in %s: %s", filepath, paste(missing, collapse = ", ")))
  }

  df <- df[df$chr == df$chr2, ]
  if (nrow(df) == 0) return(NULL)

  if (length(exclude_chroms) > 0) {
    n_before <- nrow(df)
    df <- df[!df$chr %in% exclude_chroms, ]
    n_removed <- n_before - nrow(df)
    if (n_removed > 0) {
      cat(sprintf("  Filtered %d stripes on excluded chromosomes\n", n_removed))
    }
  }
  if (nrow(df) == 0) return(NULL)

  df$direction <- case_when(
    df$pos3 >= df$pos2 ~ "3prime",
    df$pos4 <= df$pos1 ~ "5prime",
    TRUE ~ "ambiguous"
  )

  ambig_idx <- df$direction == "ambiguous"
  if (any(ambig_idx)) {
    anchor_mid <- (df$pos1[ambig_idx] + df$pos2[ambig_idx]) / 2
    span_mid <- (df$pos3[ambig_idx] + df$pos4[ambig_idx]) / 2
    df$direction[ambig_idx] <- ifelse(span_mid > anchor_mid, "3prime", "5prime")
  }

  df$anchor_center <- (df$pos1 + df$pos2) / 2
  df$source_file <- source_name

  df
}

# ==============================================================================
# FUNCTION: Convert to GRanges for anchor matching
# ==============================================================================
stripes_to_granges <- function(stripe_df) {
  if (is.null(stripe_df) || nrow(stripe_df) == 0) return(GRanges())

  GRanges(
    seqnames = stripe_df$chr,
    ranges = IRanges(start = stripe_df$pos1, end = stripe_df$pos2),
    direction = stripe_df$direction,
    pval = stripe_df$pvalue,
    stripiness = stripe_df$Stripiness,
    mean_signal = stripe_df$Mean,
    span_y1 = stripe_df$pos3,
    span_y2 = stripe_df$pos4,
    stripe_length = stripe_df$length,
    stripe_width = stripe_df$width,
    anchor_center = stripe_df$anchor_center
  )
}

# ==============================================================================
# FUNCTION: Match anchors with direction awareness (50kb tolerance)
# ==============================================================================
match_anchors <- function(gr_ctrl, gr_mut, tolerance) {
  if (length(gr_ctrl) == 0 || length(gr_mut) == 0) {
    return(data.frame(ctrl_idx = integer(0), mut_idx = integer(0),
                      same_direction = logical(0)))
  }

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
# FUNCTION: Build union set
# ==============================================================================
build_union_set <- function(ctrl_df, mut_df, tolerance) {
  gr_ctrl <- stripes_to_granges(ctrl_df)
  gr_mut <- stripes_to_granges(mut_df)

  cat(sprintf("  Control merged: %d stripes\n", length(gr_ctrl)))
  cat(sprintf("  Mutant merged:  %d stripes\n", length(gr_mut)))

  matches <- match_anchors(gr_ctrl, gr_mut, tolerance)

  shared <- matches[matches$same_direction, ]
  shared <- shared[!duplicated(shared$ctrl_idx), ]
  shared <- shared[!duplicated(shared$mut_idx), ]

  ctrl_shared <- unique(shared$ctrl_idx)
  mut_shared <- unique(shared$mut_idx)
  ctrl_only <- setdiff(seq_len(length(gr_ctrl)), ctrl_shared)
  mut_only <- setdiff(seq_len(length(gr_mut)), mut_shared)

  cat(sprintf("  Shared: %d, Control-only: %d, Mutant-only: %d\n",
              nrow(shared), length(ctrl_only), length(mut_only)))

  union_list <- list()

  if (nrow(shared) > 0) {
    for (i in seq_len(nrow(shared))) {
      ci <- shared$ctrl_idx[i]
      mi <- shared$mut_idx[i]
      union_list[[length(union_list) + 1]] <- data.frame(
        chr = as.character(seqnames(gr_ctrl)[ci]),
        pos1 = start(gr_ctrl)[ci], pos2 = end(gr_ctrl)[ci],
        chr2 = as.character(seqnames(gr_ctrl)[ci]),
        pos3 = gr_ctrl$span_y1[ci], pos4 = gr_ctrl$span_y2[ci],
        direction_type = gr_ctrl$direction[ci],
        source = "shared",
        pval_ctrl = gr_ctrl$pval[ci], pval_mut = gr_mut$pval[mi],
        stripiness_ctrl = gr_ctrl$stripiness[ci], stripiness_mut = gr_mut$stripiness[mi],
        mean_ctrl = gr_ctrl$mean_signal[ci], mean_mut = gr_mut$mean_signal[mi],
        stripe_length = gr_ctrl$stripe_length[ci],
        stripe_width = gr_ctrl$stripe_width[ci],
        anchor_center = gr_ctrl$anchor_center[ci],
        stringsAsFactors = FALSE
      )
    }
  }

  for (i in ctrl_only) {
    union_list[[length(union_list) + 1]] <- data.frame(
      chr = as.character(seqnames(gr_ctrl)[i]),
      pos1 = start(gr_ctrl)[i], pos2 = end(gr_ctrl)[i],
      chr2 = as.character(seqnames(gr_ctrl)[i]),
      pos3 = gr_ctrl$span_y1[i], pos4 = gr_ctrl$span_y2[i],
      direction_type = gr_ctrl$direction[i],
      source = "control_only",
      pval_ctrl = gr_ctrl$pval[i], pval_mut = NA_real_,
      stripiness_ctrl = gr_ctrl$stripiness[i], stripiness_mut = NA_real_,
      mean_ctrl = gr_ctrl$mean_signal[i], mean_mut = NA_real_,
      stripe_length = gr_ctrl$stripe_length[i],
      stripe_width = gr_ctrl$stripe_width[i],
      anchor_center = gr_ctrl$anchor_center[i],
      stringsAsFactors = FALSE
    )
  }

  for (i in mut_only) {
    union_list[[length(union_list) + 1]] <- data.frame(
      chr = as.character(seqnames(gr_mut)[i]),
      pos1 = start(gr_mut)[i], pos2 = end(gr_mut)[i],
      chr2 = as.character(seqnames(gr_mut)[i]),
      pos3 = gr_mut$span_y1[i], pos4 = gr_mut$span_y2[i],
      direction_type = gr_mut$direction[i],
      source = "mutant_only",
      pval_ctrl = NA_real_, pval_mut = gr_mut$pval[i],
      stripiness_ctrl = NA_real_, stripiness_mut = gr_mut$stripiness[i],
      mean_ctrl = NA_real_, mean_mut = gr_mut$mean_signal[i],
      stripe_length = gr_mut$stripe_length[i],
      stripe_width = gr_mut$stripe_width[i],
      anchor_center = gr_mut$anchor_center[i],
      stringsAsFactors = FALSE
    )
  }

  if (length(union_list) == 0) stop("No stripes in union set")

  union_df <- do.call(rbind, union_list)
  union_df <- union_df[order(union_df$chr, union_df$pos1), ]
  union_df$stripe_id <- paste0("stripe_", sprintf("%04d", seq_len(nrow(union_df))))

  union_df
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
calls_dir <- file.path(DATA_DIR, "outputs", TIMEPOINT, paste0("res_", RES_KB, "kb"), "calls")
ctrl_result <- file.path(calls_dir, "ctrl_merged", "result_filtered.tsv")
mut_result <- file.path(calls_dir, "mut_merged", "result_filtered.tsv")

cat("\n=== Loading Stripenn Results ===\n")
cat(sprintf("  ctrl: %s\n", ctrl_result))
cat(sprintf("  mut:  %s\n", mut_result))

ctrl_merged <- parse_stripenn_result(ctrl_result, "ctrl_merged")
mut_merged <- parse_stripenn_result(mut_result, "mut_merged")

if (is.null(ctrl_merged) || is.null(mut_merged)) {
  stop("Failed to load result_filtered.tsv files (empty after filtering?)")
}

cat(sprintf("\n  Loaded ctrl: %d stripes (after filtering)\n", nrow(ctrl_merged)))
cat(sprintf("  Loaded mut:  %d stripes (after filtering)\n", nrow(mut_merged)))

cat("\n=== Building Union Set ===\n")
union_df <- build_union_set(ctrl_merged, mut_merged, tolerance)
cat(sprintf("Total union: %d stripes\n", nrow(union_df)))

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

output_dir <- file.path(DATA_DIR, "outputs", TIMEPOINT, paste0("res_", RES_KB, "kb"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tsv_path <- file.path(output_dir, "union_stripes.tsv")
write.table(union_df, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Full TSV: %s (%d rows)\n", tsv_path, nrow(union_df)))

rds_path <- file.path(output_dir, "union_stripes.rds")
saveRDS(union_df, rds_path)
cat(sprintf("  RDS:      %s\n", rds_path))

bedpe_path <- file.path(output_dir, "union_stripes.bedpe")
bedpe_df <- union_df[, c("chr", "pos1", "pos2", "chr2", "pos3", "pos4")]
write.table(bedpe_df, bedpe_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  BEDPE:    %s (%d rows)\n", bedpe_path, nrow(bedpe_df)))

cat("\n=== Summary ===\n")
cat(sprintf("  Shared:       %d\n", sum(union_df$source == "shared")))
cat(sprintf("  Control-only: %d\n", sum(union_df$source == "control_only")))
cat(sprintf("  Mutant-only:  %d\n", sum(union_df$source == "mutant_only")))
cat(sprintf("  Total:        %d\n", nrow(union_df)))

dir_tab <- table(union_df$direction_type)
cat("\nDirection distribution:\n")
for (d in names(dir_tab)) {
  cat(sprintf("  %s: %d\n", d, dir_tab[d]))
}

source_dir_tab <- table(union_df$source, union_df$direction_type)
cat("\nSource x Direction:\n")
print(source_dir_tab)

cat(sprintf("\nMedian stripe length: %d bp\n", median(union_df$stripe_length)))
cat(sprintf("Median stripe width:  %d bp\n", median(union_df$stripe_width)))

cat("\n========================================\n")
cat("Stage 2 complete.\n")
cat(sprintf("End: %s\n", Sys.time()))
cat("========================================\n")
