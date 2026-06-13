#!/usr/bin/env Rscript
# loops/scripts/build_anchor_span_signal.R
#
# Build the CUT&RUN / ATAC SIGNAL JUGGERNAUT (Sheet 2) for the paper.
#
# For every FDR<0.05 significant loop, quantify mean bigwig signal over three regions:
#   anchor1 = [start1, end1]
#   anchor2 = [start2, end2]
#   span    = [end1 + 1, start2 - 1]   <-- STRICT interior between the anchors, so span
#                                            signal is independent of anchor signal
# for six marks (K119ub, K27ac, K27me3, K4me3, ATAC, CTCF), with KO-WT deltas.
#
#   {role}_{mark}_ctrl   mean WT signal
#   {role}_{mark}_mut    mean KO signal
#   {role}_{mark}_delta  mut - ctrl
#   {role}_{mark}_log2fc log2((mut+pc)/(ctrl+pc))
#   CTCF is control-only (ENCODE Bl/6 adult cerebellum, different mice) -> {role}_CTCF_ctrl only.
#
# Method reuses the established rtracklayer import.bw + coverage(weight="score") + viewMeans
# pattern from loops/scripts/preprocess_k119ub_anchor_signal.R, generalised to every mark and
# region. Coverage is padded to full chromosome length so large spans never go out of bounds.
# Hub/shared anchors legitimately repeat across loops (rows are per-loop).
#
# Output is keyed by loop_id and joins 1:1 to master_annotated_loops.tsv.
#
# Usage (from repo root, with the loops/ Bioconductor R, e.g. mariner_env):
#   Rscript loops/scripts/build_anchor_span_signal.R
#   Rscript loops/scripts/build_anchor_span_signal.R --normalize median-ratio

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

# =============================================================================
# CONFIG — paths, marks and toggles (edit here to add/drop a mark or a region)
# =============================================================================

CONFIG <- list(
  spine_tsv   = "loops/outputs/250402-late_outputs/merged_loops/merged_all_results.tsv",
  bigwig_dir  = "/Users/zakiralibhai/sdsc/bigwigs",
  out_dir     = "loops/output/master_loop_table/late",
  out_file    = "anchor_span_signal.tsv",

  # name -> {ctrl bigwig, mut bigwig}.  mut = NA  => control-only (no delta/log2fc).
  marks = list(
    K119ub = list(ctrl = "H2AK119ubCtrl.bw", mut = "H2AK119ubMut.bw"),
    K27ac  = list(ctrl = "H3K27acCtrl.bw",   mut = "H3K27acMut.bw"),
    K27me3 = list(ctrl = "H3K27me3Ctrl.bw",  mut = "H3K27me3Mut.bw"),
    K4me3  = list(ctrl = "H3K4me3Ctrl.bw",   mut = "H3K4me3Mut.bw"),
    ATAC   = list(ctrl = "ATACctrl.bw",      mut = "ATACmut.bw"),
    CTCF   = list(ctrl = "CTCF.bw",          mut = NA)               # ctrl-only ENCODE track
  ),

  fdr_threshold = 0.05,
  normalize     = "none",   # "none" (raw, default) or "median-ratio" (scale mut to ctrl per mark)
  pc_quantile   = 0.05      # adaptive pseudocount = max(quantile(nonzero, pc_quantile), 0.01)
)

# =============================================================================
# CLI OVERRIDES
# =============================================================================

parse_cli <- function(cfg) {
  args <- commandArgs(trailingOnly = TRUE)
  flag_map <- c("--spine" = "spine_tsv", "--bigwig-dir" = "bigwig_dir",
                "--out-dir" = "out_dir", "--out-file" = "out_file",
                "--fdr" = "fdr_threshold", "--normalize" = "normalize",
                "--pc-quantile" = "pc_quantile")
  numeric_keys <- c("fdr_threshold", "pc_quantile")
  i <- 1
  while (i <= length(args)) {
    key <- flag_map[[args[i]]]
    if (is.null(key) || is.na(key)) stop(sprintf("Unknown argument: %s", args[i]))
    if (i == length(args)) stop(sprintf("Missing value for %s", args[i]))
    cfg[[key]] <- if (key %in% numeric_keys) as.numeric(args[i + 1]) else args[i + 1]
    i <- i + 2
  }
  if (!cfg$normalize %in% c("none", "median-ratio"))
    stop("--normalize must be 'none' or 'median-ratio'")
  cfg
}

# =============================================================================
# HELPERS
# =============================================================================

fmt_int <- function(x) sprintf("%.0f", as.numeric(x))
key3    <- function(chr, s, e) paste(chr, fmt_int(s), fmt_int(e), sep = ":")
lookup  <- function(named_vec, keys) unname(named_vec[keys])

require_file <- function(path, what) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", what, path))
  invisible(path)
}

# Spine loader identical to build_master_loop_table.R so the two sheets align row-for-row.
load_significant_spine <- function(path, fdr) {
  require_file(path, "Spine (merged_all_results.tsv)")
  df <- read.delim(path, header = TRUE, sep = "\t", quote = "",
                   comment.char = "", stringsAsFactors = FALSE)
  need <- c("loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2", "FDR", "significant")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop(sprintf("Spine missing columns: %s", paste(miss, collapse = ", ")))
  is_sig <- as.character(df$significant) == "TRUE" & !is.na(df$FDR) & df$FDR < fdr
  sig <- df[is_sig, , drop = FALSE]
  for (cc in c("start1", "end1", "start2", "end2")) sig[[cc]] <- as.integer(sig[[cc]])
  cat(sprintf("  Spine: %d total -> %d significant (FDR < %g)\n", nrow(df), nrow(sig), fdr))
  sig
}

# Mean bigwig signal per region. Missing data counts as 0 (deepTools --missingDataAsZero);
# coverage is padded to full chromosome length so views never run out of bounds.
compute_mean_signal <- function(bw_path, regions) {
  si <- seqinfo(BigWigFile(bw_path))
  bw <- import.bw(bw_path, which = regions)
  seqlevels(bw)  <- seqlevels(si)
  seqlengths(bw) <- seqlengths(si)
  cov <- coverage(bw, weight = "score")

  res <- numeric(length(regions))
  chr_names <- as.character(seqnames(regions))
  rbc <- split(ranges(regions), seqnames(regions))
  for (chr in names(rbc)) {
    rr <- rbc[[chr]]
    if (length(rr) == 0 || !chr %in% names(cov)) next
    res[which(chr_names == chr)] <- viewMeans(Views(cov[[chr]], rr))
  }
  res
}

adaptive_pseudocount <- function(values, q) {
  nz <- values[values > 0 & !is.na(values)]
  if (length(nz) == 0) return(0.01)
  max(as.numeric(quantile(nz, q)), 0.01)
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  cfg <- parse_cli(CONFIG)
  cat("================================================================================\n")
  cat("BUILD ANCHOR / SPAN SIGNAL JUGGERNAUT\n")
  cat(sprintf("  normalize = %s\n", cfg$normalize))
  cat("================================================================================\n")
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- spine + region geometry --------------------------------------------
  cat("\n[1] Spine + regions\n")
  loops <- load_significant_spine(cfg$spine_tsv, cfg$fdr_threshold)
  if (any(loops$chr1 != loops$chr2))
    stop("Trans-chromosomal loop(s) in significant set — span undefined; investigate.")

  span_start <- loops$end1 + 1L
  span_end   <- loops$start2 - 1L
  if (any(span_start > span_end))
    stop("Degenerate span (anchors adjacent/overlapping) — investigate before quantifying.")

  # Per-loop region keys for the three roles.
  roles <- list(
    anchor1 = list(chr = loops$chr1, start = loops$start1, end = loops$end1),
    anchor2 = list(chr = loops$chr2, start = loops$start2, end = loops$end2),
    span    = list(chr = loops$chr1, start = span_start,   end = span_end)
  )
  role_keys <- lapply(roles, function(r) key3(r$chr, r$start, r$end))

  # Unique region set across all roles (hub anchors collapse to one computation).
  all_df <- do.call(rbind, lapply(roles, function(r)
    data.frame(chr = r$chr, start = r$start, end = r$end, stringsAsFactors = FALSE)))
  all_df$key <- key3(all_df$chr, all_df$start, all_df$end)
  uniq <- all_df[!duplicated(all_df$key), , drop = FALSE]
  uniq_gr <- GRanges(uniq$chr, IRanges(uniq$start, uniq$end))
  cat(sprintf("  Loops=%d  unique regions to quantify=%d\n", nrow(loops), nrow(uniq)))

  # --- per-mark signal over unique regions --------------------------------
  cat("\n[2] Quantifying bigwig signal\n")
  out <- data.frame(
    loop_id = loops$loop_id,
    chr1 = loops$chr1, start1 = loops$start1, end1 = loops$end1,
    chr2 = loops$chr2, start2 = loops$start2, end2 = loops$end2,
    span_start = span_start, span_end = span_end,
    span_length = loops$start2 - loops$end1,
    stringsAsFactors = FALSE
  )

  for (mark in names(cfg$marks)) {
    spec <- cfg$marks[[mark]]
    ctrl_path <- file.path(cfg$bigwig_dir, spec$ctrl)
    require_file(ctrl_path, sprintf("%s ctrl bigwig", mark))
    cat(sprintf("  [%s] ctrl=%s\n", mark, spec$ctrl))
    ctrl_u <- setNames(compute_mean_signal(ctrl_path, uniq_gr), uniq$key)

    has_mut <- !is.null(spec$mut) && !is.na(spec$mut)
    if (has_mut) {
      mut_path <- file.path(cfg$bigwig_dir, spec$mut)
      require_file(mut_path, sprintf("%s mut bigwig", mark))
      cat(sprintf("  [%s] mut =%s\n", mark, spec$mut))
      mut_u <- setNames(compute_mean_signal(mut_path, uniq_gr), uniq$key)

      if (cfg$normalize == "median-ratio") {
        sf <- median(ctrl_u[ctrl_u > 0]) / median(mut_u[mut_u > 0])
        mut_u <- mut_u * sf
        cat(sprintf("       median-ratio scale (ctrl/mut) = %.4f\n", sf))
      }
      pc <- adaptive_pseudocount(c(ctrl_u, mut_u), cfg$pc_quantile)
    }

    for (role in names(roles)) {
      rk <- role_keys[[role]]
      ctrl_v <- lookup(ctrl_u, rk)
      out[[sprintf("%s_%s_ctrl", role, mark)]] <- round(ctrl_v, 6)
      if (has_mut) {
        mut_v <- lookup(mut_u, rk)
        out[[sprintf("%s_%s_mut",    role, mark)]] <- round(mut_v, 6)
        out[[sprintf("%s_%s_delta",  role, mark)]] <- round(mut_v - ctrl_v, 6)
        out[[sprintf("%s_%s_log2fc", role, mark)]] <- round(log2((mut_v + pc) / (ctrl_v + pc)), 6)
      }
    }
  }

  # --- write ---------------------------------------------------------------
  out_path <- file.path(cfg$out_dir, cfg$out_file)
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n--------------------------------------------------------------------------------\n")
  cat(sprintf("  Wrote %s  (%d rows x %d cols)\n", out_path, nrow(out), ncol(out)))

  # quick non-zero summary per mark (anchor1) as a sanity check
  cat("  anchor1 non-zero fraction (ctrl):\n")
  for (mark in names(cfg$marks)) {
    col <- sprintf("anchor1_%s_ctrl", mark)
    cat(sprintf("    %-8s %.3f\n", mark, mean(out[[col]] > 0)))
  }
  cat("Done.\n")
}

if (!interactive()) main()
