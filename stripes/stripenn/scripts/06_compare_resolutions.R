#!/usr/bin/env Rscript
# stripes/stripenn/scripts/06_compare_resolutions.R
# Stage 6: Cross-resolution comparison (5kb vs 10kb) per timepoint.
# Matches stripes by anchor overlap, identifies high-confidence set,
# and generates concordance metrics and plots.
#
# Usage: Rscript 06_compare_resolutions.R <timepoint>
#   <timepoint> : 250831 | 250402

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(yaml)
  library(dplyr)
  library(svglite)
  library(VennDiagram)
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
if (length(args) < 1) {
  stop("Usage: Rscript 06_compare_resolutions.R <timepoint>")
}
TIMEPOINT <- args[1]

CODE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"

cat("\n========================================\n")
cat("Stage 6: Cross-Resolution Comparison\n")
cat(sprintf("Timepoint: %s\n", TIMEPOINT))
cat("========================================\n\n")

config_file <- file.path(CODE_DIR, "config", "stripenn_config.yaml")
config <- yaml::read_yaml(config_file)

tolerance <- config$detection$anchor_tolerance_bp
tp_base <- file.path(DATA_DIR, "outputs", TIMEPOINT)
plots_dir <- file.path(tp_base, "cross_res_plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# ==============================================================================
# LOAD BOTH RESOLUTIONS
# ==============================================================================
cat("=== Loading Results ===\n")

res_5kb_file <- file.path(tp_base, "res_5kb", "05_final_differential.tsv")
res_10kb_file <- file.path(tp_base, "res_10kb", "05_final_differential.tsv")

if (!file.exists(res_5kb_file)) stop(sprintf("5kb results not found: %s", res_5kb_file))
if (!file.exists(res_10kb_file)) stop(sprintf("10kb results not found: %s", res_10kb_file))

res_5kb <- read.delim(res_5kb_file, header = TRUE, stringsAsFactors = FALSE)
res_10kb <- read.delim(res_10kb_file, header = TRUE, stringsAsFactors = FALSE)

cat(sprintf("  5kb:  %d stripes\n", nrow(res_5kb)))
cat(sprintf("  10kb: %d stripes\n", nrow(res_10kb)))

# ==============================================================================
# BASIC STATISTICS PER RESOLUTION
# ==============================================================================
cat("\n=== Per-Resolution Statistics ===\n")

for (res_label in c("5kb", "10kb")) {
  df <- if (res_label == "5kb") res_5kb else res_10kb
  n_total <- nrow(df)
  n_sig <- sum(df$significant_FDR05, na.rm = TRUE)
  n_lost <- sum(df$direction == "lost", na.rm = TRUE)
  n_gained <- sum(df$direction == "gained", na.rm = TRUE)
  n_str <- sum(df$direction == "strengthened", na.rm = TRUE)
  n_wk <- sum(df$direction == "weakened", na.rm = TRUE)
  cat(sprintf("\n  %s: %d total, %d sig FDR<0.05\n", res_label, n_total, n_sig))
  cat(sprintf("    lost=%d, gained=%d, strengthened=%d, weakened=%d\n",
              n_lost, n_gained, n_str, n_wk))
}

# ==============================================================================
# ANCHOR-BASED MATCHING
# ==============================================================================
cat("\n=== Cross-Resolution Matching ===\n")

gr_5kb <- GRanges(
  seqnames = res_5kb$chr,
  ranges = IRanges(start = res_5kb$pos1, end = res_5kb$pos2),
  stripe_id = res_5kb$stripe_id,
  direction = res_5kb$direction,
  logFC = res_5kb$logFC,
  FDR = res_5kb$FDR,
  source = res_5kb$source,
  significant = res_5kb$significant_FDR05
)

gr_10kb <- GRanges(
  seqnames = res_10kb$chr,
  ranges = IRanges(start = res_10kb$pos1, end = res_10kb$pos2),
  stripe_id = res_10kb$stripe_id,
  direction = res_10kb$direction,
  logFC = res_10kb$logFC,
  FDR = res_10kb$FDR,
  source = res_10kb$source,
  significant = res_10kb$significant_FDR05
)

gr_5kb_ext <- resize(gr_5kb, width(gr_5kb) + 2 * tolerance, fix = "center")
gr_10kb_ext <- resize(gr_10kb, width(gr_10kb) + 2 * tolerance, fix = "center")

overlaps <- findOverlaps(gr_5kb_ext, gr_10kb_ext)

if (length(overlaps) > 0) {
  match_df <- data.frame(
    idx_5kb = queryHits(overlaps),
    idx_10kb = subjectHits(overlaps)
  )
  match_df <- match_df[!duplicated(match_df$idx_5kb), ]
  match_df <- match_df[!duplicated(match_df$idx_10kb), ]
} else {
  match_df <- data.frame(idx_5kb = integer(0), idx_10kb = integer(0))
}

n_matched <- nrow(match_df)
n_5kb_only <- nrow(res_5kb) - n_matched
n_10kb_only <- nrow(res_10kb) - n_matched

cat(sprintf("  Matched:  %d stripes\n", n_matched))
cat(sprintf("  5kb-only: %d stripes\n", n_5kb_only))
cat(sprintf("  10kb-only:%d stripes\n", n_10kb_only))

# ==============================================================================
# CONCORDANCE ANALYSIS
# ==============================================================================
cat("\n=== Concordance Analysis ===\n")

if (n_matched > 0) {
  matched_5kb <- res_5kb[match_df$idx_5kb, ]
  matched_10kb <- res_10kb[match_df$idx_10kb, ]

  both_sig <- matched_5kb$significant_FDR05 & matched_10kb$significant_FDR05
  n_both_sig <- sum(both_sig, na.rm = TRUE)

  same_direction <- matched_5kb$direction == matched_10kb$direction
  n_concordant <- sum(same_direction & both_sig, na.rm = TRUE)

  fc_cor <- cor(matched_5kb$logFC, matched_10kb$logFC, use = "complete.obs")

  cat(sprintf("  Both significant: %d / %d matched\n", n_both_sig, n_matched))
  cat(sprintf("  Direction concordant (both sig): %d / %d\n", n_concordant, n_both_sig))
  cat(sprintf("  logFC correlation: %.3f\n", fc_cor))
} else {
  n_both_sig <- 0
  n_concordant <- 0
  fc_cor <- NA
  cat("  No matched stripes — skipping concordance\n")
}

# ==============================================================================
# HIGH-CONFIDENCE SET
# ==============================================================================
cat("\n=== Building High-Confidence Set ===\n")

cross_res <- res_5kb
cross_res$in_10kb <- FALSE
cross_res$logFC_10kb <- NA_real_
cross_res$FDR_10kb <- NA_real_
cross_res$direction_10kb <- NA_character_
cross_res$resolution_support <- "5kb_only"

if (n_matched > 0) {
  cross_res$in_10kb[match_df$idx_5kb] <- TRUE
  cross_res$logFC_10kb[match_df$idx_5kb] <- matched_10kb$logFC
  cross_res$FDR_10kb[match_df$idx_5kb] <- matched_10kb$FDR
  cross_res$direction_10kb[match_df$idx_5kb] <- matched_10kb$direction

  concordant_mask <- cross_res$in_10kb &
    !is.na(cross_res$direction_10kb) &
    cross_res$direction == cross_res$direction_10kb
  cross_res$resolution_support[concordant_mask] <- "both_concordant"

  discordant_mask <- cross_res$in_10kb &
    !is.na(cross_res$direction_10kb) &
    cross_res$direction != cross_res$direction_10kb
  cross_res$resolution_support[discordant_mask] <- "both_discordant"
}

# Add 10kb-only stripes (not matched to any 5kb stripe)
matched_10kb_idx <- match_df$idx_10kb
unmatched_10kb_idx <- setdiff(seq_len(nrow(res_10kb)), matched_10kb_idx)

if (length(unmatched_10kb_idx) > 0) {
  extra_10kb <- res_10kb[unmatched_10kb_idx, ]
  extra_10kb$in_10kb <- TRUE
  extra_10kb$logFC_10kb <- extra_10kb$logFC
  extra_10kb$FDR_10kb <- extra_10kb$FDR
  extra_10kb$direction_10kb <- extra_10kb$direction
  extra_10kb$resolution_support <- "10kb_only"

  common_cols <- intersect(colnames(cross_res), colnames(extra_10kb))
  cross_res <- rbind(cross_res[, common_cols], extra_10kb[, common_cols])
}

cat(sprintf("  both_concordant: %d\n", sum(cross_res$resolution_support == "both_concordant")))
cat(sprintf("  both_discordant: %d\n", sum(cross_res$resolution_support == "both_discordant")))
cat(sprintf("  5kb_only:        %d\n", sum(cross_res$resolution_support == "5kb_only")))
cat(sprintf("  10kb_only:       %d\n", sum(cross_res$resolution_support == "10kb_only")))
cat(sprintf("  Total:           %d\n", nrow(cross_res)))

# ==============================================================================
# PLOTS
# ==============================================================================
cat("\n=== Generating Plots ===\n")

if (n_matched > 0 && !is.na(fc_cor)) {
  save_multiformat(quote({
    par(mar = c(5, 5, 4, 2))
    plot(matched_5kb$logFC, matched_10kb$logFC,
         pch = 16, cex = 0.8,
         col = adjustcolor(ifelse(both_sig, "firebrick3", "gray60"), alpha.f = 0.7),
         xlab = "logFC (5kb)", ylab = "logFC (10kb)",
         main = sprintf("logFC Correlation (%s)\nr = %.3f, n = %d matched",
                        TIMEPOINT, fc_cor, n_matched))
    abline(0, 1, lty = 2, col = "black")
    abline(h = 0, v = 0, lty = 3, col = "gray")
    legend("topleft",
           legend = c(sprintf("Both significant (%d)", n_both_sig),
                      sprintf("Not both (%d)", n_matched - n_both_sig)),
           col = c("firebrick3", "gray60"), pch = 16, bty = "n")
  }), file.path(plots_dir, "logFC_correlation"), width = 8, height = 8)
}

n_sig_5kb <- sum(res_5kb$significant_FDR05, na.rm = TRUE)
n_sig_10kb <- sum(res_10kb$significant_FDR05, na.rm = TRUE)

save_multiformat(quote({
  par(mar = c(5, 5, 4, 2))
  bardata <- matrix(c(
    sum(res_5kb$direction == "lost", na.rm = TRUE),
    sum(res_10kb$direction == "lost", na.rm = TRUE),
    sum(res_5kb$direction == "gained", na.rm = TRUE),
    sum(res_10kb$direction == "gained", na.rm = TRUE),
    sum(res_5kb$direction == "strengthened", na.rm = TRUE),
    sum(res_10kb$direction == "strengthened", na.rm = TRUE),
    sum(res_5kb$direction == "weakened", na.rm = TRUE),
    sum(res_10kb$direction == "weakened", na.rm = TRUE)
  ), nrow = 2, dimnames = list(c("5kb", "10kb"),
                                c("Lost", "Gained", "Strengthened", "Weakened")))
  barplot(bardata, beside = TRUE,
          col = c("steelblue", "coral"),
          main = sprintf("Differential Stripes by Resolution (%s)", TIMEPOINT),
          ylab = "Count", las = 1)
  legend("topright", legend = c("5kb", "10kb"),
         fill = c("steelblue", "coral"), bty = "n")
}), file.path(plots_dir, "direction_by_resolution"), width = 10, height = 7)

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================
cat("\n=== Saving Outputs ===\n")

write.table(cross_res, file.path(tp_base, "cross_res_merged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: cross_res_merged.tsv (%d stripes)\n", nrow(cross_res)))

summary_text <- sprintf("
========================================
Stage 6: Cross-Resolution Comparison (Stripenn)
Timepoint: %s
========================================

5kb stripes:  %d total, %d significant
10kb stripes: %d total, %d significant

MATCHING (anchor tolerance = %dkb):
  Matched:  %d
  5kb-only: %d
  10kb-only:%d

CONCORDANCE (matched stripes):
  Both significant:  %d / %d
  Direction concordant (both sig): %d / %d
  logFC correlation: %s

RESOLUTION SUPPORT:
  both_concordant: %d (high-confidence)
  both_discordant: %d
  5kb_only:        %d (exploratory)
  10kb_only:       %d (exploratory)
",
  TIMEPOINT,
  nrow(res_5kb), n_sig_5kb,
  nrow(res_10kb), n_sig_10kb,
  tolerance / 1000,
  n_matched, n_5kb_only, n_10kb_only,
  n_both_sig, n_matched,
  n_concordant, n_both_sig,
  ifelse(is.na(fc_cor), "N/A", sprintf("%.3f", fc_cor)),
  sum(cross_res$resolution_support == "both_concordant"),
  sum(cross_res$resolution_support == "both_discordant"),
  sum(cross_res$resolution_support == "5kb_only"),
  sum(cross_res$resolution_support == "10kb_only")
)

writeLines(summary_text, file.path(tp_base, "cross_res_summary.txt"))
cat("  Saved: cross_res_summary.txt\n")

cat("\n========================================\n")
cat("Stage 6 complete.\n")
cat(sprintf("End: %s\n", Sys.time()))
cat("========================================\n")
