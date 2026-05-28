# ML/cmpts/scripts/A3_epigenomic_validation.R
# Stage A3: Epigenomic validation of CALDER2 subcompartment calls.
#
# Computes fold-enrichment of 9 epigenomic marks per subcompartment (A.1, A.2,
# B.1, B.2) for both timepoints. Produces ctrl validation, mut validation, and
# differential heatmaps. Replicates SNIPER paper Figure 2c approach.
#
# Processes both timepoints in one invocation — BigWigs are loaded once per TP
# since the same marks are shared across timepoints.
#
# Usage:
#   Rscript A3_epigenomic_validation.R <data_root> <code_root>
#     <data_root> : HPC data directory (e.g. /expanse/.../sniper)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript A3_epigenomic_validation.R <data_root> <code_root>")
}

DATA_ROOT <- args[1]
CODE_ROOT <- args[2]

# ── Constants ──────────────────────────────────────────────────────────────────

TPS <- c("250402", "250831")
TP_LABELS <- c("250402" = "late/adult", "250831" = "early/P12")

BIGWIG_DIR <- "/expanse/lustre/projects/csd940/zalibhai/bigwigs"

LABEL_ORDER <- c("A.1", "A.2", "B.1", "B.2")
LABEL_COLORS <- c(
  "A.1" = "#e41a1c",
  "A.2" = "#ff7f00",
  "B.1" = "#4daf4a",
  "B.2" = "#377eb8"
)

MARKS <- data.table::data.table(
  mark     = c("H3K27ac", "H3K4me3", "ATAC", "H3K36me3", "RNA",
               "DNAmethylation", "H3K27me1", "H2AK119ub", "H3K27me3"),
  ctrl_bw  = c("H3K27acCtrl.bw", "H3K4me3Ctrl.bw", "ATACctrl.bw",
               "H3K36me3Ctrl.bw", "RNActrl.bw", "DNAmethylationCtrl.bw",
               "H3K27me1Ctrl.bw", "H2AK119ubCtrl.bw", "H3K27me3Ctrl.bw"),
  mut_bw   = c("H3K27acMut.bw", "H3K4me3Mut.bw", "ATACmut.bw",
               "H3K36me3Mut.bw", "RNAmut.bw", "DNAmethylationMut.bw",
               "H3K27me1Mut.bw", "H2AK119ubMut.bw", "H3K27me3Mut.bw"),
  category = c("Active", "Active", "Active", "Gene body", "Gene body",
               "Methylation", "Methylation", "Repressive", "Repressive")
)

MARK_ORDER <- MARKS$mark

OUT_DIR   <- file.path(CODE_ROOT, "outputs", "calder2")
UTIL_PATH <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

# ── Libraries ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
  library(IRanges)
  library(ggplot2)
})
source(UTIL_PATH)

# ── Header ─────────────────────────────────────────────────────────────────────

start_time <- proc.time()

cat("===========================================\n")
cat("A3: Epigenomic Validation\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("BIGWIG_DIR: %s\n", BIGWIG_DIR))
cat(sprintf("Timepoints: %s\n", paste(TPS, collapse = ", ")))
cat(sprintf("Marks:      %d marks x 2 conditions = %d BigWigs\n",
            nrow(MARKS), 2L * nrow(MARKS)))
cat(sprintf("Output dir: %s\n", OUT_DIR))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# ── Pre-flight validation ──────────────────────────────────────────────────────

cat("=== Pre-flight validation ===\n")

all_bw_paths <- c(
  file.path(BIGWIG_DIR, MARKS$ctrl_bw),
  file.path(BIGWIG_DIR, MARKS$mut_bw)
)
missing_bws <- all_bw_paths[!file.exists(all_bw_paths)]
if (length(missing_bws) > 0) {
  stop(sprintf("Missing BigWig files (%d):\n  %s",
               length(missing_bws), paste(missing_bws, collapse = "\n  ")))
}
cat(sprintf("  All %d BigWig files found.\n", length(all_bw_paths)))

for (tp in TPS) {
  labels_path <- file.path(OUT_DIR,
                           sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  if (!file.exists(labels_path))
    stop(sprintf("Missing A2 output: %s", labels_path))
  if (file.info(labels_path)$size == 0)
    stop(sprintf("Empty A2 output: %s", labels_path))
  cat(sprintf("  OK: %s (%s bytes)\n", basename(labels_path),
              format(file.info(labels_path)$size, big.mark = ",")))
}

if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("  Pre-flight passed.\n")

# ── Function definitions ───────────────────────────────────────────────────────

extract_bigwig_signal <- function(bw_path, bins_gr) {
  bw_gr     <- import.bw(bw_path, which = bins_gr)
  score_cov <- coverage(bw_gr, weight = "score")

  result  <- numeric(length(bins_gr))
  chr_vec <- as.character(seqnames(bins_gr))

  for (chr in unique(chr_vec)) {
    idx        <- which(chr_vec == chr)
    chr_ranges <- ranges(bins_gr[idx])

    if (chr %in% names(score_cov) && length(score_cov[[chr]]) > 0) {
      chr_cov <- score_cov[[chr]]
      max_end <- max(end(chr_ranges))
      if (length(chr_cov) < max_end)
        chr_cov <- c(chr_cov, Rle(0L, max_end - length(chr_cov)))
      result[idx] <- viewMeans(Views(chr_cov, chr_ranges), na.rm = TRUE)
    }
  }

  result
}

compute_enrichment <- function(signals_dt, label_col, mark_cols, mark_names) {
  rbindlist(lapply(seq_along(mark_cols), function(i) {
    col           <- mark_cols[i]
    vals          <- signals_dt[[col]]
    genome_median <- median(vals, na.rm = TRUE)

    if (genome_median == 0)
      warning(sprintf("Genome-wide median is zero for %s", col))

    by_subcomp <- signals_dt[!is.na(get(label_col)),
                             .(median_signal = median(get(col), na.rm = TRUE),
                               n_bins        = .N),
                             by = label_col]
    setnames(by_subcomp, label_col, "subcompartment")

    by_subcomp[, `:=`(mark = mark_names[i], genome_median = genome_median)]

    if (genome_median > 0) {
      by_subcomp[, fold_enrichment := median_signal / genome_median]
      by_subcomp[, log2_fold := fifelse(
        median_signal > 0,
        log2(median_signal / genome_median),
        NA_real_
      )]
    } else {
      by_subcomp[, `:=`(fold_enrichment = NA_real_, log2_fold = NA_real_)]
    }

    by_subcomp
  }))
}

compute_differential <- function(signals_dt, label_col, marks_dt) {
  rbindlist(lapply(marks_dt$mark, function(m) {
    ctrl_col <- paste0(m, "_ctrl")
    mut_col  <- paste0(m, "_mut")

    by_subcomp <- signals_dt[!is.na(get(label_col)),
                             .(ctrl_median = median(get(ctrl_col), na.rm = TRUE),
                               mut_median  = median(get(mut_col),  na.rm = TRUE),
                               n_bins      = .N),
                             by = label_col]
    setnames(by_subcomp, label_col, "subcompartment")

    by_subcomp[, `:=`(
      mark    = m,
      log2_fc = fifelse(ctrl_median > 0 & mut_median > 0,
                        log2(mut_median / ctrl_median),
                        NA_real_)
    )]

    by_subcomp
  }))
}

plot_enrichment_heatmap <- function(enrichment_dt, title, mark_order,
                                    subcomp_order, fill_col = "log2_fold",
                                    label_col = "fold_enrichment",
                                    legend_title = expression(log[2]~"fold enrichment"),
                                    subtitle = NULL) {
  dt <- copy(enrichment_dt)
  dt[, mark          := factor(mark, levels = rev(mark_order))]
  dt[, subcompartment := factor(subcompartment, levels = subcomp_order)]

  fill_vals <- dt[[fill_col]]
  max_abs   <- max(abs(fill_vals[is.finite(fill_vals)]), na.rm = TRUE)

  label_vals <- dt[[label_col]]
  dt[, cell_label := fifelse(is.finite(label_vals),
                             sprintf("%.2f", label_vals), "NA")]

  ggplot(dt, aes(x = subcompartment, y = mark)) +
    geom_tile(aes(fill = .data[[fill_col]]),
              color = "white", linewidth = 0.8) +
    geom_text(aes(label = cell_label), size = 3.5) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = legend_title,
      limits = c(-max_abs, max_abs),
      na.value = "grey80"
    ) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x   = element_text(face = "bold", size = 12),
      axis.text.y   = element_text(size = 11),
      panel.grid    = element_blank(),
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle)
}

# ── Main pipeline ──────────────────────────────────────────────────────────────

for (tp in TPS) {
  tp_start <- proc.time()
  tp_label <- sprintf("%s (%s)", tp, TP_LABELS[tp])

  cat(sprintf("\n\n###################################################\n"))
  cat(sprintf("### Timepoint: %s\n", tp_label))
  cat(sprintf("###################################################\n\n"))

  # ── 1. Load subcompartment labels ────────────────────────────────────────────

  cat("=== Loading subcompartment labels ===\n")
  labels_path <- file.path(OUT_DIR,
                           sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  labels_dt <- fread(labels_path)

  callable <- labels_dt[!is.na(ctrl_label) & !is.na(mut_label)]
  cat(sprintf("  Loaded: %d bins total, %d callable (both labels present)\n",
              nrow(labels_dt), nrow(callable)))

  # ── 2. Build GRanges ────────────────────────────────────────────────────────

  bins_gr <- GRanges(
    seqnames = callable$chr,
    ranges   = IRanges(start = callable$bin_start, end = callable$bin_end)
  )

  # ── 3. Extract BigWig signals ───────────────────────────────────────────────

  cat(sprintf("\n=== Extracting signals from %d BigWigs ===\n",
              2L * nrow(MARKS)))

  for (i in seq_len(nrow(MARKS))) {
    mark <- MARKS$mark[i]

    t0 <- proc.time()
    ctrl_path <- file.path(BIGWIG_DIR, MARKS$ctrl_bw[i])
    callable[, (paste0(mark, "_ctrl")) :=
               extract_bigwig_signal(ctrl_path, bins_gr)]
    t_ctrl <- (proc.time() - t0)[3]

    t0 <- proc.time()
    mut_path <- file.path(BIGWIG_DIR, MARKS$mut_bw[i])
    callable[, (paste0(mark, "_mut")) :=
               extract_bigwig_signal(mut_path, bins_gr)]
    t_mut <- (proc.time() - t0)[3]

    cat(sprintf("  [%2d/%d] %-18s ctrl=%.1fs  mut=%.1fs\n",
                i, nrow(MARKS), mark, t_ctrl, t_mut))
  }

  # ── 4. Save per-bin signal table ────────────────────────────────────────────

  out_signals <- file.path(OUT_DIR, sprintf("%s_bin_signals.tsv", tp))
  fwrite(callable, out_signals, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("\n  Written: %s (%d rows x %d cols)\n",
              basename(out_signals), nrow(callable), ncol(callable)))

  # ── 5. Compute enrichment ───────────────────────────────────────────────────

  cat("\n=== Computing enrichment ===\n")

  ctrl_mark_cols <- paste0(MARKS$mark, "_ctrl")
  mut_mark_cols  <- paste0(MARKS$mark, "_mut")

  enrichment_ctrl <- compute_enrichment(
    callable, "ctrl_label", ctrl_mark_cols, MARKS$mark)
  enrichment_ctrl[, condition := "ctrl"]
  cat(sprintf("  Ctrl validation: %d values (label=ctrl_label, signal=ctrl)\n",
              nrow(enrichment_ctrl)))

  enrichment_mut <- compute_enrichment(
    callable, "mut_label", mut_mark_cols, MARKS$mark)
  enrichment_mut[, condition := "mut"]
  cat(sprintf("  Mut  validation: %d values (label=mut_label, signal=mut)\n",
              nrow(enrichment_mut)))

  diff_dt <- compute_differential(callable, "ctrl_label", MARKS)
  cat(sprintf("  Differential:    %d values (label=ctrl_label, signal=mut/ctrl)\n",
              nrow(diff_dt)))

  # ── 6. Save enrichment tables ───────────────────────────────────────────────

  all_enrichment <- rbind(enrichment_ctrl, enrichment_mut, fill = TRUE)
  all_enrichment[, tp := tp]
  out_enrich <- file.path(OUT_DIR, sprintf("%s_enrichment_matrix.tsv", tp))
  fwrite(all_enrichment, out_enrich, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s\n", basename(out_enrich)))

  diff_dt[, tp := tp]
  out_diff <- file.path(OUT_DIR, sprintf("%s_differential_matrix.tsv", tp))
  fwrite(diff_dt, out_diff, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s\n", basename(out_diff)))

  # ── 7. Generate heatmaps ────────────────────────────────────────────────────

  cat("\n=== Generating heatmaps ===\n")

  n_callable <- nrow(callable)

  # 7a. Ctrl validation
  p_ctrl <- plot_enrichment_heatmap(
    enrichment_ctrl,
    title      = sprintf("Ctrl Subcompartment Enrichment: %s", tp_label),
    subtitle   = sprintf("Ctrl labels x Ctrl signal (%s callable bins)",
                         format(n_callable, big.mark = ",")),
    mark_order = MARK_ORDER, subcomp_order = LABEL_ORDER
  )
  save_multiformat_ggplot(p_ctrl,
    file.path(OUT_DIR, sprintf("%s_enrichment_heatmap_ctrl", tp)),
    width = 7, height = 8)

  # 7b. Mut validation
  p_mut <- plot_enrichment_heatmap(
    enrichment_mut,
    title      = sprintf("Mut Subcompartment Enrichment: %s", tp_label),
    subtitle   = sprintf("Mut labels x Mut signal (%s callable bins)",
                         format(n_callable, big.mark = ",")),
    mark_order = MARK_ORDER, subcomp_order = LABEL_ORDER
  )
  save_multiformat_ggplot(p_mut,
    file.path(OUT_DIR, sprintf("%s_enrichment_heatmap_mut", tp)),
    width = 7, height = 8)

  # 7c. Differential (mut vs ctrl signal within ctrl-defined subcompartments)
  p_diff <- plot_enrichment_heatmap(
    diff_dt,
    title        = sprintf("Mut vs Ctrl Signal Change: %s", tp_label),
    subtitle     = "Ctrl labels x log2(mut/ctrl signal)",
    mark_order   = MARK_ORDER,
    subcomp_order = LABEL_ORDER,
    fill_col     = "log2_fc",
    label_col    = "log2_fc",
    legend_title = expression(log[2]~"(mut / ctrl)")
  )
  save_multiformat_ggplot(p_diff,
    file.path(OUT_DIR, sprintf("%s_enrichment_heatmap_diff", tp)),
    width = 7, height = 8)

  # 7d. Combined 3-panel figure
  combined_long <- rbind(
    enrichment_ctrl[, .(mark, subcompartment,
                        fill_value  = log2_fold,
                        label_value = fold_enrichment,
                        panel       = "Ctrl Enrichment")],
    enrichment_mut[, .(mark, subcompartment,
                       fill_value  = log2_fold,
                       label_value = fold_enrichment,
                       panel       = "Mut Enrichment")],
    diff_dt[, .(mark, subcompartment,
                fill_value  = log2_fc,
                label_value = log2_fc,
                panel       = "Differential (log2 FC)")]
  )
  combined_long[, mark          := factor(mark, levels = rev(MARK_ORDER))]
  combined_long[, subcompartment := factor(subcompartment, levels = LABEL_ORDER)]
  combined_long[, panel := factor(panel, levels = c("Ctrl Enrichment",
                                                     "Mut Enrichment",
                                                     "Differential (log2 FC)"))]

  max_abs_combined <- max(
    abs(combined_long$fill_value[is.finite(combined_long$fill_value)]),
    na.rm = TRUE)
  combined_long[, cell_label := fifelse(is.finite(label_value),
                                        sprintf("%.2f", label_value), "NA")]

  p_combined <- ggplot(combined_long,
                       aes(x = subcompartment, y = mark, fill = fill_value)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = cell_label), size = 2.8) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = expression(log[2]~"value"),
      limits = c(-max_abs_combined, max_abs_combined),
      na.value = "grey80"
    ) +
    facet_wrap(~ panel, nrow = 1) +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid  = element_blank(),
      strip.text  = element_text(face = "bold", size = 11),
      plot.title  = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(x = NULL, y = NULL,
         title = sprintf("Epigenomic Subcompartment Validation: %s", tp_label))

  save_multiformat_ggplot(p_combined,
    file.path(OUT_DIR, sprintf("%s_enrichment_combined", tp)),
    width = 18, height = 8)

  # ── 8. Verification ─────────────────────────────────────────────────────────

  cat("\n=== Verification ===\n")
  any_warn <- FALSE

  cat("  Ctrl fold-enrichment matrix:\n")
  ctrl_wide <- dcast(enrichment_ctrl, mark ~ subcompartment,
                     value.var = "fold_enrichment")
  ctrl_wide[, mark := factor(mark, levels = MARK_ORDER)]
  ctrl_wide <- ctrl_wide[order(mark)]
  print(ctrl_wide, row.names = FALSE)
  cat("\n")

  # H3K27ac monotonic decrease A.1 > A.2 > B.1 > B.2
  k27ac_ctrl <- enrichment_ctrl[mark == "H3K27ac"]
  k27ac_ctrl <- k27ac_ctrl[match(LABEL_ORDER, subcompartment)]
  k27ac_vals <- k27ac_ctrl$fold_enrichment
  if (all(is.finite(k27ac_vals)) && all(diff(k27ac_vals) < 0)) {
    cat(sprintf("  H3K27ac gradient OK: %s\n",
                paste(sprintf("%.2f", k27ac_vals), collapse = " > ")))
  } else {
    warning(sprintf("H3K27ac ctrl NOT monotonically decreasing: %s",
                    paste(sprintf("%.2f", k27ac_vals), collapse = ", ")))
    any_warn <- TRUE
  }

  # H3K27me3 enriched in B.1 > 1.5x
  k27me3_b1 <- enrichment_ctrl[mark == "H3K27me3" &
                                  subcompartment == "B.1", fold_enrichment]
  if (length(k27me3_b1) == 1 && is.finite(k27me3_b1)) {
    if (k27me3_b1 >= 1.5) {
      cat(sprintf("  H3K27me3 B.1 enrichment: %.2f (OK, >= 1.5)\n", k27me3_b1))
    } else {
      warning(sprintf("H3K27me3 B.1 fold=%.2f (expected >= 1.5)", k27me3_b1))
      any_warn <- TRUE
    }
  }

  # H2AK119ub enriched in B.1
  k119_b1 <- enrichment_ctrl[mark == "H2AK119ub" &
                               subcompartment == "B.1", fold_enrichment]
  if (length(k119_b1) == 1 && is.finite(k119_b1))
    cat(sprintf("  H2AK119ub B.1 enrichment: %.2f\n", k119_b1))

  # NaN/Inf check
  all_vals <- c(enrichment_ctrl$fold_enrichment,
                enrichment_mut$fold_enrichment,
                diff_dt$log2_fc)
  n_bad <- sum(!is.finite(all_vals) & !is.na(all_vals))
  if (n_bad > 0) {
    warning(sprintf("%d NaN/Inf values in enrichment results", n_bad))
    any_warn <- TRUE
  } else {
    cat("  No NaN/Inf values detected.\n")
  }

  # Output file check
  expected_outputs <- c(
    sprintf("%s_bin_signals.tsv", tp),
    sprintf("%s_enrichment_matrix.tsv", tp),
    sprintf("%s_differential_matrix.tsv", tp)
  )
  for (f in expected_outputs) {
    fpath <- file.path(OUT_DIR, f)
    if (!file.exists(fpath)) {
      warning(sprintf("Missing output: %s", f))
      any_warn <- TRUE
    }
  }

  if (!any_warn) cat("  All verification checks passed.\n")

  tp_elapsed <- (proc.time() - tp_start)[3]
  cat(sprintf("\n  Timepoint %s complete: %.1f min\n", tp, tp_elapsed / 60))
}

# ── Final summary ──────────────────────────────────────────────────────────────

total_elapsed <- (proc.time() - start_time)[3]

cat("\n===========================================\n")
cat("A3 COMPLETE\n")
cat(sprintf("Total runtime: %.1f min\n", total_elapsed / 60))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
