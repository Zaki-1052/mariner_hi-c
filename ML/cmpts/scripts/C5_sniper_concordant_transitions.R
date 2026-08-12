# ML/cmpts/scripts/C5_sniper_concordant_transitions.R
# Stage C5: SNIPER differential analysis + concordant transition heatmap.
#
# Three parts:
#   Part 1 — SNIPER-based differential subcompartment analysis (A2 equivalent)
#   Part 2 — Epigenomic validation of SNIPER calls (A3 equivalent, no BigWigs)
#   Part 3 — Concordant transition heatmap: bins where SNIPER and CALDER2 agree
#
# Processes both timepoints in one invocation.
#
# Usage:
#   Rscript C5_sniper_concordant_transitions.R <data_root> <code_root>
#     <data_root> : HPC data directory (kept for CLI consistency)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript C5_sniper_concordant_transitions.R <data_root> <code_root>")
}

DATA_ROOT <- args[1]
CODE_ROOT <- args[2]

# ── Constants ──────────────────────────────────────────────────────────────────

TPS <- c("250402", "250831")
TP_LABELS <- c("250402" = "late/adult", "250831" = "early/P12")
TP_DIRS   <- c("250402" = "late", "250831" = "early")

LABEL_ORDER <- c("A.1", "A.2", "B.1", "B.2")
LABEL_COLORS <- c(
  "A.1"    = "#e41a1c",
  "A.2"    = "#ff7f00",
  "B.1"    = "#4daf4a",
  "B.2"    = "#377eb8",
  "change" = "#cccccc"
)

SNIPER_TO_CALDER <- c(
  "A1" = "A.1", "A2" = "A.2", "B1" = "B.1", "B2" = "B.2"
)

MARK_NAMES <- c("H3K27ac", "H3K4me3", "ATAC", "H3K36me3", "RNA",
                "DNAmethylation", "H3K27me1", "H2AK119ub", "H3K27me3")

CALDER2_BASE <- file.path(CODE_ROOT, "outputs", "calder2")
SNIPER_BASE  <- file.path(CODE_ROOT, "outputs", "sniper")
OUT_DIR      <- file.path(CODE_ROOT, "outputs", "integration", "sniper_concordant")
TSV_DIR      <- file.path(OUT_DIR, "tsv")
PLOT_DIR     <- file.path(OUT_DIR, "plots")
UTIL_PATH    <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

# ── Libraries ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(scales)
  library(patchwork)
})
source(UTIL_PATH)

# ── Header ─────────────────────────────────────────────────────────────────────

dir.create(TSV_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("===========================================\n")
cat("C5: SNIPER Concordant Transitions\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("Output dir: %s\n", OUT_DIR))
cat(sprintf("Timepoints: %s\n", paste(TPS, collapse = ", ")))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# ── Pre-flight validation ──────────────────────────────────────────────────────

cat("=== Pre-flight validation ===\n")
for (tp in TPS) {
  tp_dir <- TP_DIRS[tp]
  files_needed <- list(
    sniper_ctrl = file.path(SNIPER_BASE, "predictions", tp, "ctrl_merged", "predictions.bed"),
    sniper_mut  = file.path(SNIPER_BASE, "predictions", tp, "mut_merged",  "predictions.bed"),
    calder_labels = file.path(CALDER2_BASE, tp_dir, sprintf("%s_subcompartment_labels_100kb.tsv", tp)),
    bin_signals   = file.path(CALDER2_BASE, tp_dir, sprintf("%s_bin_signals.tsv", tp)),
    concordance   = file.path(SNIPER_BASE, "tsvs", sprintf("%s_transition_concordance.tsv", tp))
  )
  for (nm in names(files_needed)) {
    path <- files_needed[[nm]]
    if (!file.exists(path)) stop(sprintf("Missing: %s (%s)", nm, path))
    if (file.info(path)$size == 0) stop(sprintf("Empty: %s (%s)", nm, path))
    cat(sprintf("  OK: %s/%s (%s bytes)\n", tp, nm,
                format(file.info(path)$size, big.mark = ",")))
  }
}
if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
cat("  Pre-flight passed.\n")

# ── Function definitions ───────────────────────────────────────────────────────

load_sniper_bed <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t",
              col.names = c("chr", "start", "end", "label", "score", "strand",
                            "thickStart", "thickEnd", "itemRgb"))
  dt[, calder_bin_start := start + 1L]
  dt[, sniper_label := SNIPER_TO_CALDER[label]]
  unknown <- dt[is.na(sniper_label)]
  if (nrow(unknown) > 0) {
    stop(sprintf("Unknown SNIPER labels in %s: %s",
                 basename(path), paste(unique(unknown$label), collapse = ", ")))
  }
  dt[, .(chr, calder_bin_start, sniper_label)]
}

compute_transition_matrix <- function(dt, from_col, to_col) {
  valid  <- dt[!is.na(get(from_col)) & !is.na(get(to_col))]
  counts <- valid[, .N, by = c(from_col, to_col)]
  mat    <- matrix(0L, nrow = 4, ncol = 4, dimnames = list(LABEL_ORDER, LABEL_ORDER))
  for (i in seq_len(nrow(counts))) {
    r  <- as.character(counts[[from_col]][i])
    cc <- as.character(counts[[to_col]][i])
    mat[r, cc] <- counts$N[i]
  }
  mat
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

plot_transition_heatmap <- function(trans_mat, title, subtitle) {
  mat_long <- as.data.table(as.table(trans_mat))
  names(mat_long) <- c("ctrl_label", "mut_label", "count")
  mat_long[, log10_count := log10(count + 1)]
  mat_long[, on_diagonal := (as.character(ctrl_label) == as.character(mut_label))]
  mat_long[, ctrl_label := factor(ctrl_label, levels = rev(LABEL_ORDER))]
  mat_long[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]

  ggplot(mat_long, aes(x = mut_label, y = ctrl_label, fill = log10_count)) +
    geom_tile(aes(color = on_diagonal), linewidth = 1.2) +
    geom_text(aes(label = format(count, big.mark = ",")), fontface = "bold", size = 5) +
    scale_fill_gradient(low = "white", high = "#2166ac", name = "log10(count+1)") +
    scale_color_manual(values = c("TRUE" = "#d7301f", "FALSE" = "grey85"), guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = "Mutant label", y = "Control label", title = title, subtitle = subtitle) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          axis.text  = element_text(face = "bold", size = 12))
}

plot_transition_sankey <- function(alluv_dt, title, subtitle) {
  dt <- copy(alluv_dt)
  dt[, flow_type := fifelse(as.character(ctrl_label) == as.character(mut_label),
                            as.character(ctrl_label), "change")]
  dt[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
  dt[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]

  ggplot(dt, aes(axis1 = ctrl_label, axis2 = mut_label, y = N)) +
    geom_alluvium(aes(fill = flow_type), width = 1/4, alpha = 0.75, knot.pos = 0.4) +
    geom_stratum(width = 1/4, fill = "white", color = "grey40", linewidth = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 4, fontface = "bold") +
    scale_fill_manual(values = LABEL_COLORS, name = "Flow") +
    scale_x_discrete(limits = c("Control", "Mutant"), expand = c(0.05, 0.05)) +
    labs(y = "Number of 100kb bins", title = title, subtitle = subtitle) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
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

plot_confirmation_heatmap <- function(rate_mat, n_mat, tp_label) {
  rate_long <- as.data.table(as.table(rate_mat))
  n_long    <- as.data.table(as.table(n_mat))
  names(rate_long) <- c("ctrl_label", "mut_label", "rate")
  names(n_long)    <- c("ctrl_label", "mut_label", "count")

  dt <- merge(rate_long, n_long, by = c("ctrl_label", "mut_label"))
  dt[, on_diagonal := (as.character(ctrl_label) == as.character(mut_label))]
  dt[, cell_label := fifelse(
    on_diagonal, "",
    fifelse(is.finite(rate) & count > 0,
            sprintf("%d\n(%.0f%%)", as.integer(count), 100 * rate),
            "0"))]
  dt[, ctrl_label := factor(ctrl_label, levels = rev(LABEL_ORDER))]
  dt[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]
  dt[on_diagonal == TRUE, rate := NA_real_]

  ggplot(dt, aes(x = mut_label, y = ctrl_label, fill = rate)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = cell_label), fontface = "bold", size = 4, lineheight = 0.85) +
    scale_fill_gradient(low = "white", high = "#1a9850", name = "Confirmation\nrate",
                        limits = c(0, 1), labels = percent,
                        na.value = "grey90") +
    scale_x_discrete(position = "top") +
    labs(x = "Mutant label", y = "Control label",
         title = sprintf("SNIPER Confirmation of CALDER2 Transitions: %s", tp_label),
         subtitle = "Fraction of CALDER2 transitions confirmed by SNIPER") +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          axis.text  = element_text(face = "bold", size = 12))
}

# ── Main pipeline ──────────────────────────────────────────────────────────────

summary_rows <- list()

for (tp in TPS) {
  tp_dir   <- TP_DIRS[tp]
  tp_label <- sprintf("%s (%s)", tp, TP_LABELS[tp])

  cat(sprintf("\n\n###################################################\n"))
  cat(sprintf("### Timepoint: %s\n", tp_label))
  cat(sprintf("###################################################\n\n"))

  # ══════════════════════════════════════════════════════════════════════════════
  # Part 1: SNIPER Differential Analysis (A2-equivalent)
  # ══════════════════════════════════════════════════════════════════════════════

  cat("=== Part 1: SNIPER Differential ===\n")

  sniper_ctrl_path <- file.path(SNIPER_BASE, "predictions", tp, "ctrl_merged", "predictions.bed")
  sniper_mut_path  <- file.path(SNIPER_BASE, "predictions", tp, "mut_merged",  "predictions.bed")

  sniper_ctrl <- load_sniper_bed(sniper_ctrl_path)
  sniper_mut  <- load_sniper_bed(sniper_mut_path)
  cat(sprintf("  Loaded SNIPER ctrl: %d bins\n", nrow(sniper_ctrl)))
  cat(sprintf("  Loaded SNIPER mut:  %d bins\n", nrow(sniper_mut)))

  sniper_diff <- merge(
    sniper_ctrl[, .(chr, calder_bin_start, sniper_ctrl = sniper_label)],
    sniper_mut[,  .(chr, calder_bin_start, sniper_mut  = sniper_label)],
    by = c("chr", "calder_bin_start")
  )
  setnames(sniper_diff, "calder_bin_start", "bin_start")
  sniper_diff[, sniper_changed := (sniper_ctrl != sniper_mut)]

  n_sniper_total   <- nrow(sniper_diff)
  n_sniper_changed <- sum(sniper_diff$sniper_changed)
  pct_sniper       <- 100 * n_sniper_changed / n_sniper_total

  cat(sprintf("  Joined bins:    %d\n", n_sniper_total))
  cat(sprintf("  Changed bins:   %d (%.1f%%)\n", n_sniper_changed, pct_sniper))

  sniper_mat <- compute_transition_matrix(sniper_diff, "sniper_ctrl", "sniper_mut")
  cat(sprintf("  Matrix total:   %d\n", sum(sniper_mat)))

  sniper_chisq <- chisq.test(sniper_mat, simulate.p.value = FALSE)
  sniper_cramer <- sqrt(sniper_chisq$statistic /
                        (sum(sniper_mat) * (min(nrow(sniper_mat), ncol(sniper_mat)) - 1)))

  cat(sprintf("  Chi-squared: X2=%.1f, df=%d, p=%.3e\n",
              sniper_chisq$statistic, sniper_chisq$parameter, sniper_chisq$p.value))
  cat(sprintf("  Cramer's V:  %.4f\n", unname(sniper_cramer)))

  # Write SNIPER differential TSVs
  fwrite(sniper_diff, file.path(TSV_DIR, sprintf("%s_sniper_differential.tsv", tp)),
         sep = "\t", quote = FALSE, na = "NA")

  write.table(sniper_mat, file.path(TSV_DIR, sprintf("%s_sniper_transition_matrix.tsv", tp)),
              sep = "\t", quote = FALSE, col.names = NA)

  sniper_summary <- as.data.table(as.table(sniper_mat))
  names(sniper_summary) <- c("ctrl_label", "mut_label", "count")
  sniper_summary[, pct_of_total := round(100 * count / n_sniper_total, 2)]
  sniper_summary[, is_diagonal := (as.character(ctrl_label) == as.character(mut_label))]
  fwrite(sniper_summary[order(ctrl_label, mut_label)],
         file.path(TSV_DIR, sprintf("%s_sniper_transition_summary.tsv", tp)),
         sep = "\t", quote = FALSE)
  cat(sprintf("  Written: 3 SNIPER differential TSVs\n"))

  # SNIPER figures
  cat("\n  Generating SNIPER figures...\n")

  p_sniper_heat <- plot_transition_heatmap(
    sniper_mat,
    title    = sprintf("SNIPER Subcompartment Transitions: %s", tp_label),
    subtitle = sprintf("X²=%.1f, p=%.2e, %.1f%% bins changed (%d/%d)",
                       sniper_chisq$statistic, sniper_chisq$p.value,
                       pct_sniper, n_sniper_changed, n_sniper_total))
  save_multiformat_ggplot(p_sniper_heat,
    file.path(PLOT_DIR, sprintf("c5_sniper_transition_heatmap_%s", tp_dir)),
    width = 7, height = 6)

  sniper_alluv <- sniper_diff[, .N, by = .(ctrl_label = sniper_ctrl, mut_label = sniper_mut)]
  p_sniper_sankey <- plot_transition_sankey(
    sniper_alluv,
    title    = sprintf("SNIPER Subcompartment Flow: %s", tp_label),
    subtitle = sprintf("%.1f%% of bins change subcompartment", pct_sniper))
  save_multiformat_ggplot(p_sniper_sankey,
    file.path(PLOT_DIR, sprintf("c5_sniper_transition_sankey_%s", tp_dir)),
    width = 8, height = 7)

  # ══════════════════════════════════════════════════════════════════════════════
  # Part 2: SNIPER Epigenomic Validation (A3-equivalent, no BigWigs)
  # ══════════════════════════════════════════════════════════════════════════════

  cat("\n=== Part 2: SNIPER Epigenomic Validation ===\n")

  bin_signals_path <- file.path(CALDER2_BASE, tp_dir, sprintf("%s_bin_signals.tsv", tp))
  bin_signals <- fread(bin_signals_path)
  cat(sprintf("  Loaded bin_signals: %d rows x %d cols\n",
              nrow(bin_signals), ncol(bin_signals)))

  ctrl_mark_cols <- paste0(MARK_NAMES, "_ctrl")
  mut_mark_cols  <- paste0(MARK_NAMES, "_mut")

  # Join SNIPER ctrl labels to bin_signals
  sniper_ctrl_signals <- merge(
    bin_signals,
    sniper_ctrl[, .(chr, bin_start = calder_bin_start, sniper_ctrl_label = sniper_label)],
    by = c("chr", "bin_start")
  )
  cat(sprintf("  SNIPER ctrl x bin_signals join: %d rows\n", nrow(sniper_ctrl_signals)))

  enrichment_ctrl <- compute_enrichment(
    sniper_ctrl_signals, "sniper_ctrl_label", ctrl_mark_cols, MARK_NAMES)
  enrichment_ctrl[, condition := "ctrl"]
  cat(sprintf("  Ctrl enrichment: %d values\n", nrow(enrichment_ctrl)))

  # Join SNIPER mut labels to bin_signals
  sniper_mut_signals <- merge(
    bin_signals,
    sniper_mut[, .(chr, bin_start = calder_bin_start, sniper_mut_label = sniper_label)],
    by = c("chr", "bin_start")
  )
  cat(sprintf("  SNIPER mut x bin_signals join: %d rows\n", nrow(sniper_mut_signals)))

  enrichment_mut <- compute_enrichment(
    sniper_mut_signals, "sniper_mut_label", mut_mark_cols, MARK_NAMES)
  enrichment_mut[, condition := "mut"]
  cat(sprintf("  Mut enrichment:  %d values\n", nrow(enrichment_mut)))

  all_enrichment <- rbind(enrichment_ctrl, enrichment_mut, fill = TRUE)
  all_enrichment[, tp := tp]
  fwrite(all_enrichment,
         file.path(TSV_DIR, sprintf("%s_sniper_enrichment.tsv", tp)),
         sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s_sniper_enrichment.tsv\n", tp))

  # Enrichment figures
  cat("  Generating enrichment heatmaps...\n")

  p_enrich_ctrl <- plot_enrichment_heatmap(
    enrichment_ctrl,
    title      = sprintf("SNIPER Ctrl Enrichment: %s", tp_label),
    subtitle   = sprintf("SNIPER ctrl labels x Ctrl signal (%s bins)",
                         format(nrow(sniper_ctrl_signals), big.mark = ",")),
    mark_order = MARK_NAMES, subcomp_order = LABEL_ORDER)
  save_multiformat_ggplot(p_enrich_ctrl,
    file.path(PLOT_DIR, sprintf("c5_sniper_enrichment_ctrl_%s", tp_dir)),
    width = 7, height = 8)

  p_enrich_mut <- plot_enrichment_heatmap(
    enrichment_mut,
    title      = sprintf("SNIPER Mut Enrichment: %s", tp_label),
    subtitle   = sprintf("SNIPER mut labels x Mut signal (%s bins)",
                         format(nrow(sniper_mut_signals), big.mark = ",")),
    mark_order = MARK_NAMES, subcomp_order = LABEL_ORDER)
  save_multiformat_ggplot(p_enrich_mut,
    file.path(PLOT_DIR, sprintf("c5_sniper_enrichment_mut_%s", tp_dir)),
    width = 7, height = 8)

  # ══════════════════════════════════════════════════════════════════════════════
  # Part 3: Concordant Transitions (key deliverable)
  # ══════════════════════════════════════════════════════════════════════════════

  cat("\n=== Part 3: Concordant Transitions ===\n")

  conc_path <- file.path(SNIPER_BASE, "tsvs",
                         sprintf("%s_transition_concordance.tsv", tp))
  conc_dt <- fread(conc_path)
  cat(sprintf("  Loaded concordance: %d rows\n", nrow(conc_dt)))
  cat(sprintf("  Agreement breakdown:\n"))
  for (ag in c("Both_stable", "Both_change_agree", "Both_change_disagree",
               "SNIPER_only", "CALDER2_only")) {
    n <- sum(conc_dt$agreement == ag)
    cat(sprintf("    %-25s %5d (%.1f%%)\n", ag, n, 100 * n / nrow(conc_dt)))
  }

  # ── 3a. Concordant bins (Both_change_agree) ─────────────────────────────────

  concordant <- conc_dt[agreement == "Both_change_agree"]
  n_concordant <- nrow(concordant)
  cat(sprintf("\n  Concordant bins: %d\n", n_concordant))

  concordant_mat <- compute_transition_matrix(concordant, "calder_ctrl", "calder_mut")
  cat(sprintf("  Concordant matrix total: %d\n", sum(concordant_mat)))

  concordant_chisq <- chisq.test(concordant_mat, simulate.p.value = TRUE, B = 10000)
  cat(sprintf("  Chi-squared (MC): X2=%.1f, p=%.3e\n",
              concordant_chisq$statistic, concordant_chisq$p.value))

  # ── 3b. Discordant bins (Both_change_disagree) ──────────────────────────────

  discordant <- conc_dt[agreement == "Both_change_disagree"]
  n_discordant <- nrow(discordant)
  cat(sprintf("  Discordant bins: %d\n", n_discordant))

  discordant_mat_calder <- compute_transition_matrix(discordant, "calder_ctrl", "calder_mut")
  discordant_mat_sniper <- compute_transition_matrix(discordant, "sniper_ctrl", "sniper_mut")

  # ── 3c. Confirmation rates ──────────────────────────────────────────────────

  calder_labels_path <- file.path(CALDER2_BASE, tp_dir,
                                  sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  calder_labels <- fread(calder_labels_path)
  calder_changed <- calder_labels[label_changed == TRUE &
                                  !is.na(ctrl_label) & !is.na(mut_label)]
  calder_full_mat <- compute_transition_matrix(calder_labels, "ctrl_label", "mut_label")

  n_calder_changed <- nrow(calder_changed)
  overall_confirmation <- n_concordant / n_calder_changed
  cat(sprintf("  CALDER2 changed bins: %d\n", n_calder_changed))
  cat(sprintf("  Overall confirmation rate: %d/%d = %.1f%%\n",
              n_concordant, n_calder_changed, 100 * overall_confirmation))

  calder_off_diag <- compute_transition_matrix(calder_changed, "ctrl_label", "mut_label")

  rate_mat <- matrix(NA_real_, nrow = 4, ncol = 4,
                     dimnames = list(LABEL_ORDER, LABEL_ORDER))
  n_conc_mat <- matrix(0L, nrow = 4, ncol = 4,
                       dimnames = list(LABEL_ORDER, LABEL_ORDER))
  for (r in LABEL_ORDER) {
    for (cc in LABEL_ORDER) {
      if (r != cc) {
        denom <- calder_off_diag[r, cc]
        numer <- concordant_mat[r, cc]
        n_conc_mat[r, cc] <- numer
        rate_mat[r, cc] <- if (denom > 0) numer / denom else NA_real_
      }
    }
  }

  # ── Write Part 3 TSVs ──────────────────────────────────────────────────────

  fwrite(concordant, file.path(TSV_DIR, sprintf("%s_concordant_transitions.tsv", tp)),
         sep = "\t", quote = FALSE)

  write.table(concordant_mat,
              file.path(TSV_DIR, sprintf("%s_concordant_transition_matrix.tsv", tp)),
              sep = "\t", quote = FALSE, col.names = NA)

  fwrite(discordant, file.path(TSV_DIR, sprintf("%s_discordant_transitions.tsv", tp)),
         sep = "\t", quote = FALSE)

  rate_long <- as.data.table(as.table(rate_mat))
  names(rate_long) <- c("ctrl_label", "mut_label", "confirmation_rate")
  n_long <- as.data.table(as.table(n_conc_mat))
  names(n_long) <- c("ctrl_label", "mut_label", "n_concordant")
  calder_off_long <- as.data.table(as.table(calder_off_diag))
  names(calder_off_long) <- c("ctrl_label", "mut_label", "n_calder_changed")
  confirm_dt <- merge(merge(rate_long, n_long, by = c("ctrl_label", "mut_label")),
                      calder_off_long, by = c("ctrl_label", "mut_label"))
  confirm_dt <- confirm_dt[as.character(ctrl_label) != as.character(mut_label)]
  fwrite(confirm_dt[order(ctrl_label, mut_label)],
         file.path(TSV_DIR, sprintf("%s_confirmation_rates.tsv", tp)),
         sep = "\t", quote = FALSE, na = "NA")

  summary_dt <- data.table(
    timepoint       = tp,
    tp_label        = TP_LABELS[tp],
    n_sniper_bins   = n_sniper_total,
    n_sniper_changed = n_sniper_changed,
    pct_sniper_changed = round(pct_sniper, 2),
    n_calder_changed = n_calder_changed,
    n_concordant     = n_concordant,
    n_discordant     = n_discordant,
    confirmation_rate = round(overall_confirmation, 4),
    sniper_chisq_stat = round(unname(sniper_chisq$statistic), 1),
    sniper_cramer_v   = round(unname(sniper_cramer), 4),
    concordant_chisq_p = signif(concordant_chisq$p.value, 3)
  )
  fwrite(summary_dt,
         file.path(TSV_DIR, sprintf("%s_concordant_summary.tsv", tp)),
         sep = "\t", quote = FALSE)
  summary_rows[[tp]] <- summary_dt

  cat(sprintf("  Written: 5 concordant TSVs\n"))

  # ── Part 3 figures ──────────────────────────────────────────────────────────

  cat("\n  Generating concordant figures...\n")

  # 3.1 Concordant transition heatmap — the key deliverable
  pct_conc <- 100 * n_concordant / n_sniper_total
  p_conc_heat <- plot_transition_heatmap(
    concordant_mat,
    title    = sprintf("Concordant Transitions (SNIPER+CALDER2): %s", tp_label),
    subtitle = sprintf("%d bins where both tools agree (%.1f%% of SNIPER bins)",
                       n_concordant, pct_conc))
  save_multiformat_ggplot(p_conc_heat,
    file.path(PLOT_DIR, sprintf("c5_concordant_heatmap_%s", tp_dir)),
    width = 7, height = 6)

  # 3.2 Concordant vs discordant side-by-side
  p_disc_heat <- plot_transition_heatmap(
    discordant_mat_calder,
    title    = sprintf("Discordant Transitions (CALDER2 view): %s", tp_label),
    subtitle = sprintf("%d bins where tools disagree on transition", n_discordant))

  p_conc_vs_disc <- p_conc_heat + p_disc_heat +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = sprintf("Concordant vs Discordant Transitions: %s", tp_label),
      theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))
  save_multiformat_ggplot(p_conc_vs_disc,
    file.path(PLOT_DIR, sprintf("c5_concordant_vs_discordant_%s", tp_dir)),
    width = 14, height = 6)

  # 3.3 Three-way comparison: CALDER2 | SNIPER | Concordant
  p_calder_heat <- plot_transition_heatmap(
    calder_full_mat,
    title    = "CALDER2",
    subtitle = sprintf("%d total bins", sum(calder_full_mat)))
  p_sniper_3way <- plot_transition_heatmap(
    sniper_mat,
    title    = "SNIPER",
    subtitle = sprintf("%d total bins", sum(sniper_mat)))
  p_conc_3way <- plot_transition_heatmap(
    concordant_mat,
    title    = "Concordant",
    subtitle = sprintf("%d confirmed bins", sum(concordant_mat)))

  p_three_way <- p_calder_heat + p_sniper_3way + p_conc_3way +
    plot_layout(ncol = 3) +
    plot_annotation(
      title = sprintf("Transition Matrix Comparison: %s", tp_label),
      theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))
  save_multiformat_ggplot(p_three_way,
    file.path(PLOT_DIR, sprintf("c5_three_way_comparison_%s", tp_dir)),
    width = 21, height = 6)

  # 3.4 Confirmation rate heatmap
  p_confirm <- plot_confirmation_heatmap(rate_mat, n_conc_mat, tp_label)
  save_multiformat_ggplot(p_confirm,
    file.path(PLOT_DIR, sprintf("c5_confirmation_rates_%s", tp_dir)),
    width = 7, height = 6)

  # 3.5 Concordant transition sankey
  conc_alluv <- concordant[, .N, by = .(ctrl_label = calder_ctrl,
                                         mut_label = calder_mut)]
  p_conc_sankey <- plot_transition_sankey(
    conc_alluv,
    title    = sprintf("Concordant Flow (SNIPER+CALDER2): %s", tp_label),
    subtitle = sprintf("%d high-confidence transitions", n_concordant))
  save_multiformat_ggplot(p_conc_sankey,
    file.path(PLOT_DIR, sprintf("c5_concordant_sankey_%s", tp_dir)),
    width = 8, height = 7)

  # ── Verification ────────────────────────────────────────────────────────────

  cat("\n=== Verification ===\n")
  any_warn <- FALSE

  # Check 1: SNIPER matrix sum
  if (sum(sniper_mat) < 19000 || sum(sniper_mat) > 21000) {
    warning(sprintf("SNIPER matrix sum %d outside [19000, 21000]", sum(sniper_mat)))
    any_warn <- TRUE
  }

  # Check 2: Concordant matrix sum
  if (sum(concordant_mat) < 400 || sum(concordant_mat) > 1200) {
    warning(sprintf("Concordant matrix sum %d outside [400, 1200]", sum(concordant_mat)))
    any_warn <- TRUE
  }

  # Check 3: H3K27ac gradient on SNIPER ctrl labels
  k27ac_sniper <- enrichment_ctrl[mark == "H3K27ac"]
  k27ac_sniper <- k27ac_sniper[match(LABEL_ORDER, subcompartment)]
  k27ac_vals <- k27ac_sniper$fold_enrichment
  if (all(is.finite(k27ac_vals)) && all(diff(k27ac_vals) < 0)) {
    cat(sprintf("  H3K27ac SNIPER gradient OK: %s\n",
                paste(sprintf("%.2f", k27ac_vals), collapse = " > ")))
  } else {
    warning(sprintf("H3K27ac SNIPER ctrl NOT monotonically decreasing: %s",
                    paste(sprintf("%.2f", k27ac_vals), collapse = ", ")))
    any_warn <- TRUE
  }

  # Check 4: Concordant bins are a subset of CALDER2 changed bins
  conc_keys <- concordant[, paste(chr, bin_start)]
  calder_changed_keys <- calder_changed[, paste(chr, bin_start)]
  n_not_in_calder <- sum(!conc_keys %in% calder_changed_keys)
  if (n_not_in_calder > 0) {
    warning(sprintf("%d concordant bins NOT in CALDER2 changed set", n_not_in_calder))
    any_warn <- TRUE
  } else {
    cat(sprintf("  All %d concordant bins confirmed in CALDER2 changed set.\n",
                n_concordant))
  }

  # Check 5: Confirmation rate sanity
  if (overall_confirmation < 0.05 || overall_confirmation > 0.50) {
    warning(sprintf("Confirmation rate %.1f%% outside [5%%, 50%%]",
                    100 * overall_confirmation))
    any_warn <- TRUE
  } else {
    cat(sprintf("  Confirmation rate: %.1f%% (plausible)\n",
                100 * overall_confirmation))
  }

  # Check 6: TSV existence
  expected_tsvs <- c(
    sprintf("%s_sniper_differential.tsv", tp),
    sprintf("%s_sniper_transition_matrix.tsv", tp),
    sprintf("%s_sniper_transition_summary.tsv", tp),
    sprintf("%s_sniper_enrichment.tsv", tp),
    sprintf("%s_concordant_transitions.tsv", tp),
    sprintf("%s_concordant_transition_matrix.tsv", tp),
    sprintf("%s_discordant_transitions.tsv", tp),
    sprintf("%s_confirmation_rates.tsv", tp),
    sprintf("%s_concordant_summary.tsv", tp)
  )
  for (f in expected_tsvs) {
    fpath <- file.path(TSV_DIR, f)
    if (!file.exists(fpath) || file.info(fpath)$size == 0) {
      warning(sprintf("Missing or empty: %s", f))
      any_warn <- TRUE
    }
  }

  # Check 7: Figure directory counts
  expected_plots <- c(
    sprintf("c5_sniper_transition_heatmap_%s", tp_dir),
    sprintf("c5_sniper_transition_sankey_%s", tp_dir),
    sprintf("c5_sniper_enrichment_ctrl_%s", tp_dir),
    sprintf("c5_sniper_enrichment_mut_%s", tp_dir),
    sprintf("c5_concordant_heatmap_%s", tp_dir),
    sprintf("c5_concordant_vs_discordant_%s", tp_dir),
    sprintf("c5_three_way_comparison_%s", tp_dir),
    sprintf("c5_confirmation_rates_%s", tp_dir),
    sprintf("c5_concordant_sankey_%s", tp_dir)
  )
  for (d in expected_plots) {
    dpath <- file.path(PLOT_DIR, d)
    if (!dir.exists(dpath) || length(list.files(dpath)) < 4) {
      warning(sprintf("Missing or incomplete figure dir: %s", d))
      any_warn <- TRUE
    }
  }

  if (!any_warn) cat("  All verification checks passed.\n")

  cat(sprintf("\n  Timepoint %s complete.\n", tp))
}

# ── Cross-timepoint summary ───────────────────────────────────────────────────

combined_summary <- rbindlist(summary_rows)
fwrite(combined_summary,
       file.path(TSV_DIR, "c5_combined_summary.tsv"),
       sep = "\t", quote = FALSE)
cat(sprintf("\n  Written: c5_combined_summary.tsv\n"))

cat("\n===========================================\n")
cat("C5 COMPLETE\n")
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
