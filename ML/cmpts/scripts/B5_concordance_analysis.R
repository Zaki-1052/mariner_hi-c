# ML/cmpts/scripts/B5_concordance_analysis.R
# Stage B5: Concordance analysis between SNIPER and CALDER2 subcompartment calls.
#
# Quantifies training fidelity (ctrl — SNIPER trained on CALDER2 labels) and
# generalization (mut — unseen by SNIPER during training). Computes confusion
# matrices, Cohen's kappa, per-class precision/recall/F1, and differential
# transition concordance.
#
# Processes both timepoints in one invocation.
#
# Usage:
#   Rscript B5_concordance_analysis.R <data_root> <code_root>
#     <data_root> : HPC data directory (e.g. /expanse/.../sniper)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript B5_concordance_analysis.R <data_root> <code_root>")
}

DATA_ROOT <- args[1]
CODE_ROOT <- args[2]

# ── Constants ──────────────────────────────────────────────────────────────────

TPS <- c("250402", "250831")
TP_LABELS <- c("250402" = "late/adult", "250831" = "early/P12")
TP_DIRS   <- c("250402" = "late", "250831" = "early")

SAMPLES     <- c("ctrl_merged", "mut_merged")
LABEL_ORDER <- c("A.1", "A.2", "B.1", "B.2")
LABEL_COLORS <- c(
  "A.1" = "#e41a1c",
  "A.2" = "#ff7f00",
  "B.1" = "#4daf4a",
  "B.2" = "#377eb8"
)

SNIPER_TO_CALDER <- c(
  "A1" = "A.1", "A2" = "A.2", "B1" = "B.1", "B2" = "B.2"
)

CONDITION_LABELS <- c("ctrl_merged" = "ctrl_label", "mut_merged" = "mut_label")

OUT_BASE  <- file.path(CODE_ROOT, "outputs", "sniper")
UTIL_PATH <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

AGREEMENT_COLORS <- c(
  "Both_stable"          = "#4daf4a",
  "Both_change_agree"    = "#377eb8",
  "Both_change_disagree" = "#e41a1c",
  "SNIPER_only"          = "#ff7f00",
  "CALDER2_only"         = "#984ea3"
)

# ── Libraries ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(scales)
})
source(UTIL_PATH)

# ── Header ─────────────────────────────────────────────────────────────────────

cat("===========================================\n")
cat("B5: SNIPER-CALDER2 Concordance Analysis\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("Output dir: %s\n", OUT_BASE))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

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

load_calder2_labels <- function(path, label_col) {
  dt <- fread(path)
  if (!label_col %in% names(dt)) {
    stop(sprintf("Column '%s' not found in %s", label_col, basename(path)))
  }
  setnames(dt, label_col, "calder_label")
  dt <- dt[!is.na(calder_label) & calder_label != "NA"]
  dt[, .(chr, bin_start, calder_label)]
}

join_predictions <- function(sniper_dt, calder_dt) {
  merged <- merge(sniper_dt, calder_dt,
                  by.x = c("chr", "calder_bin_start"),
                  by.y = c("chr", "bin_start"),
                  all.x = TRUE)
  n_total  <- nrow(merged)
  n_na     <- sum(is.na(merged$calder_label))
  merged   <- merged[!is.na(calder_label)]
  n_joined <- nrow(merged)
  cat(sprintf("    SNIPER bins: %d, matched: %d, unmatched/NA: %d\n",
              n_total, n_joined, n_na))
  merged[, concordant := (sniper_label == calder_label)]
  merged
}

compute_confusion_matrix <- function(joined_dt) {
  counts <- joined_dt[, .N, by = .(sniper_label, calder_label)]
  mat <- matrix(0L, nrow = 4, ncol = 4,
                dimnames = list(SNIPER = LABEL_ORDER, CALDER2 = LABEL_ORDER))
  for (i in seq_len(nrow(counts))) {
    r  <- counts$sniper_label[i]
    cc <- counts$calder_label[i]
    mat[r, cc] <- counts$N[i]
  }
  mat
}

compute_cohens_kappa <- function(conf_mat) {
  n  <- sum(conf_mat)
  p_o <- sum(diag(conf_mat)) / n
  row_m <- rowSums(conf_mat)
  col_m <- colSums(conf_mat)
  p_e   <- sum(row_m * col_m) / n^2
  kappa <- (p_o - p_e) / (1 - p_e)
  list(kappa = kappa, p_observed = p_o, p_expected = p_e, n = n)
}

compute_per_class_metrics <- function(conf_mat) {
  classes <- rownames(conf_mat)
  metrics <- lapply(classes, function(cls) {
    tp <- conf_mat[cls, cls]
    row_sum <- sum(conf_mat[cls, ])
    col_sum <- sum(conf_mat[, cls])
    precision <- if (row_sum > 0) tp / row_sum else NA_real_
    recall    <- if (col_sum > 0) tp / col_sum else NA_real_
    f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA_real_
    }
    data.table(class = cls, precision = precision, recall = recall, f1 = f1,
               support_sniper = row_sum, support_calder2 = col_sum)
  })
  rbindlist(metrics)
}

compute_transition_concordance <- function(sniper_ctrl_dt, sniper_mut_dt, calder_labels_path) {
  calder_full <- fread(calder_labels_path)
  calder_both <- calder_full[!is.na(ctrl_label) & ctrl_label != "NA" &
                             !is.na(mut_label)  & mut_label  != "NA",
                             .(chr, bin_start, calder_ctrl = ctrl_label,
                               calder_mut = mut_label)]

  sniper_both <- merge(
    sniper_ctrl_dt[, .(chr, calder_bin_start, sniper_ctrl = sniper_label)],
    sniper_mut_dt[,  .(chr, calder_bin_start, sniper_mut  = sniper_label)],
    by = c("chr", "calder_bin_start")
  )

  quad <- merge(sniper_both, calder_both,
                by.x = c("chr", "calder_bin_start"),
                by.y = c("chr", "bin_start"))
  cat(sprintf("    Bins with all 4 labels: %d\n", nrow(quad)))

  quad[, sniper_changed := (sniper_ctrl != sniper_mut)]
  quad[, calder_changed := (calder_ctrl != calder_mut)]
  quad[, same_transition := (sniper_ctrl == calder_ctrl & sniper_mut == calder_mut)]

  quad[, agreement := fifelse(
    !sniper_changed & !calder_changed, "Both_stable",
    fifelse(sniper_changed & calder_changed & same_transition, "Both_change_agree",
    fifelse(sniper_changed & calder_changed & !same_transition, "Both_change_disagree",
    fifelse(sniper_changed & !calder_changed, "SNIPER_only",
            "CALDER2_only"))))]

  quad
}

plot_confusion_heatmap <- function(conf_mat, title, subtitle) {
  mat_long <- as.data.table(as.table(conf_mat))
  names(mat_long) <- c("sniper_label", "calder_label", "count")
  mat_long[, log10_count := log10(count + 1)]
  mat_long[, on_diagonal := (as.character(sniper_label) == as.character(calder_label))]
  mat_long[, sniper_label := factor(sniper_label, levels = rev(LABEL_ORDER))]
  mat_long[, calder_label := factor(calder_label, levels = LABEL_ORDER)]

  ggplot(mat_long, aes(x = calder_label, y = sniper_label, fill = log10_count)) +
    geom_tile(aes(color = on_diagonal), linewidth = 1.2) +
    geom_text(aes(label = format(count, big.mark = ",")), fontface = "bold", size = 5) +
    scale_fill_gradient(low = "white", high = "#2166ac", name = expression(log[10](n+1))) +
    scale_color_manual(values = c("TRUE" = "#d7301f", "FALSE" = "grey85"), guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = "CALDER2 label", y = "SNIPER label",
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          axis.text  = element_text(face = "bold", size = 12))
}

plot_discordant_alluvial <- function(joined_ctrl, joined_mut, tp_label) {
  make_alluv_dt <- function(joined_dt, cond_label) {
    disc <- joined_dt[concordant == FALSE]
    if (nrow(disc) == 0) return(data.table())
    counts <- disc[, .N, by = .(sniper_label, calder_label)]
    counts[, condition := cond_label]
    counts
  }

  alluv_dt <- rbind(
    make_alluv_dt(joined_ctrl, "Ctrl"),
    make_alluv_dt(joined_mut,  "Mut")
  )
  if (nrow(alluv_dt) == 0) return(NULL)

  n_disc_ctrl <- sum(joined_ctrl$concordant == FALSE)
  n_disc_mut  <- sum(joined_mut$concordant == FALSE)
  pct_ctrl <- 100 * n_disc_ctrl / nrow(joined_ctrl)
  pct_mut  <- 100 * n_disc_mut  / nrow(joined_mut)

  alluv_dt[, sniper_label := factor(sniper_label, levels = LABEL_ORDER)]
  alluv_dt[, calder_label := factor(calder_label, levels = LABEL_ORDER)]
  alluv_dt[, condition := factor(condition, levels = c("Ctrl", "Mut"))]

  ggplot(alluv_dt, aes(axis1 = sniper_label, axis2 = calder_label, y = N)) +
    geom_alluvium(aes(fill = sniper_label), width = 1/4, alpha = 0.75, knot.pos = 0.4) +
    geom_stratum(width = 1/4, fill = "white", color = "grey40", linewidth = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = LABEL_COLORS, name = "SNIPER label") +
    scale_x_discrete(limits = c("SNIPER", "CALDER2"), expand = c(0.05, 0.05)) +
    facet_wrap(~condition, scales = "free_y") +
    labs(y = "Number of discordant 100kb bins",
         title = sprintf("Discordant Bins: SNIPER vs CALDER2 (%s)", tp_label),
         subtitle = sprintf("Ctrl: %d bins (%.1f%%)  |  Mut: %d bins (%.1f%%)",
                            n_disc_ctrl, pct_ctrl, n_disc_mut, pct_mut)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12))
}

plot_per_class_concordance <- function(metrics_ctrl, metrics_mut, tp_label) {
  metrics_ctrl[, condition := "Ctrl (training)"]
  metrics_mut[,  condition := "Mut (unseen)"]
  combined <- rbind(metrics_ctrl, metrics_mut)

  long <- melt(combined, id.vars = c("class", "condition", "support_sniper", "support_calder2"),
               measure.vars = c("precision", "recall", "f1"),
               variable.name = "metric", value.name = "value")
  long[, metric := factor(metric, levels = c("precision", "recall", "f1"),
                          labels = c("Precision", "Recall", "F1"))]
  long[, class := factor(class, levels = LABEL_ORDER)]
  long[, condition := factor(condition, levels = c("Ctrl (training)", "Mut (unseen)"))]

  ggplot(long, aes(x = class, y = value, fill = condition)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.2f", value)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
    facet_wrap(~metric) +
    scale_fill_manual(values = c("Ctrl (training)" = "#377eb8", "Mut (unseen)" = "#e41a1c"),
                      name = "Condition") +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
    labs(x = "Subcompartment", y = "Score",
         title = sprintf("Per-Class Concordance: %s", tp_label),
         subtitle = "SNIPER vs CALDER2 (SNIPER = predicted, CALDER2 = reference)") +
    theme_minimal(base_size = 13) +
    theme(strip.text = element_text(face = "bold", size = 12))
}

plot_transition_agreement <- function(trans_dt, tp_label) {
  summary_dt <- trans_dt[, .N, by = agreement]
  summary_dt[, pct := 100 * N / sum(N)]
  summary_dt[, agreement := factor(agreement,
    levels = c("Both_stable", "Both_change_agree", "Both_change_disagree",
               "SNIPER_only", "CALDER2_only"))]
  summary_dt <- summary_dt[order(agreement)]

  ggplot(summary_dt, aes(x = agreement, y = N, fill = agreement)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(N, big.mark = ","), pct)),
              vjust = -0.3, size = 3.5, lineheight = 0.85) +
    scale_fill_manual(values = AGREEMENT_COLORS) +
    scale_x_discrete(labels = c(
      "Both_stable"          = "Both\nstable",
      "Both_change_agree"    = "Both change\n(agree)",
      "Both_change_disagree" = "Both change\n(disagree)",
      "SNIPER_only"          = "SNIPER\nonly",
      "CALDER2_only"         = "CALDER2\nonly"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Transition agreement category", y = "Number of 100kb bins",
         title = sprintf("Ctrl-to-Mut Transition Agreement: %s", tp_label),
         subtitle = "Do SNIPER and CALDER2 agree on which bins change subcompartment?") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank())
}

# ── Main loop ──────────────────────────────────────────────────────────────────

combined_summary <- list()

for (tp in TPS) {

  tp_label <- sprintf("%s (%s)", tp, TP_LABELS[tp])
  tp_dir   <- TP_DIRS[tp]
  OUT_DIR  <- OUT_BASE

  cat(sprintf("\n###############################################\n"))
  cat(sprintf("### Timepoint: %s\n", tp_label))
  cat(sprintf("###############################################\n\n"))

  # ── Pre-flight validation ──────────────────────────────────────────────────

  cat("=== Pre-flight validation ===\n")

  sniper_ctrl_path <- file.path(CODE_ROOT, "outputs", "sniper", "predictions",
                                tp, "ctrl_merged", "predictions.bed")
  sniper_mut_path  <- file.path(CODE_ROOT, "outputs", "sniper", "predictions",
                                tp, "mut_merged", "predictions.bed")
  calder_labels_path <- file.path(CODE_ROOT, "outputs", "calder2", tp_dir,
                                  sprintf("%s_subcompartment_labels_100kb.tsv", tp))

  input_paths <- c(sniper_ctrl = sniper_ctrl_path,
                   sniper_mut  = sniper_mut_path,
                   calder2     = calder_labels_path)
  for (nm in names(input_paths)) {
    p <- input_paths[nm]
    if (!file.exists(p)) stop(sprintf("Missing input: %s (%s)", nm, p))
    if (file.info(p)$size == 0) stop(sprintf("Empty input: %s (%s)", nm, p))
    cat(sprintf("  OK: %s (%s bytes)\n", nm, format(file.info(p)$size, big.mark = ",")))
  }

  # ── Load data ──────────────────────────────────────────────────────────────

  cat("\n=== Loading data ===\n")

  cat("  Loading SNIPER ctrl predictions...\n")
  sniper_ctrl <- load_sniper_bed(sniper_ctrl_path)
  cat(sprintf("    %d bins\n", nrow(sniper_ctrl)))

  cat("  Loading SNIPER mut predictions...\n")
  sniper_mut <- load_sniper_bed(sniper_mut_path)
  cat(sprintf("    %d bins\n", nrow(sniper_mut)))

  cat("  Loading CALDER2 ctrl labels...\n")
  calder_ctrl <- load_calder2_labels(calder_labels_path, "ctrl_label")
  cat(sprintf("    %d callable bins\n", nrow(calder_ctrl)))

  cat("  Loading CALDER2 mut labels...\n")
  calder_mut <- load_calder2_labels(calder_labels_path, "mut_label")
  cat(sprintf("    %d callable bins\n", nrow(calder_mut)))

  # ── Per-condition concordance ──────────────────────────────────────────────

  tp_results <- list()

  for (cond in c("ctrl", "mut")) {
    cat(sprintf("\n=== Concordance: %s ===\n", toupper(cond)))

    sniper_dt <- if (cond == "ctrl") sniper_ctrl else sniper_mut
    calder_dt <- if (cond == "ctrl") calder_ctrl else calder_mut

    cat("  Joining predictions...\n")
    joined <- join_predictions(sniper_dt, calder_dt)

    conf_mat <- compute_confusion_matrix(joined)
    kappa    <- compute_cohens_kappa(conf_mat)
    metrics  <- compute_per_class_metrics(conf_mat)
    accuracy <- sum(diag(conf_mat)) / sum(conf_mat)

    cat(sprintf("  Accuracy: %.1f%% (%d/%d)\n",
                100 * accuracy, sum(diag(conf_mat)), sum(conf_mat)))
    cat(sprintf("  Cohen's kappa: %.4f (p_o=%.4f, p_e=%.4f)\n",
                kappa$kappa, kappa$p_observed, kappa$p_expected))
    cat("  Per-class metrics:\n")
    for (i in seq_len(nrow(metrics))) {
      cat(sprintf("    %-4s  P=%.3f  R=%.3f  F1=%.3f  (SNIPER=%d, CALDER2=%d)\n",
                  metrics$class[i], metrics$precision[i], metrics$recall[i],
                  metrics$f1[i], metrics$support_sniper[i], metrics$support_calder2[i]))
    }

    tp_results[[cond]] <- list(
      joined   = joined,
      conf_mat = conf_mat,
      kappa    = kappa,
      metrics  = metrics,
      accuracy = accuracy
    )

    # Write per-bin concordance TSV
    out_conc <- file.path(OUT_DIR, sprintf("%s_concordance_%s.tsv", tp, cond))
    fwrite(joined[, .(chr, bin_start = calder_bin_start, sniper_label, calder_label, concordant)],
           out_conc, sep = "\t", quote = FALSE)
    cat(sprintf("  Written: %s\n", basename(out_conc)))

    # Write confusion matrix TSV
    out_cm <- file.path(OUT_DIR, sprintf("%s_confusion_matrix_%s.tsv", tp, cond))
    write.table(conf_mat, out_cm, sep = "\t", quote = FALSE, col.names = NA)
    cat(sprintf("  Written: %s\n", basename(out_cm)))
  }

  # Write per-class metrics (both conditions stacked)
  metrics_ctrl <- copy(tp_results$ctrl$metrics)
  metrics_mut  <- copy(tp_results$mut$metrics)
  metrics_ctrl[, condition := "ctrl"]
  metrics_mut[,  condition := "mut"]
  all_metrics <- rbind(metrics_ctrl, metrics_mut)
  out_metrics <- file.path(OUT_DIR, sprintf("%s_per_class_metrics.tsv", tp))
  fwrite(all_metrics, out_metrics, sep = "\t", quote = FALSE)
  cat(sprintf("\n  Written: %s\n", basename(out_metrics)))

  # ── Differential transition concordance ────────────────────────────────────

  cat("\n=== Differential transition concordance ===\n")
  trans_dt <- compute_transition_concordance(sniper_ctrl, sniper_mut,
                                             calder_labels_path)

  trans_summary <- trans_dt[, .N, by = agreement]
  trans_summary[, pct := round(100 * N / sum(N), 2)]
  trans_summary <- trans_summary[order(match(agreement,
    c("Both_stable", "Both_change_agree", "Both_change_disagree",
      "SNIPER_only", "CALDER2_only")))]

  cat("  Agreement categories:\n")
  for (i in seq_len(nrow(trans_summary))) {
    cat(sprintf("    %-25s %6d (%5.1f%%)\n",
                trans_summary$agreement[i], trans_summary$N[i], trans_summary$pct[i]))
  }

  out_trans <- file.path(OUT_DIR, sprintf("%s_transition_concordance.tsv", tp))
  fwrite(trans_dt[, .(chr, bin_start = calder_bin_start,
                       sniper_ctrl, sniper_mut, calder_ctrl, calder_mut,
                       sniper_changed, calder_changed, agreement)],
         out_trans, sep = "\t", quote = FALSE)
  cat(sprintf("  Written: %s\n", basename(out_trans)))

  out_trans_summ <- file.path(OUT_DIR, sprintf("%s_transition_agreement_summary.tsv", tp))
  fwrite(trans_summary, out_trans_summ, sep = "\t", quote = FALSE)
  cat(sprintf("  Written: %s\n", basename(out_trans_summ)))

  # ── Figures ──────────────────────────────────────────────────────────────

  cat("\n=== Generating figures ===\n")

  # 1. Confusion heatmap — ctrl
  p_cm_ctrl <- plot_confusion_heatmap(
    tp_results$ctrl$conf_mat,
    sprintf("SNIPER vs CALDER2 Concordance (Ctrl): %s", tp_label),
    sprintf("Training data | Accuracy=%.1f%%, kappa=%.3f",
            100 * tp_results$ctrl$accuracy, tp_results$ctrl$kappa$kappa))
  save_multiformat_ggplot(p_cm_ctrl,
    file.path(OUT_DIR, sprintf("%s_confusion_heatmap_ctrl", tp)),
    width = 7, height = 6)

  # 2. Confusion heatmap — mut
  p_cm_mut <- plot_confusion_heatmap(
    tp_results$mut$conf_mat,
    sprintf("SNIPER vs CALDER2 Concordance (Mut): %s", tp_label),
    sprintf("Unseen data | Accuracy=%.1f%%, kappa=%.3f",
            100 * tp_results$mut$accuracy, tp_results$mut$kappa$kappa))
  save_multiformat_ggplot(p_cm_mut,
    file.path(OUT_DIR, sprintf("%s_confusion_heatmap_mut", tp)),
    width = 7, height = 6)

  # 3. Discordant alluvial
  p_disc <- plot_discordant_alluvial(
    tp_results$ctrl$joined, tp_results$mut$joined, tp_label)
  if (!is.null(p_disc)) {
    save_multiformat_ggplot(p_disc,
      file.path(OUT_DIR, sprintf("%s_discordant_alluvial", tp)),
      width = 10, height = 7)
  }

  # 4. Per-class concordance bar chart
  metrics_ctrl_plot <- copy(tp_results$ctrl$metrics)
  metrics_mut_plot  <- copy(tp_results$mut$metrics)
  p_class <- plot_per_class_concordance(metrics_ctrl_plot, metrics_mut_plot, tp_label)
  save_multiformat_ggplot(p_class,
    file.path(OUT_DIR, sprintf("%s_per_class_concordance", tp)),
    width = 9, height = 6)

  # 5. Transition agreement bar chart
  p_trans <- plot_transition_agreement(trans_dt, tp_label)
  save_multiformat_ggplot(p_trans,
    file.path(OUT_DIR, sprintf("%s_transition_agreement", tp)),
    width = 8, height = 6)

  # ── Build summary row for combined output ────────────────────────────────

  for (cond in c("ctrl", "mut")) {
    r <- tp_results[[cond]]
    m <- r$metrics
    combined_summary[[length(combined_summary) + 1]] <- data.table(
      timepoint    = tp,
      tp_label     = TP_LABELS[tp],
      condition    = cond,
      n_bins       = r$kappa$n,
      accuracy     = round(r$accuracy, 4),
      kappa        = round(r$kappa$kappa, 4),
      f1_A1        = round(m[class == "A.1", f1], 4),
      f1_A2        = round(m[class == "A.2", f1], 4),
      f1_B1        = round(m[class == "B.1", f1], 4),
      f1_B2        = round(m[class == "B.2", f1], 4)
    )
  }

  # ── Verification ─────────────────────────────────────────────────────────

  cat("\n=== Verification ===\n")
  any_warn <- FALSE

  # Check 1: join coverage
  for (cond in c("ctrl", "mut")) {
    sniper_n <- if (cond == "ctrl") nrow(sniper_ctrl) else nrow(sniper_mut)
    joined_n <- nrow(tp_results[[cond]]$joined)
    coverage <- joined_n / sniper_n
    cat(sprintf("  %s join coverage: %.1f%% (%d/%d)\n",
                cond, 100 * coverage, joined_n, sniper_n))
    if (coverage < 0.95) {
      warning(sprintf("Low join coverage for %s %s: %.1f%%", tp, cond, 100 * coverage))
      any_warn <- TRUE
    }
  }

  # Check 2: ctrl accuracy
  if (tp_results$ctrl$accuracy < 0.70) {
    warning(sprintf("Low ctrl accuracy for %s: %.1f%%", tp, 100 * tp_results$ctrl$accuracy))
    any_warn <- TRUE
  }

  # Check 3: ctrl kappa
  if (tp_results$ctrl$kappa$kappa < 0.50) {
    warning(sprintf("Low ctrl kappa for %s: %.3f", tp, tp_results$ctrl$kappa$kappa))
    any_warn <- TRUE
  }

  # Check 4: no empty classes
  for (cond in c("ctrl", "mut")) {
    cm <- tp_results[[cond]]$conf_mat
    empty_rows <- LABEL_ORDER[rowSums(cm) == 0]
    empty_cols <- LABEL_ORDER[colSums(cm) == 0]
    if (length(empty_rows) > 0) {
      warning(sprintf("Empty SNIPER classes in %s %s: %s",
                      tp, cond, paste(empty_rows, collapse = ", ")))
      any_warn <- TRUE
    }
    if (length(empty_cols) > 0) {
      warning(sprintf("Empty CALDER2 classes in %s %s: %s",
                      tp, cond, paste(empty_cols, collapse = ", ")))
      any_warn <- TRUE
    }
  }

  # Check 5: mut kappa
  if (tp_results$mut$kappa$kappa < 0.30) {
    warning(sprintf("Low mut kappa for %s: %.3f (poor generalization)",
                    tp, tp_results$mut$kappa$kappa))
    any_warn <- TRUE
  }

  # Check 6: Both_stable dominance
  stable_pct <- trans_summary[agreement == "Both_stable", pct]
  if (length(stable_pct) > 0 && stable_pct < 80) {
    warning(sprintf("Both_stable < 80%% for %s: %.1f%%", tp, stable_pct))
    any_warn <- TRUE
  }

  if (!any_warn) cat("  All checks passed.\n")
}

# ── Combined summary ──────────────────────────────────────────────────────────

cat("\n===========================================\n")
cat("Combined Summary\n")
cat("===========================================\n")

summary_dt <- rbindlist(combined_summary)
out_combined <- file.path(OUT_BASE, "combined_concordance_summary.tsv")
fwrite(summary_dt, out_combined, sep = "\t", quote = FALSE)
cat(sprintf("Written: %s\n\n", basename(out_combined)))

print(summary_dt)

cat("\n===========================================\n")
cat(sprintf("B5 COMPLETE\n"))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
