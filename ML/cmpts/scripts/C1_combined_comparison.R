# ML/cmpts/scripts/C1_combined_comparison.R
# Stage C1: Combined subcompartment comparison across conditions and timepoints.
#
# Produces genome-wide karyotype tracks, combined transition heatmaps,
# developmental comparison of transition rates, and cross-timepoint
# subcompartment stability analysis. Processes both timepoints in one
# invocation.
#
# Usage:
#   Rscript C1_combined_comparison.R <data_root> <code_root>
#     <data_root> : HPC data directory (kept for CLI consistency)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript C1_combined_comparison.R <data_root> <code_root>")
}

DATA_ROOT <- args[1]
CODE_ROOT <- args[2]

# -- Constants ----------------------------------------------------------------

TPS <- c("250402", "250831")
TP_LABELS <- c("250402" = "late/adult", "250831" = "early/P12")
TP_DIRS   <- c("250402" = "late", "250831" = "early")

BIN_SIZE    <- 100000L
LABEL_ORDER <- c("A.1", "A.2", "B.1", "B.2")
LABEL_COLORS <- c(
  "A.1" = "#e41a1c",
  "A.2" = "#ff7f00",
  "B.1" = "#4daf4a",
  "B.2" = "#377eb8"
)

MM10_CHROM_SIZES <- c(
  chr1  = 195471971L, chr2  = 182113224L, chr3  = 160039680L, chr4  = 156508116L,
  chr5  = 151834684L, chr6  = 149736546L, chr7  = 145441459L, chr8  = 129401213L,
  chr9  = 124595110L, chr10 = 130694993L, chr11 = 122082543L, chr12 = 120129022L,
  chr13 = 120421639L, chr14 = 124902244L, chr15 = 104043685L, chr16 = 98207768L,
  chr17 = 94987271L,  chr18 = 90702639L,  chr19 = 61431566L
)

CHROM_ORDER <- paste0("chr", 1:19)

OUT_DIR   <- file.path(CODE_ROOT, "outputs", "calder2", "combined")
UTIL_PATH <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

TRANSITION_CATEGORIES <- c(
  "A_to_B"   = "#d01c8b",
  "B_to_A"   = "#4dac26",
  "Within_A" = "#ff7f00",
  "Within_B" = "#377eb8"
)

CONDITION_COLORS <- c(
  "Ctrl developmental"  = "#1b9e77",
  "Mut developmental"   = "#d95f02",
  "Early genotype"      = "#7570b3",
  "Late genotype"       = "#e7298a"
)

# -- Libraries ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
source(UTIL_PATH)

# -- Header -------------------------------------------------------------------

start_time <- proc.time()

cat("===========================================\n")
cat("C1: Combined Subcompartment Comparison\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("Output dir: %s\n", OUT_DIR))
cat(sprintf("Timepoints: %s\n", paste(TPS, collapse = ", ")))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# -- Pre-flight validation ----------------------------------------------------

cat("=== Pre-flight validation ===\n")
input_paths <- list()
for (tp in TPS) {
  path <- file.path(CODE_ROOT, "outputs", "calder2", TP_DIRS[tp],
                    sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  if (!file.exists(path)) stop(sprintf("Missing A2 output: %s", path))
  if (file.info(path)$size == 0) stop(sprintf("Empty A2 output: %s", path))
  input_paths[[tp]] <- path
  cat(sprintf("  OK: %s (%s bytes)\n", basename(path),
              format(file.info(path)$size, big.mark = ",")))
}

if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -- Function definitions -----------------------------------------------------

compute_transition_matrix <- function(dt, from_col, to_col) {
  callable <- dt[!is.na(get(from_col)) & !is.na(get(to_col))]
  callable[, from := factor(get(from_col), levels = LABEL_ORDER)]
  callable[, to   := factor(get(to_col),   levels = LABEL_ORDER)]
  mat <- table(callable$from, callable$to)
  dimnames(mat) <- list(from = LABEL_ORDER, to = LABEL_ORDER)
  as.matrix(mat)
}

classify_transition <- function(from_label, to_label) {
  from_A <- grepl("^A", from_label)
  to_A   <- grepl("^A", to_label)
  same   <- (from_label == to_label)
  fifelse(same, "Stable",
    fifelse(from_A & !to_A, "A_to_B",
      fifelse(!from_A & to_A, "B_to_A",
        fifelse(from_A & to_A, "Within_A", "Within_B"))))
}

matrix_to_long <- function(mat, tp = NA_character_, tp_label = NA_character_) {
  dt <- as.data.table(as.table(mat))
  names(dt) <- c("from_label", "to_label", "count")
  total <- sum(dt$count)
  dt[, `:=`(
    rate         = count / total,
    pct_of_total = round(100 * count / total, 2),
    is_diagonal  = (as.character(from_label) == as.character(to_label)),
    category     = classify_transition(as.character(from_label), as.character(to_label))
  )]
  if (!is.na(tp)) dt[, tp := tp]
  if (!is.na(tp_label)) dt[, tp_label := tp_label]
  dt
}

build_karyotype_data <- function(labels_list) {
  karyotype_dt <- rbindlist(lapply(TPS, function(tp) {
    dt <- labels_list[[tp]]
    rbind(
      dt[, .(chr, bin_start, bin_end, label = ctrl_label,
             condition = "Ctrl", tp = tp, tp_label = TP_LABELS[tp])],
      dt[, .(chr, bin_start, bin_end, label = mut_label,
             condition = "Mut",  tp = tp, tp_label = TP_LABELS[tp])]
    )
  }))

  karyotype_dt[, chr := factor(chr, levels = rev(CHROM_ORDER))]
  karyotype_dt[, chr_idx := as.integer(chr)]
  karyotype_dt[, y_center := chr_idx + fifelse(condition == "Ctrl", 0.15, -0.15)]
  karyotype_dt[, ymin := y_center - 0.12]
  karyotype_dt[, ymax := y_center + 0.12]
  karyotype_dt[, xmin := (bin_start - 1) / 1e6]
  karyotype_dt[, xmax := bin_end / 1e6]
  karyotype_dt[, tp_label := factor(tp_label, levels = TP_LABELS)]

  karyotype_dt
}

plot_karyotype_track <- function(karyotype_dt) {
  max_mb <- max(MM10_CHROM_SIZES) / 1e6

  ggplot(karyotype_dt) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = label), color = NA) +
    scale_fill_manual(
      values = LABEL_COLORS, na.value = "#e0e0e0",
      name = "Subcompartment",
      breaks = LABEL_ORDER
    ) +
    facet_wrap(~ tp_label, ncol = 1, strip.position = "top") +
    scale_y_continuous(
      breaks = seq_along(CHROM_ORDER),
      labels = rev(CHROM_ORDER),
      expand = expansion(add = 0.4)
    ) +
    scale_x_continuous(
      name = "Position (Mb)",
      limits = c(0, max_mb),
      expand = expansion(mult = 0.01)
    ) +
    labs(y = NULL,
         title = "Genome-wide Subcompartment Map",
         subtitle = "Ctrl (upper) vs Mut (lower) per chromosome") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.text          = element_text(face = "bold", size = 12),
      axis.text.y         = element_text(size = 9),
      legend.position     = "bottom"
    )
}

plot_combined_heatmap <- function(mat_long) {
  mat_long[, from_label := factor(from_label, levels = rev(LABEL_ORDER))]
  mat_long[, to_label   := factor(to_label,   levels = LABEL_ORDER)]
  mat_long[, tp_label   := factor(tp_label, levels = TP_LABELS)]

  ggplot(mat_long, aes(x = to_label, y = from_label, fill = log10(count + 1))) +
    geom_tile(aes(color = is_diagonal), linewidth = 1.2) +
    geom_text(aes(label = format(count, big.mark = ",")),
              fontface = "bold", size = 4.5) +
    scale_fill_gradient(low = "white", high = "#2166ac", name = "log10(count+1)") +
    scale_color_manual(values = c("TRUE" = "#d7301f", "FALSE" = "grey85"),
                       guide = "none") +
    scale_x_discrete(position = "top") +
    facet_wrap(~ tp_label, nrow = 1) +
    labs(x = "Mutant label", y = "Control label",
         title = "Subcompartment Transitions: Ctrl to Mut") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_text(face = "bold", size = 12),
      strip.text   = element_text(face = "bold", size = 13)
    )
}

plot_developmental_comparison <- function(rates_combined) {
  off_diag <- rates_combined[is_diagonal == FALSE]
  off_diag[, transition := sprintf("%s -> %s", from_label, to_label)]
  off_diag[, tp_label := factor(tp_label, levels = TP_LABELS)]
  off_diag[, category := factor(category,
    levels = c("A_to_B", "B_to_A", "Within_A", "Within_B"))]

  ggplot(off_diag, aes(x = transition, y = pct_of_total, fill = tp_label)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    facet_wrap(~ category, scales = "free_x", nrow = 1) +
    scale_fill_manual(
      values = c("late/adult" = "#e7298a", "early/P12" = "#7570b3"),
      name = "Timepoint"
    ) +
    labs(x = NULL, y = "% of callable bins",
         title = "Transition Rates: Early vs Late",
         subtitle = "Per-transition comparison between timepoints") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
      strip.text   = element_text(face = "bold"),
      legend.position = "top"
    )
}

plot_stability_summary <- function(stability_dt) {
  stability_dt[, comparison := factor(comparison, levels = comparison)]

  ggplot(stability_dt, aes(x = comparison, y = pct_changed, fill = comparison)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%", pct_changed)),
              vjust = -0.5, fontface = "bold", size = 4.5) +
    scale_fill_manual(values = unname(CONDITION_COLORS)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, NA)) +
    labs(x = NULL, y = "% bins changed",
         title = "Subcompartment Dynamics",
         subtitle = "Genotype effect (ctrl vs mut) vs developmental change (early vs late)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 11),
      panel.grid.major.x = element_blank()
    )
}

plot_developmental_heatmap <- function(dev_long) {
  dev_long[, from_label := factor(from_label, levels = rev(LABEL_ORDER))]
  dev_long[, to_label   := factor(to_label,   levels = LABEL_ORDER)]
  dev_long[, condition   := factor(condition, levels = c("Ctrl", "Mut"))]

  ggplot(dev_long, aes(x = to_label, y = from_label, fill = log10(count + 1))) +
    geom_tile(aes(color = is_diagonal), linewidth = 1.2) +
    geom_text(aes(label = format(count, big.mark = ",")),
              fontface = "bold", size = 4.5) +
    scale_fill_gradient(low = "white", high = "#2166ac", name = "log10(count+1)") +
    scale_color_manual(values = c("TRUE" = "#d7301f", "FALSE" = "grey85"),
                       guide = "none") +
    scale_x_discrete(position = "top") +
    facet_wrap(~ condition, nrow = 1) +
    labs(x = "Late label", y = "Early label",
         title = "Developmental Transitions: Early (P12) to Late (Adult)") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_text(face = "bold", size = 12),
      strip.text   = element_text(face = "bold", size = 13)
    )
}

# -- Load data ----------------------------------------------------------------

cat("=== Loading A2 label files ===\n")
labels <- list()
for (tp in TPS) {
  labels[[tp]] <- fread(input_paths[[tp]])
  cat(sprintf("  %s: %d bins (%d callable)\n", tp, nrow(labels[[tp]]),
              labels[[tp]][!is.na(ctrl_label) & !is.na(mut_label), .N]))
}

# -- Phase 1: Within-timepoint transitions ------------------------------------

cat("\n=== Phase 1: Within-timepoint transitions ===\n")

trans_matrices <- list()
trans_stats    <- list()
rates_all      <- list()

for (tp in TPS) {
  mat <- compute_transition_matrix(labels[[tp]], "ctrl_label", "mut_label")
  trans_matrices[[tp]] <- mat

  total <- sum(mat)
  n_changed <- total - sum(diag(mat))
  pct_changed <- 100 * n_changed / total

  chisq_result <- suppressWarnings(chisq.test(mat))
  n <- sum(mat)
  k <- min(nrow(mat), ncol(mat))
  cramer_v <- sqrt(chisq_result$statistic / (n * (k - 1)))

  trans_stats[[tp]] <- list(
    total       = total,
    n_changed   = n_changed,
    pct_changed = pct_changed,
    chisq       = chisq_result,
    cramer_v    = unname(cramer_v)
  )

  rates_all[[tp]] <- matrix_to_long(mat, tp = tp, tp_label = TP_LABELS[tp])

  cat(sprintf("  %s (%s):\n", tp, TP_LABELS[tp]))
  cat(sprintf("    Callable: %d | Changed: %d (%.1f%%)\n",
              total, n_changed, pct_changed))
  cat(sprintf("    X2=%.1f, p=%.2e, V=%.3f\n",
              chisq_result$statistic, chisq_result$p.value, unname(cramer_v)))
}

rates_combined <- rbindlist(rates_all)

# -- Phase 2: Cross-timepoint developmental join ------------------------------

cat("\n=== Phase 2: Cross-timepoint developmental join ===\n")

dev_join <- merge(
  labels[["250831"]][, .(chr, bin_start,
                         ctrl_early = ctrl_label, mut_early = mut_label,
                         rank_ctrl_early = continous_rank_ctrl,
                         rank_mut_early  = continous_rank_mut)],
  labels[["250402"]][, .(chr, bin_start,
                         ctrl_late = ctrl_label, mut_late = mut_label,
                         rank_ctrl_late = continous_rank_ctrl,
                         rank_mut_late  = continous_rank_mut)],
  by = c("chr", "bin_start"),
  all = TRUE
)

n_ctrl_callable <- dev_join[!is.na(ctrl_early) & !is.na(ctrl_late), .N]
n_mut_callable  <- dev_join[!is.na(mut_early)  & !is.na(mut_late),  .N]
cat(sprintf("  Total bins: %d\n", nrow(dev_join)))
cat(sprintf("  Ctrl callable both TPs: %d\n", n_ctrl_callable))
cat(sprintf("  Mut callable both TPs:  %d\n", n_mut_callable))

dev_join[, ctrl_dev_changed := (!is.na(ctrl_early) & !is.na(ctrl_late) &
                                 ctrl_early != ctrl_late)]
dev_join[, mut_dev_changed  := (!is.na(mut_early)  & !is.na(mut_late) &
                                 mut_early != mut_late)]

ctrl_dev_pct <- 100 * sum(dev_join$ctrl_dev_changed, na.rm = TRUE) / n_ctrl_callable
mut_dev_pct  <- 100 * sum(dev_join$mut_dev_changed,  na.rm = TRUE) / n_mut_callable

cat(sprintf("  Ctrl developmental change: %d bins (%.1f%%)\n",
            sum(dev_join$ctrl_dev_changed, na.rm = TRUE), ctrl_dev_pct))
cat(sprintf("  Mut developmental change:  %d bins (%.1f%%)\n",
            sum(dev_join$mut_dev_changed, na.rm = TRUE), mut_dev_pct))

# -- Phase 3: Developmental transition matrices -------------------------------

cat("\n=== Phase 3: Developmental transition matrices ===\n")

dev_mat_ctrl <- compute_transition_matrix(dev_join, "ctrl_early", "ctrl_late")
dev_mat_mut  <- compute_transition_matrix(dev_join, "mut_early",  "mut_late")

dev_chisq_ctrl <- suppressWarnings(chisq.test(dev_mat_ctrl))
dev_chisq_mut  <- suppressWarnings(chisq.test(dev_mat_mut))

cat(sprintf("  Ctrl early->late: X2=%.1f, p=%.2e\n",
            dev_chisq_ctrl$statistic, dev_chisq_ctrl$p.value))
cat(sprintf("  Mut early->late:  X2=%.1f, p=%.2e\n",
            dev_chisq_mut$statistic, dev_chisq_mut$p.value))

dev_long_ctrl <- matrix_to_long(dev_mat_ctrl)
dev_long_ctrl[, condition := "Ctrl"]
dev_long_mut  <- matrix_to_long(dev_mat_mut)
dev_long_mut[, condition := "Mut"]
dev_long_all  <- rbindlist(list(dev_long_ctrl, dev_long_mut))

# -- Phase 4: Stability summary ----------------------------------------------

cat("\n=== Phase 4: Stability summary ===\n")

stability_dt <- data.table(
  comparison  = c("Ctrl developmental", "Mut developmental",
                  "Early genotype", "Late genotype"),
  n_callable  = c(n_ctrl_callable, n_mut_callable,
                  trans_stats[["250831"]]$total, trans_stats[["250402"]]$total),
  n_changed   = c(sum(dev_join$ctrl_dev_changed, na.rm = TRUE),
                  sum(dev_join$mut_dev_changed, na.rm = TRUE),
                  trans_stats[["250831"]]$n_changed,
                  trans_stats[["250402"]]$n_changed),
  axis        = c("Developmental", "Developmental", "Genotype", "Genotype")
)
stability_dt[, pct_changed := 100 * n_changed / n_callable]

cat("  Change rates:\n")
for (i in seq_len(nrow(stability_dt))) {
  cat(sprintf("    %s: %.1f%% (%d / %d)\n",
              stability_dt$comparison[i], stability_dt$pct_changed[i],
              stability_dt$n_changed[i], stability_dt$n_callable[i]))
}

homogeneity_test <- suppressWarnings(
  chisq.test(rbind(as.vector(dev_mat_ctrl), as.vector(dev_mat_mut)))
)
cat(sprintf("  Ctrl vs Mut developmental homogeneity: X2=%.1f, p=%.2e\n",
            homogeneity_test$statistic, homogeneity_test$p.value))

# -- Write TSVs ---------------------------------------------------------------

cat("\n=== Writing TSV outputs ===\n")

out_rates <- file.path(OUT_DIR, "combined_transition_rates.tsv")
fwrite(rates_combined, out_rates, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_rates), nrow(rates_combined)))

dev_out <- dev_join[, .(chr, bin_start,
                        ctrl_early, ctrl_late, mut_early, mut_late,
                        ctrl_dev_changed, mut_dev_changed)]
out_stability <- file.path(OUT_DIR, "developmental_stability.tsv")
fwrite(dev_out, out_stability, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_stability), nrow(dev_out)))

out_dev_mat <- file.path(OUT_DIR, "developmental_transition_matrix.tsv")
fwrite(dev_long_all, out_dev_mat, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_dev_mat), nrow(dev_long_all)))

# -- Figures ------------------------------------------------------------------

cat("\n=== Generating figures ===\n")

# Figure 1: Genome-wide karyotype track
cat("  Figure 1: Karyotype track...\n")
karyotype_dt <- build_karyotype_data(labels)
p_karyotype <- plot_karyotype_track(karyotype_dt)
save_multiformat_ggplot(p_karyotype,
  file.path(OUT_DIR, "combined_karyotype_track"),
  width = 16, height = 12)

# Figure 2: Combined transition heatmap (early | late)
cat("  Figure 2: Combined transition heatmap...\n")
heatmap_long <- rbindlist(lapply(TPS, function(tp) {
  dt <- as.data.table(as.table(trans_matrices[[tp]]))
  names(dt) <- c("from_label", "to_label", "count")
  dt[, `:=`(tp_label = TP_LABELS[tp],
            is_diagonal = (as.character(from_label) == as.character(to_label)))]
  dt
}))
p_heatmap <- plot_combined_heatmap(heatmap_long)
save_multiformat_ggplot(p_heatmap,
  file.path(OUT_DIR, "combined_transition_heatmap"),
  width = 12, height = 6)

# Figure 3: Developmental comparison (early vs late transition rates)
cat("  Figure 3: Developmental comparison...\n")
p_dev_compare <- plot_developmental_comparison(rates_combined)
save_multiformat_ggplot(p_dev_compare,
  file.path(OUT_DIR, "combined_developmental_comparison"),
  width = 12, height = 7)

# Figure 4a: Stability barplot
cat("  Figure 4a: Stability barplot...\n")
p_stability <- plot_stability_summary(stability_dt)
save_multiformat_ggplot(p_stability,
  file.path(OUT_DIR, "combined_stability_barplot"),
  width = 7, height = 5.5)

# Figure 4b: Developmental transition heatmap (ctrl | mut, early -> late)
cat("  Figure 4b: Developmental transition heatmap...\n")
p_dev_heatmap <- plot_developmental_heatmap(dev_long_all)
save_multiformat_ggplot(p_dev_heatmap,
  file.path(OUT_DIR, "combined_developmental_heatmap"),
  width = 12, height = 6)

# -- Verification -------------------------------------------------------------

cat("\n=== Verification ===\n")
any_warn <- FALSE

for (tp in TPS) {
  total <- trans_stats[[tp]]$total
  cat(sprintf("  %s matrix total: %d\n", tp, total))
  if (total < 20000 || total > 25000) {
    warning(sprintf("Unexpected bin count for %s: %d", tp, total))
    any_warn <- TRUE
  }
}

if (n_ctrl_callable < 23000) {
  warning(sprintf("Low cross-TP ctrl coverage: %d", n_ctrl_callable))
  any_warn <- TRUE
}
if (n_mut_callable < 23000) {
  warning(sprintf("Low cross-TP mut coverage: %d", n_mut_callable))
  any_warn <- TRUE
}
cat(sprintf("  Cross-TP coverage: ctrl=%d, mut=%d\n",
            n_ctrl_callable, n_mut_callable))

for (rate_name in c("ctrl_dev_pct", "mut_dev_pct")) {
  rate_val <- get(rate_name)
  if (rate_val < 5 || rate_val > 50) {
    warning(sprintf("Unusual developmental change rate %s: %.1f%%",
                    rate_name, rate_val))
    any_warn <- TRUE
  }
}
cat(sprintf("  Developmental rates: ctrl=%.1f%%, mut=%.1f%%\n",
            ctrl_dev_pct, mut_dev_pct))

a2_expected <- c("250402" = 15.3, "250831" = 18.3)
for (tp in TPS) {
  observed <- trans_stats[[tp]]$pct_changed
  expected <- a2_expected[tp]
  if (abs(observed - expected) > 1.0) {
    warning(sprintf("Within-TP %s change rate %.1f%% differs from A2 (%.1f%%)",
                    tp, observed, expected))
    any_warn <- TRUE
  }
}
cat(sprintf("  Within-TP rates: 250402=%.1f%% (A2: 15.3%%), 250831=%.1f%% (A2: 18.3%%)\n",
            trans_stats[["250402"]]$pct_changed,
            trans_stats[["250831"]]$pct_changed))

expected_tsvs <- c("combined_transition_rates.tsv",
                    "developmental_stability.tsv",
                    "developmental_transition_matrix.tsv")
expected_figs <- c("combined_karyotype_track",
                   "combined_transition_heatmap",
                   "combined_developmental_comparison",
                   "combined_stability_barplot",
                   "combined_developmental_heatmap")

for (f in expected_tsvs) {
  fpath <- file.path(OUT_DIR, f)
  if (!file.exists(fpath) || file.info(fpath)$size == 0) {
    warning(sprintf("Missing or empty output: %s", f))
    any_warn <- TRUE
  }
}
for (fig in expected_figs) {
  fig_dir <- file.path(OUT_DIR, fig)
  if (!dir.exists(fig_dir)) {
    warning(sprintf("Missing figure directory: %s", fig))
    any_warn <- TRUE
  } else {
    n_files <- length(list.files(fig_dir))
    if (n_files < 4) {
      warning(sprintf("Incomplete figure directory %s: %d files (expected 4)", fig, n_files))
      any_warn <- TRUE
    }
  }
}

if (!any_warn) cat("  All checks passed.\n")

elapsed <- (proc.time() - start_time)["elapsed"]
cat(sprintf("\n===========================================\n"))
cat(sprintf("C1 COMPLETE\n"))
cat(sprintf("Runtime: %.1f seconds\n", elapsed))
cat(sprintf("Outputs: %d TSVs + %d figure sets\n",
            length(expected_tsvs), length(expected_figs)))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
