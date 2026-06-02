# ML/cmpts/scripts/C4_homer_decomposition.R
# Stage C4: HOMER A/B x subcompartment decomposition.
#
# Decomposes HOMER's coarse A/B compartment transitions into subcompartment-
# resolution switches using the A4 100kb handoff file joined natively to
# CALDER2 labels. Adds PC1 effect-size statistics, genome-wide "iceberg"
# fractions, cross-timepoint comparison, SNIPER validation overlay, and
# per-chromosome breakdown.
#
# Usage:
#   Rscript C4_homer_decomposition.R <data_root> <code_root>
#     <data_root> : HPC data directory (kept for CLI consistency)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript C4_homer_decomposition.R <data_root> <code_root>")
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
  "A.1"    = "#e41a1c",
  "A.2"    = "#ff7f00",
  "B.1"    = "#4daf4a",
  "B.2"    = "#377eb8",
  "change" = "#cccccc"
)

AUTOSOMES <- paste0("chr", 1:19)
CHROM_ORDER <- AUTOSOMES

TRANSITION_LEVELS <- c("True_A_to_B", "True_B_to_A",
                        "Within_A_shift", "Within_B_shift",
                        "Stable", "Uncallable")
TRANSITION_COLORS <- c(
  True_A_to_B    = "#4dac26",
  True_B_to_A    = "#d01c8b",
  Within_A_shift = "#ff7f00",
  Within_B_shift = "#377eb8",
  Stable         = "#999999",
  Uncallable     = "#e0e0e0"
)

DIRECTION_LABELS <- c(
  "A_to_B_in_Mutant" = "HOMER: A→B",
  "B_to_A_in_Mutant" = "HOMER: B→A"
)

ICEBERG_ORDER <- c("Non_significant", "Sig_stable",
                    "Sig_within_shift", "Sig_true_flip", "Sig_uncallable")
ICEBERG_LABELS <- c(
  Non_significant  = "Not significant",
  Sig_stable       = "Significant + stable subcompartment",
  Sig_within_shift = "Significant + within-compartment shift",
  Sig_true_flip    = "Significant + true A/B flip",
  Sig_uncallable   = "Significant + uncallable"
)
ICEBERG_COLORS <- c(
  Non_significant  = "#eeeeee",
  Sig_stable       = "#999999",
  Sig_within_shift = "#ff7f00",
  Sig_true_flip    = "#e41a1c",
  Sig_uncallable   = "#e0e0e0"
)

AGREEMENT_COLORS <- c(
  Both_stable        = "#4daf4a",
  Both_change_agree  = "#377eb8",
  Both_change_disagree = "#e41a1c",
  SNIPER_only        = "#ff7f00",
  CALDER2_only       = "#984ea3"
)

CALDER2_BASE  <- file.path(CODE_ROOT, "outputs", "calder2")
INTEG_BASE    <- file.path(CODE_ROOT, "outputs", "integration")
SNIPER_TSV    <- file.path(CODE_ROOT, "outputs", "sniper", "tsvs")
OUT_DIR       <- file.path(CODE_ROOT, "outputs", "integration", "homer_decomposition")
TSV_DIR       <- file.path(OUT_DIR, "tsv")
PLOT_DIR      <- file.path(OUT_DIR, "plots")
UTIL_PATH     <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

# -- Libraries ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggalluvial)
  library(patchwork)
})
source(UTIL_PATH)

# -- Helper functions ---------------------------------------------------------

classify_transition <- function(ctrl_label, mut_label) {
  ctrl_chr <- substr(as.character(ctrl_label), 1, 1)
  mut_chr  <- substr(as.character(mut_label), 1, 1)

  result <- rep("Uncallable", length(ctrl_label))
  callable <- !is.na(ctrl_label) & !is.na(mut_label)
  same     <- callable & as.character(ctrl_label) == as.character(mut_label)
  a_to_b   <- callable & ctrl_chr == "A" & mut_chr == "B"
  b_to_a   <- callable & ctrl_chr == "B" & mut_chr == "A"
  within_a <- callable & !same & ctrl_chr == "A" & mut_chr == "A"
  within_b <- callable & !same & ctrl_chr == "B" & mut_chr == "B"

  result[same]     <- "Stable"
  result[a_to_b]   <- "True_A_to_B"
  result[b_to_a]   <- "True_B_to_A"
  result[within_a] <- "Within_A_shift"
  result[within_b] <- "Within_B_shift"
  result
}

assign_iceberg_category <- function(any_sig, transition_category) {
  result <- rep("Non_significant", length(any_sig))
  sig <- !is.na(any_sig) & any_sig == TRUE
  result[sig & transition_category == "Stable"]       <- "Sig_stable"
  result[sig & transition_category %in% c("Within_A_shift", "Within_B_shift")] <- "Sig_within_shift"
  result[sig & transition_category %in% c("True_A_to_B", "True_B_to_A")]      <- "Sig_true_flip"
  result[sig & transition_category == "Uncallable"]   <- "Sig_uncallable"
  result
}

sig_label <- function(p) {
  ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**",
      ifelse(p < 0.05, "*", "ns")))
}

# -- Header -------------------------------------------------------------------

start_time <- proc.time()

cat("===========================================\n")
cat("C4: HOMER A/B x Subcompartment Decomposition\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("Output dir: %s\n", OUT_DIR))
cat(sprintf("Timepoints: %s\n", paste(TPS, collapse = ", ")))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# -- Pre-flight validation ----------------------------------------------------

cat("=== Pre-flight validation ===\n")

homer_paths  <- list()
calder_paths <- list()
sniper_paths <- list()
SNIPER_AVAILABLE <- c()

for (tp in TPS) {
  hp <- file.path(INTEG_BASE, TP_DIRS[tp],
                  sprintf("%s_homer_100kb_aggregated.tsv", tp))
  if (!file.exists(hp)) stop(sprintf("Missing HOMER 100kb aggregated: %s", hp))
  if (file.info(hp)$size == 0) stop(sprintf("Empty HOMER 100kb aggregated: %s", hp))
  homer_paths[[tp]] <- hp
  cat(sprintf("  OK: HOMER 100kb %s (%s bytes)\n", tp,
              format(file.info(hp)$size, big.mark = ",")))

  cp <- file.path(CALDER2_BASE, TP_DIRS[tp],
                  sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  if (!file.exists(cp)) stop(sprintf("Missing CALDER2 labels: %s", cp))
  if (file.info(cp)$size == 0) stop(sprintf("Empty CALDER2 labels: %s", cp))
  calder_paths[[tp]] <- cp
  cat(sprintf("  OK: CALDER2 %s (%s bytes)\n", tp,
              format(file.info(cp)$size, big.mark = ",")))

  sp <- file.path(SNIPER_TSV, sprintf("%s_transition_concordance.tsv", tp))
  if (file.exists(sp) && file.info(sp)$size > 0) {
    sniper_paths[[tp]] <- sp
    SNIPER_AVAILABLE[tp] <- TRUE
    cat(sprintf("  OK: SNIPER concordance %s (%s bytes)\n", tp,
                format(file.info(sp)$size, big.mark = ",")))
  } else {
    SNIPER_AVAILABLE[tp] <- FALSE
    cat(sprintf("  [SKIP] SNIPER concordance not available for %s\n", tp))
  }
}

weakening_path <- file.path(INTEG_BASE, "combined_weakening_decomposition.tsv")
if (!file.exists(weakening_path)) {
  cat("  [NOTE] combined_weakening_decomposition.tsv not found; will compute rates from join\n")
  weakening_path <- NULL
} else {
  cat(sprintf("  OK: A4 weakening decomposition (%s bytes)\n",
              format(file.info(weakening_path)$size, big.mark = ",")))
}

if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
dir.create(TSV_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("  Pre-flight passed.\n\n")

# =============================================================================
# Phase 1: Load & Native 100kb Join
# =============================================================================

cat("=== Phase 1: Load & native 100kb join ===\n")

joined <- list()
n_tsvs <- 0L
n_figs <- 0L

for (tp in TPS) {
  cat(sprintf("\n  --- %s (%s) ---\n", tp, TP_LABELS[tp]))

  homer <- fread(homer_paths[[tp]])
  calder <- fread(calder_paths[[tp]])
  cat(sprintf("  HOMER 100kb: %d rows\n", nrow(homer)))
  cat(sprintf("  CALDER2:     %d rows (%d callable)\n",
              nrow(calder), calder[!is.na(ctrl_label) & !is.na(mut_label), .N]))

  j <- merge(homer, calder[, .(chr, bin_start, ctrl_label, mut_label,
                                continous_rank_ctrl, continous_rank_mut, label_changed)],
             by.x = c("Chr", "calder_bin_start"),
             by.y = c("chr", "bin_start"),
             all.x = TRUE)

  j[, transition_category := classify_transition(ctrl_label, mut_label)]
  j[, transition_category := factor(transition_category, levels = TRANSITION_LEVELS)]
  j[, iceberg_category := assign_iceberg_category(any_sig, as.character(transition_category))]
  j[, iceberg_category := factor(iceberg_category, levels = ICEBERG_ORDER)]

  n_callable  <- j[!is.na(ctrl_label), .N]
  n_sig       <- j[any_sig == TRUE, .N]
  n_true_flip <- j[iceberg_category == "Sig_true_flip", .N]
  n_within    <- j[iceberg_category == "Sig_within_shift", .N]
  n_stable    <- j[iceberg_category == "Sig_stable", .N]
  cat(sprintf("  Joined: %d rows (%d callable, %d uncallable)\n",
              nrow(j), n_callable, nrow(j) - n_callable))
  cat(sprintf("  Significant 100kb bins: %d (%.1f%% of total)\n",
              n_sig, 100 * n_sig / nrow(j)))
  cat(sprintf("    True flip:    %d (%.1f%% of sig)\n",
              n_true_flip, 100 * n_true_flip / max(1, n_sig)))
  cat(sprintf("    Within-shift: %d (%.1f%% of sig)\n",
              n_within, 100 * n_within / max(1, n_sig)))
  cat(sprintf("    Stable:       %d (%.1f%% of sig)\n",
              n_stable, 100 * n_stable / max(1, n_sig)))

  joined[[tp]] <- j

  out_file <- file.path(TSV_DIR, sprintf("%s_native_join_100kb.tsv", tp))
  fwrite(j, out_file, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s (%d rows)\n", basename(out_file), nrow(j)))
  n_tsvs <- n_tsvs + 1L
}
cat("\n")

# =============================================================================
# Phase 2: PC1 Effect-Size Distributions & Wilcoxon Tests
# =============================================================================

cat("=== Phase 2: PC1 effect-size distributions ===\n")

effect_summaries <- list()
wilcox_results   <- list()

for (tp in TPS) {
  cat(sprintf("\n  --- %s (%s) ---\n", tp, TP_LABELS[tp]))
  sig_bins <- joined[[tp]][any_sig == TRUE & transition_category != "Uncallable"]
  sig_bins[, abs_diff := abs(mean_Difference)]

  sig_bins[, effect_group := fifelse(
    transition_category %in% c("True_A_to_B", "True_B_to_A"), "True_flip",
    fifelse(transition_category %in% c("Within_A_shift", "Within_B_shift"), "Within_shift",
            "Stable"))]
  sig_bins[, effect_group := factor(effect_group,
    levels = c("True_flip", "Within_shift", "Stable"))]

  eff_summ <- sig_bins[, .(
    n            = .N,
    median_absDiff = median(abs_diff),
    mean_absDiff   = mean(abs_diff),
    sd_absDiff     = sd(abs_diff),
    q25_absDiff    = quantile(abs_diff, 0.25),
    q75_absDiff    = quantile(abs_diff, 0.75)
  ), by = effect_group]
  eff_summ[, tp := tp]
  eff_summ[, tp_label := TP_LABELS[tp]]
  effect_summaries[[tp]] <- eff_summ
  cat(sprintf("  Effect size summary:\n"))
  for (i in seq_len(nrow(eff_summ))) {
    cat(sprintf("    %s: n=%d, median|dPC1|=%.3f, mean=%.3f\n",
                eff_summ$effect_group[i], eff_summ$n[i],
                eff_summ$median_absDiff[i], eff_summ$mean_absDiff[i]))
  }

  groups <- c("True_flip", "Within_shift", "Stable")
  pairs  <- list(c("True_flip", "Within_shift"),
                 c("True_flip", "Stable"),
                 c("Within_shift", "Stable"))
  wt_rows <- list()
  for (pair in pairs) {
    g1_vals <- sig_bins[effect_group == pair[1], abs_diff]
    g2_vals <- sig_bins[effect_group == pair[2], abs_diff]
    if (length(g1_vals) >= 3 && length(g2_vals) >= 3) {
      wt <- wilcox.test(g1_vals, g2_vals, alternative = "two.sided")
      wt_rows[[length(wt_rows) + 1]] <- data.table(
        tp = tp, tp_label = TP_LABELS[tp],
        group1 = pair[1], group2 = pair[2],
        n1 = length(g1_vals), n2 = length(g2_vals),
        W_statistic = wt$statistic, p_value = wt$p.value
      )
      cat(sprintf("    Wilcoxon %s vs %s: W=%s, p=%.2e\n",
                  pair[1], pair[2],
                  format(wt$statistic, big.mark = ","), wt$p.value))
    }
  }
  wt_dt <- rbindlist(wt_rows)
  wt_dt[, p_adj := p.adjust(p_value, method = "BH")]
  wilcox_results[[tp]] <- wt_dt

  out_eff <- file.path(TSV_DIR, sprintf("%s_effect_size_summary.tsv", tp))
  fwrite(eff_summ, out_eff, sep = "\t", quote = FALSE, na = "NA")
  n_tsvs <- n_tsvs + 1L

  out_wt <- file.path(TSV_DIR, sprintf("%s_wilcoxon_tests.tsv", tp))
  fwrite(wt_dt, out_wt, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s, %s\n", basename(out_eff), basename(out_wt)))
  n_tsvs <- n_tsvs + 1L
}
cat("\n")

# =============================================================================
# Phase 3: Genome-Wide Fractions ("The Iceberg")
# =============================================================================

cat("=== Phase 3: Genome-wide iceberg fractions ===\n")

iceberg_rows <- list()
for (tp in TPS) {
  ib <- joined[[tp]][, .N, by = iceberg_category]
  ib[, pct := round(100 * N / sum(N), 2)]
  ib[, tp := tp]
  ib[, tp_label := TP_LABELS[tp]]
  iceberg_rows[[tp]] <- ib
  cat(sprintf("  %s (%s):\n", tp, TP_LABELS[tp]))
  for (i in seq_len(nrow(ib))) {
    cat(sprintf("    %s: %s (%.1f%%)\n", ib$iceberg_category[i],
                format(ib$N[i], big.mark = ","), ib$pct[i]))
  }
}
iceberg_dt <- rbindlist(iceberg_rows)

out_ice <- file.path(TSV_DIR, "iceberg_summary.tsv")
fwrite(iceberg_dt, out_ice, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n\n", basename(out_ice), nrow(iceberg_dt)))
n_tsvs <- n_tsvs + 1L

# =============================================================================
# Phase 4: Cross-Timepoint Comparison
# =============================================================================

cat("=== Phase 4: Cross-timepoint comparison ===\n")

# 4a: Decomposition rates per direction (from A4 or recomputed)
decomp_rows <- list()
for (tp in TPS) {
  sig_j <- joined[[tp]][any_sig == TRUE & !is.na(dominant_homer_direction) &
                         dominant_homer_direction != ""]
  decomp <- sig_j[transition_category != "Uncallable",
                   .N, by = .(dominant_homer_direction, transition_category)]
  decomp[, pct_of_direction := round(100 * N / sum(N), 2),
         by = dominant_homer_direction]
  decomp[, tp := tp]
  decomp[, tp_label := TP_LABELS[tp]]
  decomp_rows[[tp]] <- decomp
}
decomp_dt <- rbindlist(decomp_rows)
decomp_dt[, transition_category := factor(transition_category, levels = TRANSITION_LEVELS)]

out_decomp <- file.path(TSV_DIR, "decomposition_rates_combined.tsv")
fwrite(decomp_dt, out_decomp, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_decomp), nrow(decomp_dt)))
n_tsvs <- n_tsvs + 1L

# 4b: Bins significant in BOTH timepoints
sig_late  <- joined[["250402"]][any_sig == TRUE,
  .(Chr, calder_bin_start, transition_category, ctrl_label, mut_label, mean_Difference)]
sig_early <- joined[["250831"]][any_sig == TRUE,
  .(Chr, calder_bin_start, transition_category, ctrl_label, mut_label, mean_Difference)]

both_sig <- merge(sig_late, sig_early,
                  by = c("Chr", "calder_bin_start"),
                  suffixes = c("_late", "_early"))
both_sig[, concordant := as.character(transition_category_late) ==
                         as.character(transition_category_early)]
n_both <- nrow(both_sig)
n_concord <- both_sig[concordant == TRUE, .N]
cat(sprintf("  Bins significant in both TPs: %d\n", n_both))
cat(sprintf("  Concordant transition category: %d (%.1f%%)\n",
            n_concord, 100 * n_concord / max(1, n_both)))

out_both <- file.path(TSV_DIR, "cross_timepoint_sig_bins.tsv")
fwrite(both_sig, out_both, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n\n", basename(out_both), nrow(both_sig)))
n_tsvs <- n_tsvs + 1L

# =============================================================================
# Phase 5: SNIPER Validation Overlay
# =============================================================================

cat("=== Phase 5: SNIPER validation overlay ===\n")

sniper_validation <- list()
for (tp in TPS) {
  if (!SNIPER_AVAILABLE[tp]) {
    cat(sprintf("  [SKIP] SNIPER data not available for %s\n", tp))
    next
  }
  cat(sprintf("  --- %s (%s) ---\n", tp, TP_LABELS[tp]))

  sniper <- fread(sniper_paths[[tp]])
  sig_j <- joined[[tp]][any_sig == TRUE & transition_category != "Uncallable"]
  sig_j[, bin_start_1based := calder_bin_start]

  sv <- merge(sig_j, sniper[, .(chr, bin_start, sniper_changed, calder_changed, agreement)],
              by.x = c("Chr", "bin_start_1based"),
              by.y = c("chr", "bin_start"),
              all.x = FALSE)

  cat(sprintf("  Matched: %d / %d sig callable bins (%.1f%%)\n",
              nrow(sv), nrow(sig_j), 100 * nrow(sv) / max(1, nrow(sig_j))))

  sv[, effect_group := fifelse(
    transition_category %in% c("True_A_to_B", "True_B_to_A"), "True_flip",
    fifelse(transition_category %in% c("Within_A_shift", "Within_B_shift"), "Within_shift",
            "Stable"))]

  val_summary <- sv[, .(
    N = .N,
    n_both_stable     = sum(agreement == "Both_stable"),
    n_both_agree      = sum(agreement == "Both_change_agree"),
    n_both_disagree   = sum(agreement == "Both_change_disagree"),
    n_sniper_only     = sum(agreement == "SNIPER_only"),
    n_calder_only     = sum(agreement == "CALDER2_only")
  ), by = effect_group]
  val_summary[, pct_sniper_confirms := round(100 * (n_both_agree + n_both_stable) / N, 1)]
  val_summary[, tp := tp]
  val_summary[, tp_label := TP_LABELS[tp]]

  sniper_validation[[tp]] <- val_summary
  for (i in seq_len(nrow(val_summary))) {
    cat(sprintf("    %s: n=%d, SNIPER confirms=%.1f%%\n",
                val_summary$effect_group[i], val_summary$N[i],
                val_summary$pct_sniper_confirms[i]))
  }

  out_sv <- file.path(TSV_DIR, sprintf("%s_sniper_homer_validation.tsv", tp))
  fwrite(val_summary, out_sv, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s\n", basename(out_sv)))
  n_tsvs <- n_tsvs + 1L

  sniper_validation[[tp]] <- list(summary = val_summary, detail = sv)
}
cat("\n")

# =============================================================================
# Phase 6: Per-Chromosome Breakdown
# =============================================================================

cat("=== Phase 6: Per-chromosome breakdown ===\n")

chrom_rows <- list()
for (tp in TPS) {
  cb <- joined[[tp]][!is.na(ctrl_label), .(
    n_total      = .N,
    n_any_sig    = sum(any_sig, na.rm = TRUE),
    n_true_flip  = sum(iceberg_category == "Sig_true_flip", na.rm = TRUE),
    n_within_shift = sum(iceberg_category == "Sig_within_shift", na.rm = TRUE),
    n_stable_sig = sum(iceberg_category == "Sig_stable", na.rm = TRUE)
  ), by = Chr]
  cb[, pct_sig := round(100 * n_any_sig / n_total, 1)]
  cb[, pct_true_flip_of_sig := round(100 * n_true_flip / pmax(1L, n_any_sig), 1)]
  cb[, tp := tp]
  cb[, tp_label := TP_LABELS[tp]]
  cb[, Chr := factor(Chr, levels = CHROM_ORDER)]
  chrom_rows[[tp]] <- cb
}
chrom_dt <- rbindlist(chrom_rows)

out_chrom <- file.path(TSV_DIR, "chromosome_breakdown.tsv")
fwrite(chrom_dt, out_chrom, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_chrom), nrow(chrom_dt)))
n_tsvs <- n_tsvs + 1L

# Summary table
summary_rows <- list()
for (tp in TPS) {
  j <- joined[[tp]]
  n_sig <- j[any_sig == TRUE, .N]
  summary_rows[[tp]] <- data.table(
    tp = tp, tp_label = TP_LABELS[tp],
    total_bins           = nrow(j),
    callable_bins        = j[!is.na(ctrl_label), .N],
    sig_bins             = n_sig,
    pct_sig              = round(100 * n_sig / nrow(j), 1),
    n_true_flip          = j[iceberg_category == "Sig_true_flip", .N],
    n_within_shift       = j[iceberg_category == "Sig_within_shift", .N],
    n_sig_stable         = j[iceberg_category == "Sig_stable", .N],
    pct_true_flip_of_sig = round(100 * j[iceberg_category == "Sig_true_flip", .N] /
                                  max(1, n_sig), 1),
    median_absDiff_flip  = effect_summaries[[tp]][effect_group == "True_flip", median_absDiff],
    median_absDiff_shift = effect_summaries[[tp]][effect_group == "Within_shift", median_absDiff],
    n_both_sig_bins      = n_both,
    sniper_available     = SNIPER_AVAILABLE[tp]
  )
}
summary_dt <- rbindlist(summary_rows)
out_summ <- file.path(TSV_DIR, "c4_combined_summary.tsv")
fwrite(summary_dt, out_summ, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n\n", basename(out_summ), nrow(summary_dt)))
n_tsvs <- n_tsvs + 1L

# =============================================================================
# Figures
# =============================================================================

cat("=== Generating figures ===\n")

theme_c4 <- theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

# --- Figure 1: Iceberg stacked bar -------------------------------------------

cat("  Fig 1: Iceberg genome-wide fractions\n")

ice_plot <- copy(iceberg_dt)
ice_plot[, iceberg_label := factor(ICEBERG_LABELS[as.character(iceberg_category)],
  levels = rev(ICEBERG_LABELS[ICEBERG_ORDER]))]
ice_plot[, tp_label := factor(tp_label, levels = TP_LABELS)]

p_iceberg <- ggplot(ice_plot[iceberg_category != "Sig_uncallable"],
                    aes(x = tp_label, y = pct, fill = iceberg_label)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(pct >= 3, sprintf("%.1f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = setNames(
    rev(ICEBERG_COLORS[ICEBERG_ORDER[ICEBERG_ORDER != "Sig_uncallable"]]),
    rev(ICEBERG_LABELS[ICEBERG_ORDER[ICEBERG_ORDER != "Sig_uncallable"]])),
    name = "Category") +
  labs(x = "Timepoint", y = "% of 100kb bins",
       title = "Genome-Wide Decomposition of HOMER Compartment Changes",
       subtitle = "Most bins show no significant change; true A/B flips are rare") +
  theme_c4 +
  coord_flip()

save_multiformat_ggplot(p_iceberg,
  file.path(PLOT_DIR, "c4_iceberg_stacked"), width = 10, height = 5)
n_figs <- n_figs + 1L

# --- Figure 2: Effect-size violin + box --------------------------------------

cat("  Fig 2: PC1 effect-size distributions\n")

violin_data <- list()
for (tp in TPS) {
  d <- joined[[tp]][any_sig == TRUE & transition_category != "Uncallable"]
  d[, abs_diff := abs(mean_Difference)]
  d[, effect_group := fifelse(
    transition_category %in% c("True_A_to_B", "True_B_to_A"), "True_flip",
    fifelse(transition_category %in% c("Within_A_shift", "Within_B_shift"), "Within_shift",
            "Stable"))]
  d[, effect_group := factor(effect_group, levels = c("True_flip", "Within_shift", "Stable"))]
  d[, tp_label := TP_LABELS[tp]]
  violin_data[[tp]] <- d[, .(abs_diff, effect_group, tp_label)]
}
violin_dt <- rbindlist(violin_data)
violin_dt[, tp_label := factor(tp_label, levels = TP_LABELS)]

eff_colors <- c(True_flip = "#e41a1c", Within_shift = "#ff7f00", Stable = "#999999")

p_violin <- ggplot(violin_dt, aes(x = effect_group, y = abs_diff, fill = effect_group)) +
  geom_violin(scale = "width", alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
  facet_wrap(~ tp_label, ncol = 2) +
  scale_fill_manual(values = eff_colors, guide = "none") +
  scale_y_continuous(breaks = seq(0, 2, 0.2)) +
  labs(x = "Transition category", y = "|Mean ΔPC1| (100kb bins)",
       title = "True Compartment Flips Have Larger PC1 Changes",
       subtitle = "Significant HOMER 100kb bins stratified by CALDER2 subcompartment transition") +
  theme_c4

wt_all <- rbindlist(wilcox_results)
anno_data <- list()
for (tp in TPS) {
  wt <- wilcox_results[[tp]]
  for (i in seq_len(nrow(wt))) {
    ymax <- violin_dt[tp_label == TP_LABELS[tp], max(abs_diff, na.rm = TRUE)]
    anno_data[[length(anno_data) + 1]] <- data.table(
      tp_label = TP_LABELS[tp],
      group1 = wt$group1[i], group2 = wt$group2[i],
      label = sig_label(wt$p_adj[i]),
      y = ymax * (0.85 + 0.08 * i)
    )
  }
}
anno_dt <- rbindlist(anno_data)
anno_dt[, tp_label := factor(tp_label, levels = TP_LABELS)]

p_violin <- p_violin +
  geom_text(data = anno_dt,
            aes(x = 1.5, y = y, label = sprintf("%s vs %s: %s", group1, group2, label)),
            inherit.aes = FALSE, size = 3, hjust = 0.5)

save_multiformat_ggplot(p_violin,
  file.path(PLOT_DIR, "c4_effect_size_violin"), width = 11, height = 7)
n_figs <- n_figs + 1L

# --- Figure 3: Decomposition comparison stacked bar --------------------------

cat("  Fig 3: Decomposition comparison (early vs late)\n")

decomp_plot <- copy(decomp_dt)
decomp_plot[, direction_label := DIRECTION_LABELS[dominant_homer_direction]]
decomp_plot[, direction_label := factor(direction_label, levels = DIRECTION_LABELS)]
decomp_plot[, tp_label := factor(tp_label, levels = TP_LABELS)]
decomp_plot[, transition_category := factor(transition_category, levels = rev(TRANSITION_LEVELS))]

p_decomp <- ggplot(decomp_plot,
                    aes(x = direction_label, y = pct_of_direction, fill = transition_category)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(pct_of_direction >= 5,
                                sprintf("%.0f%%", pct_of_direction), "")),
            position = position_stack(vjust = 0.5), size = 3.2, fontface = "bold") +
  facet_wrap(~ tp_label, ncol = 2) +
  scale_fill_manual(values = TRANSITION_COLORS, name = "Transition type",
                    breaks = TRANSITION_LEVELS) +
  labs(x = "HOMER direction", y = "% of significant bins",
       title = "Subcompartment Decomposition of HOMER A/B Transitions",
       subtitle = "True flips vs within-compartment shifts vs stable bins") +
  theme_c4

save_multiformat_ggplot(p_decomp,
  file.path(PLOT_DIR, "c4_decomposition_comparison"), width = 12, height = 7)
n_figs <- n_figs + 1L

# --- Figure 4: Cross-TP alluvial ---------------------------------------------

cat("  Fig 4: Cross-timepoint alluvial\n")

if (n_both >= 10) {
  both_alluv <- both_sig[, .N,
    by = .(transition_category_early, transition_category_late)]
  both_alluv[, transition_category_early := factor(
    as.character(transition_category_early), levels = TRANSITION_LEVELS)]
  both_alluv[, transition_category_late := factor(
    as.character(transition_category_late), levels = TRANSITION_LEVELS)]

  p_alluv <- ggplot(both_alluv,
    aes(axis1 = transition_category_early, axis2 = transition_category_late,
        y = N)) +
    geom_alluvium(aes(fill = transition_category_late), alpha = 0.7) +
    geom_stratum(width = 0.3, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Early (P12)", "Late (Adult)"), expand = c(0.15, 0.05)) +
    scale_fill_manual(values = TRANSITION_COLORS, name = "Late classification") +
    labs(y = "Number of 100kb bins",
         title = sprintf("Transition Category Stability Across Timepoints (n=%d bins)", n_both),
         subtitle = "Bins significant in BOTH early and late") +
    theme_c4

  save_multiformat_ggplot(p_alluv,
    file.path(PLOT_DIR, "c4_cross_tp_alluvial"), width = 10, height = 8)
  n_figs <- n_figs + 1L
} else {
  cat("    [SKIP] Too few bins significant in both TPs for alluvial\n")
}

# --- Figure 5: Per-chromosome heatmap ----------------------------------------

cat("  Fig 5: Per-chromosome heatmap\n")

chrom_plot <- copy(chrom_dt)
chrom_plot[, Chr := factor(Chr, levels = rev(CHROM_ORDER))]
chrom_plot[, tp_label := factor(tp_label, levels = TP_LABELS)]

p_chrom <- ggplot(chrom_plot, aes(x = tp_label, y = Chr,
                                   fill = pct_true_flip_of_sig)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d", n_true_flip)), size = 3) +
  scale_fill_gradient2(low = "white", mid = "#fee0d2", high = "#e41a1c",
                       midpoint = 15, name = "% true flips\n(of sig bins)") +
  labs(x = "Timepoint", y = "Chromosome",
       title = "Chromosomal Heterogeneity in True Compartment Flips",
       subtitle = "Count of true-flip 100kb bins shown in each cell") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())

save_multiformat_ggplot(p_chrom,
  file.path(PLOT_DIR, "c4_chromosome_heatmap"), width = 7, height = 10)
n_figs <- n_figs + 1L

# --- Figure 6: SNIPER validation bar -----------------------------------------

cat("  Fig 6: SNIPER validation\n")

any_sniper <- any(SNIPER_AVAILABLE)
if (any_sniper) {
  sv_rows <- list()
  for (tp in TPS) {
    if (!SNIPER_AVAILABLE[tp]) next
    sv <- sniper_validation[[tp]]$detail
    sv_long <- sv[, .(N = .N), by = .(effect_group, agreement)]
    sv_long[, pct := round(100 * N / sum(N), 1), by = effect_group]
    sv_long[, tp := tp]
    sv_long[, tp_label := TP_LABELS[tp]]
    sv_rows[[length(sv_rows) + 1]] <- sv_long
  }
  sv_plot <- rbindlist(sv_rows)
  sv_plot[, effect_group := factor(effect_group,
    levels = c("True_flip", "Within_shift", "Stable"))]
  sv_plot[, agreement := factor(agreement,
    levels = names(AGREEMENT_COLORS))]
  sv_plot[, tp_label := factor(tp_label, levels = TP_LABELS)]

  p_sniper <- ggplot(sv_plot, aes(x = effect_group, y = pct, fill = agreement)) +
    geom_col(position = "stack", width = 0.6) +
    geom_text(aes(label = ifelse(pct >= 5, sprintf("%.0f%%", pct), "")),
              position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
    facet_wrap(~ tp_label, ncol = 2) +
    scale_fill_manual(values = AGREEMENT_COLORS, name = "SNIPER-CALDER2\nagreement") +
    labs(x = "Transition category (HOMER x CALDER2)", y = "% of matched bins",
         title = "SNIPER Validation of Subcompartment Transitions",
         subtitle = "For significant HOMER bins: does SNIPER independently confirm the transition?") +
    theme_c4

  save_multiformat_ggplot(p_sniper,
    file.path(PLOT_DIR, "c4_sniper_validation"), width = 12, height = 7)
  n_figs <- n_figs + 1L
} else {
  cat("    [SKIP] No SNIPER data available\n")
}

# --- Figure 7: Publication panel (patchwork) ----------------------------------

cat("  Fig 7: Publication panel\n")

p_pub <- (p_iceberg | p_decomp) /
         (p_violin  | p_chrom) +
  plot_annotation(
    title = "HOMER A/B Compartment Transitions Decomposed by Subcompartment",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

save_multiformat_ggplot(p_pub,
  file.path(PLOT_DIR, "c4_publication_panel"), width = 20, height = 16)
n_figs <- n_figs + 1L

cat(sprintf("  Generated %d figure sets\n\n", n_figs))

# -- Verification --------------------------------------------------------------

cat("=== Verification ===\n")
any_warn <- FALSE

for (tp in TPS) {
  n_joined <- nrow(joined[[tp]])
  n_homer  <- nrow(fread(homer_paths[[tp]], select = 1L))
  if (n_joined != n_homer) {
    warning(sprintf("Row count mismatch %s: joined=%d, homer=%d", tp, n_joined, n_homer))
    any_warn <- TRUE
  } else {
    cat(sprintf("  [PASS] %s: join preserved %d rows\n", tp, n_joined))
  }
}

for (tp in TPS) {
  n_flip <- joined[[tp]][iceberg_category == "Sig_true_flip", .N]
  if (n_flip == 0) {
    warning(sprintf("No true flip bins for %s", tp))
    any_warn <- TRUE
  } else {
    cat(sprintf("  [PASS] %s: %d true-flip bins\n", tp, n_flip))
  }
}

for (tp in TPS) {
  wt <- wilcox_results[[tp]]
  if (nrow(wt) < 3) {
    warning(sprintf("Incomplete Wilcoxon tests for %s: %d rows", tp, nrow(wt)))
    any_warn <- TRUE
  } else {
    cat(sprintf("  [PASS] %s: %d Wilcoxon tests computed\n", tp, nrow(wt)))
  }
}

expected_tsvs <- c("iceberg_summary.tsv", "decomposition_rates_combined.tsv",
                    "cross_timepoint_sig_bins.tsv", "chromosome_breakdown.tsv",
                    "c4_combined_summary.tsv")
for (tp in TPS) {
  expected_tsvs <- c(expected_tsvs,
    sprintf("%s_native_join_100kb.tsv", tp),
    sprintf("%s_effect_size_summary.tsv", tp),
    sprintf("%s_wilcoxon_tests.tsv", tp))
  if (SNIPER_AVAILABLE[tp]) {
    expected_tsvs <- c(expected_tsvs,
      sprintf("%s_sniper_homer_validation.tsv", tp))
  }
}
for (f in expected_tsvs) {
  fpath <- file.path(TSV_DIR, f)
  if (!file.exists(fpath) || file.info(fpath)$size == 0) {
    warning(sprintf("Missing or empty output: %s", f))
    any_warn <- TRUE
  }
}

expected_figs <- c("c4_iceberg_stacked", "c4_effect_size_violin",
                    "c4_decomposition_comparison", "c4_chromosome_heatmap",
                    "c4_publication_panel")
if (n_both >= 10) expected_figs <- c(expected_figs, "c4_cross_tp_alluvial")
if (any_sniper) expected_figs <- c(expected_figs, "c4_sniper_validation")

for (fig in expected_figs) {
  fig_dir <- file.path(PLOT_DIR, fig)
  if (!dir.exists(fig_dir)) {
    warning(sprintf("Missing figure directory: %s", fig))
    any_warn <- TRUE
  } else {
    n_files <- length(list.files(fig_dir))
    if (n_files < 4) {
      warning(sprintf("Incomplete figure directory %s: %d files", fig, n_files))
      any_warn <- TRUE
    }
  }
}

if (!any_warn) cat("  All checks passed.\n")

elapsed <- (proc.time() - start_time)["elapsed"]
cat(sprintf("\n===========================================\n"))
cat(sprintf("C4 COMPLETE\n"))
cat(sprintf("Runtime: %.1f seconds\n", elapsed))
cat(sprintf("Outputs: %d TSVs + %d figure sets\n", n_tsvs, n_figs))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
