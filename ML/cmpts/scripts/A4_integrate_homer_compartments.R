# ML/cmpts/scripts/A4_integrate_homer_compartments.R
# Stage A4: Integration of HOMER A/B compartment calls with CALDER2 subcompartments.
#
# Cross-references significant HOMER differential bins (25kb) with CALDER2
# subcompartment labels (100kb) to decompose coarse A/B transitions into
# subcompartment-level events: true flips, within-compartment shifts, or
# subcompartment-stable bins.
#
# Processes both timepoints in a single invocation (like A3).
#
# Usage:
#   Rscript A4_integrate_homer_compartments.R <data_root> <code_root>
#     <data_root> : HPC data directory (e.g. /expanse/.../sniper)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript A4_integrate_homer_compartments.R <data_root> <code_root>")
}

DATA_ROOT <- args[1]
CODE_ROOT <- args[2]

# ── Constants ──────────────────────────────────────────────────────────────────

TPS <- c("250402", "250831")
TP_LABELS <- c("250402" = "late/adult", "250831" = "early/P12")

LABEL_ORDER <- c("A.1", "A.2", "B.1", "B.2")
LABEL_COLORS <- c(
  "A.1"    = "#e41a1c",
  "A.2"    = "#ff7f00",
  "B.1"    = "#4daf4a",
  "B.2"    = "#377eb8",
  "change" = "#cccccc"
)

MM10_CHROM_SIZES <- c(
  chr1  = 195471971L, chr2  = 182113224L, chr3  = 160039680L, chr4  = 156508116L,
  chr5  = 151834684L, chr6  = 149736546L, chr7  = 145441459L, chr8  = 129401213L,
  chr9  = 124595110L, chr10 = 130694993L, chr11 = 122082543L, chr12 = 120129022L,
  chr13 = 120421639L, chr14 = 124902244L, chr15 = 104043685L, chr16 = 98207768L,
  chr17 = 94987271L,  chr18 = 90702639L,  chr19 = 61431566L
)
AUTOSOMES <- names(MM10_CHROM_SIZES)

HOMER_SIG_FDR  <- 0.05
HOMER_SIG_DIFF <- 0.30

REPO_ROOT <- dirname(dirname(CODE_ROOT))
HOMER_PATHS <- c(
  "250402" = file.path(REPO_ROOT, "tads", "tad-pc-analysis", "output",
                       "compartment_analysis", "compartment_all_annotated.tsv"),
  "250831" = file.path(REPO_ROOT, "tads", "tad-pc-analysis", "output",
                       "compartment_analysis_early", "compartment_all_annotated.tsv")
)

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

TP_DIRS     <- c("250402" = "late", "250831" = "early")
CALDER2_BASE <- file.path(CODE_ROOT, "outputs", "calder2")
OUT_BASE     <- file.path(CODE_ROOT, "outputs", "integration")
UTIL_PATH   <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

# ── Libraries ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(scales)
})
source(UTIL_PATH)

# ── Header ─────────────────────────────────────────────────────────────────────

start_time <- proc.time()

cat("===========================================\n")
cat("A4: HOMER A/B x CALDER2 Subcompartment Integration\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("REPO_ROOT:  %s\n", REPO_ROOT))
cat(sprintf("Timepoints: %s\n", paste(TPS, collapse = ", ")))
cat(sprintf("Output base: %s\n", OUT_BASE))
cat(sprintf("Sig thresh: FDR<%.2f, |Diff|>%.2f\n", HOMER_SIG_FDR, HOMER_SIG_DIFF))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# ── Pre-flight validation ──────────────────────────────────────────────────────

cat("=== Pre-flight validation ===\n")

for (tp in TPS) {
  hp <- HOMER_PATHS[tp]
  if (!file.exists(hp)) stop(sprintf("Missing HOMER file: %s", hp))
  if (file.info(hp)$size == 0) stop(sprintf("Empty HOMER file: %s", hp))
  cat(sprintf("  OK: HOMER %s (%s bytes)\n", tp,
              format(file.info(hp)$size, big.mark = ",")))
}

for (tp in TPS) {
  lp <- file.path(CALDER2_BASE, TP_DIRS[tp],
                   sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  if (!file.exists(lp)) stop(sprintf("Missing CALDER2 labels: %s", lp))
  if (file.info(lp)$size == 0) stop(sprintf("Empty CALDER2 labels: %s", lp))
  cat(sprintf("  OK: CALDER2 %s (%s bytes)\n", tp,
              format(file.info(lp)$size, big.mark = ",")))
}

if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
cat("  Pre-flight passed.\n")

# ── Function definitions ───────────────────────────────────────────────────────

load_homer <- function(path, autosomes) {
  dt <- fread(path)
  n_raw <- nrow(dt)
  dt <- dt[Chr %in% autosomes]
  n_auto <- nrow(dt)
  cat(sprintf("  Loaded HOMER: %s total, %s autosomal (%d on sex/other chroms dropped)\n",
              format(n_raw, big.mark = ","),
              format(n_auto, big.mark = ","),
              n_raw - n_auto))

  if ("Annotation" %in% names(dt) && !all(is.na(dt$Annotation))) {
    is_intergenic   <- startsWith(dt$Annotation, "Intergenic")
    genic_mean      <- mean(dt$ctrl_avg_PC1[!is_intergenic], na.rm = TRUE)
    intergenic_mean <- mean(dt$ctrl_avg_PC1[is_intergenic],  na.rm = TRUE)
    cat(sprintf("  Polarity check: genic mean=%+.4f, intergenic mean=%+.4f\n",
                genic_mean, intergenic_mean))

    if (intergenic_mean > genic_mean) {
      cat("  ** PC1 polarity FLIPPED — correcting Difference and direction.\n")
      dt[, ctrl_avg_PC1 := -ctrl_avg_PC1]
      dt[, mut_avg_PC1  := -mut_avg_PC1]
      dt[, Difference   := -Difference]
      dt[, direction := fifelse(Difference > 0,
                                "B_to_A_in_Mutant", "A_to_B_in_Mutant")]
      dt[, direction_label := fifelse(Difference > 0,
                                      "More Active (B->A)", "More Inactive (A->B)")]
    } else {
      cat("  Polarity OK (genic > intergenic).\n")
    }
  }

  dt[, is_significant := (adj_pvalue < HOMER_SIG_FDR & abs(Difference) > HOMER_SIG_DIFF)]
  dt[, calder_bin_start := floor(Start / 100000L) * 100000L + 1L]
  dt
}

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

build_crosstab <- function(sig_dt) {
  ct <- sig_dt[!is.na(ctrl_label) & !is.na(mut_label),
               .N, by = .(direction, ctrl_label, mut_label)]
  ct[, transition_category := classify_transition(ctrl_label, mut_label)]
  ct[, pct_of_direction := round(100 * N / sum(N), 2), by = direction]
  ct[order(direction, ctrl_label, mut_label)]
}

aggregate_homer_100kb <- function(homer_dt) {
  agg <- homer_dt[, .(
    mean_ctrl_PC1  = mean(ctrl_avg_PC1, na.rm = TRUE),
    mean_mut_PC1   = mean(mut_avg_PC1, na.rm = TRUE),
    mean_Difference = mean(Difference, na.rm = TRUE),
    n_25kb_bins    = .N,
    n_sig_bins     = sum(is_significant),
    min_adj_pvalue = min(adj_pvalue, na.rm = TRUE)
  ), by = .(Chr, calder_bin_start)]

  agg[, any_sig := (n_sig_bins > 0L)]

  sig_dirs <- homer_dt[is_significant == TRUE,
                       .N, by = .(Chr, calder_bin_start, direction)]
  top_dir <- sig_dirs[sig_dirs[, .I[which.max(N)], by = .(Chr, calder_bin_start)]$V1]
  agg <- merge(agg, top_dir[, .(Chr, calder_bin_start,
                                 dominant_homer_direction = direction)],
               by = c("Chr", "calder_bin_start"), all.x = TRUE)

  agg[order(Chr, calder_bin_start)]
}

plot_sankey <- function(sig_dt, tp, tp_label, n_sig) {
  alluv_dt <- sig_dt[!is.na(ctrl_label) & !is.na(mut_label),
                      .N, by = .(direction, ctrl_label, mut_label)]
  alluv_dt[, direction_lbl := DIRECTION_LABELS[direction]]
  alluv_dt[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
  alluv_dt[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]

  alluv_dt[, flow_fill := fifelse(
    as.character(ctrl_label) == as.character(mut_label),
    as.character(ctrl_label), "change"
  )]

  p <- ggplot(alluv_dt,
              aes(axis1 = direction_lbl, axis2 = ctrl_label,
                  axis3 = mut_label, y = N)) +
    geom_alluvium(aes(fill = flow_fill),
                  width = 1/5, alpha = 0.75, knot.pos = 0.4) +
    geom_stratum(width = 1/5, fill = "white", color = "grey40", linewidth = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = LABEL_COLORS, name = "Flow") +
    scale_x_discrete(limits = c("HOMER call", "CALDER2 ctrl", "CALDER2 mut"),
                     expand = c(0.1, 0.1)) +
    labs(y = "Number of 25kb bins",
         title = sprintf("HOMER → CALDER2 Decomposition: %s", tp_label),
         subtitle = sprintf("Significant HOMER bins (FDR<%.2f, |Δ|>%.2f): n=%s",
                            HOMER_SIG_FDR, HOMER_SIG_DIFF,
                            format(n_sig, big.mark = ","))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank())
  p
}

plot_dotplot <- function(crosstab_dt, tp, tp_label) {
  ct <- copy(crosstab_dt)
  ct[, direction_lbl := DIRECTION_LABELS[direction]]
  ct[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
  ct[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]
  ct[, transition_category := factor(transition_category, levels = TRANSITION_LEVELS)]

  full_grid <- CJ(direction_lbl = unique(ct$direction_lbl),
                   ctrl_label = factor(LABEL_ORDER, levels = LABEL_ORDER),
                   mut_label  = factor(LABEL_ORDER, levels = LABEL_ORDER))
  ct <- merge(full_grid, ct, by = c("direction_lbl", "ctrl_label", "mut_label"),
              all.x = TRUE)

  p <- ggplot(ct[!is.na(N) & N > 0],
              aes(x = ctrl_label, y = mut_label, size = N,
                  color = transition_category)) +
    geom_point(alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey60", linewidth = 0.5) +
    facet_wrap(~ direction_lbl, nrow = 1) +
    scale_size_area(max_size = 14, name = "Bin count",
                    labels = label_comma()) +
    scale_color_manual(values = TRANSITION_COLORS, name = "Transition type",
                       drop = FALSE) +
    labs(x = "CALDER2 control label", y = "CALDER2 mutant label",
         title = sprintf("Subcompartment Transitions: %s", tp_label),
         subtitle = "Significant HOMER bins decomposed by subcompartment") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 12))
  p
}

# ── Main processing loop ──────────────────────────────────────────────────────

decomp_list <- list()

for (tp in TPS) {

  tp_start <- proc.time()
  tp_label <- sprintf("%s (%s)", tp, TP_LABELS[tp])
  CALDER2_DIR <- file.path(CALDER2_BASE, TP_DIRS[tp])
  OUT_DIR     <- file.path(OUT_BASE, TP_DIRS[tp])
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  cat("\n===========================================\n")
  cat(sprintf("  Timepoint: %s\n", tp_label))
  cat("===========================================\n")

  # --- Load CALDER2 labels ---
  cat("\n=== Load CALDER2 labels ===\n")
  calder_path <- file.path(CALDER2_DIR,
                           sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  calder <- fread(calder_path)
  calder[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
  calder[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]
  calder <- calder[chr %in% AUTOSOMES]
  n_calder <- nrow(calder)
  n_callable <- sum(!is.na(calder$ctrl_label) & !is.na(calder$mut_label))
  cat(sprintf("  CALDER2: %s 100kb bins, %s callable\n",
              format(n_calder, big.mark = ","),
              format(n_callable, big.mark = ",")))

  # --- Load HOMER ---
  cat("\n=== Load HOMER data ===\n")
  homer <- load_homer(HOMER_PATHS[tp], AUTOSOMES)
  n_total <- nrow(homer)
  n_sig   <- sum(homer$is_significant)
  cat(sprintf("  Significant bins: %s (FDR<%.2f, |Diff|>%.2f)\n",
              format(n_sig, big.mark = ","), HOMER_SIG_FDR, HOMER_SIG_DIFF))

  n_sig_a2b <- sum(homer$is_significant & homer$direction == "A_to_B_in_Mutant")
  n_sig_b2a <- sum(homer$is_significant & homer$direction == "B_to_A_in_Mutant")
  cat(sprintf("    A→B: %s,  B→A: %s\n",
              format(n_sig_a2b, big.mark = ","),
              format(n_sig_b2a, big.mark = ",")))

  # --- Aggregate HOMER to 100kb (for C4) ---
  cat("\n=== Aggregate HOMER to 100kb ===\n")
  homer_100kb <- aggregate_homer_100kb(homer)
  out_agg <- file.path(OUT_DIR, sprintf("%s_homer_100kb_aggregated.tsv", tp))
  fwrite(homer_100kb, out_agg, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s (%s rows)\n", basename(out_agg),
              format(nrow(homer_100kb), big.mark = ",")))

  # --- Join HOMER 25kb -> CALDER2 100kb ---
  cat("\n=== Join HOMER x CALDER2 ===\n")
  joined <- merge(
    homer[, .(Chr, Start, End, Difference, adj_pvalue,
              ctrl_avg_PC1, mut_avg_PC1, direction,
              is_significant, calder_bin_start)],
    calder[, .(chr, bin_start, bin_end, ctrl_label, mut_label,
               continous_rank_ctrl, continous_rank_mut)],
    by.x = c("Chr", "calder_bin_start"),
    by.y = c("chr", "bin_start"),
    all.x = TRUE
  )

  joined[, transition_category := classify_transition(ctrl_label, mut_label)]
  joined[, transition_category := factor(transition_category,
                                         levels = TRANSITION_LEVELS)]

  n_joined <- nrow(joined)
  n_with_label <- sum(!is.na(joined$ctrl_label))
  cat(sprintf("  Joined: %s rows (%s with CALDER2 labels, %s without)\n",
              format(n_joined, big.mark = ","),
              format(n_with_label, big.mark = ","),
              format(n_joined - n_with_label, big.mark = ",")))

  # --- Save full join ---
  out_joined <- file.path(OUT_DIR, sprintf("%s_homer_calder2_joined.tsv", tp))
  fwrite(joined[, .(Chr, Start, End, Difference, adj_pvalue,
                     ctrl_avg_PC1, mut_avg_PC1, direction,
                     is_significant, ctrl_label, mut_label,
                     transition_category)],
         out_joined, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s\n", basename(out_joined)))

  # --- Cross-tabulate significant bins ---
  cat("\n=== Cross-tabulate significant bins ===\n")
  sig_joined <- joined[is_significant == TRUE]
  n_sig_callable <- sum(!is.na(sig_joined$ctrl_label) &
                        !is.na(sig_joined$mut_label))
  n_sig_uncallable <- nrow(sig_joined) - n_sig_callable
  cat(sprintf("  Significant bins: %s total, %s callable, %s uncallable\n",
              format(nrow(sig_joined), big.mark = ","),
              format(n_sig_callable, big.mark = ","),
              format(n_sig_uncallable, big.mark = ",")))

  crosstab <- build_crosstab(sig_joined)
  out_crosstab <- file.path(OUT_DIR,
                            sprintf("%s_homer_calder2_crosstab.tsv", tp))
  fwrite(crosstab, out_crosstab, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written: %s (%d rows)\n", basename(out_crosstab),
              nrow(crosstab)))

  # --- Weakening decomposition ---
  cat("\n=== Weakening decomposition ===\n")
  decomp <- sig_joined[, .N, by = .(direction, transition_category)]
  decomp[, pct_of_direction := round(100 * N / sum(N), 2), by = direction]
  decomp[, tp := tp]
  decomp[, tp_label := TP_LABELS[tp]]

  for (dir in unique(decomp$direction)) {
    dir_lbl <- DIRECTION_LABELS[dir]
    dir_total <- sum(decomp[direction == dir, N])
    cat(sprintf("\n  %s (n=%s):\n", dir_lbl, format(dir_total, big.mark = ",")))
    sub <- decomp[direction == dir][order(-N)]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("    %-16s  %5s (%5.1f%%)\n",
                  sub$transition_category[i],
                  format(sub$N[i], big.mark = ","),
                  sub$pct_of_direction[i]))
    }
  }

  out_decomp <- file.path(OUT_DIR,
                           sprintf("%s_homer_weakening_decomposition.tsv", tp))
  fwrite(decomp, out_decomp, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("\n  Written: %s\n", basename(out_decomp)))
  decomp_list[[tp]] <- decomp

  # --- Figures ---
  cat("\n=== Generating figures ===\n")

  p_sankey <- plot_sankey(sig_joined, tp, tp_label, n_sig)
  save_multiformat_ggplot(p_sankey,
    file.path(OUT_DIR, sprintf("%s_homer_calder2_sankey", tp)),
    width = 10, height = 8)

  p_dot <- plot_dotplot(crosstab, tp, tp_label)
  save_multiformat_ggplot(p_dot,
    file.path(OUT_DIR, sprintf("%s_homer_calder2_dotplot", tp)),
    width = 10, height = 6)

  # --- Per-timepoint verification ---
  cat("\n=== Verification ===\n")
  any_warn_tp <- FALSE

  if (n_joined != n_total) {
    warning(sprintf("Row count mismatch: joined=%d, homer=%d", n_joined, n_total))
    any_warn_tp <- TRUE
  }

  if (sum(joined$is_significant) != n_sig) {
    warning(sprintf("Sig count mismatch: joined=%d, expected=%d",
                    sum(joined$is_significant), n_sig))
    any_warn_tp <- TRUE
  }

  pct_uncallable <- 100 * n_sig_uncallable / n_sig
  cat(sprintf("  Uncallable: %d of %d sig bins (%.1f%%)\n",
              n_sig_uncallable, n_sig, pct_uncallable))
  if (pct_uncallable > 30) {
    warning(sprintf("High uncallable fraction: %.1f%%", pct_uncallable))
    any_warn_tp <- TRUE
  }

  decomp_sum <- sum(decomp$N)
  if (decomp_sum != n_sig) {
    warning(sprintf("Decomposition sum=%d != n_sig=%d", decomp_sum, n_sig))
    any_warn_tp <- TRUE
  }

  n_true_flips <- sum(sig_joined$transition_category %in%
                      c("True_A_to_B", "True_B_to_A"), na.rm = TRUE)
  cat(sprintf("  True compartment flips: %s of %s sig bins (%.1f%%)\n",
              format(n_true_flips, big.mark = ","),
              format(n_sig, big.mark = ","),
              100 * n_true_flips / n_sig))
  if (n_true_flips == 0) {
    warning("No true compartment flips detected")
    any_warn_tp <- TRUE
  }

  expected_files <- c(
    sprintf("%s_homer_calder2_joined.tsv", tp),
    sprintf("%s_homer_100kb_aggregated.tsv", tp),
    sprintf("%s_homer_calder2_crosstab.tsv", tp),
    sprintf("%s_homer_weakening_decomposition.tsv", tp)
  )
  for (f in expected_files) {
    fpath <- file.path(OUT_DIR, f)
    if (!file.exists(fpath)) {
      warning(sprintf("Missing output: %s", fpath))
      any_warn_tp <- TRUE
    }
  }

  if (!any_warn_tp) cat("  All checks passed.\n")

  tp_elapsed <- (proc.time() - tp_start)[3]
  cat(sprintf("\n  Timepoint %s complete: %.1f min\n", tp, tp_elapsed / 60))
}

# ── Cross-timepoint combined summary ──────────────────────────────────────────

cat("\n=== Cross-timepoint summary ===\n")
combined <- rbindlist(decomp_list)
out_combined <- file.path(OUT_BASE, "combined_weakening_decomposition.tsv")
fwrite(combined, out_combined, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_combined), nrow(combined)))

# ── Final verification ─────────────────────────────────────────────────────────

cat("\n=== Final verification ===\n")
any_warn <- FALSE

for (tp in TPS) {
  tp_out <- file.path(OUT_BASE, TP_DIRS[tp])
  per_tp_expected <- c(
    sprintf("%s_homer_calder2_joined.tsv", tp),
    sprintf("%s_homer_100kb_aggregated.tsv", tp),
    sprintf("%s_homer_calder2_crosstab.tsv", tp),
    sprintf("%s_homer_weakening_decomposition.tsv", tp)
  )
  for (f in per_tp_expected) {
    fpath <- file.path(tp_out, f)
    if (!file.exists(fpath)) {
      warning(sprintf("Missing output: %s", f))
      any_warn <- TRUE
    } else if (file.info(fpath)$size == 0) {
      warning(sprintf("Empty output: %s", f))
      any_warn <- TRUE
    }
  }
}
combined_path <- file.path(OUT_BASE, "combined_weakening_decomposition.tsv")
if (!file.exists(combined_path)) {
  warning("Missing combined_weakening_decomposition.tsv")
  any_warn <- TRUE
}

if (!any_warn) cat("  All checks passed.\n")

total_elapsed <- (proc.time() - start_time)[3]
cat("\n===========================================\n")
cat("A4 COMPLETE\n")
cat(sprintf("Total runtime: %.1f min\n", total_elapsed / 60))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
