# ML/cmpts/scripts/C3_loops_stripes_integration.R
# Stage C3: Loop cluster and stripe x subcompartment integration.
#
# Maps loop anchors (Popay cluster pipeline, 38,948 loops x 6 clusters) and
# stripe anchors (Stripenn) to CALDER2 100kb subcompartment labels (A.1/A.2/
# B.1/B.2). Tests whether gained Polycomb loops (clust5) sit in B.1
# (facultative heterochromatin). Loops mapped against both timepoints.
#
# Usage:
#   Rscript C3_loops_stripes_integration.R <data_root> <code_root>
#     <data_root> : HPC data directory (kept for CLI consistency)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript C3_loops_stripes_integration.R <data_root> <code_root>")
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

CLUSTER_ORDER <- paste0("clust", 1:6)
CLUSTER_LABELS <- c(
  "clust1" = "Unchanged (high)",
  "clust2" = "Unchanged (mod)",
  "clust3" = "Moderate loss",
  "clust4" = "Moderate gain",
  "clust5" = "Strong gain",
  "clust6" = "Strong loss"
)
CLUSTER_COLORS <- c(
  "clust1" = "#1b9e77", "clust2" = "#d95f02", "clust3" = "#7570b3",
  "clust4" = "#e7298a", "clust5" = "#984ea3", "clust6" = "#a65628"
)

DIRECTION_ORDER <- c("gained", "lost", "unchanged")
DIRECTION_COLORS <- c(
  "gained"    = "#d01c8b",
  "lost"      = "#4dac26",
  "unchanged" = "#999999"
)

REPO_ROOT <- dirname(dirname(CODE_ROOT))
LOOP_PATH <- file.path(REPO_ROOT, "cluster", "outputs", "bap1_late",
                        "cluster3", "k-6", "data", "combined-clusters.txt")
STRIPE_PATHS <- c(
  "250402" = file.path(REPO_ROOT, "stripes", "stripenn", "outputs",
                        "250402", "cross_res_merged.tsv"),
  "250831" = file.path(REPO_ROOT, "stripes", "stripenn", "outputs",
                        "250831", "cross_res_merged.tsv")
)

CALDER2_BASE <- file.path(CODE_ROOT, "outputs", "calder2")
OUT_DIR      <- file.path(CODE_ROOT, "outputs", "integration", "loops_stripes")
TSV_DIR      <- file.path(OUT_DIR, "tsv")
PLOT_DIR     <- file.path(OUT_DIR, "plots")
UTIL_PATH    <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

# -- Libraries ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggalluvial)
})
source(UTIL_PATH)

# -- Helper functions ---------------------------------------------------------

map_to_bin_0based <- function(pos) {
  floor(pos / BIN_SIZE) * BIN_SIZE + 1L
}

map_to_bin_1based <- function(pos) {
  floor((pos - 1) / BIN_SIZE) * BIN_SIZE + 1L
}

compute_enrichment <- function(dt, group_col, label_col) {
  dt <- copy(dt)
  total <- sum(dt$N)
  dt[, group_total := sum(N), by = group_col]
  dt[, label_total := sum(N), by = label_col]
  dt[, expected := (group_total * label_total) / total]
  dt[, obs_exp := fifelse(expected > 0, N / expected, NA_real_)]
  dt[, log2_obs_exp := fifelse(!is.na(obs_exp) & obs_exp > 0,
                                log2(obs_exp), NA_real_)]
  dt
}

fisher_enrichment <- function(anchor_dt, group_val, label_val,
                               group_col = "GROUP", label_col = "ctrl_label") {
  dt <- anchor_dt[!is.na(get(label_col))]
  in_group  <- dt[[group_col]] == group_val
  in_label  <- dt[[label_col]] == label_val
  mat <- matrix(c(
    sum( in_group &  in_label),
    sum( in_group & !in_label),
    sum(!in_group &  in_label),
    sum(!in_group & !in_label)
  ), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat)
  data.table(
    group = group_val, label = label_val, condition = label_col,
    a = mat[1, 1], b = mat[1, 2], c = mat[2, 1], d = mat[2, 2],
    odds_ratio = as.numeric(ft$estimate),
    p_value = ft$p.value,
    ci_low = ft$conf.int[1], ci_high = ft$conf.int[2]
  )
}

classify_transition <- function(ctrl, mut) {
  from_A <- grepl("^A", ctrl)
  to_A   <- grepl("^A", mut)
  same   <- (ctrl == mut)
  fifelse(same, "Stable",
    fifelse(from_A & !to_A, "A_to_B",
      fifelse(!from_A & to_A, "B_to_A",
        fifelse(from_A & to_A, "Within_A", "Within_B"))))
}

# -- Header -------------------------------------------------------------------

start_time <- proc.time()

cat("===========================================\n")
cat("C3: Loop/Stripe x Subcompartment Integration\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("REPO_ROOT:  %s\n", REPO_ROOT))
cat(sprintf("Output dir: %s\n", OUT_DIR))
cat(sprintf("Timepoints: %s\n", paste(TPS, collapse = ", ")))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# -- Pre-flight validation ----------------------------------------------------

cat("=== Pre-flight validation ===\n")
calder_paths <- list()
for (tp in TPS) {
  path <- file.path(CALDER2_BASE, TP_DIRS[tp],
                    sprintf("%s_subcompartment_labels_100kb.tsv", tp))
  if (!file.exists(path)) stop(sprintf("Missing CALDER2 output: %s", path))
  if (file.info(path)$size == 0) stop(sprintf("Empty CALDER2 output: %s", path))
  calder_paths[[tp]] <- path
  cat(sprintf("  OK: %s (%s bytes)\n", basename(path),
              format(file.info(path)$size, big.mark = ",")))
}

if (!file.exists(LOOP_PATH)) stop(sprintf("Missing loop clusters: %s", LOOP_PATH))
cat(sprintf("  OK: %s (%s bytes)\n", basename(LOOP_PATH),
            format(file.info(LOOP_PATH)$size, big.mark = ",")))

for (tp in TPS) {
  sp <- STRIPE_PATHS[tp]
  if (!file.exists(sp)) stop(sprintf("Missing stripe data: %s", sp))
  cat(sprintf("  OK: %s (%s bytes)\n", basename(sp),
              format(file.info(sp)$size, big.mark = ",")))
}

if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
dir.create(TSV_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("\n")

# =============================================================================
# Phase 1: Load data
# =============================================================================

cat("=== Phase 1: Loading data ===\n")

labels <- list()
for (tp in TPS) {
  labels[[tp]] <- fread(calder_paths[[tp]])
  n_callable <- labels[[tp]][!is.na(ctrl_label) & !is.na(mut_label), .N]
  cat(sprintf("  CALDER2 %s (%s): %d bins, %d callable\n",
              tp, TP_LABELS[tp], nrow(labels[[tp]]), n_callable))
}

loops <- fread(LOOP_PATH)
cat(sprintf("  Loops: %d rows, %d clusters\n", nrow(loops), uniqueN(loops$GROUP)))
stopifnot(all(loops$chr1 == loops$chr2))

stripes <- list()
for (tp in TPS) {
  stripes[[tp]] <- fread(STRIPE_PATHS[tp])
  cat(sprintf("  Stripes %s: %d rows\n", tp, nrow(stripes[[tp]])))
}
cat("\n")

# =============================================================================
# Phase 2: Loop anchor x subcompartment mapping (both timepoints)
# =============================================================================

cat("=== Phase 2: Loop anchor x subcompartment mapping ===\n")

loops[, anchor1_mid := (x1 + x2) / 2]
loops[, anchor2_mid := (y1 + y2) / 2]
loops[, anchor1_bin := map_to_bin_0based(anchor1_mid)]
loops[, anchor2_bin := map_to_bin_0based(anchor2_mid)]

loop_enrich_all <- list()
loop_crosstab_all <- list()
anchor_labeled_all <- list()
fisher_results <- list()

for (tp in TPS) {
  cat(sprintf("\n  --- %s (%s) ---\n", tp, TP_LABELS[tp]))

  cdt <- labels[[tp]][, .(chr, bin_start, ctrl_label, mut_label)]

  tmp <- copy(loops)
  tmp <- merge(tmp, cdt, by.x = c("chr1", "anchor1_bin"),
               by.y = c("chr", "bin_start"), all.x = TRUE)
  setnames(tmp, c("ctrl_label", "mut_label"), c("a1_ctrl", "a1_mut"))
  tmp <- merge(tmp, cdt, by.x = c("chr2", "anchor2_bin"),
               by.y = c("chr", "bin_start"), all.x = TRUE)
  setnames(tmp, c("ctrl_label", "mut_label"), c("a2_ctrl", "a2_mut"))

  anchor_long <- rbind(
    tmp[, .(GROUP, chr = chr1, mid = anchor1_mid, bin_start = anchor1_bin,
            ctrl_label = a1_ctrl, mut_label = a1_mut, anchor_side = "anchor1")],
    tmp[, .(GROUP, chr = chr2, mid = anchor2_mid, bin_start = anchor2_bin,
            ctrl_label = a2_ctrl, mut_label = a2_mut, anchor_side = "anchor2")]
  )

  n_total  <- nrow(anchor_long)
  callable <- anchor_long[!is.na(ctrl_label)]
  n_call   <- nrow(callable)
  cat(sprintf("  Total anchors: %d, callable: %d (%.1f%%)\n",
              n_total, n_call, 100 * n_call / n_total))

  # Cross-tabulation and enrichment (ctrl)
  ctrl_counts <- callable[, .N, by = .(GROUP, ctrl_label)]
  ctrl_enrich <- compute_enrichment(ctrl_counts, "GROUP", "ctrl_label")

  ctrl_tab <- dcast(ctrl_counts, GROUP ~ ctrl_label, value.var = "N", fill = 0L)
  ctrl_mat <- as.matrix(ctrl_tab[, LABEL_ORDER, with = FALSE])
  rownames(ctrl_mat) <- ctrl_tab$GROUP
  ctrl_chisq <- chisq.test(ctrl_mat)
  cat(sprintf("  Ctrl chi-squared: X2=%.1f, p=%.2e, df=%d\n",
              ctrl_chisq$statistic, ctrl_chisq$p.value, ctrl_chisq$parameter))

  # Cross-tabulation and enrichment (mut)
  mut_counts <- callable[, .N, by = .(GROUP, mut_label)]
  mut_enrich <- compute_enrichment(mut_counts, "GROUP", "mut_label")

  mut_tab <- dcast(mut_counts, GROUP ~ mut_label, value.var = "N", fill = 0L)
  mut_mat <- as.matrix(mut_tab[, LABEL_ORDER, with = FALSE])
  rownames(mut_mat) <- mut_tab$GROUP
  mut_chisq <- chisq.test(mut_mat)
  cat(sprintf("  Mut  chi-squared: X2=%.1f, p=%.2e, df=%d\n",
              mut_chisq$statistic, mut_chisq$p.value, mut_chisq$parameter))

  # Fisher's exact for key tests
  ft_list <- list()
  for (cond in c("ctrl_label", "mut_label")) {
    for (cl in c("clust5", "clust6")) {
      for (lb in LABEL_ORDER) {
        ft_list[[length(ft_list) + 1]] <- fisher_enrichment(
          callable, cl, lb, "GROUP", cond)
      }
    }
  }
  ft_dt <- rbindlist(ft_list)
  ft_dt[, timepoint := tp]
  fisher_results[[tp]] <- ft_dt

  # Report key tests
  key <- ft_dt[group == "clust5" & label == "B.1" & condition == "ctrl_label"]
  cat(sprintf("  KEY: clust5 in B.1 (ctrl): OR=%.2f, p=%.2e [%.2f, %.2f]\n",
              key$odds_ratio, key$p_value, key$ci_low, key$ci_high))
  key_mut <- ft_dt[group == "clust5" & label == "B.1" & condition == "mut_label"]
  cat(sprintf("  KEY: clust5 in B.1 (mut):  OR=%.2f, p=%.2e [%.2f, %.2f]\n",
              key_mut$odds_ratio, key_mut$p_value, key_mut$ci_low, key_mut$ci_high))
  key_c6 <- ft_dt[group == "clust6" & label == "A.1" & condition == "ctrl_label"]
  cat(sprintf("  clust6 in A.1 (ctrl): OR=%.2f, p=%.2e\n",
              key_c6$odds_ratio, key_c6$p_value))

  # Per-cluster subcompartment distribution
  cat("  Per-cluster ctrl distribution (%):\n")
  for (cl in CLUSTER_ORDER) {
    cl_n <- callable[GROUP == cl, .N]
    pcts <- sapply(LABEL_ORDER, function(lb) {
      round(100 * callable[GROUP == cl & ctrl_label == lb, .N] / cl_n, 1)
    })
    cat(sprintf("    %s (n=%d): A.1=%.1f  A.2=%.1f  B.1=%.1f  B.2=%.1f\n",
                cl, cl_n, pcts[1], pcts[2], pcts[3], pcts[4]))
  }

  # Store
  ctrl_enrich[, timepoint := tp]
  ctrl_enrich[, condition := "ctrl"]
  mut_enrich[, timepoint := tp]
  mut_enrich[, condition := "mut"]
  loop_enrich_all[[paste0(tp, "_ctrl")]] <- ctrl_enrich
  loop_enrich_all[[paste0(tp, "_mut")]]  <- mut_enrich

  loop_crosstab_all[[paste0(tp, "_ctrl")]] <- ctrl_tab
  loop_crosstab_all[[paste0(tp, "_mut")]]  <- mut_tab

  callable[, timepoint := tp]
  anchor_labeled_all[[tp]] <- callable
}

cat("\n")

# =============================================================================
# Phase 3: Transition analysis at loop anchors
# =============================================================================

cat("=== Phase 3: Transition analysis at loop anchors ===\n")

trans_results <- list()
for (tp in TPS) {
  a <- anchor_labeled_all[[tp]]
  a_both <- a[!is.na(ctrl_label) & !is.na(mut_label)]
  a_both[, transition := classify_transition(
    as.character(ctrl_label), as.character(mut_label))]

  trans_summ <- a_both[, .N, by = .(GROUP, transition)]
  trans_summ[, pct := round(100 * N / sum(N), 1), by = GROUP]
  trans_summ[, timepoint := tp]
  trans_results[[tp]] <- trans_summ

  cat(sprintf("\n  %s (%s):\n", tp, TP_LABELS[tp]))
  for (cl in c("clust5", "clust6")) {
    cat(sprintf("    %s:\n", cl))
    cl_dt <- trans_summ[GROUP == cl][order(-N)]
    for (i in seq_len(nrow(cl_dt))) {
      cat(sprintf("      %s: %d (%.1f%%)\n",
                  cl_dt$transition[i], cl_dt$N[i], cl_dt$pct[i]))
    }
  }

  anchor_labeled_all[[tp]] <- a_both
}

cat("\n")

# =============================================================================
# Phase 4: Stripe x subcompartment mapping
# =============================================================================

cat("=== Phase 4: Stripe x subcompartment mapping ===\n")

stripe_enrich_all <- list()
stripe_body_all   <- list()
stripe_fisher_all <- list()

for (tp in TPS) {
  cat(sprintf("\n  --- %s (%s) ---\n", tp, TP_LABELS[tp]))
  sdt <- copy(stripes[[tp]])

  sdt[, anchor_bin := map_to_bin_1based(anchor_center)]
  sdt <- merge(sdt, labels[[tp]][, .(chr, bin_start, ctrl_label, mut_label)],
               by.x = c("chr", "anchor_bin"), by.y = c("chr", "bin_start"),
               all.x = TRUE)

  sdt[direction %in% c("weakened", "strengthened"), direction := "unchanged"]
  sdt[, direction := factor(direction, levels = DIRECTION_ORDER)]

  s_callable <- sdt[!is.na(ctrl_label)]
  cat(sprintf("  Callable: %d / %d (%.1f%%)\n",
              nrow(s_callable), nrow(sdt), 100 * nrow(s_callable) / nrow(sdt)))

  # Enrichment
  s_counts <- s_callable[, .N, by = .(direction, ctrl_label)]
  s_enrich <- compute_enrichment(s_counts, "direction", "ctrl_label")
  s_enrich[, timepoint := tp]
  stripe_enrich_all[[tp]] <- s_enrich

  s_tab <- dcast(s_counts, direction ~ ctrl_label, value.var = "N", fill = 0L)
  s_mat <- as.matrix(s_tab[, intersect(LABEL_ORDER, names(s_tab)), with = FALSE])
  rownames(s_mat) <- as.character(s_tab$direction)
  s_chisq <- suppressWarnings(chisq.test(s_mat))
  cat(sprintf("  Chi-squared: X2=%.1f, p=%.2e\n",
              s_chisq$statistic, s_chisq$p.value))

  # Fisher's: gained in B.1
  ft_g <- fisher_enrichment(s_callable, "gained", "B.1", "direction", "ctrl_label")
  ft_g[, timepoint := tp]
  stripe_fisher_all[[tp]] <- ft_g
  cat(sprintf("  Gained in B.1: OR=%.2f, p=%.2e\n", ft_g$odds_ratio, ft_g$p_value))

  # Stripe body composition
  body_sdt <- s_callable[!is.na(pos3) & !is.na(pos4) & pos4 > pos3]
  if (nrow(body_sdt) > 0) {
    calder_key <- labels[[tp]][, .(chr, bin_start, body_label = ctrl_label)]

    body_list <- lapply(seq_len(nrow(body_sdt)), function(i) {
      row <- body_sdt[i]
      from <- map_to_bin_1based(row$pos3)
      to   <- map_to_bin_1based(row$pos4)
      if (from > to) return(NULL)
      bins <- seq(from, to, by = BIN_SIZE)
      data.table(stripe_id = row$stripe_id,
                 direction = as.character(row$direction),
                 chr = row$chr, bin_start = bins)
    })
    body_bins_dt <- rbindlist(body_list[!vapply(body_list, is.null, logical(1))])

    body_bins_dt <- merge(body_bins_dt, calder_key,
                          by = c("chr", "bin_start"), all.x = TRUE)

    body_summary <- body_bins_dt[!is.na(body_label), .(
      body_bins = .N,
      frac_A1 = round(sum(body_label == "A.1") / .N, 4),
      frac_A2 = round(sum(body_label == "A.2") / .N, 4),
      frac_B1 = round(sum(body_label == "B.1") / .N, 4),
      frac_B2 = round(sum(body_label == "B.2") / .N, 4)
    ), by = .(stripe_id, direction)]
    body_summary[, timepoint := tp]
    stripe_body_all[[tp]] <- body_summary

    cat(sprintf("  Body composition: %d stripes analyzed, median bins/stripe: %d\n",
                nrow(body_summary), median(body_summary$body_bins)))

    for (dir in DIRECTION_ORDER) {
      dir_body <- body_summary[as.character(direction) == dir]
      if (nrow(dir_body) > 0) {
        cat(sprintf("    %s (n=%d): mean B.1=%.1f%%, B.2=%.1f%%\n",
                    dir, nrow(dir_body),
                    100 * mean(dir_body$frac_B1),
                    100 * mean(dir_body$frac_B2)))
      }
    }
  }
}

cat("\n")

# =============================================================================
# Phase 5: Write TSVs
# =============================================================================

cat("=== Phase 5: Writing TSV outputs ===\n")

n_tsvs <- 0L

for (tp in TPS) {
  # Per-anchor labels
  fwrite(anchor_labeled_all[[tp]],
         file.path(TSV_DIR, sprintf("%s_loop_anchor_labels.tsv", tp)),
         sep = "\t", quote = FALSE, na = "NA")
  n_tsvs <- n_tsvs + 1L

  # Crosstab ctrl
  fwrite(loop_crosstab_all[[paste0(tp, "_ctrl")]],
         file.path(TSV_DIR, sprintf("%s_loop_subcompartment_crosstab_ctrl.tsv", tp)),
         sep = "\t", quote = FALSE)
  n_tsvs <- n_tsvs + 1L

  # Crosstab mut
  fwrite(loop_crosstab_all[[paste0(tp, "_mut")]],
         file.path(TSV_DIR, sprintf("%s_loop_subcompartment_crosstab_mut.tsv", tp)),
         sep = "\t", quote = FALSE)
  n_tsvs <- n_tsvs + 1L

  # Transition summary
  fwrite(trans_results[[tp]],
         file.path(TSV_DIR, sprintf("%s_loop_transition_summary.tsv", tp)),
         sep = "\t", quote = FALSE)
  n_tsvs <- n_tsvs + 1L

  # Stripe enrichment
  fwrite(stripe_enrich_all[[tp]],
         file.path(TSV_DIR, sprintf("%s_stripe_subcompartment_enrichment.tsv", tp)),
         sep = "\t", quote = FALSE, na = "NA")
  n_tsvs <- n_tsvs + 1L

  # Stripe body
  if (!is.null(stripe_body_all[[tp]])) {
    fwrite(stripe_body_all[[tp]],
           file.path(TSV_DIR, sprintf("%s_stripe_body_composition.tsv", tp)),
           sep = "\t", quote = FALSE)
    n_tsvs <- n_tsvs + 1L
  }
}

# Combined enrichment (all timepoints, both conditions)
all_enrich <- rbindlist(loop_enrich_all, fill = TRUE)
fwrite(all_enrich,
       file.path(TSV_DIR, "loop_subcompartment_enrichment.tsv"),
       sep = "\t", quote = FALSE, na = "NA")
n_tsvs <- n_tsvs + 1L

# Combined Fisher's tests
all_fisher <- rbindlist(fisher_results, fill = TRUE)
fwrite(all_fisher,
       file.path(TSV_DIR, "loop_enrichment_tests.tsv"),
       sep = "\t", quote = FALSE)
n_tsvs <- n_tsvs + 1L

# Combined summary
summary_rows <- list()
for (tp in TPS) {
  ft <- fisher_results[[tp]]
  c5b1_ctrl <- ft[group == "clust5" & label == "B.1" & condition == "ctrl_label"]
  c5b1_mut  <- ft[group == "clust5" & label == "B.1" & condition == "mut_label"]
  c6a1_ctrl <- ft[group == "clust6" & label == "A.1" & condition == "ctrl_label"]
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    timepoint = tp, analysis = "loop",
    test = "clust5_B.1_ctrl", OR = c5b1_ctrl$odds_ratio, p = c5b1_ctrl$p_value)
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    timepoint = tp, analysis = "loop",
    test = "clust5_B.1_mut", OR = c5b1_mut$odds_ratio, p = c5b1_mut$p_value)
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    timepoint = tp, analysis = "loop",
    test = "clust6_A.1_ctrl", OR = c6a1_ctrl$odds_ratio, p = c6a1_ctrl$p_value)
  s_ft <- stripe_fisher_all[[tp]]
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    timepoint = tp, analysis = "stripe",
    test = "gained_B.1_ctrl", OR = s_ft$odds_ratio, p = s_ft$p_value)
}
fwrite(rbindlist(summary_rows),
       file.path(TSV_DIR, "combined_summary.tsv"),
       sep = "\t", quote = FALSE)
n_tsvs <- n_tsvs + 1L

cat(sprintf("  Wrote %d TSV files\n\n", n_tsvs))

# =============================================================================
# Phase 6: Figures
# =============================================================================

cat("=== Phase 6: Generating figures ===\n")
n_figs <- 0L

# --- Figure 1: Balloon plot (loop cluster x subcompartment enrichment) -------

cat("  Fig 1: Loop subcompartment enrichment balloon plot\n")

enrich_ctrl <- rbindlist(loop_enrich_all[grepl("_ctrl$", names(loop_enrich_all))])
enrich_ctrl[, GROUP := factor(GROUP, levels = rev(CLUSTER_ORDER))]
enrich_ctrl[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
enrich_ctrl[, tp_label := TP_LABELS[timepoint]]
enrich_ctrl[, tp_label := factor(tp_label, levels = TP_LABELS)]

p_balloon <- ggplot(enrich_ctrl, aes(x = ctrl_label, y = GROUP)) +
  geom_point(aes(size = N, fill = log2_obs_exp),
             shape = 21, color = "grey30", stroke = 0.3) +
  facet_wrap(~ tp_label, ncol = 2) +
  scale_size_continuous(name = "Anchor count", range = c(2, 14),
                        breaks = pretty_breaks(4)) +
  scale_fill_gradient2(name = expression(log[2]~"(obs/exp)"),
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1.5, 1.5), oob = squish) +
  scale_y_discrete(labels = function(x) {
    paste0(x, "\n", CLUSTER_LABELS[x])
  }) +
  labs(x = "Subcompartment (ctrl)", y = NULL,
       title = "Loop Cluster x Subcompartment Enrichment") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

save_multiformat_ggplot(p_balloon,
  file.path(PLOT_DIR, "c3_loop_subcompartment_enrichment"),
  width = 10, height = 8)
n_figs <- n_figs + 1L

# --- Figure 2: Stacked bar (subcompartment % per cluster, ctrl vs mut) -------

cat("  Fig 2: Subcompartment stacked bar per cluster\n")

a_late <- anchor_labeled_all[["250402"]]
bar_data <- rbind(
  a_late[!is.na(ctrl_label), .(GROUP, label = ctrl_label, condition = "Ctrl")],
  a_late[!is.na(mut_label),  .(GROUP, label = mut_label,  condition = "Mut")]
)
bar_data[, label := factor(label, levels = LABEL_ORDER)]
bar_data[, GROUP := factor(GROUP, levels = CLUSTER_ORDER)]
bar_data[, condition := factor(condition, levels = c("Ctrl", "Mut"))]
bar_summary <- bar_data[, .N, by = .(GROUP, label, condition)]
bar_summary[, pct := 100 * N / sum(N), by = .(GROUP, condition)]

p_stacked <- ggplot(bar_summary, aes(x = condition, y = pct, fill = label)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~ GROUP, nrow = 1, labeller = labeller(GROUP = CLUSTER_LABELS)) +
  scale_fill_manual(values = LABEL_COLORS, name = "Subcompartment") +
  labs(x = NULL, y = "% of anchors",
       title = "Subcompartment Distribution per Cluster (Late/Adult)",
       subtitle = "Ctrl vs Mut conditions") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 9))

save_multiformat_ggplot(p_stacked,
  file.path(PLOT_DIR, "c3_loop_subcompartment_stacked"),
  width = 14, height = 6)
n_figs <- n_figs + 1L

# --- Figure 3: Heatmap (obs/exp with values) --------------------------------

cat("  Fig 3: Enrichment heatmap\n")

heat_data <- rbindlist(loop_enrich_all[grepl("_ctrl$", names(loop_enrich_all))])
heat_data[, GROUP := factor(GROUP, levels = rev(CLUSTER_ORDER))]
heat_data[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
heat_data[, tp_label := factor(TP_LABELS[timepoint], levels = TP_LABELS)]

p_heatmap <- ggplot(heat_data, aes(x = ctrl_label, y = GROUP)) +
  geom_tile(aes(fill = log2_obs_exp), color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", obs_exp)), size = 3.2, color = "black") +
  facet_wrap(~ tp_label, ncol = 2) +
  scale_fill_gradient2(name = expression(log[2]~"(obs/exp)"),
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0) +
  scale_y_discrete(labels = function(x) {
    paste0(x, " ", CLUSTER_LABELS[x])
  }) +
  labs(x = "Subcompartment (ctrl)", y = NULL,
       title = "Loop Anchor Enrichment in Subcompartments") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

save_multiformat_ggplot(p_heatmap,
  file.path(PLOT_DIR, "c3_loop_subcompartment_heatmap"),
  width = 10, height = 8)
n_figs <- n_figs + 1L

# --- Figure 4: Stripe enrichment dot plot ------------------------------------

cat("  Fig 4: Stripe subcompartment enrichment\n")

s_enrich <- rbindlist(stripe_enrich_all)
s_enrich[, direction := factor(direction, levels = rev(DIRECTION_ORDER))]
s_enrich[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
s_enrich[, tp_label := factor(TP_LABELS[timepoint], levels = TP_LABELS)]

p_stripe <- ggplot(s_enrich, aes(x = ctrl_label, y = direction)) +
  geom_point(aes(size = N, fill = log2_obs_exp),
             shape = 21, color = "grey30", stroke = 0.3) +
  facet_wrap(~ tp_label, ncol = 2) +
  scale_size_continuous(name = "Stripe count", range = c(2, 12),
                        breaks = pretty_breaks(3)) +
  scale_fill_gradient2(name = expression(log[2]~"(obs/exp)"),
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.5, 0.5), oob = squish) +
  labs(x = "Subcompartment (ctrl)", y = "Stripe direction",
       title = "Stripe Anchor x Subcompartment Enrichment") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

save_multiformat_ggplot(p_stripe,
  file.path(PLOT_DIR, "c3_stripe_subcompartment_enrichment"),
  width = 10, height = 5)
n_figs <- n_figs + 1L

# --- Figure 5: Alluvial (clust5/clust6 transitions, late) -------------------

cat("  Fig 5: clust5/clust6 transition alluvial\n")

a_late_both <- anchor_labeled_all[["250402"]][
  GROUP %in% c("clust5", "clust6") & !is.na(ctrl_label) & !is.na(mut_label)]

alluv_dt <- a_late_both[, .N, by = .(GROUP, ctrl_label, mut_label)]
alluv_dt[, flow_type := fifelse(as.character(ctrl_label) == as.character(mut_label),
                                 as.character(ctrl_label), "change")]
alluv_dt[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
alluv_dt[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]
alluv_dt[, GROUP := factor(GROUP, levels = c("clust5", "clust6"),
                            labels = c("clust5 (Strong gain)", "clust6 (Strong loss)"))]

p_alluv <- ggplot(alluv_dt, aes(axis1 = ctrl_label, axis2 = mut_label, y = N)) +
  geom_alluvium(aes(fill = flow_type), width = 1/4, alpha = 0.75, knot.pos = 0.4) +
  geom_stratum(width = 1/4, fill = "white", color = "grey40", linewidth = 0.5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5, fontface = "bold") +
  facet_wrap(~ GROUP, scales = "free_y") +
  scale_fill_manual(values = LABEL_COLORS, name = "Flow") +
  scale_x_discrete(limits = c("Ctrl", "Mut"), expand = c(0.05, 0.05)) +
  labs(y = "Number of loop anchors",
       title = "Subcompartment Transitions at Differential Loop Anchors (Late)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11))

save_multiformat_ggplot(p_alluv,
  file.path(PLOT_DIR, "c3_clust5_clust6_transitions"),
  width = 12, height = 7)
n_figs <- n_figs + 1L

# --- Figure 6: Developmental comparison (clust5 early vs late) ---------------

cat("  Fig 6: Developmental comparison for clust5\n")

dev_data <- list()
for (tp in TPS) {
  a <- anchor_labeled_all[[tp]]
  c5 <- a[GROUP == "clust5" & !is.na(ctrl_label)]
  dist <- c5[, .N, by = ctrl_label]
  dist[, pct := 100 * N / sum(N)]
  dist[, timepoint := tp]
  dist[, tp_label := TP_LABELS[tp]]
  dev_data[[tp]] <- dist
}
dev_dt <- rbindlist(dev_data)
dev_dt[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
dev_dt[, tp_label := factor(tp_label, levels = TP_LABELS)]

p_dev <- ggplot(dev_dt, aes(x = tp_label, y = pct, fill = ctrl_label)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = LABEL_COLORS, name = "Subcompartment") +
  labs(x = "Timepoint", y = "% of clust5 anchors (ctrl)",
       title = "Clust5 (Gained Polycomb Loops): Subcompartment at Early vs Late",
       subtitle = "Same loop anchors mapped against each timepoint's CALDER2 labels") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

save_multiformat_ggplot(p_dev,
  file.path(PLOT_DIR, "c3_developmental_comparison"),
  width = 8, height = 6)
n_figs <- n_figs + 1L

cat(sprintf("  Generated %d figure sets\n\n", n_figs))

# -- Verification -------------------------------------------------------------

cat("=== Verification ===\n")
any_warn <- FALSE

n_anchors <- nrow(loops) * 2L
cat(sprintf("  Expected anchors: %d\n", n_anchors))

for (tp in TPS) {
  a <- anchor_labeled_all[[tp]]
  callable_frac <- nrow(a) / n_anchors
  cat(sprintf("  %s callable: %d / %d (%.1f%%)\n",
              tp, nrow(a), n_anchors, 100 * callable_frac))
  if (callable_frac < 0.85) {
    warning(sprintf("Low callable fraction for %s: %.1f%%", tp, 100 * callable_frac))
    any_warn <- TRUE
  }
}

n_c5 <- anchor_labeled_all[["250402"]][GROUP == "clust5", .N]
n_c6 <- anchor_labeled_all[["250402"]][GROUP == "clust6", .N]
cat(sprintf("  clust5 anchors (late): %d, clust6: %d\n", n_c5, n_c6))

ft_late <- fisher_results[["250402"]]
key_late <- ft_late[group == "clust5" & label == "B.1" & condition == "ctrl_label"]
if (key_late$p_value > 0.05) {
  cat("  [NOTE] clust5 NOT significantly enriched in B.1 (ctrl, late)\n")
} else {
  cat(sprintf("  [PASS] clust5 enriched in B.1 (ctrl, late): OR=%.2f, p=%.2e\n",
              key_late$odds_ratio, key_late$p_value))
}

expected_tsvs <- c("loop_subcompartment_enrichment.tsv",
                    "loop_enrichment_tests.tsv",
                    "combined_summary.tsv")
for (tp in TPS) {
  expected_tsvs <- c(expected_tsvs,
    sprintf("%s_loop_anchor_labels.tsv", tp),
    sprintf("%s_loop_subcompartment_crosstab_ctrl.tsv", tp),
    sprintf("%s_loop_subcompartment_crosstab_mut.tsv", tp),
    sprintf("%s_loop_transition_summary.tsv", tp),
    sprintf("%s_stripe_subcompartment_enrichment.tsv", tp))
}

for (f in expected_tsvs) {
  fpath <- file.path(TSV_DIR, f)
  if (!file.exists(fpath) || file.info(fpath)$size == 0) {
    warning(sprintf("Missing or empty output: %s", f))
    any_warn <- TRUE
  }
}

expected_figs <- c("c3_loop_subcompartment_enrichment",
                    "c3_loop_subcompartment_stacked",
                    "c3_loop_subcompartment_heatmap",
                    "c3_stripe_subcompartment_enrichment",
                    "c3_clust5_clust6_transitions",
                    "c3_developmental_comparison")
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
cat(sprintf("C3 COMPLETE\n"))
cat(sprintf("Runtime: %.1f seconds\n", elapsed))
cat(sprintf("Outputs: %d TSVs + %d figure sets\n", n_tsvs, n_figs))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
