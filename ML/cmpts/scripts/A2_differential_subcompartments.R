# ML/cmpts/scripts/A2_differential_subcompartments.R
# Stage A2: Differential subcompartment analysis for a single timepoint.
#
# Bins variable-width CALDER2 segments to uniform 100kb resolution via
# plurality vote, then computes transition matrices and chi-squared tests
# comparing ctrl vs mut subcompartment assignments.
#
# Usage:
#   Rscript --vanilla A2_differential_subcompartments.R <timepoint> <data_root> <code_root>
#     <timepoint> : 250402 | 250831
#     <data_root> : HPC data directory (kept for CLI consistency with A1)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript A2_differential_subcompartments.R <timepoint> <data_root> <code_root>")
}

TP        <- args[1]
DATA_ROOT <- args[2]
CODE_ROOT <- args[3]

# ── Constants ──────────────────────────────────────────────────────────────────

SAMPLES     <- c("ctrl_merged", "mut_merged")
BIN_SIZE    <- 100000L
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

TP_LABELS <- c("250402" = "late/adult", "250831" = "early/P12")

OUT_DIR   <- file.path(CODE_ROOT, "outputs", "calder2")
UTIL_PATH <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

# ── Libraries ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
})
source(UTIL_PATH)

# ── Header ─────────────────────────────────────────────────────────────────────

cat("===========================================\n")
cat("A2: Differential Subcompartment Analysis\n")
cat("===========================================\n")
cat(sprintf("Timepoint:  %s (%s)\n", TP, TP_LABELS[TP]))
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("Output dir: %s\n", OUT_DIR))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# ── Pre-flight validation ──────────────────────────────────────────────────────

cat("=== Pre-flight validation ===\n")
input_files <- list()
for (samp in SAMPLES) {
  path <- file.path(CODE_ROOT, "outputs", "calder2", TP, samp,
                    "sub_compartments", "all_sub_compartments.tsv")
  if (!file.exists(path)) stop(sprintf("Missing A1 output: %s", path))
  if (file.info(path)$size == 0) stop(sprintf("Empty A1 output: %s", path))
  input_files[[samp]] <- path
  cat(sprintf("  OK: %s/%s (%s bytes)\n", TP, samp,
              format(file.info(path)$size, big.mark = ",")))
}

if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Function definitions ───────────────────────────────────────────────────────

truncate_to_depth2 <- function(comp_name_vec) {
  sub("^([AB]\\.[12]).*", "\\1", comp_name_vec)
}

build_100kb_grid <- function(chrom_sizes) {
  grid_list <- lapply(names(chrom_sizes), function(chr) {
    sz     <- chrom_sizes[[chr]]
    n_bins <- ceiling(sz / BIN_SIZE)
    data.table(
      chr       = chr,
      bin_index = seq_len(n_bins),
      bin_start = (seq_len(n_bins) - 1L) * BIN_SIZE + 1L,
      bin_end   = pmin(seq_len(n_bins) * BIN_SIZE, sz)
    )
  })
  rbindlist(grid_list)
}

bin_plurality_vote <- function(comp_dt, grid_dt) {
  dt <- copy(comp_dt)

  dt[, bin_first := ceiling(pos_start / BIN_SIZE)]
  dt[, bin_last  := ceiling(pos_end   / BIN_SIZE)]
  dt[bin_first == 0L, bin_first := 1L]
  dt[, n_bins := bin_last - bin_first + 1L]

  n_bins_vec <- dt$n_bins
  expanded   <- dt[rep(seq_len(.N), n_bins_vec)]
  expanded[, bin_index := dt$bin_first[rep(seq_len(nrow(dt)), n_bins_vec)] +
                          sequence(n_bins_vec) - 1L]

  expanded[, seg_bin_start := (bin_index - 1L) * BIN_SIZE + 1L]
  expanded[, seg_bin_end   := bin_index * BIN_SIZE]
  expanded[, overlap := pmin(pos_end, seg_bin_end) - pmax(pos_start, seg_bin_start) + 1L]
  expanded <- expanded[overlap > 0L]

  by_label <- expanded[, .(total_overlap = sum(overlap),
                           cr_sum        = sum(continous_rank * overlap)),
                       by = .(chr, bin_index, label_d2)]

  result <- by_label[by_label[, .I[which.max(total_overlap)], by = .(chr, bin_index)]$V1]
  result[, continous_rank_wt := cr_sum / total_overlap]
  result[, label_d2 := factor(label_d2, levels = LABEL_ORDER)]

  result <- merge(grid_dt[, .(chr, bin_index, bin_start, bin_end)],
                  result[, .(chr, bin_index, label_d2, continous_rank_wt)],
                  by = c("chr", "bin_index"), all.x = TRUE)
  result[order(chr, bin_index)]
}

compute_transition_matrix <- function(labeled_dt) {
  valid  <- labeled_dt[!is.na(ctrl_label) & !is.na(mut_label)]
  counts <- valid[, .N, by = .(ctrl_label, mut_label)]
  mat    <- matrix(0L, nrow = 4, ncol = 4, dimnames = list(LABEL_ORDER, LABEL_ORDER))
  for (i in seq_len(nrow(counts))) {
    r <- as.character(counts$ctrl_label[i])
    cc <- as.character(counts$mut_label[i])
    mat[r, cc] <- counts$N[i]
  }
  mat
}

# ── Build reference 100kb grid ─────────────────────────────────────────────────

cat("\n=== Building 100kb reference grid ===\n")
grid_dt <- build_100kb_grid(MM10_CHROM_SIZES)
cat(sprintf("  Total bins: %d across %d chromosomes\n",
            nrow(grid_dt), length(MM10_CHROM_SIZES)))

# ── Load CALDER2 files for this timepoint ─────────────────────────────────────

cat("\n=== Loading CALDER2 outputs ===\n")
calder_data <- list()
for (samp in SAMPLES) {
  path <- input_files[[samp]]
  comp <- fread(path)

  expected_cols <- c("chr", "pos_start", "pos_end", "comp_name",
                     "comp_rank", "continous_rank")
  if (!all(expected_cols %in% names(comp))) {
    stop(sprintf("Unexpected columns in %s: %s",
                 path, paste(names(comp), collapse = ", ")))
  }

  comp[, label_d2 := truncate_to_depth2(comp_name)]
  unknown_labels <- setdiff(unique(comp$label_d2), LABEL_ORDER)
  if (length(unknown_labels) > 0) {
    warning(sprintf("Unexpected depth-2 labels in %s/%s: %s",
                    TP, samp, paste(unknown_labels, collapse = ", ")))
  }

  cat(sprintf("  %-25s %7d rows  chroms=%d\n",
              samp, nrow(comp), length(unique(comp$chr))))
  calder_data[[samp]] <- comp[, .(chr, pos_start, pos_end, label_d2, continous_rank)]
}

# ── Plurality-vote binning ─────────────────────────────────────────────────────

cat("\n=== Binning to 100kb grid ===\n")

cat("  Binning ctrl to 100kb grid...\n")
ctrl_binned <- bin_plurality_vote(calder_data[["ctrl_merged"]], grid_dt)
setnames(ctrl_binned,
         c("label_d2", "continous_rank_wt"),
         c("ctrl_label", "continous_rank_ctrl"))

cat("  Binning mut to 100kb grid...\n")
mut_binned <- bin_plurality_vote(calder_data[["mut_merged"]], grid_dt)
setnames(mut_binned,
         c("label_d2", "continous_rank_wt"),
         c("mut_label", "continous_rank_mut"))

# ── Join ctrl and mut ──────────────────────────────────────────────────────────

labeled <- merge(
  ctrl_binned[, .(chr, bin_index, bin_start, bin_end,
                  ctrl_label, continous_rank_ctrl)],
  mut_binned[, .(chr, bin_index, mut_label, continous_rank_mut)],
  by = c("chr", "bin_index"), all = TRUE
)
labeled[, label_changed := (!is.na(ctrl_label) & !is.na(mut_label) &
                              as.character(ctrl_label) != as.character(mut_label))]

n_total    <- nrow(labeled)
n_callable <- sum(!is.na(labeled$ctrl_label) & !is.na(labeled$mut_label))
n_changed  <- sum(labeled$label_changed, na.rm = TRUE)
pct_changed <- 100 * n_changed / n_callable

cat(sprintf("\n  Total bins:    %d\n", n_total))
cat(sprintf("  Callable bins: %d (both ctrl and mut labeled)\n", n_callable))
cat(sprintf("  Changed bins:  %d (%.1f%% of callable)\n", n_changed, pct_changed))

# ── Write labels TSV ───────────────────────────────────────────────────────────

out_labels <- file.path(OUT_DIR, sprintf("%s_subcompartment_labels_100kb.tsv", TP))
fwrite(labeled[, .(chr, bin_start, bin_end, ctrl_label, mut_label,
                   continous_rank_ctrl, continous_rank_mut, label_changed)],
       out_labels, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s\n", basename(out_labels)))

# ── Transition matrix ──────────────────────────────────────────────────────────

cat("\n=== Transition analysis ===\n")
trans_mat <- compute_transition_matrix(labeled)
cat(sprintf("  Transition matrix total: %d (expected ~24,639)\n", sum(trans_mat)))

# ── Chi-squared test ───────────────────────────────────────────────────────────

chisq_result <- chisq.test(trans_mat, simulate.p.value = FALSE)
cramer_v     <- sqrt(chisq_result$statistic / (sum(trans_mat) * (min(nrow(trans_mat), ncol(trans_mat)) - 1)))

cat(sprintf("  Chi-squared: X2=%.1f, df=%d, p=%.3e\n",
            chisq_result$statistic, chisq_result$parameter, chisq_result$p.value))
cat(sprintf("  Cramer's V:  %.4f\n", unname(cramer_v)))
cat(sprintf("  %% bins changed: %.1f%%\n", pct_changed))

# ── Write transition matrix ───────────────────────────────────────────────────

out_mat <- file.path(OUT_DIR, sprintf("%s_transition_matrix.tsv", TP))
write.table(trans_mat, out_mat, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("  Written: %s\n", basename(out_mat)))

# ── Write transition summary ──────────────────────────────────────────────────

trans_long <- as.data.table(as.table(trans_mat))
names(trans_long) <- c("ctrl_label", "mut_label", "count")
trans_long[, pct_of_total := round(100 * count / n_callable, 2)]
trans_long[, is_diagonal := (as.character(ctrl_label) == as.character(mut_label))]
trans_long <- trans_long[order(ctrl_label, mut_label)]
out_summary <- file.path(OUT_DIR, sprintf("%s_transition_summary.tsv", TP))
fwrite(trans_long, out_summary, sep = "\t", quote = FALSE)
cat(sprintf("  Written: %s\n", basename(out_summary)))

# ── Figures ────────────────────────────────────────────────────────────────────

cat("\n=== Generating figures ===\n")

tp_label <- sprintf("%s (%s)", TP, TP_LABELS[TP])

# 1. Transition heatmap
mat_long <- as.data.table(as.table(trans_mat))
names(mat_long) <- c("ctrl_label", "mut_label", "count")
mat_long[, log10_count := log10(count + 1)]
mat_long[, on_diagonal := (as.character(ctrl_label) == as.character(mut_label))]
mat_long[, ctrl_label := factor(ctrl_label, levels = rev(LABEL_ORDER))]
mat_long[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]

p_heatmap <- ggplot(mat_long, aes(x = mut_label, y = ctrl_label, fill = log10_count)) +
  geom_tile(aes(color = on_diagonal), linewidth = 1.2) +
  geom_text(aes(label = format(count, big.mark = ",")), fontface = "bold", size = 5) +
  scale_fill_gradient(low = "white", high = "#2166ac", name = "log10(count+1)") +
  scale_color_manual(values = c("TRUE" = "#d7301f", "FALSE" = "grey85"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = "Mutant label", y = "Control label",
       title = sprintf("Subcompartment Transitions: %s", tp_label),
       subtitle = sprintf("X²=%.1f, p=%.2e, %.1f%% bins changed",
                          chisq_result$statistic, chisq_result$p.value, pct_changed)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text  = element_text(face = "bold", size = 12))

save_multiformat_ggplot(p_heatmap,
  file.path(OUT_DIR, sprintf("%s_transition_heatmap", TP)),
  width = 7, height = 6)

# 2. Sankey / alluvial diagram
alluv_dt <- labeled[!is.na(ctrl_label) & !is.na(mut_label),
                     .N, by = .(ctrl_label, mut_label)]
alluv_dt[, flow_type := fifelse(as.character(ctrl_label) == as.character(mut_label),
                                 as.character(ctrl_label), "change")]
alluv_dt[, ctrl_label := factor(ctrl_label, levels = LABEL_ORDER)]
alluv_dt[, mut_label  := factor(mut_label,  levels = LABEL_ORDER)]

p_sankey <- ggplot(alluv_dt, aes(axis1 = ctrl_label, axis2 = mut_label, y = N)) +
  geom_alluvium(aes(fill = flow_type), width = 1/4, alpha = 0.75, knot.pos = 0.4) +
  geom_stratum(width = 1/4, fill = "white", color = "grey40", linewidth = 0.5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = LABEL_COLORS, name = "Flow") +
  scale_x_discrete(limits = c("Control", "Mutant"), expand = c(0.05, 0.05)) +
  labs(y = "Number of 100kb bins",
       title = sprintf("Subcompartment Flow: %s", tp_label),
       subtitle = sprintf("%.1f%% of bins change subcompartment", pct_changed)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

save_multiformat_ggplot(p_sankey,
  file.path(OUT_DIR, sprintf("%s_transition_sankey", TP)),
  width = 8, height = 7)

# 3. Genome percentage barplot
pct_long <- rbind(
  labeled[!is.na(ctrl_label), .(label = ctrl_label, condition = "Control")],
  labeled[!is.na(mut_label),  .(label = mut_label,  condition = "Mutant")]
)
pct_long[, condition := factor(condition, levels = c("Control", "Mutant"))]
pct_long[, label := factor(label, levels = LABEL_ORDER)]
pct_summary <- pct_long[, .N, by = .(condition, label)]
pct_summary[, pct := 100 * N / sum(N), by = condition]

p_pct <- ggplot(pct_summary, aes(x = condition, y = pct, fill = label)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = LABEL_COLORS[LABEL_ORDER], name = "Subcompartment") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  labs(x = NULL, y = "% of 100kb bins",
       title = sprintf("Subcompartment Fractions: %s", tp_label)) +
  theme_minimal(base_size = 14)

save_multiformat_ggplot(p_pct,
  file.path(OUT_DIR, sprintf("%s_subcompartment_genome_pct", TP)),
  width = 5, height = 7)

# ── Verification ───────────────────────────────────────────────────────────────

cat("\n=== Verification ===\n")
any_warn <- FALSE

total <- sum(trans_mat)
cat(sprintf("  Matrix total: %d (expected ~24,639)\n", total))
if (total < 20000 || total > 25000) {
  warning(sprintf("Unexpected bin count for %s: %d", TP, total))
  any_warn <- TRUE
}

if (TP == "250402" && chisq_result$p.value >= 0.05) {
  warning(sprintf("Late timepoint %s chi-squared p=%.4f (expected p<0.05)",
                  TP, chisq_result$p.value))
  any_warn <- TRUE
}

col_totals <- colSums(trans_mat)
for (lbl in LABEL_ORDER) {
  pct <- 100 * col_totals[lbl] / total
  if (pct < 5 || pct > 60) {
    warning(sprintf("Unexpected %s fraction in mut: %.1f%%", lbl, pct))
    any_warn <- TRUE
  }
}

expected_files <- c(
  sprintf("%s_subcompartment_labels_100kb.tsv", TP),
  sprintf("%s_transition_matrix.tsv", TP),
  sprintf("%s_transition_summary.tsv", TP)
)
for (f in expected_files) {
  fpath <- file.path(OUT_DIR, f)
  if (!file.exists(fpath)) {
    warning(sprintf("Missing output: %s", fpath))
    any_warn <- TRUE
  }
}

if (!any_warn) cat("  All checks passed.\n")

cat("\n===========================================\n")
cat(sprintf("A2 COMPLETE: %s\n", TP))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
