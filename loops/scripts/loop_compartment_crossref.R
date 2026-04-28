# scripts/loop_compartment_crossref.R
# Loop-Compartment/TAD Cross-Reference Analysis
#
# Tests whether gained/lost differential loops localize to specific
# compartment types (A vs B) and compartment shift regions, providing
# independent validation of loop rewriting in BAP1-KO.

# ============================================================================
# SECTION 0: CONFIGURATION
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

source("scripts/utils/multi_format_output.R")

# --- Input paths ---
LOOPS_FILE <- "outputs/250402-late_outputs/loop_annotation_extended/late/extended_characterized_loops.tsv"
DIFFPC_FILE <- "tads/tad-pc-analysis/inputs/late/diffPC/diffcompartments.txt"
ASHIFT_FILE <- "tads/tad-pc-analysis/inputs/late/diffPC/regions.Up_mut_vs_ctrl.regions.txt"
BSHIFT_FILE <- "tads/tad-pc-analysis/inputs/late/diffPC/regions.Down_mut_vs_ctrl.regions.txt"
DIFFTAD_FILE <- "tads/tad-pc-analysis/inputs/late/diffTAD/Bap1.diff.tad.txt"

# --- Output directory ---
OUTPUT_DIR <- "output/loop_compartment_crossref"
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")

dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Color palette ---
COLORS <- list(
  gained = "#E74C3C",
  lost = "#3498DB",
  a_shift = "#2ECC71",
  b_shift = "#9B59B6",
  polycomb = "#8E44AD",
  non_polycomb = "#BDC3C7"
)

cat("=================================================================\n")
cat("LOOP-COMPARTMENT/TAD CROSS-REFERENCE ANALYSIS\n")
cat("=================================================================\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Output: %s\n\n", OUTPUT_DIR))

# ============================================================================
# SECTION 1: LOAD LOOP DATA
# ============================================================================

cat("\n--- Section 1: Loading loop data ---\n")

stopifnot(file.exists(LOOPS_FILE))
loops <- read_tsv(LOOPS_FILE, show_col_types = FALSE)

required_cols <- c("loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2",
                   "logFC", "FDR", "direction", "loop_distance",
                   "anchor1_type_extended", "anchor2_type_extended")
missing <- setdiff(required_cols, colnames(loops))
if (length(missing) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
}

gained <- loops %>% filter(direction == "up_in_mutant")
lost <- loops %>% filter(direction == "down_in_mutant")

cat(sprintf("  Total loops: %d\n", nrow(loops)))
cat(sprintf("  Gained (up_in_mutant): %d\n", nrow(gained)))
cat(sprintf("  Lost (down_in_mutant): %d\n", nrow(lost)))

# Helper to create GRanges from loop anchors
make_anchor_gr <- function(df, anchor = 1) {
  if (anchor == 1) {
    GRanges(seqnames = df$chr1,
            ranges = IRanges(start = df$start1, end = df$end1),
            loop_id = df$loop_id)
  } else {
    GRanges(seqnames = df$chr2,
            ranges = IRanges(start = df$start2, end = df$end2),
            loop_id = df$loop_id)
  }
}

gained_a1 <- make_anchor_gr(gained, 1)
gained_a2 <- make_anchor_gr(gained, 2)
lost_a1 <- make_anchor_gr(lost, 1)
lost_a2 <- make_anchor_gr(lost, 2)

gained_anchors <- c(gained_a1, gained_a2)
lost_anchors <- c(lost_a1, lost_a2)

cat(sprintf("  Gained anchors: %d (2 per loop)\n", length(gained_anchors)))
cat(sprintf("  Lost anchors: %d (2 per loop)\n", length(lost_anchors)))

# ============================================================================
# SECTION 2: LOAD HOMER COMPARTMENT DATA
# ============================================================================

cat("\n--- Section 2: Loading HOMER compartment data ---\n")

stopifnot(file.exists(DIFFPC_FILE))
pc_raw <- read.table(DIFFPC_FILE, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE, quote = "")

cat(sprintf("  Loaded %d compartment bins\n", nrow(pc_raw)))

# Parse by column position (HOMER header is messy with full paths)
# Col 2: Chr, Col 3: Start, Col 4: End
# Cols 20-22: ctrl PC1 values, Cols 23-25: mut PC1 values
# Col 26: Difference, Col 27: p-value, Col 28: adj. p-value
pc_df <- data.frame(
  chr = pc_raw[[2]],
  start = as.numeric(pc_raw[[3]]),
  end = as.numeric(pc_raw[[4]]),
  ctrl_M1 = as.numeric(pc_raw[[20]]),
  ctrl_M2 = as.numeric(pc_raw[[21]]),
  ctrl_M3 = as.numeric(pc_raw[[22]]),
  mut_M1 = as.numeric(pc_raw[[23]]),
  mut_M2 = as.numeric(pc_raw[[24]]),
  mut_M3 = as.numeric(pc_raw[[25]]),
  pc1_diff = as.numeric(pc_raw[[26]]),
  pvalue = as.numeric(pc_raw[[27]]),
  fdr = as.numeric(pc_raw[[28]]),
  stringsAsFactors = FALSE
)

pc_df$mean_ctrl_pc1 <- rowMeans(pc_df[, c("ctrl_M1", "ctrl_M2", "ctrl_M3")], na.rm = TRUE)
pc_df$mean_mut_pc1 <- rowMeans(pc_df[, c("mut_M1", "mut_M2", "mut_M3")], na.rm = TRUE)

# Clean NaN from rowMeans (all-NA rows)
pc_df$mean_ctrl_pc1[is.nan(pc_df$mean_ctrl_pc1)] <- NA
pc_df$mean_mut_pc1[is.nan(pc_df$mean_mut_pc1)] <- NA

pc_df <- pc_df[!is.na(pc_df$pc1_diff) & !is.na(pc_df$chr), ]
cat(sprintf("  After cleanup: %d bins\n", nrow(pc_df)))
cat(sprintf("  PC1 difference range: [%.3f, %.3f]\n",
            min(pc_df$pc1_diff), max(pc_df$pc1_diff)))

pc_gr <- GRanges(seqnames = pc_df$chr,
                 ranges = IRanges(start = pc_df$start, end = pc_df$end))
mcols(pc_gr) <- pc_df[, c("mean_ctrl_pc1", "mean_mut_pc1", "pc1_diff", "pvalue", "fdr")]

# Load pre-computed shift regions
cat("  Loading shift regions...\n")
stopifnot(file.exists(ASHIFT_FILE))
stopifnot(file.exists(BSHIFT_FILE))

ashift_raw <- read.table(ASHIFT_FILE, sep = "\t", header = TRUE,
                         comment.char = "", check.names = FALSE, quote = "")
bshift_raw <- read.table(BSHIFT_FILE, sep = "\t", header = TRUE,
                         comment.char = "", check.names = FALSE, quote = "")

# Col 2: chr, Col 3: start, Col 4: end
ashift_gr <- GRanges(seqnames = ashift_raw[[2]],
                     ranges = IRanges(start = as.numeric(ashift_raw[[3]]),
                                     end = as.numeric(ashift_raw[[4]])))
bshift_gr <- GRanges(seqnames = bshift_raw[[2]],
                     ranges = IRanges(start = as.numeric(bshift_raw[[3]]),
                                     end = as.numeric(bshift_raw[[4]])))

cat(sprintf("  A-shift regions (more active in mut): %d\n", length(ashift_gr)))
cat(sprintf("  B-shift regions (more inactive in mut): %d\n", length(bshift_gr)))

# ============================================================================
# SECTION 3: LOAD HOMER DIFFERENTIAL TAD DATA
# ============================================================================

cat("\n--- Section 3: Loading HOMER differential TAD data ---\n")

stopifnot(file.exists(DIFFTAD_FILE))
tad_raw <- read.table(DIFFTAD_FILE, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, comment.char = "",
                      check.names = FALSE, quote = "")

cat(sprintf("  Loaded %d TADs\n", nrow(tad_raw)))

# Col 2: chr1, Col 3: start1, Col 4: end1
# Cols 11-13: ctrl IR, Cols 14-16: mut IR
# Col 17: Difference, Col 18: p-value, Col 19: adj. p-value
tad_df <- data.frame(
  chr = tad_raw[[2]],
  start = as.numeric(tad_raw[[3]]),
  end = as.numeric(tad_raw[[4]]),
  ir_diff = as.numeric(tad_raw[[17]]),
  pvalue = as.numeric(tad_raw[[18]]),
  fdr = as.numeric(tad_raw[[19]]),
  stringsAsFactors = FALSE
)

tad_df <- tad_df[!is.na(tad_df$ir_diff) & !is.na(tad_df$chr), ]

tad_gr <- GRanges(seqnames = tad_df$chr,
                  ranges = IRanges(start = tad_df$start, end = tad_df$end))
mcols(tad_gr) <- tad_df[, c("ir_diff", "pvalue", "fdr")]

cat(sprintf("  After cleanup: %d TADs\n", nrow(tad_df)))
cat(sprintf("  TAD IR difference range: [%.3f, %.3f]\n",
            min(tad_df$ir_diff), max(tad_df$ir_diff)))

# ============================================================================
# SECTION 4: Q1 — Do loop anchors overlap compartment shift regions?
# ============================================================================

cat("\n--- Section 4: Q1 — Loop anchor overlap with compartment shift regions ---\n")

compute_shift_overlap <- function(anchor_gr, ashift_gr, bshift_gr) {
  a_overlap <- countOverlaps(anchor_gr, ashift_gr) > 0
  b_overlap <- countOverlaps(anchor_gr, bshift_gr) > 0
  data.frame(
    n_anchors = length(anchor_gr),
    n_a_shift = sum(a_overlap),
    n_b_shift = sum(b_overlap),
    n_neither = sum(!a_overlap & !b_overlap),
    pct_a_shift = 100 * sum(a_overlap) / length(anchor_gr),
    pct_b_shift = 100 * sum(b_overlap) / length(anchor_gr),
    pct_neither = 100 * sum(!a_overlap & !b_overlap) / length(anchor_gr)
  )
}

gained_shift <- compute_shift_overlap(gained_anchors, ashift_gr, bshift_gr)
lost_shift <- compute_shift_overlap(lost_anchors, ashift_gr, bshift_gr)

cat("  Gained loop anchors:\n")
cat(sprintf("    A-shift: %d (%.1f%%)\n", gained_shift$n_a_shift, gained_shift$pct_a_shift))
cat(sprintf("    B-shift: %d (%.1f%%)\n", gained_shift$n_b_shift, gained_shift$pct_b_shift))
cat(sprintf("    Neither: %d (%.1f%%)\n", gained_shift$n_neither, gained_shift$pct_neither))
cat("  Lost loop anchors:\n")
cat(sprintf("    A-shift: %d (%.1f%%)\n", lost_shift$n_a_shift, lost_shift$pct_a_shift))
cat(sprintf("    B-shift: %d (%.1f%%)\n", lost_shift$n_b_shift, lost_shift$pct_b_shift))
cat(sprintf("    Neither: %d (%.1f%%)\n", lost_shift$n_neither, lost_shift$pct_neither))

# Fisher's exact test: per-loop (either anchor) to avoid double-counting
gained_any_bshift <- countOverlaps(gained_a1, bshift_gr) > 0 |
                     countOverlaps(gained_a2, bshift_gr) > 0
lost_any_bshift <- countOverlaps(lost_a1, bshift_gr) > 0 |
                   countOverlaps(lost_a2, bshift_gr) > 0

fisher_bshift <- matrix(
  c(sum(gained_any_bshift), sum(lost_any_bshift),
    sum(!gained_any_bshift), sum(!lost_any_bshift)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("B-shift overlap", "No B-shift"),
                  c("Gained", "Lost"))
)
fisher_bshift_result <- fisher.test(fisher_bshift)

cat(sprintf("\n  Fisher's exact test (B-shift enrichment, gained vs lost):\n"))
cat(sprintf("    OR = %.2f, p = %.2e\n",
            fisher_bshift_result$estimate, fisher_bshift_result$p.value))

# Also test A-shift
gained_any_ashift <- countOverlaps(gained_a1, ashift_gr) > 0 |
                     countOverlaps(gained_a2, ashift_gr) > 0
lost_any_ashift <- countOverlaps(lost_a1, ashift_gr) > 0 |
                   countOverlaps(lost_a2, ashift_gr) > 0

fisher_ashift <- matrix(
  c(sum(gained_any_ashift), sum(lost_any_ashift),
    sum(!gained_any_ashift), sum(!lost_any_ashift)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("A-shift overlap", "No A-shift"),
                  c("Gained", "Lost"))
)
fisher_ashift_result <- fisher.test(fisher_ashift)

cat(sprintf("  Fisher's exact test (A-shift enrichment, gained vs lost):\n"))
cat(sprintf("    OR = %.2f, p = %.2e\n",
            fisher_ashift_result$estimate, fisher_ashift_result$p.value))

# --- Visualization: Grouped bar plot ---
shift_plot_df <- data.frame(
  direction = rep(c("Gained", "Lost"), each = 3),
  region = rep(c("A-shift\n(more active)", "B-shift\n(more inactive)", "Neither"), 2),
  pct = c(gained_shift$pct_a_shift, gained_shift$pct_b_shift, gained_shift$pct_neither,
          lost_shift$pct_a_shift, lost_shift$pct_b_shift, lost_shift$pct_neither),
  count = c(gained_shift$n_a_shift, gained_shift$n_b_shift, gained_shift$n_neither,
            lost_shift$n_a_shift, lost_shift$n_b_shift, lost_shift$n_neither)
)
shift_plot_df$direction <- factor(shift_plot_df$direction, levels = c("Gained", "Lost"))
shift_plot_df$region <- factor(shift_plot_df$region,
                               levels = c("A-shift\n(more active)",
                                          "B-shift\n(more inactive)", "Neither"))

p_shift <- ggplot(shift_plot_df, aes(x = region, y = pct, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = count), position = position_dodge(0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Gained" = COLORS$gained, "Lost" = COLORS$lost),
                    name = "Loop Direction") +
  labs(title = "Loop Anchor Overlap with Compartment Shift Regions",
       subtitle = sprintf("Fisher B-shift: OR=%.2f, p=%.2e | A-shift: OR=%.2f, p=%.2e",
                          fisher_bshift_result$estimate, fisher_bshift_result$p.value,
                          fisher_ashift_result$estimate, fisher_ashift_result$p.value),
       x = "Compartment Shift Region",
       y = "% of Anchors") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10))

save_multiformat_ggplot(p_shift,
                        file.path(PLOTS_DIR, "q1_shift_overlap_barplot"),
                        width = 9, height = 7)

# ============================================================================
# SECTION 5: Q2 — Mean PC1 value and PC1 difference at loop anchors
# ============================================================================

cat("\n--- Section 5: Q2 — PC1 values at loop anchors ---\n")

# Extract PC1 values at both anchors of each loop
extract_pc1_at_anchors <- function(loops_df, pc_gr) {
  a1_gr <- make_anchor_gr(loops_df, 1)
  a2_gr <- make_anchor_gr(loops_df, 2)

  ov1 <- findOverlaps(a1_gr, pc_gr)
  ov2 <- findOverlaps(a2_gr, pc_gr)

  # Average PC1 across overlapping bins per anchor
  a1_pc1_ctrl <- tapply(pc_gr$mean_ctrl_pc1[subjectHits(ov1)], queryHits(ov1), mean)
  a1_pc1_mut  <- tapply(pc_gr$mean_mut_pc1[subjectHits(ov1)],  queryHits(ov1), mean)
  a1_pc1_diff <- tapply(pc_gr$pc1_diff[subjectHits(ov1)],      queryHits(ov1), mean)

  a2_pc1_ctrl <- tapply(pc_gr$mean_ctrl_pc1[subjectHits(ov2)], queryHits(ov2), mean)
  a2_pc1_mut  <- tapply(pc_gr$mean_mut_pc1[subjectHits(ov2)],  queryHits(ov2), mean)
  a2_pc1_diff <- tapply(pc_gr$pc1_diff[subjectHits(ov2)],      queryHits(ov2), mean)

  n <- nrow(loops_df)
  result <- data.frame(
    anchor1_pc1_ctrl = rep(NA_real_, n),
    anchor1_pc1_mut  = rep(NA_real_, n),
    anchor1_pc1_diff = rep(NA_real_, n),
    anchor2_pc1_ctrl = rep(NA_real_, n),
    anchor2_pc1_mut  = rep(NA_real_, n),
    anchor2_pc1_diff = rep(NA_real_, n)
  )

  result$anchor1_pc1_ctrl[as.integer(names(a1_pc1_ctrl))] <- a1_pc1_ctrl
  result$anchor1_pc1_mut[as.integer(names(a1_pc1_mut))]   <- a1_pc1_mut
  result$anchor1_pc1_diff[as.integer(names(a1_pc1_diff))] <- a1_pc1_diff
  result$anchor2_pc1_ctrl[as.integer(names(a2_pc1_ctrl))] <- a2_pc1_ctrl
  result$anchor2_pc1_mut[as.integer(names(a2_pc1_mut))]   <- a2_pc1_mut
  result$anchor2_pc1_diff[as.integer(names(a2_pc1_diff))] <- a2_pc1_diff

  # Average across both anchors (NaN → NA for all-NA cases)
  result$mean_pc1_ctrl <- rowMeans(result[, c("anchor1_pc1_ctrl", "anchor2_pc1_ctrl")], na.rm = TRUE)
  result$mean_pc1_mut  <- rowMeans(result[, c("anchor1_pc1_mut", "anchor2_pc1_mut")],   na.rm = TRUE)
  result$mean_pc1_diff <- rowMeans(result[, c("anchor1_pc1_diff", "anchor2_pc1_diff")], na.rm = TRUE)

  result$mean_pc1_ctrl[is.nan(result$mean_pc1_ctrl)] <- NA
  result$mean_pc1_mut[is.nan(result$mean_pc1_mut)]   <- NA
  result$mean_pc1_diff[is.nan(result$mean_pc1_diff)] <- NA

  return(result)
}

gained_pc1 <- extract_pc1_at_anchors(gained, pc_gr)
lost_pc1   <- extract_pc1_at_anchors(lost, pc_gr)

cat(sprintf("  Gained loops with PC1 data: %d / %d (%.1f%%)\n",
            sum(!is.na(gained_pc1$mean_pc1_diff)), nrow(gained),
            100 * sum(!is.na(gained_pc1$mean_pc1_diff)) / nrow(gained)))
cat(sprintf("  Lost loops with PC1 data: %d / %d (%.1f%%)\n",
            sum(!is.na(lost_pc1$mean_pc1_diff)), nrow(lost),
            100 * sum(!is.na(lost_pc1$mean_pc1_diff)) / nrow(lost)))

# Wilcoxon tests
wilcox_pc1_diff <- wilcox.test(
  gained_pc1$mean_pc1_diff[!is.na(gained_pc1$mean_pc1_diff)],
  lost_pc1$mean_pc1_diff[!is.na(lost_pc1$mean_pc1_diff)]
)

wilcox_ctrl_pc1 <- wilcox.test(
  gained_pc1$mean_pc1_ctrl[!is.na(gained_pc1$mean_pc1_ctrl)],
  lost_pc1$mean_pc1_ctrl[!is.na(lost_pc1$mean_pc1_ctrl)]
)

cat(sprintf("\n  Wilcoxon test — PC1 difference (gained vs lost):\n"))
cat(sprintf("    Gained median PC1_diff: %.4f\n",
            median(gained_pc1$mean_pc1_diff, na.rm = TRUE)))
cat(sprintf("    Lost median PC1_diff: %.4f\n",
            median(lost_pc1$mean_pc1_diff, na.rm = TRUE)))
cat(sprintf("    W = %.0f, p = %.2e\n",
            wilcox_pc1_diff$statistic, wilcox_pc1_diff$p.value))

cat(sprintf("\n  Wilcoxon test — Baseline ctrl PC1 (gained vs lost):\n"))
cat(sprintf("    Gained median ctrl_PC1: %.4f\n",
            median(gained_pc1$mean_pc1_ctrl, na.rm = TRUE)))
cat(sprintf("    Lost median ctrl_PC1: %.4f\n",
            median(lost_pc1$mean_pc1_ctrl, na.rm = TRUE)))
cat(sprintf("    W = %.0f, p = %.2e\n",
            wilcox_ctrl_pc1$statistic, wilcox_ctrl_pc1$p.value))

# --- Visualization 5a: Violin — PC1 difference ---
violin_df <- data.frame(
  direction = c(rep("Gained", nrow(gained_pc1)), rep("Lost", nrow(lost_pc1))),
  pc1_diff = c(gained_pc1$mean_pc1_diff, lost_pc1$mean_pc1_diff)
) %>% filter(!is.na(pc1_diff))
violin_df$direction <- factor(violin_df$direction, levels = c("Gained", "Lost"))

p_pc1diff <- ggplot(violin_df, aes(x = direction, y = pc1_diff, fill = direction)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Gained" = COLORS$gained, "Lost" = COLORS$lost)) +
  labs(title = "PC1 Difference (mut - ctrl) at Loop Anchors",
       subtitle = sprintf("Wilcoxon p = %.2e | Positive = A-shift, Negative = B-shift",
                          wilcox_pc1_diff$p.value),
       x = "Loop Direction", y = "Mean PC1 Difference (mut - ctrl)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

save_multiformat_ggplot(p_pc1diff,
                        file.path(PLOTS_DIR, "q2_pc1_diff_violin"),
                        width = 7, height = 7)

# --- Visualization 5b: Violin — baseline ctrl PC1 ---
violin_ctrl_df <- data.frame(
  direction = c(rep("Gained", nrow(gained_pc1)), rep("Lost", nrow(lost_pc1))),
  ctrl_pc1 = c(gained_pc1$mean_pc1_ctrl, lost_pc1$mean_pc1_ctrl)
) %>% filter(!is.na(ctrl_pc1))
violin_ctrl_df$direction <- factor(violin_ctrl_df$direction, levels = c("Gained", "Lost"))

p_ctrl_pc1 <- ggplot(violin_ctrl_df, aes(x = direction, y = ctrl_pc1, fill = direction)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Gained" = COLORS$gained, "Lost" = COLORS$lost)) +
  labs(title = "Baseline (ctrl) PC1 at Loop Anchors",
       subtitle = sprintf("Wilcoxon p = %.2e | A compartment > 0, B compartment < 0",
                          wilcox_ctrl_pc1$p.value),
       x = "Loop Direction", y = "Mean Control PC1") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

save_multiformat_ggplot(p_ctrl_pc1,
                        file.path(PLOTS_DIR, "q2_ctrl_pc1_violin"),
                        width = 7, height = 7)

# --- Visualization 5c: Scatter — distance vs PC1 difference ---
scatter_df <- data.frame(
  direction = c(rep("Gained", nrow(gained)), rep("Lost", nrow(lost))),
  loop_distance = c(gained$loop_distance, lost$loop_distance),
  pc1_diff = c(gained_pc1$mean_pc1_diff, lost_pc1$mean_pc1_diff)
) %>% filter(!is.na(pc1_diff))
scatter_df$direction <- factor(scatter_df$direction, levels = c("Gained", "Lost"))

p_scatter <- ggplot(scatter_df, aes(x = loop_distance / 1e6, y = pc1_diff, color = direction)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Gained" = COLORS$gained, "Lost" = COLORS$lost),
                     name = "Direction") +
  scale_x_log10() +
  labs(title = "Loop Distance vs PC1 Difference at Anchors",
       x = "Loop Distance (Mb, log scale)",
       y = "Mean PC1 Difference (mut - ctrl)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 14))

save_multiformat_ggplot(p_scatter,
                        file.path(PLOTS_DIR, "q2_distance_vs_pc1diff_scatter"),
                        width = 9, height = 7)

# ============================================================================
# SECTION 6: Q3 — TAD compaction changes at loop locations
# ============================================================================

cat("\n--- Section 6: Q3 — TAD inclusion ratio changes at loop locations ---\n")

# Find the TAD containing each loop's midpoint
extract_tad_ir <- function(loops_df, tad_gr) {
  mid_chr <- loops_df$chr1
  mid_pos <- as.integer((loops_df$start1 + loops_df$end2) / 2)
  mid_gr <- GRanges(seqnames = mid_chr,
                    ranges = IRanges(start = mid_pos, width = 1))

  ov <- findOverlaps(mid_gr, tad_gr)

  ir_vals <- tapply(tad_gr$ir_diff[subjectHits(ov)], queryHits(ov), mean)

  result <- rep(NA_real_, nrow(loops_df))
  result[as.integer(names(ir_vals))] <- ir_vals
  return(result)
}

gained_ir <- extract_tad_ir(gained, tad_gr)
lost_ir   <- extract_tad_ir(lost, tad_gr)

cat(sprintf("  Gained loops with TAD data: %d / %d (%.1f%%)\n",
            sum(!is.na(gained_ir)), nrow(gained),
            100 * sum(!is.na(gained_ir)) / nrow(gained)))
cat(sprintf("  Lost loops with TAD data: %d / %d (%.1f%%)\n",
            sum(!is.na(lost_ir)), nrow(lost),
            100 * sum(!is.na(lost_ir)) / nrow(lost)))

wilcox_ir <- wilcox.test(gained_ir[!is.na(gained_ir)],
                         lost_ir[!is.na(lost_ir)])

cat(sprintf("\n  Wilcoxon test — TAD IR difference (gained vs lost):\n"))
cat(sprintf("    Gained median IR_diff: %.4f\n", median(gained_ir, na.rm = TRUE)))
cat(sprintf("    Lost median IR_diff: %.4f\n", median(lost_ir, na.rm = TRUE)))
cat(sprintf("    W = %.0f, p = %.2e\n", wilcox_ir$statistic, wilcox_ir$p.value))

# --- Visualization: Box plot ---
ir_df <- data.frame(
  direction = c(rep("Gained", length(gained_ir)), rep("Lost", length(lost_ir))),
  ir_diff = c(gained_ir, lost_ir)
) %>% filter(!is.na(ir_diff))
ir_df$direction <- factor(ir_df$direction, levels = c("Gained", "Lost"))

p_ir <- ggplot(ir_df, aes(x = direction, y = ir_diff, fill = direction)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Gained" = COLORS$gained, "Lost" = COLORS$lost)) +
  labs(title = "TAD Inclusion Ratio Change at Loop Locations",
       subtitle = sprintf("Wilcoxon p = %.2e | Positive = more insulated in mut",
                          wilcox_ir$p.value),
       x = "Loop Direction",
       y = "TAD IR Difference (mut - ctrl)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

save_multiformat_ggplot(p_ir,
                        file.path(PLOTS_DIR, "q3_tad_ir_boxplot"),
                        width = 7, height = 7)

# ============================================================================
# SECTION 7: Q4 — Distance x compartment interaction
# ============================================================================

cat("\n--- Section 7: Q4 — Distance x compartment interaction ---\n")

all_dist_pc1 <- data.frame(
  direction = c(rep("Gained", nrow(gained)), rep("Lost", nrow(lost))),
  loop_distance = c(gained$loop_distance, lost$loop_distance),
  pc1_diff = c(gained_pc1$mean_pc1_diff, lost_pc1$mean_pc1_diff),
  ctrl_pc1 = c(gained_pc1$mean_pc1_ctrl, lost_pc1$mean_pc1_ctrl)
) %>%
  filter(!is.na(pc1_diff)) %>%
  mutate(
    distance_bin = case_when(
      loop_distance < 250000 ~ "Short (<250kb)",
      loop_distance <= 1000000 ~ "Medium (250kb-1Mb)",
      TRUE ~ "Long (>1Mb)"
    ),
    distance_bin = factor(distance_bin,
                          levels = c("Short (<250kb)", "Medium (250kb-1Mb)", "Long (>1Mb)"))
  )
all_dist_pc1$direction <- factor(all_dist_pc1$direction, levels = c("Gained", "Lost"))

dist_summary <- all_dist_pc1 %>%
  group_by(distance_bin, direction) %>%
  summarise(
    n = n(),
    median_pc1_diff = median(pc1_diff, na.rm = TRUE),
    mean_pc1_diff = mean(pc1_diff, na.rm = TRUE),
    median_ctrl_pc1 = median(ctrl_pc1, na.rm = TRUE),
    .groups = "drop"
  )

cat("  Distance bin x direction summary:\n")
print(as.data.frame(dist_summary), row.names = FALSE)

# Wilcoxon within each distance bin
cat("\n  Wilcoxon tests within distance bins:\n")
dist_wilcox_results <- list()
for (bin in levels(all_dist_pc1$distance_bin)) {
  g <- all_dist_pc1 %>% filter(distance_bin == bin, direction == "Gained") %>% pull(pc1_diff)
  l <- all_dist_pc1 %>% filter(distance_bin == bin, direction == "Lost") %>% pull(pc1_diff)
  if (length(g) >= 3 && length(l) >= 3) {
    wt <- wilcox.test(g, l)
    dist_wilcox_results[[bin]] <- wt
    cat(sprintf("    %s: W=%.0f, p=%.2e (n_gained=%d, n_lost=%d)\n",
                bin, wt$statistic, wt$p.value, length(g), length(l)))
  } else {
    cat(sprintf("    %s: insufficient data (n_gained=%d, n_lost=%d)\n",
                bin, length(g), length(l)))
  }
}

# --- Visualization: Faceted violin ---
p_dist_pc1 <- ggplot(all_dist_pc1, aes(x = direction, y = pc1_diff, fill = direction)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~distance_bin, scales = "free_x") +
  scale_fill_manual(values = c("Gained" = COLORS$gained, "Lost" = COLORS$lost)) +
  labs(title = "PC1 Difference by Distance Bin and Loop Direction",
       x = "Loop Direction",
       y = "Mean PC1 Difference (mut - ctrl)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14))

save_multiformat_ggplot(p_dist_pc1,
                        file.path(PLOTS_DIR, "q4_distance_pc1_violin"),
                        width = 12, height = 7)

# ============================================================================
# SECTION 8: POLYCOMB-SPECIFIC COMPARTMENT ANALYSIS
# ============================================================================

cat("\n--- Section 8: Polycomb-specific compartment analysis ---\n")

# Identify gained loops with at least one Polycomb anchor
gained_has_polycomb <- grepl("Polycomb", gained$anchor1_type_extended) |
                       grepl("Polycomb", gained$anchor2_type_extended)

n_poly <- sum(gained_has_polycomb)
n_nonpoly <- sum(!gained_has_polycomb)
cat(sprintf("  Gained loops with Polycomb anchor: %d\n", n_poly))
cat(sprintf("  Gained loops without Polycomb anchor: %d\n", n_nonpoly))

# Compare baseline ctrl PC1
poly_ctrl_pc1 <- gained_pc1$mean_pc1_ctrl[gained_has_polycomb]
nonpoly_ctrl_pc1 <- gained_pc1$mean_pc1_ctrl[!gained_has_polycomb]
poly_ctrl_pc1 <- poly_ctrl_pc1[!is.na(poly_ctrl_pc1)]
nonpoly_ctrl_pc1 <- nonpoly_ctrl_pc1[!is.na(nonpoly_ctrl_pc1)]

wilcox_poly_ctrl <- wilcox.test(poly_ctrl_pc1, nonpoly_ctrl_pc1)

cat(sprintf("\n  Wilcoxon test — Baseline ctrl PC1 (Polycomb vs non-Polycomb gained):\n"))
cat(sprintf("    Polycomb median ctrl_PC1: %.4f (n=%d)\n",
            median(poly_ctrl_pc1), length(poly_ctrl_pc1)))
cat(sprintf("    Non-Polycomb median ctrl_PC1: %.4f (n=%d)\n",
            median(nonpoly_ctrl_pc1), length(nonpoly_ctrl_pc1)))
cat(sprintf("    W = %.0f, p = %.2e\n",
            wilcox_poly_ctrl$statistic, wilcox_poly_ctrl$p.value))

# Compare PC1 difference
poly_pc1_diff <- gained_pc1$mean_pc1_diff[gained_has_polycomb]
nonpoly_pc1_diff <- gained_pc1$mean_pc1_diff[!gained_has_polycomb]
poly_pc1_diff <- poly_pc1_diff[!is.na(poly_pc1_diff)]
nonpoly_pc1_diff <- nonpoly_pc1_diff[!is.na(nonpoly_pc1_diff)]

wilcox_poly_diff <- wilcox.test(poly_pc1_diff, nonpoly_pc1_diff)

cat(sprintf("\n  Wilcoxon test — PC1 difference (Polycomb vs non-Polycomb gained):\n"))
cat(sprintf("    Polycomb median PC1_diff: %.4f (n=%d)\n",
            median(poly_pc1_diff), length(poly_pc1_diff)))
cat(sprintf("    Non-Polycomb median PC1_diff: %.4f (n=%d)\n",
            median(nonpoly_pc1_diff), length(nonpoly_pc1_diff)))
cat(sprintf("    W = %.0f, p = %.2e\n",
            wilcox_poly_diff$statistic, wilcox_poly_diff$p.value))

# --- Visualization: Violin plot ---
poly_violin_df <- data.frame(
  category = c(rep("Polycomb", length(poly_ctrl_pc1)),
               rep("Non-Polycomb", length(nonpoly_ctrl_pc1))),
  ctrl_pc1 = c(poly_ctrl_pc1, nonpoly_ctrl_pc1)
)
poly_violin_df$category <- factor(poly_violin_df$category,
                                  levels = c("Polycomb", "Non-Polycomb"))

p_poly <- ggplot(poly_violin_df, aes(x = category, y = ctrl_pc1, fill = category)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Polycomb" = COLORS$polycomb,
                               "Non-Polycomb" = COLORS$non_polycomb)) +
  labs(title = "Baseline ctrl PC1 at Gained Loop Anchors",
       subtitle = sprintf("Wilcoxon p = %.2e | Polycomb vs Non-Polycomb",
                          wilcox_poly_ctrl$p.value),
       x = "Anchor Classification",
       y = "Mean Control PC1 (A > 0, B < 0)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

save_multiformat_ggplot(p_poly,
                        file.path(PLOTS_DIR, "q5_polycomb_pc1_violin"),
                        width = 7, height = 7)

# ============================================================================
# SECTION 9: SUMMARY & OUTPUT
# ============================================================================

cat("\n--- Section 9: Writing outputs ---\n")

# --- Annotated loop table ---
all_pc1 <- extract_pc1_at_anchors(loops, pc_gr)
loops_annotated <- cbind(loops, all_pc1)
loops_annotated$tad_ir_diff <- extract_tad_ir(loops, tad_gr)

all_a1 <- make_anchor_gr(loops, 1)
all_a2 <- make_anchor_gr(loops, 2)
loops_annotated$overlaps_Ashift <- countOverlaps(all_a1, ashift_gr) > 0 |
                                   countOverlaps(all_a2, ashift_gr) > 0
loops_annotated$overlaps_Bshift <- countOverlaps(all_a1, bshift_gr) > 0 |
                                   countOverlaps(all_a2, bshift_gr) > 0

write_tsv(loops_annotated, file.path(TABLES_DIR, "loop_compartment_annotated.tsv"))
cat(sprintf("  Saved: loop_compartment_annotated.tsv (%d rows x %d cols)\n",
            nrow(loops_annotated), ncol(loops_annotated)))

# --- Summary report ---
report_lines <- c(
  "=================================================================",
  "LOOP-COMPARTMENT/TAD CROSS-REFERENCE ANALYSIS",
  "=================================================================",
  "",
  sprintf("Date: %s", Sys.time()),
  "",
  "INPUT DATA:",
  sprintf("  Loops: %s", LOOPS_FILE),
  sprintf("    Gained: %d | Lost: %d | Total: %d",
          nrow(gained), nrow(lost), nrow(loops)),
  sprintf("  Compartment bins: %s (%d bins)", DIFFPC_FILE, nrow(pc_df)),
  sprintf("  A-shift regions: %d", length(ashift_gr)),
  sprintf("  B-shift regions: %d", length(bshift_gr)),
  sprintf("  Differential TADs: %d", nrow(tad_df)),
  "",
  "SIGN CONVENTION:",
  "  Positive PC1 difference = mut_PC1 > ctrl_PC1 = A-shift (more active in mutant)",
  "  Negative PC1 difference = mut_PC1 < ctrl_PC1 = B-shift (more inactive in mutant)",
  "",
  "=================================================================",
  "Q1: LOOP ANCHOR OVERLAP WITH COMPARTMENT SHIFT REGIONS",
  "=================================================================",
  "",
  "Gained loop anchors:",
  sprintf("  A-shift overlap: %d / %d (%.1f%%)",
          gained_shift$n_a_shift, gained_shift$n_anchors, gained_shift$pct_a_shift),
  sprintf("  B-shift overlap: %d / %d (%.1f%%)",
          gained_shift$n_b_shift, gained_shift$n_anchors, gained_shift$pct_b_shift),
  sprintf("  Neither: %d / %d (%.1f%%)",
          gained_shift$n_neither, gained_shift$n_anchors, gained_shift$pct_neither),
  "",
  "Lost loop anchors:",
  sprintf("  A-shift overlap: %d / %d (%.1f%%)",
          lost_shift$n_a_shift, lost_shift$n_anchors, lost_shift$pct_a_shift),
  sprintf("  B-shift overlap: %d / %d (%.1f%%)",
          lost_shift$n_b_shift, lost_shift$n_anchors, lost_shift$pct_b_shift),
  sprintf("  Neither: %d / %d (%.1f%%)",
          lost_shift$n_neither, lost_shift$n_anchors, lost_shift$pct_neither),
  "",
  "Fisher's exact tests (per-loop, either-anchor):",
  sprintf("  B-shift enrichment (gained vs lost): OR=%.3f, p=%.2e",
          fisher_bshift_result$estimate, fisher_bshift_result$p.value),
  sprintf("    95%% CI: [%.3f, %.3f]",
          fisher_bshift_result$conf.int[1], fisher_bshift_result$conf.int[2]),
  sprintf("    Gained loops in B-shift: %d / %d (%.1f%%)",
          sum(gained_any_bshift), nrow(gained),
          100 * sum(gained_any_bshift) / nrow(gained)),
  sprintf("    Lost loops in B-shift: %d / %d (%.1f%%)",
          sum(lost_any_bshift), nrow(lost),
          100 * sum(lost_any_bshift) / nrow(lost)),
  "",
  sprintf("  A-shift enrichment (gained vs lost): OR=%.3f, p=%.2e",
          fisher_ashift_result$estimate, fisher_ashift_result$p.value),
  sprintf("    95%% CI: [%.3f, %.3f]",
          fisher_ashift_result$conf.int[1], fisher_ashift_result$conf.int[2]),
  sprintf("    Gained loops in A-shift: %d / %d (%.1f%%)",
          sum(gained_any_ashift), nrow(gained),
          100 * sum(gained_any_ashift) / nrow(gained)),
  sprintf("    Lost loops in A-shift: %d / %d (%.1f%%)",
          sum(lost_any_ashift), nrow(lost),
          100 * sum(lost_any_ashift) / nrow(lost)),
  "",
  "=================================================================",
  "Q2: PC1 VALUES AT LOOP ANCHORS",
  "=================================================================",
  "",
  "PC1 difference (mut - ctrl) at anchors:",
  sprintf("  Gained: median=%.4f, mean=%.4f (n=%d)",
          median(gained_pc1$mean_pc1_diff, na.rm = TRUE),
          mean(gained_pc1$mean_pc1_diff, na.rm = TRUE),
          sum(!is.na(gained_pc1$mean_pc1_diff))),
  sprintf("  Lost: median=%.4f, mean=%.4f (n=%d)",
          median(lost_pc1$mean_pc1_diff, na.rm = TRUE),
          mean(lost_pc1$mean_pc1_diff, na.rm = TRUE),
          sum(!is.na(lost_pc1$mean_pc1_diff))),
  sprintf("  Wilcoxon: W=%.0f, p=%.2e",
          wilcox_pc1_diff$statistic, wilcox_pc1_diff$p.value),
  "",
  "Baseline ctrl PC1 at anchors:",
  sprintf("  Gained: median=%.4f, mean=%.4f (n=%d)",
          median(gained_pc1$mean_pc1_ctrl, na.rm = TRUE),
          mean(gained_pc1$mean_pc1_ctrl, na.rm = TRUE),
          sum(!is.na(gained_pc1$mean_pc1_ctrl))),
  sprintf("  Lost: median=%.4f, mean=%.4f (n=%d)",
          median(lost_pc1$mean_pc1_ctrl, na.rm = TRUE),
          mean(lost_pc1$mean_pc1_ctrl, na.rm = TRUE),
          sum(!is.na(lost_pc1$mean_pc1_ctrl))),
  sprintf("  Wilcoxon: W=%.0f, p=%.2e",
          wilcox_ctrl_pc1$statistic, wilcox_ctrl_pc1$p.value),
  "",
  "=================================================================",
  "Q3: TAD INCLUSION RATIO CHANGES AT LOOP LOCATIONS",
  "=================================================================",
  "",
  sprintf("  Gained: median IR_diff=%.4f (n=%d)",
          median(gained_ir, na.rm = TRUE), sum(!is.na(gained_ir))),
  sprintf("  Lost: median IR_diff=%.4f (n=%d)",
          median(lost_ir, na.rm = TRUE), sum(!is.na(lost_ir))),
  sprintf("  Wilcoxon: W=%.0f, p=%.2e",
          wilcox_ir$statistic, wilcox_ir$p.value),
  "",
  "=================================================================",
  "Q4: DISTANCE x COMPARTMENT INTERACTION",
  "=================================================================",
  ""
)

for (i in seq_len(nrow(dist_summary))) {
  row <- dist_summary[i, ]
  report_lines <- c(report_lines,
    sprintf("  %s | %s: n=%d, median_PC1_diff=%.4f, median_ctrl_PC1=%.4f",
            row$distance_bin, row$direction, row$n,
            row$median_pc1_diff, row$median_ctrl_pc1))
}

# Add within-bin Wilcoxon results
report_lines <- c(report_lines, "")
for (bin in names(dist_wilcox_results)) {
  wt <- dist_wilcox_results[[bin]]
  report_lines <- c(report_lines,
    sprintf("  %s Wilcoxon: W=%.0f, p=%.2e", bin, wt$statistic, wt$p.value))
}

report_lines <- c(report_lines,
  "",
  "=================================================================",
  "POLYCOMB-SPECIFIC COMPARTMENT ANALYSIS (Gained loops only)",
  "=================================================================",
  "",
  sprintf("  Polycomb-anchor gained loops: %d", n_poly),
  sprintf("  Non-Polycomb gained loops: %d", n_nonpoly),
  "",
  "Baseline ctrl PC1:",
  sprintf("  Polycomb: median=%.4f (n=%d)",
          median(poly_ctrl_pc1), length(poly_ctrl_pc1)),
  sprintf("  Non-Polycomb: median=%.4f (n=%d)",
          median(nonpoly_ctrl_pc1), length(nonpoly_ctrl_pc1)),
  sprintf("  Wilcoxon: W=%.0f, p=%.2e",
          wilcox_poly_ctrl$statistic, wilcox_poly_ctrl$p.value),
  "",
  "PC1 difference (mut - ctrl):",
  sprintf("  Polycomb: median=%.4f (n=%d)",
          median(poly_pc1_diff), length(poly_pc1_diff)),
  sprintf("  Non-Polycomb: median=%.4f (n=%d)",
          median(nonpoly_pc1_diff), length(nonpoly_pc1_diff)),
  sprintf("  Wilcoxon: W=%.0f, p=%.2e",
          wilcox_poly_diff$statistic, wilcox_poly_diff$p.value),
  "",
  "=================================================================",
  "INTERPRETATION",
  "=================================================================",
  ""
)

# Conditional interpretation
gained_median_ctrl <- median(gained_pc1$mean_pc1_ctrl, na.rm = TRUE)
lost_median_ctrl <- median(lost_pc1$mean_pc1_ctrl, na.rm = TRUE)

if (wilcox_ctrl_pc1$p.value < 0.05 && gained_median_ctrl < lost_median_ctrl) {
  report_lines <- c(report_lines,
    "Gained loops are located in regions with LOWER baseline PC1 (B-compartment,",
    "heterochromatic) compared to lost loops (A-compartment, active).",
    "This independently validates the loop rewriting model.")
} else if (wilcox_ctrl_pc1$p.value < 0.05) {
  report_lines <- c(report_lines,
    "Significant difference in baseline PC1, but direction is unexpected.",
    "Gained loops are NOT in lower-PC1 regions than lost loops.")
} else {
  report_lines <- c(report_lines,
    "No significant difference in baseline PC1 between gained and lost loops.")
}

if (fisher_bshift_result$p.value < 0.05 && fisher_bshift_result$estimate > 1) {
  report_lines <- c(report_lines,
    "",
    "Gained loops are ENRICHED at B-shift (heterochromatin-shifting) regions,",
    "consistent with BAP1-KO driving Polycomb-mediated chromatin compaction.")
}

if (fisher_ashift_result$p.value < 0.05 && fisher_ashift_result$estimate < 1) {
  report_lines <- c(report_lines,
    "",
    "Gained loops are DEPLETED at A-shift regions relative to lost loops,",
    "further supporting compartment-specific loop rewriting.")
}

if (wilcox_poly_ctrl$p.value < 0.05 &&
    median(poly_ctrl_pc1) < median(nonpoly_ctrl_pc1)) {
  report_lines <- c(report_lines,
    "",
    "Polycomb-anchor gained loops sit in DEEPER B-compartment than",
    "non-Polycomb gained loops, confirming the Polycomb-compaction axis.")
}

report_lines <- c(report_lines,
  "",
  "=================================================================",
  "OUTPUT FILES",
  "=================================================================",
  "",
  "tables/loop_compartment_annotated.tsv - Full loop table with PC1/TAD annotations",
  "",
  "plots/q1_shift_overlap_barplot/ - Shift region overlap (gained vs lost)",
  "plots/q2_pc1_diff_violin/ - PC1 difference at loop anchors",
  "plots/q2_ctrl_pc1_violin/ - Baseline ctrl PC1 at loop anchors",
  "plots/q2_distance_vs_pc1diff_scatter/ - Distance vs PC1 difference",
  "plots/q3_tad_ir_boxplot/ - TAD IR change at loop locations",
  "plots/q4_distance_pc1_violin/ - Distance bin x direction interaction",
  "plots/q5_polycomb_pc1_violin/ - Polycomb vs non-Polycomb ctrl PC1",
  "",
  "=================================================================",
  "END REPORT",
  "================================================================="
)

writeLines(report_lines, file.path(OUTPUT_DIR, "summary_report.txt"))
cat("  Saved: summary_report.txt\n")

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("=================================================================\n")
