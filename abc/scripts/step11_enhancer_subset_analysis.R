# abc/scripts/step11_enhancer_subset_analysis.R
# Stratified analysis of enhancer phenotypic classes across loop metrics,
# ABC scores, and contact profiles.
#
# Tests whether K119ub-only enhancers show contact/loop changes despite no
# K27ac/ATAC changes — the key biological question.
#
# Inputs (all relative to abc/ working directory):
#   enhancer_subsets/enhancer_classes_*.tsv    — 4 enhancer class files
#   ../25042-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe
#   ../characterized_loops.tsv                — differential loops with annotations
#   results/delta_abc_all_pairs.tsv           — 180K E-G pairs
#   results/delta_abc_with_rnaseq.tsv         — E-G pairs + RNA-seq
#   results/gene_level_summary.tsv            — per-gene summary with DE results
#   reference/mm10_tss.bed                    — TSS regions
#
# Outputs:
#   results/enhancer_subset_analysis/         — TSVs + plots (PDF + PNG)
#
# Usage:
#   cd abc && Rscript scripts/step11_enhancer_subset_analysis.R

# =============================================================================
# CONFIGURATION
# =============================================================================

ENHANCER_FILES <- c(
  Activity_Lost = "enhancer_subsets/enhancer_classes_activity_lost.tsv",
  K119ub_Only   = "enhancer_subsets/enhancer_classes_k119ub_only.tsv",
  Activity_Gain = "enhancer_subsets/enhancer_classes_activity_gain.tsv",
  Stable        = "enhancer_subsets/enhancer_classes_stable.tsv"
)

LOOPS_FILE       <- "../outputs/250402-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe"
CHAR_LOOPS_FILE  <- "characterized_loops.tsv"
ABC_PAIRS_FILE   <- "results/delta_abc_all_pairs.tsv"
ABC_RNASEQ_FILE  <- "results/delta_abc_with_rnaseq.tsv"
GENE_SUMMARY_FILE <- "results/gene_level_summary.tsv"
TSS_FILE         <- "reference/mm10_tss.bed"

OUTPUT_DIR <- "results/enhancer_subset_analysis"

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

cat("================================================================================\n")
cat("STEP 11: ENHANCER SUBSET STRATIFIED ANALYSIS\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(patchwork)
  library(scales)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating input files...\n")
for (nm in names(ENHANCER_FILES)) {
  stopifnot(file.exists(ENHANCER_FILES[[nm]]))
  cat(sprintf("  [OK] %s: %s\n", nm, ENHANCER_FILES[[nm]]))
}
stopifnot(file.exists(LOOPS_FILE))
cat(sprintf("  [OK] Loops: %s\n", LOOPS_FILE))
stopifnot(file.exists(CHAR_LOOPS_FILE))
cat(sprintf("  [OK] Characterized loops: %s\n", CHAR_LOOPS_FILE))
stopifnot(file.exists(ABC_PAIRS_FILE))
cat(sprintf("  [OK] ABC pairs: %s\n", ABC_PAIRS_FILE))
stopifnot(file.exists(ABC_RNASEQ_FILE))
cat(sprintf("  [OK] ABC + RNA-seq: %s\n", ABC_RNASEQ_FILE))
stopifnot(file.exists(GENE_SUMMARY_FILE))
cat(sprintf("  [OK] Gene summary: %s\n", GENE_SUMMARY_FILE))
stopifnot(file.exists(TSS_FILE))
cat(sprintf("  [OK] TSS: %s\n", TSS_FILE))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nOutput directory: %s\n\n", OUTPUT_DIR))

# =============================================================================
# HELPERS
# =============================================================================

# Publication theme
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Multi-format save (PDF + SVG + JPG in subfolders)
UTIL_PATH <- "../scripts/utils/multi_format_output.R"
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

# Format p-values for plot annotations
fmt_p <- function(p) {
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

# --- Enhancer classes ---
enh_list <- lapply(names(ENHANCER_FILES), function(cls) {
  df <- read.table(ENHANCER_FILES[[cls]], sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)
  df$enhancer_class <- cls
  df
})
enh_all <- do.call(rbind, enh_list)
enh_all$enhancer_class <- factor(enh_all$enhancer_class, levels = CLASS_ORDER)
cat(sprintf("  Enhancers: %d total across %d classes\n", nrow(enh_all),
            length(CLASS_ORDER)))
for (cls in CLASS_ORDER) {
  cat(sprintf("    %s: %d\n", cls, sum(enh_all$enhancer_class == cls)))
}

# --- Loops ---
loops <- read.table(LOOPS_FILE, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
cat(sprintf("  Loops: %d\n", nrow(loops)))

# --- Characterized loops (differential, with anchor annotations) ---
char_loops <- read.table(CHAR_LOOPS_FILE, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE)
cat(sprintf("  Characterized loops: %d\n", nrow(char_loops)))

# --- ABC pairs ---
abc <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d\n", nrow(abc)))

# --- ABC + RNA-seq ---
abc_rna <- read.table(ABC_RNASEQ_FILE, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE)
cat(sprintf("  ABC + RNA-seq pairs: %d\n", nrow(abc_rna)))

# --- Gene-level summary ---
gene_summary <- read.table(GENE_SUMMARY_FILE, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE)
cat(sprintf("  Gene summary: %d genes\n", nrow(gene_summary)))

# --- TSS ---
tss <- read.table(TSS_FILE, sep = "\t", header = TRUE, comment.char = "",
                  stringsAsFactors = FALSE)
# Fix header: first column is "#chr" -> "chr"
colnames(tss)[1] <- "chr"
cat(sprintf("  TSS regions: %d\n", nrow(tss)))

cat("\n")

# =============================================================================
# BUILD GENOMIC RANGES
# =============================================================================

cat("Building GRanges objects...\n")

# Enhancer GRanges (with class metadata)
enh_gr <- GRanges(
  seqnames = enh_all$chr,
  ranges = IRanges(start = enh_all$start, end = enh_all$end),
  enhancer_class = enh_all$enhancer_class,
  enh_idx = seq_len(nrow(enh_all))
)

# Loop anchor1 and anchor2 GRanges
loop_anchor1_gr <- GRanges(
  seqnames = loops$chr1,
  ranges = IRanges(start = loops$start1, end = loops$end1),
  loop_idx = seq_len(nrow(loops))
)
loop_anchor2_gr <- GRanges(
  seqnames = loops$chr2,
  ranges = IRanges(start = loops$start2, end = loops$end2),
  loop_idx = seq_len(nrow(loops))
)

# TSS GRanges
tss_gr <- GRanges(
  seqnames = tss$chr,
  ranges = IRanges(start = tss$start, end = tss$end),
  gene_name = tss$name
)

# ABC enhancer GRanges (unique enhancer positions)
abc_enh_gr <- GRanges(
  seqnames = abc$chr,
  ranges = IRanges(start = abc$start, end = abc$end)
)

cat("GRanges built.\n\n")


# #############################################################################
# PART A: LOOP ANCHORING (Metrics 1 + 2)
# #############################################################################

cat("=== PART A: Loop Anchoring ===\n\n")

# Find overlaps: enhancer ↔ loop anchors (either anchor)
hits_a1 <- findOverlaps(enh_gr, loop_anchor1_gr, type = "any")
hits_a2 <- findOverlaps(enh_gr, loop_anchor2_gr, type = "any")

cat(sprintf("  Enhancer-anchor1 overlaps: %d\n", length(hits_a1)))
cat(sprintf("  Enhancer-anchor2 overlaps: %d\n", length(hits_a2)))

# Combine: for each enhancer, collect all overlapping loop indices
enh_loop_map <- data.frame(
  enh_idx = c(queryHits(hits_a1), queryHits(hits_a2)),
  loop_idx = c(subjectHits(hits_a1), subjectHits(hits_a2))
)
# Deduplicate (an enhancer may overlap both anchors of the same loop)
enh_loop_map <- unique(enh_loop_map)

cat(sprintf("  Unique enhancer-loop pairs: %d\n", nrow(enh_loop_map)))
cat(sprintf("  Enhancers with >=1 loop: %d / %d (%.1f%%)\n",
            length(unique(enh_loop_map$enh_idx)), nrow(enh_all),
            100 * length(unique(enh_loop_map$enh_idx)) / nrow(enh_all)))

# Attach loop distance (midpoint-to-midpoint distance between anchors)
enh_loop_map$loop_distance <- abs(
  (loops$start1[enh_loop_map$loop_idx] + loops$end1[enh_loop_map$loop_idx]) / 2 -
  (loops$start2[enh_loop_map$loop_idx] + loops$end2[enh_loop_map$loop_idx]) / 2
)
enh_loop_map$loop_logFC <- loops$logFC[enh_loop_map$loop_idx]
enh_loop_map$loop_FDR <- loops$FDR[enh_loop_map$loop_idx]
enh_loop_map$loop_direction <- loops$direction[enh_loop_map$loop_idx]
enh_loop_map$enhancer_class <- enh_all$enhancer_class[enh_loop_map$enh_idx]

# Per-enhancer: count loops and mean distance
enh_loop_agg <- aggregate(
  loop_idx ~ enh_idx,
  data = enh_loop_map,
  FUN = function(x) length(unique(x))
)
colnames(enh_loop_agg)[2] <- "n_loops"

enh_dist_agg <- aggregate(
  loop_distance ~ enh_idx,
  data = enh_loop_map,
  FUN = mean
)
colnames(enh_dist_agg)[2] <- "mean_loop_distance"

# Merge per-enhancer metrics
enh_metrics <- data.frame(
  enh_idx = seq_len(nrow(enh_all)),
  chr = enh_all$chr,
  start = enh_all$start,
  end = enh_all$end,
  enhancer_class = enh_all$enhancer_class
)
enh_metrics <- merge(enh_metrics, enh_loop_agg, by = "enh_idx", all.x = TRUE)
enh_metrics <- merge(enh_metrics, enh_dist_agg, by = "enh_idx", all.x = TRUE)
enh_metrics$n_loops[is.na(enh_metrics$n_loops)] <- 0

# Per-class summary
cat("\n  Per-class loop summary:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_metrics[enh_metrics$enhancer_class == cls, ]
  cat(sprintf("    %s (n=%d): median loops=%d, mean loops=%.2f, zero-loop=%.1f%%\n",
              cls, nrow(sub), median(sub$n_loops), mean(sub$n_loops),
              100 * mean(sub$n_loops == 0)))
}

# --- Stats: Kruskal-Wallis ---
kw_loops <- kruskal.test(n_loops ~ enhancer_class, data = enh_metrics)
cat(sprintf("\n  Kruskal-Wallis (loops per enhancer): chi2=%.1f, p=%s\n",
            kw_loops$statistic, fmt_p(kw_loops$p.value)))

# Pairwise Wilcoxon for loop counts
pw_loops <- pairwise.wilcox.test(enh_metrics$n_loops, enh_metrics$enhancer_class,
                                 p.adjust.method = "holm")
cat("  Pairwise Wilcoxon (loops per enhancer):\n")
print(pw_loops$p.value)

# Distance stats (only enhancers with loops)
enh_with_loops <- enh_metrics[enh_metrics$n_loops > 0, ]
kw_dist <- kruskal.test(mean_loop_distance ~ enhancer_class, data = enh_with_loops)
cat(sprintf("\n  Kruskal-Wallis (loop distance): chi2=%.1f, p=%s\n",
            kw_dist$statistic, fmt_p(kw_dist$p.value)))

# --- Plot 01: Boxplot of loops per enhancer ---
cat("\nGenerating Part A plots...\n")
p01 <- ggplot(enh_metrics, aes(x = enhancer_class, y = n_loops,
                                fill = enhancer_class)) +
  geom_boxplot(outlier.size = 0.3, notch = TRUE) +
  scale_fill_manual(values = CLASS_COLORS) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 2, 5, 10, 20, 50)) +
  labs(
    x = "Enhancer class", y = "Loops per enhancer (pseudo-log scale)",
    title = "Loop anchoring by enhancer class",
    subtitle = sprintf("Kruskal-Wallis %s", fmt_p(kw_loops$p.value))
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p01, "01_loops_per_enhancer_boxplot")

# --- Plot 02: Violin of loop distance ---
p02 <- ggplot(enh_with_loops, aes(x = enhancer_class, y = mean_loop_distance / 1e6,
                                   fill = enhancer_class)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  scale_fill_manual(values = CLASS_COLORS) +
  labs(
    x = "Enhancer class", y = "Mean loop distance (Mb)",
    title = "Loop distance by enhancer class",
    subtitle = sprintf("Kruskal-Wallis %s", fmt_p(kw_dist$p.value))
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p02, "02_loop_distance_violin")


# #############################################################################
# PART B: PROMOTER LOOPS & GENE EXPRESSION (Metric 3)
# #############################################################################

cat("\n=== PART B: Promoter Loops & Gene Expression ===\n\n")

# For each enhancer-loop pair, check whether the OTHER anchor overlaps a TSS
# First: identify which anchor the enhancer overlaps
hits_a1_df <- data.frame(enh_idx = queryHits(hits_a1), loop_idx = subjectHits(hits_a1),
                         enh_on_anchor = 1L)
hits_a2_df <- data.frame(enh_idx = queryHits(hits_a2), loop_idx = subjectHits(hits_a2),
                         enh_on_anchor = 2L)
enh_anchor_info <- rbind(hits_a1_df, hits_a2_df)

# The "other" anchor for each row
# If enhancer on anchor1, check anchor2 for TSS (and vice versa)
other_anchor_gr <- ifelse(enh_anchor_info$enh_on_anchor == 1L, 2L, 1L)

# Build GRanges for "other" anchors
other_chr <- ifelse(other_anchor_gr == 1L,
                    loops$chr1[enh_anchor_info$loop_idx],
                    loops$chr2[enh_anchor_info$loop_idx])
other_start <- ifelse(other_anchor_gr == 1L,
                      loops$start1[enh_anchor_info$loop_idx],
                      loops$start2[enh_anchor_info$loop_idx])
other_end <- ifelse(other_anchor_gr == 1L,
                    loops$end1[enh_anchor_info$loop_idx],
                    loops$end2[enh_anchor_info$loop_idx])
other_gr <- GRanges(seqnames = other_chr,
                    ranges = IRanges(start = other_start, end = other_end))

# Overlap other anchor with TSS
tss_hits <- findOverlaps(other_gr, tss_gr, type = "any")
cat(sprintf("  Other-anchor ↔ TSS overlaps: %d\n", length(tss_hits)))

# Mark promoter-connected rows
enh_anchor_info$is_promoter_loop <- FALSE
enh_anchor_info$gene_name <- NA_character_
promoter_rows <- queryHits(tss_hits)
enh_anchor_info$is_promoter_loop[promoter_rows] <- TRUE
# Take first overlapping gene if multiple
first_gene <- tapply(subjectHits(tss_hits), queryHits(tss_hits), function(x) x[1])
enh_anchor_info$gene_name[as.integer(names(first_gene))] <-
  tss$name[as.integer(first_gene)]

# Attach enhancer class
enh_anchor_info$enhancer_class <- enh_all$enhancer_class[enh_anchor_info$enh_idx]

# Deduplicate to unique enhancer-loop-gene tuples
promoter_links <- enh_anchor_info[enh_anchor_info$is_promoter_loop, ]
promoter_links <- unique(promoter_links[, c("enh_idx", "loop_idx",
                                             "enhancer_class", "gene_name")])
cat(sprintf("  Unique enhancer-loop-gene promoter links: %d\n", nrow(promoter_links)))

# Count promoter loops per enhancer
prom_count <- aggregate(loop_idx ~ enh_idx, data = promoter_links,
                        FUN = function(x) length(unique(x)))
colnames(prom_count)[2] <- "n_promoter_loops"
enh_metrics <- merge(enh_metrics, prom_count, by = "enh_idx", all.x = TRUE)
enh_metrics$n_promoter_loops[is.na(enh_metrics$n_promoter_loops)] <- 0

# Promoter loop proportion per class
cat("\n  Promoter loop proportions:\n")
prom_prop_stats <- list()
for (cls in CLASS_ORDER) {
  sub <- enh_metrics[enh_metrics$enhancer_class == cls, ]
  n_with_prom <- sum(sub$n_promoter_loops > 0)
  n_total <- nrow(sub)
  cat(sprintf("    %s: %d / %d (%.1f%%) have promoter loops\n",
              cls, n_with_prom, n_total, 100 * n_with_prom / n_total))
  prom_prop_stats[[cls]] <- c(with_prom = n_with_prom, without_prom = n_total - n_with_prom)
}

# Fisher's exact: each class vs Stable
cat("\n  Fisher's exact (vs Stable):\n")
stable_vec <- prom_prop_stats[["Stable"]]
for (cls in setdiff(CLASS_ORDER, "Stable")) {
  cls_vec <- prom_prop_stats[[cls]]
  mat <- matrix(c(cls_vec, stable_vec), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat)
  cat(sprintf("    %s vs Stable: OR=%.2f, %s\n", cls, ft$estimate, fmt_p(ft$p.value)))
}

# Merge promoter-linked genes with RNA-seq DE
promoter_gene_de <- merge(
  promoter_links[, c("enh_idx", "enhancer_class", "gene_name")],
  gene_summary[, c("TargetGene", "log2FC", "padj")],
  by.x = "gene_name", by.y = "TargetGene"
)
# Deduplicate to unique enhancer-gene pairs
promoter_gene_de <- unique(promoter_gene_de[, c("enh_idx", "enhancer_class",
                                                  "gene_name", "log2FC", "padj")])
cat(sprintf("\n  Promoter-linked genes with DE data: %d unique enhancer-gene pairs\n",
            nrow(promoter_gene_de)))

# Gene logFC by class
cat("  Gene logFC by class (promoter-connected):\n")
for (cls in CLASS_ORDER) {
  sub <- promoter_gene_de[promoter_gene_de$enhancer_class == cls &
                            !is.na(promoter_gene_de$log2FC), ]
  if (nrow(sub) > 0) {
    cat(sprintf("    %s: n=%d, median logFC=%.3f, mean logFC=%.3f\n",
                cls, nrow(sub), median(sub$log2FC), mean(sub$log2FC)))
  } else {
    cat(sprintf("    %s: no genes with DE data\n", cls))
  }
}

# Wilcoxon: K119ub_Only vs Stable gene logFC
k_genes <- promoter_gene_de$log2FC[promoter_gene_de$enhancer_class == "K119ub_Only" &
                                     !is.na(promoter_gene_de$log2FC)]
s_genes <- promoter_gene_de$log2FC[promoter_gene_de$enhancer_class == "Stable" &
                                     !is.na(promoter_gene_de$log2FC)]
if (length(k_genes) >= 5 && length(s_genes) >= 5) {
  wt_gene <- wilcox.test(k_genes, s_genes)
  cat(sprintf("\n  Wilcoxon K119ub_Only vs Stable gene logFC: %s\n",
              fmt_p(wt_gene$p.value)))
}

# --- Plot 03: Stacked bar of promoter vs non-promoter ---
prom_bar_data <- data.frame(
  enhancer_class = rep(CLASS_ORDER, each = 2),
  loop_type = rep(c("Promoter", "Non-promoter"), times = length(CLASS_ORDER)),
  count = NA_real_
)
for (i in seq_along(CLASS_ORDER)) {
  cls <- CLASS_ORDER[i]
  sub <- enh_metrics[enh_metrics$enhancer_class == cls & enh_metrics$n_loops > 0, ]
  n_prom <- sum(sub$n_promoter_loops > 0)
  n_nonprom <- sum(sub$n_promoter_loops == 0)
  prom_bar_data$count[(i-1)*2 + 1] <- n_prom
  prom_bar_data$count[(i-1)*2 + 2] <- n_nonprom
}
prom_bar_data$enhancer_class <- factor(prom_bar_data$enhancer_class, levels = CLASS_ORDER)

p03 <- ggplot(prom_bar_data, aes(x = enhancer_class, y = count, fill = loop_type)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c(Promoter = "#E41A1C", `Non-promoter` = "grey70")) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Enhancer class", y = "Proportion of loop-anchored enhancers",
    fill = "Loop type",
    title = "Promoter vs non-promoter loops by enhancer class"
  ) +
  theme_pub
save_plot(p03, "03_promoter_loop_proportion")

# --- Plot 04: Boxplot of gene logFC by class ---
plot_gene_de <- promoter_gene_de[!is.na(promoter_gene_de$log2FC), ]
plot_gene_de$enhancer_class <- factor(plot_gene_de$enhancer_class, levels = CLASS_ORDER)

if (nrow(plot_gene_de) > 0) {
  p04 <- ggplot(plot_gene_de, aes(x = enhancer_class, y = log2FC,
                                   fill = enhancer_class)) +
    geom_boxplot(outlier.size = 0.3, notch = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = CLASS_COLORS) +
    labs(
      x = "Enhancer class", y = "RNA-seq log2FC (KO/WT)",
      title = "Gene expression change at promoter-connected genes",
      subtitle = "Genes connected to enhancers via loops with TSS overlap"
    ) +
    theme_pub +
    theme(legend.position = "none")
  save_plot(p04, "04_gene_logfc_by_class")
}


# #############################################################################
# PART C: LOOP STRENGTH logFC (Metric 4)
# #############################################################################

cat("\n=== PART C: Loop Strength logFC ===\n\n")

# Use enh_loop_map which has loop_logFC and enhancer_class
loop_logfc_data <- enh_loop_map[!is.na(enh_loop_map$loop_logFC), ]
loop_logfc_data$enhancer_class <- factor(loop_logfc_data$enhancer_class,
                                          levels = CLASS_ORDER)

cat("  Loop logFC per class:\n")
for (cls in CLASS_ORDER) {
  sub <- loop_logfc_data[loop_logfc_data$enhancer_class == cls, ]
  sig <- sub[!is.na(sub$loop_FDR) & sub$loop_FDR < 0.05, ]
  cat(sprintf("    %s: n=%d loops, median logFC=%.4f, sig (FDR<0.05)=%d (%.1f%%)\n",
              cls, nrow(sub), median(sub$loop_logFC), nrow(sig),
              100 * nrow(sig) / max(nrow(sub), 1)))
}

# KEY HYPOTHESIS TEST: K119ub_Only median logFC != 0
k_logfc <- loop_logfc_data$loop_logFC[loop_logfc_data$enhancer_class == "K119ub_Only"]
if (length(k_logfc) >= 5) {
  wt_k119 <- wilcox.test(k_logfc, mu = 0)
  cat(sprintf("\n  One-sample Wilcoxon K119ub_Only logFC vs 0: %s (median=%.4f)\n",
              fmt_p(wt_k119$p.value), median(k_logfc)))
  # Effect size: rank-biserial r = Z / sqrt(N)
  z_k119 <- qnorm(wt_k119$p.value / 2) * sign(median(k_logfc))
  r_k119 <- z_k119 / sqrt(length(k_logfc))
  cat(sprintf("  Effect size (rank-biserial r): %.3f\n", r_k119))
}

# K119ub_Only vs Stable
s_logfc <- loop_logfc_data$loop_logFC[loop_logfc_data$enhancer_class == "Stable"]
if (length(k_logfc) >= 5 && length(s_logfc) >= 5) {
  wt_ks <- wilcox.test(k_logfc, s_logfc)
  cat(sprintf("  Wilcoxon K119ub_Only vs Stable: %s\n", fmt_p(wt_ks$p.value)))
}

# Kruskal-Wallis across all 4
kw_logfc <- kruskal.test(loop_logFC ~ enhancer_class, data = loop_logfc_data)
cat(sprintf("  Kruskal-Wallis (loop logFC): chi2=%.1f, %s\n",
            kw_logfc$statistic, fmt_p(kw_logfc$p.value)))

pw_logfc <- pairwise.wilcox.test(loop_logfc_data$loop_logFC,
                                  loop_logfc_data$enhancer_class,
                                  p.adjust.method = "holm")
cat("  Pairwise Wilcoxon (loop logFC):\n")
print(pw_logfc$p.value)

# --- Plot 05: Violin + box of loop logFC ---
n_per_class <- as.data.frame(table(loop_logfc_data$enhancer_class))
colnames(n_per_class) <- c("enhancer_class", "n")
p05 <- ggplot(loop_logfc_data, aes(x = enhancer_class, y = loop_logFC,
                                    fill = enhancer_class)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(data = n_per_class, aes(x = enhancer_class, y = Inf,
            label = paste0("n=", format(n, big.mark = ","))),
            vjust = 1.5, size = 3.2, inherit.aes = FALSE) +
  scale_fill_manual(values = CLASS_COLORS) +
  labs(
    x = "Enhancer class", y = "Loop logFC (KO/WT)",
    title = "Loop strength change by enhancer class",
    subtitle = sprintf("Kruskal-Wallis %s", fmt_p(kw_logfc$p.value))
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p05, "05_loop_logfc_violin")

# --- Plot 06: Overlaid density curves ---
p06 <- ggplot(loop_logfc_data, aes(x = loop_logFC, color = enhancer_class,
                                    fill = enhancer_class)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = CLASS_COLORS) +
  scale_fill_manual(values = CLASS_COLORS) +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  labs(
    x = "Loop logFC (KO/WT)", y = "Density",
    color = "Enhancer class", fill = "Enhancer class",
    title = "Loop logFC density by enhancer class"
  ) +
  theme_pub
save_plot(p06, "06_loop_logfc_density")


# #############################################################################
# PART D: ABC SCORES (Metric 5)
# #############################################################################

cat("\n=== PART D: ABC Scores ===\n\n")

# Join enhancer classes → ABC pairs via findOverlaps
hits_abc <- findOverlaps(enh_gr, abc_enh_gr, type = "any")
cat(sprintf("  Enhancer ↔ ABC overlaps: %d\n", length(hits_abc)))

# Build matched dataframe
abc_matched <- data.frame(
  enh_idx = queryHits(hits_abc),
  abc_row = subjectHits(hits_abc)
)
abc_matched$enhancer_class <- enh_all$enhancer_class[abc_matched$enh_idx]

# Exclude self-promoters (class == "promoter" in ABC data)
abc_matched$abc_class <- abc$class[abc_matched$abc_row]
n_before <- nrow(abc_matched)
abc_matched <- abc_matched[abc_matched$abc_class != "promoter", ]
cat(sprintf("  After removing self-promoters: %d (removed %d)\n",
            nrow(abc_matched), n_before - nrow(abc_matched)))

# Attach ABC metrics
abc_matched$delta_ABC <- abc$delta_ABC[abc_matched$abc_row]
abc_matched$delta_unnorm <- abc$delta_unnorm[abc_matched$abc_row]
abc_matched$ABC_WT <- abc$ABC.Score_WT[abc_matched$abc_row]
abc_matched$ABC_KO <- abc$ABC.Score_KO[abc_matched$abc_row]
abc_matched$hic_contact_WT <- abc$hic_contact_WT[abc_matched$abc_row]
abc_matched$hic_contact_KO <- abc$hic_contact_KO[abc_matched$abc_row]
abc_matched$activity_WT <- abc$activity_base_WT[abc_matched$abc_row]
abc_matched$activity_KO <- abc$activity_base_KO[abc_matched$abc_row]
abc_matched$distance <- abc$distance[abc_matched$abc_row]
abc_matched$target_gene <- abc$TargetGene[abc_matched$abc_row]

# Classify ABC pairs as gained/lost/unchanged
abc_matched$abc_direction <- "unchanged"
abc_matched$abc_direction[abc_matched$delta_ABC > 0.01] <- "gained"
abc_matched$abc_direction[abc_matched$delta_ABC < -0.01] <- "lost"

# Match rate
n_enh_with_abc <- length(unique(abc_matched$enh_idx))
cat(sprintf("  Enhancers matched to ABC: %d / %d (%.1f%%)\n",
            n_enh_with_abc, nrow(enh_all),
            100 * n_enh_with_abc / nrow(enh_all)))

# Per-enhancer ABC aggregation
abc_enh_agg <- aggregate(
  cbind(delta_ABC, delta_unnorm) ~ enh_idx,
  data = abc_matched,
  FUN = mean
)
colnames(abc_enh_agg)[2:3] <- c("mean_delta_abc", "mean_delta_unnorm")

abc_count_agg <- aggregate(abc_row ~ enh_idx, data = abc_matched,
                           FUN = length)
colnames(abc_count_agg)[2] <- "n_abc_pairs"

abc_gained <- aggregate(abc_direction ~ enh_idx, data = abc_matched,
                        FUN = function(x) sum(x == "gained"))
colnames(abc_gained)[2] <- "n_gained"

abc_lost <- aggregate(abc_direction ~ enh_idx, data = abc_matched,
                      FUN = function(x) sum(x == "lost"))
colnames(abc_lost)[2] <- "n_lost"

# Merge all ABC aggregates
enh_abc_metrics <- merge(abc_enh_agg, abc_count_agg, by = "enh_idx")
enh_abc_metrics <- merge(enh_abc_metrics, abc_gained, by = "enh_idx")
enh_abc_metrics <- merge(enh_abc_metrics, abc_lost, by = "enh_idx")
enh_abc_metrics$enhancer_class <- enh_all$enhancer_class[enh_abc_metrics$enh_idx]
enh_abc_metrics$chr <- enh_all$chr[enh_abc_metrics$enh_idx]
enh_abc_metrics$start <- enh_all$start[enh_abc_metrics$enh_idx]
enh_abc_metrics$end <- enh_all$end[enh_abc_metrics$enh_idx]

cat("\n  Per-class ABC summary:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_abc_metrics[enh_abc_metrics$enhancer_class == cls, ]
  cat(sprintf("    %s: n=%d enhancers, median delta_ABC=%.5f, median delta_unnorm=%.5f\n",
              cls, nrow(sub), median(sub$mean_delta_abc), median(sub$mean_delta_unnorm)))
  cat(sprintf("      mean n_abc_pairs=%.1f, mean gained=%.2f, mean lost=%.2f\n",
              mean(sub$n_abc_pairs), mean(sub$n_gained), mean(sub$n_lost)))
}

# Stats
kw_abc <- kruskal.test(mean_delta_abc ~ enhancer_class, data = enh_abc_metrics)
cat(sprintf("\n  Kruskal-Wallis (mean delta_ABC): chi2=%.1f, %s\n",
            kw_abc$statistic, fmt_p(kw_abc$p.value)))

kw_unnorm <- kruskal.test(mean_delta_unnorm ~ enhancer_class, data = enh_abc_metrics)
cat(sprintf("  Kruskal-Wallis (mean delta_unnorm): chi2=%.1f, %s\n",
            kw_unnorm$statistic, fmt_p(kw_unnorm$p.value)))

pw_abc <- pairwise.wilcox.test(enh_abc_metrics$mean_delta_abc,
                                enh_abc_metrics$enhancer_class,
                                p.adjust.method = "holm")
cat("  Pairwise Wilcoxon (mean delta_ABC):\n")
print(pw_abc$p.value)

pw_unnorm <- pairwise.wilcox.test(enh_abc_metrics$mean_delta_unnorm,
                                   enh_abc_metrics$enhancer_class,
                                   p.adjust.method = "holm")
cat("  Pairwise Wilcoxon (mean delta_unnorm):\n")
print(pw_unnorm$p.value)

# --- Plot 07: Boxplot of mean delta_ABC ---
enh_abc_metrics$enhancer_class <- factor(enh_abc_metrics$enhancer_class,
                                          levels = CLASS_ORDER)
p07 <- ggplot(enh_abc_metrics, aes(x = enhancer_class, y = mean_delta_abc,
                                    fill = enhancer_class)) +
  geom_boxplot(outlier.size = 0.3, notch = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = CLASS_COLORS) +
  labs(
    x = "Enhancer class", y = "Mean delta-ABC per enhancer",
    title = "ABC score change by enhancer class",
    subtitle = sprintf("Kruskal-Wallis %s", fmt_p(kw_abc$p.value))
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p07, "07_delta_abc_boxplot")

# --- Plot 08: Boxplot of mean delta_unnorm ---
p08 <- ggplot(enh_abc_metrics, aes(x = enhancer_class, y = mean_delta_unnorm,
                                    fill = enhancer_class)) +
  geom_boxplot(outlier.size = 0.3, notch = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = CLASS_COLORS) +
  labs(
    x = "Enhancer class", y = "Mean delta(Activity x Contact) per enhancer",
    title = "Unnormalized ABC change by enhancer class",
    subtitle = sprintf("Kruskal-Wallis %s", fmt_p(kw_unnorm$p.value))
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p08, "08_delta_unnorm_boxplot")

# --- Plot 09: Stacked bar of gained/lost/unchanged ---
abc_dir_summary <- data.frame(
  enhancer_class = character(),
  direction = character(),
  count = integer(),
  stringsAsFactors = FALSE
)
for (cls in CLASS_ORDER) {
  sub <- abc_matched[abc_matched$enhancer_class == cls, ]
  for (d in c("gained", "unchanged", "lost")) {
    abc_dir_summary <- rbind(abc_dir_summary, data.frame(
      enhancer_class = cls,
      direction = d,
      count = sum(sub$abc_direction == d)
    ))
  }
}
abc_dir_summary$enhancer_class <- factor(abc_dir_summary$enhancer_class,
                                          levels = CLASS_ORDER)
abc_dir_summary$direction <- factor(abc_dir_summary$direction,
                                     levels = c("gained", "unchanged", "lost"))

p09 <- ggplot(abc_dir_summary, aes(x = enhancer_class, y = count,
                                    fill = direction)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c(gained = "#B2182B", unchanged = "grey70",
                                lost = "#2166AC")) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Enhancer class", y = "Proportion of ABC pairs",
    fill = "ABC change\n(|dABC| > 0.01)",
    title = "ABC pair directionality by enhancer class"
  ) +
  theme_pub
save_plot(p09, "09_abc_direction_stacked_bar")


# #############################################################################
# PART E: CONTACT DECAY PROFILE (Metric 6 — bias check)
# #############################################################################

cat("\n=== PART E: Contact Decay Profile ===\n\n")

# Use the abc_matched dataframe which has distance and contact values
contact_data <- abc_matched[, c("enhancer_class", "distance",
                                 "hic_contact_WT", "hic_contact_KO")]
contact_data <- contact_data[!is.na(contact_data$distance), ]

# Bin distance
contact_data$dist_bin <- cut(
  contact_data$distance,
  breaks = c(0, 50e3, 100e3, 200e3, 500e3, 1e6, 5e6),
  labels = c("<50kb", "50-100kb", "100-200kb", "200-500kb", "500kb-1Mb", "1-5Mb"),
  include.lowest = TRUE
)
contact_data <- contact_data[!is.na(contact_data$dist_bin), ]
contact_data$delta_contact <- contact_data$hic_contact_KO - contact_data$hic_contact_WT

# Aggregate per bin × class
decay_agg <- aggregate(
  cbind(hic_contact_WT, hic_contact_KO, delta_contact) ~ dist_bin + enhancer_class,
  data = contact_data,
  FUN = mean
)
decay_agg$enhancer_class <- factor(decay_agg$enhancer_class, levels = CLASS_ORDER)

cat("  Contact decay summary:\n")
print(decay_agg[order(decay_agg$enhancer_class, decay_agg$dist_bin), ], row.names = FALSE)

# --- Plot 10: Line plot of mean WT contact by distance ---
p10 <- ggplot(decay_agg, aes(x = dist_bin, y = hic_contact_WT,
                              color = enhancer_class, group = enhancer_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = CLASS_COLORS) +
  scale_y_log10(labels = comma_format()) +
  labs(
    x = "Genomic distance", y = "Mean WT Hi-C contact (log scale)",
    color = "Enhancer class",
    title = "Contact frequency decay by enhancer class (WT)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p10, "10_contact_decay_wt")

# --- Plot 11: Line plot of mean delta_contact by distance ---
p11 <- ggplot(decay_agg, aes(x = dist_bin, y = delta_contact,
                              color = enhancer_class, group = enhancer_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = CLASS_COLORS) +
  labs(
    x = "Genomic distance", y = "Mean delta Hi-C contact (KO - WT)",
    color = "Enhancer class",
    title = "Contact change by distance and enhancer class"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p11, "11_contact_decay_delta")


# #############################################################################
# PART F: SUMMARY INTEGRATION
# #############################################################################

cat("\n=== PART F: Summary Integration ===\n\n")

# --- Class-level summary table ---
class_summary <- data.frame(enhancer_class = CLASS_ORDER)
for (i in seq_along(CLASS_ORDER)) {
  cls <- CLASS_ORDER[i]

  # Enhancer counts
  sub_enh <- enh_metrics[enh_metrics$enhancer_class == cls, ]
  class_summary$n_enhancers[i] <- nrow(sub_enh)

  # Loop metrics
  class_summary$median_loops[i] <- median(sub_enh$n_loops)
  class_summary$mean_loops[i] <- round(mean(sub_enh$n_loops), 2)
  class_summary$frac_with_loops[i] <- round(mean(sub_enh$n_loops > 0), 3)

  # Loop distance (only enhancers with loops)
  sub_dist <- sub_enh[sub_enh$n_loops > 0, ]
  class_summary$median_loop_distance[i] <- if (nrow(sub_dist) > 0) {
    round(median(sub_dist$mean_loop_distance, na.rm = TRUE))
  } else { NA_real_ }

  # Promoter loops
  class_summary$frac_promoter_loops[i] <- round(mean(sub_enh$n_promoter_loops > 0), 3)

  # Loop logFC
  sub_logfc <- loop_logfc_data[loop_logfc_data$enhancer_class == cls, ]
  class_summary$median_loop_logfc[i] <- if (nrow(sub_logfc) > 0) {
    round(median(sub_logfc$loop_logFC), 4)
  } else { NA_real_ }
  class_summary$frac_sig_loops[i] <- if (nrow(sub_logfc) > 0) {
    round(mean(!is.na(sub_logfc$loop_FDR) & sub_logfc$loop_FDR < 0.05), 3)
  } else { NA_real_ }

  # ABC metrics
  sub_abc <- enh_abc_metrics[enh_abc_metrics$enhancer_class == cls, ]
  class_summary$n_with_abc[i] <- nrow(sub_abc)
  class_summary$median_delta_abc[i] <- if (nrow(sub_abc) > 0) {
    round(median(sub_abc$mean_delta_abc), 5)
  } else { NA_real_ }
  class_summary$median_delta_unnorm[i] <- if (nrow(sub_abc) > 0) {
    round(median(sub_abc$mean_delta_unnorm), 5)
  } else { NA_real_ }
}

cat("Class-level summary:\n")
print(class_summary, row.names = FALSE)

# --- Collect all statistical tests ---
stat_tests <- data.frame(
  test_name = character(),
  statistic = numeric(),
  p_value = numeric(),
  comparison = character(),
  stringsAsFactors = FALSE
)

# Add KW tests
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "KW_loops_per_enhancer", statistic = kw_loops$statistic,
  p_value = kw_loops$p.value, comparison = "all_4_classes"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "KW_loop_distance", statistic = kw_dist$statistic,
  p_value = kw_dist$p.value, comparison = "all_4_classes"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "KW_loop_logFC", statistic = kw_logfc$statistic,
  p_value = kw_logfc$p.value, comparison = "all_4_classes"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "KW_delta_ABC", statistic = kw_abc$statistic,
  p_value = kw_abc$p.value, comparison = "all_4_classes"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "KW_delta_unnorm", statistic = kw_unnorm$statistic,
  p_value = kw_unnorm$p.value, comparison = "all_4_classes"))

# Add key hypothesis test
if (exists("wt_k119")) {
  stat_tests <- rbind(stat_tests, data.frame(
    test_name = "Wilcoxon_K119ub_logFC_vs_zero", statistic = wt_k119$statistic,
    p_value = wt_k119$p.value, comparison = "K119ub_Only_vs_0"))
}
if (exists("wt_ks")) {
  stat_tests <- rbind(stat_tests, data.frame(
    test_name = "Wilcoxon_K119ub_vs_Stable_logFC", statistic = wt_ks$statistic,
    p_value = wt_ks$p.value, comparison = "K119ub_Only_vs_Stable"))
}

# Add pairwise Wilcoxon from loop logFC
for (r in seq_len(nrow(pw_logfc$p.value))) {
  for (c in seq_len(ncol(pw_logfc$p.value))) {
    pv <- pw_logfc$p.value[r, c]
    if (!is.na(pv)) {
      stat_tests <- rbind(stat_tests, data.frame(
        test_name = "PW_Wilcoxon_loop_logFC",
        statistic = NA_real_,
        p_value = pv,
        comparison = paste(rownames(pw_logfc$p.value)[r], "vs",
                           colnames(pw_logfc$p.value)[c])
      ))
    }
  }
}

# Add pairwise Wilcoxon from ABC
for (r in seq_len(nrow(pw_abc$p.value))) {
  for (c in seq_len(ncol(pw_abc$p.value))) {
    pv <- pw_abc$p.value[r, c]
    if (!is.na(pv)) {
      stat_tests <- rbind(stat_tests, data.frame(
        test_name = "PW_Wilcoxon_delta_ABC",
        statistic = NA_real_,
        p_value = pv,
        comparison = paste(rownames(pw_abc$p.value)[r], "vs",
                           colnames(pw_abc$p.value)[c])
      ))
    }
  }
}

cat(sprintf("\n  Total statistical tests collected: %d\n", nrow(stat_tests)))

# --- Plot 12: Multi-panel patchwork ---
cat("\nGenerating summary patchwork...\n")
p12 <- (p05 + p08) / (p06 + p11) +
  plot_annotation(
    title = "Enhancer class stratified analysis: key metrics",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
save_plot(p12, "12_summary_patchwork", w = 14, h = 12)


# #############################################################################
# PART G: K119ub_Only STRATIFIED BY K119ub MAGNITUDE (Plots 13-15)
# #############################################################################

cat("\n=== PART G: K119ub_Only Stratified by K119ub Magnitude ===\n\n")

TERTILE_COLORS <- c(T1_low = "#FDBB84", T2_mid = "#E34A33", T3_high = "#7F0000")

# K119ub_Only enhancer indices in enh_all
k119_idx <- which(enh_all$enhancer_class == "K119ub_Only")
k119_enh <- enh_all[k119_idx, ]
cat(sprintf("  K119ub_Only enhancers: %d\n", nrow(k119_enh)))

# Assign tertiles by log2FC_K119ub
breaks <- quantile(k119_enh$log2FC_K119ub, probs = c(0, 1/3, 2/3, 1))
k119_enh$tertile <- cut(k119_enh$log2FC_K119ub,
                         breaks = breaks,
                         labels = c("T1_low", "T2_mid", "T3_high"),
                         include.lowest = TRUE)

cat("  Tertile breakpoints:\n")
cat(sprintf("    T1_low:  [%.3f, %.3f]\n", breaks[1], breaks[2]))
cat(sprintf("    T2_mid:  (%.3f, %.3f]\n", breaks[2], breaks[3]))
cat(sprintf("    T3_high: (%.3f, %.3f]\n", breaks[3], breaks[4]))
for (trt in c("T1_low", "T2_mid", "T3_high")) {
  cat(sprintf("    %s: n=%d\n", trt, sum(k119_enh$tertile == trt)))
}

# Build tertile lookup: enh_idx -> tertile + log2FC_K119ub
k119_tertile <- data.frame(
  enh_idx = k119_idx,
  tertile = k119_enh$tertile,
  log2FC_K119ub = k119_enh$log2FC_K119ub
)

# Merge tertile -> enh_loop_map (for loop logFC)
k119_loop <- merge(
  enh_loop_map[enh_loop_map$enhancer_class == "K119ub_Only", ],
  k119_tertile[, c("enh_idx", "tertile", "log2FC_K119ub")],
  by = "enh_idx"
)
cat(sprintf("\n  K119ub_Only enhancer-loop pairs with tertile: %d\n", nrow(k119_loop)))

# Merge tertile -> enh_abc_metrics (for delta_ABC)
k119_abc <- merge(
  enh_abc_metrics[enh_abc_metrics$enhancer_class == "K119ub_Only", ],
  k119_tertile[, c("enh_idx", "tertile", "log2FC_K119ub")],
  by = "enh_idx"
)
cat(sprintf("  K119ub_Only enhancers with ABC + tertile: %d\n", nrow(k119_abc)))

# Tertile rank for trend tests
k119_loop$tertile_rank <- as.integer(factor(k119_loop$tertile,
                                             levels = c("T1_low", "T2_mid", "T3_high")))
k119_abc$tertile_rank <- as.integer(factor(k119_abc$tertile,
                                            levels = c("T1_low", "T2_mid", "T3_high")))

# Spearman: tertile rank vs metrics
cor_tert_loop <- cor.test(k119_loop$tertile_rank, k119_loop$loop_logFC,
                           method = "spearman")
cat(sprintf("\n  Spearman (tertile rank vs loop logFC): rho=%.4f, %s\n",
            cor_tert_loop$estimate, fmt_p(cor_tert_loop$p.value)))

cor_tert_abc <- cor.test(k119_abc$tertile_rank, k119_abc$mean_delta_abc,
                          method = "spearman")
cat(sprintf("  Spearman (tertile rank vs mean delta_ABC): rho=%.4f, %s\n",
            cor_tert_abc$estimate, fmt_p(cor_tert_abc$p.value)))

# Continuous Spearman for scatter plot
cor_cont <- cor.test(k119_loop$log2FC_K119ub, k119_loop$loop_logFC,
                      method = "spearman")
cat(sprintf("  Spearman (continuous K119ub vs loop logFC): rho=%.4f, %s\n",
            cor_cont$estimate, fmt_p(cor_cont$p.value)))

# Kruskal-Wallis across 3 tertiles
kw_tert_loop <- kruskal.test(loop_logFC ~ tertile, data = k119_loop)
cat(sprintf("\n  Kruskal-Wallis (loop logFC by tertile): chi2=%.1f, %s\n",
            kw_tert_loop$statistic, fmt_p(kw_tert_loop$p.value)))

kw_tert_abc <- kruskal.test(mean_delta_abc ~ tertile, data = k119_abc)
cat(sprintf("  Kruskal-Wallis (mean delta_ABC by tertile): chi2=%.1f, %s\n",
            kw_tert_abc$statistic, fmt_p(kw_tert_abc$p.value)))

# Wilcoxon T1 vs T3
wt_t1_t3_loop <- wilcox.test(
  k119_loop$loop_logFC[k119_loop$tertile == "T1_low"],
  k119_loop$loop_logFC[k119_loop$tertile == "T3_high"]
)
cat(sprintf("  Wilcoxon T1_low vs T3_high (loop logFC): %s\n",
            fmt_p(wt_t1_t3_loop$p.value)))

wt_t1_t3_abc <- wilcox.test(
  k119_abc$mean_delta_abc[k119_abc$tertile == "T1_low"],
  k119_abc$mean_delta_abc[k119_abc$tertile == "T3_high"]
)
cat(sprintf("  Wilcoxon T1_low vs T3_high (mean delta_ABC): %s\n",
            fmt_p(wt_t1_t3_abc$p.value)))

# Per-tertile summary
cat("\n  Per-tertile loop logFC:\n")
for (trt in c("T1_low", "T2_mid", "T3_high")) {
  sub <- k119_loop[k119_loop$tertile == trt, ]
  cat(sprintf("    %s: n=%d, median logFC=%.4f\n",
              trt, nrow(sub), median(sub$loop_logFC)))
}
cat("  Per-tertile mean delta_ABC:\n")
for (trt in c("T1_low", "T2_mid", "T3_high")) {
  sub <- k119_abc[k119_abc$tertile == trt, ]
  cat(sprintf("    %s: n=%d, median delta_ABC=%.5f\n",
              trt, nrow(sub), median(sub$mean_delta_abc)))
}

# --- Plot 13: Violin+box of loop logFC by tertile ---
cat("\nGenerating Part G plots...\n")
k119_loop$tertile <- factor(k119_loop$tertile,
                              levels = c("T1_low", "T2_mid", "T3_high"))

p13 <- ggplot(k119_loop, aes(x = tertile, y = loop_logFC, fill = tertile)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = TERTILE_COLORS) +
  labs(
    x = "K119ub gain tertile", y = "Loop logFC (KO/WT)",
    title = "Loop strength change by K119ub magnitude (K119ub_Only)",
    subtitle = sprintf("Spearman rho=%.3f, %s",
                        cor_tert_loop$estimate, fmt_p(cor_tert_loop$p.value))
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p13, "13_k119ub_tertile_loop_logfc", w = 6, h = 6)

# --- Plot 14: Boxplot of mean delta_ABC by tertile ---
k119_abc$tertile <- factor(k119_abc$tertile,
                             levels = c("T1_low", "T2_mid", "T3_high"))

p14 <- ggplot(k119_abc, aes(x = tertile, y = mean_delta_abc, fill = tertile)) +
  geom_boxplot(outlier.size = 0.3, notch = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = TERTILE_COLORS) +
  labs(
    x = "K119ub gain tertile", y = "Mean delta-ABC per enhancer",
    title = "ABC score change by K119ub magnitude (K119ub_Only)",
    subtitle = sprintf("Spearman rho=%.3f, %s",
                        cor_tert_abc$estimate, fmt_p(cor_tert_abc$p.value))
  ) +
  theme_pub +
  theme(legend.position = "none")
save_plot(p14, "14_k119ub_tertile_delta_abc", w = 6, h = 6)

# --- Plot 15: Continuous scatter: log2FC_K119ub vs loop_logFC ---
p15 <- ggplot(k119_loop, aes(x = log2FC_K119ub, y = loop_logFC)) +
  geom_point(alpha = 0.15, size = 0.6, color = "#B2182B") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "log2FC K119ub (KO/WT)", y = "Loop logFC (KO/WT)",
    title = "K119ub gain vs loop strength change (K119ub_Only enhancers)",
    subtitle = sprintf("Spearman rho=%.3f, %s",
                        cor_cont$estimate, fmt_p(cor_cont$p.value))
  ) +
  theme_pub
save_plot(p15, "15_k119ub_vs_loop_logfc_scatter", w = 7, h = 6)


# #############################################################################
# PART H: RNA-seq CONCORDANCE FOR TARGET GENES (Plots 16-18)
# #############################################################################

cat("\n=== PART H: RNA-seq Concordance for Target Genes ===\n\n")

# Build join key: abc_matched -> abc coordinates -> abc_rna
h_data <- data.frame(
  enh_idx = abc_matched$enh_idx,
  enhancer_class = abc_matched$enhancer_class,
  delta_ABC = abc_matched$delta_ABC,
  chr = abc$chr[abc_matched$abc_row],
  start = abc$start[abc_matched$abc_row],
  end = abc$end[abc_matched$abc_row],
  target_gene = abc$TargetGene[abc_matched$abc_row],
  stringsAsFactors = FALSE
)

# Merge with abc_rna on 4-column key
abc_rna_merge <- merge(
  h_data,
  abc_rna[, c("chr", "start", "end", "TargetGene", "log2FC", "padj")],
  by.x = c("chr", "start", "end", "target_gene"),
  by.y = c("chr", "start", "end", "TargetGene")
)
cat(sprintf("  ABC-RNA-seq merged pairs: %d\n", nrow(abc_rna_merge)))

# Filter for concordance: non-zero delta_ABC, non-NA and non-zero log2FC
conc_data <- abc_rna_merge[abs(abc_rna_merge$delta_ABC) > 0 &
                             !is.na(abc_rna_merge$log2FC) &
                             abc_rna_merge$log2FC != 0, ]
cat(sprintf("  After filtering zeros/NAs: %d\n", nrow(conc_data)))

# Concordance: delta_ABC and gene log2FC have same sign
conc_data$concordant <- sign(conc_data$delta_ABC) == sign(conc_data$log2FC)
conc_data$enhancer_class <- factor(conc_data$enhancer_class, levels = CLASS_ORDER)

# Per-class concordance
cat("\n  Per-class concordance:\n")
conc_by_class <- data.frame(
  enhancer_class = CLASS_ORDER,
  n_pairs = NA_integer_,
  n_concordant = NA_integer_,
  concordance_rate = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_along(CLASS_ORDER)) {
  cls <- CLASS_ORDER[i]
  sub <- conc_data[conc_data$enhancer_class == cls, ]
  conc_by_class$n_pairs[i] <- nrow(sub)
  conc_by_class$n_concordant[i] <- sum(sub$concordant)
  conc_by_class$concordance_rate[i] <- if (nrow(sub) > 0) {
    mean(sub$concordant)
  } else NA_real_
  cat(sprintf("    %s: %d / %d (%.1f%%)\n", cls, conc_by_class$n_concordant[i],
              conc_by_class$n_pairs[i],
              100 * conc_by_class$concordance_rate[i]))
}

# Binomial test: K119ub_Only concordance vs 50%
k119_conc <- conc_data[conc_data$enhancer_class == "K119ub_Only", ]
binom_k119 <- binom.test(sum(k119_conc$concordant), nrow(k119_conc), p = 0.5)
cat(sprintf("\n  Binomial K119ub_Only concordance vs 50%%: %s (rate=%.1f%%)\n",
            fmt_p(binom_k119$p.value), 100 * binom_k119$estimate))

# Fisher's exact: K119ub_Only concordance vs each other class
cat("  Fisher's exact (K119ub_Only concordance vs each class):\n")
k119_row <- conc_by_class[conc_by_class$enhancer_class == "K119ub_Only", ]
fisher_h_results <- list()
for (cls in setdiff(CLASS_ORDER, "K119ub_Only")) {
  other_row <- conc_by_class[conc_by_class$enhancer_class == cls, ]
  mat <- matrix(c(
    k119_row$n_concordant, k119_row$n_pairs - k119_row$n_concordant,
    other_row$n_concordant, other_row$n_pairs - other_row$n_concordant
  ), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat)
  fisher_h_results[[cls]] <- ft
  cat(sprintf("    K119ub_Only vs %s: OR=%.2f, %s\n",
              cls, ft$estimate, fmt_p(ft$p.value)))
}

# DE fraction among K119ub_Only target genes
k119_rna <- abc_rna_merge[abc_rna_merge$enhancer_class == "K119ub_Only" &
                            !is.na(abc_rna_merge$padj), ]
k119_de_genes <- unique(k119_rna$target_gene[k119_rna$padj < 0.05 &
                                               abs(k119_rna$log2FC) > 0.5])
k119_all_genes <- unique(k119_rna$target_gene)
cat(sprintf("\n  K119ub_Only target genes: %d total, %d DE (padj<0.05, |log2FC|>0.5) = %.1f%%\n",
            length(k119_all_genes), length(k119_de_genes),
            100 * length(k119_de_genes) / max(length(k119_all_genes), 1)))

# --- Plot 16: Stacked bar of concordance rate by class ---
cat("\nGenerating Part H plots...\n")
conc_bar <- data.frame(
  enhancer_class = rep(CLASS_ORDER, each = 2),
  type = rep(c("Concordant", "Discordant"), times = length(CLASS_ORDER)),
  count = NA_integer_,
  stringsAsFactors = FALSE
)
for (i in seq_along(CLASS_ORDER)) {
  cls <- CLASS_ORDER[i]
  row_i <- conc_by_class[conc_by_class$enhancer_class == cls, ]
  conc_bar$count[(i-1)*2 + 1] <- row_i$n_concordant
  conc_bar$count[(i-1)*2 + 2] <- row_i$n_pairs - row_i$n_concordant
}
conc_bar$enhancer_class <- factor(conc_bar$enhancer_class, levels = CLASS_ORDER)
conc_bar$type <- factor(conc_bar$type, levels = c("Concordant", "Discordant"))

p16 <- ggplot(conc_bar, aes(x = enhancer_class, y = count, fill = type)) +
  geom_col(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c(Concordant = "#4DAF4A", Discordant = "#E41A1C")) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Enhancer class", y = "Proportion of E-G pairs",
    fill = "Concordance",
    title = "ABC-RNA concordance by enhancer class",
    subtitle = sprintf("K119ub_Only: %.1f%% concordant, %s",
                        100 * binom_k119$estimate, fmt_p(binom_k119$p.value))
  ) +
  theme_pub
save_plot(p16, "16_abc_rnaseq_concordance_bar", w = 7, h = 6)

# --- Plot 17: Scatter delta_ABC vs gene log2FC colored by class ---
# Downsample Stable to 3K for visibility
set.seed(42)
scatter_data <- abc_rna_merge[!is.na(abc_rna_merge$log2FC), ]
scatter_data$enhancer_class <- factor(scatter_data$enhancer_class, levels = CLASS_ORDER)
stable_rows <- which(scatter_data$enhancer_class == "Stable")
if (length(stable_rows) > 3000) {
  keep_stable <- sample(stable_rows, 3000)
  drop_stable <- setdiff(stable_rows, keep_stable)
  scatter_data <- scatter_data[-drop_stable, ]
}
cat(sprintf("  Scatter plot data: %d pairs (Stable downsampled to <=3K)\n",
            nrow(scatter_data)))

p17 <- ggplot(scatter_data, aes(x = delta_ABC, y = log2FC, color = enhancer_class)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = CLASS_COLORS) +
  labs(
    x = "delta-ABC (KO - WT)", y = "Gene log2FC (KO/WT)",
    color = "Enhancer class",
    title = "ABC change vs gene expression change",
    subtitle = "Each point = one enhancer-gene pair"
  ) +
  theme_pub
save_plot(p17, "17_delta_abc_vs_gene_logfc_scatter", w = 8, h = 7)

# --- Plot 18: Histogram of K119ub_Only target gene log2FC ---
k119_gene_logfc <- unique(k119_rna[, c("target_gene", "log2FC", "padj")])
k119_gene_logfc <- k119_gene_logfc[!is.na(k119_gene_logfc$log2FC), ]
k119_gene_logfc$is_de <- !is.na(k119_gene_logfc$padj) &
                          k119_gene_logfc$padj < 0.05 &
                          abs(k119_gene_logfc$log2FC) > 0.5

p18 <- ggplot(k119_gene_logfc, aes(x = log2FC, fill = is_de)) +
  geom_histogram(bins = 50, color = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B2182B"),
                     labels = c("Not DE", "DE (padj<0.05, |log2FC|>0.5)")) +
  labs(
    x = "Gene log2FC (KO/WT)", y = "Count",
    fill = "",
    title = "K119ub_Only target gene expression changes",
    subtitle = sprintf("%d unique genes, %d DE (%.1f%%)",
                        nrow(k119_gene_logfc), sum(k119_gene_logfc$is_de),
                        100 * mean(k119_gene_logfc$is_de))
  ) +
  theme_pub
save_plot(p18, "18_k119ub_target_gene_logfc_hist", w = 7, h = 5)


# #############################################################################
# PART I: MOTIF ENRICHMENT BED EXPORT (No Plots)
# #############################################################################

cat("\n=== PART I: Motif Enrichment BED Export ===\n\n")

bed_dir <- file.path(OUTPUT_DIR, "beds_for_homer")
dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)

# Export 4 class-level BED files (BED6: chr, start, end, name, score, strand)
cat("  Exporting class-level BED files:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_all[enh_all$enhancer_class == cls, c("chr", "start", "end")]
  bed <- data.frame(
    chr = sub$chr, start = sub$start, end = sub$end,
    name = cls, score = 0, strand = "."
  )
  bed_path <- file.path(bed_dir, sprintf("%s.bed", cls))
  write.table(bed, bed_path, sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  cat(sprintf("    %s: %d regions\n", basename(bed_path), nrow(bed)))
}

# Export 3 K119ub_Only tertile BED files
cat("  Exporting K119ub tertile BED files:\n")
for (trt in c("T1_low", "T2_mid", "T3_high")) {
  sub_idx <- k119_tertile$enh_idx[k119_tertile$tertile == trt]
  sub <- enh_all[sub_idx, c("chr", "start", "end")]
  bed <- data.frame(
    chr = sub$chr, start = sub$start, end = sub$end,
    name = trt, score = 0, strand = "."
  )
  bed_path <- file.path(bed_dir, sprintf("K119ub_Only_%s.bed", trt))
  write.table(bed, bed_path, sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  cat(sprintf("    %s: %d regions\n", basename(bed_path), nrow(bed)))
}

# Print HOMER commands
cat("\n  HOMER findMotifsGenome.pl commands:\n")
cat("  # K119ub_Only vs Stable background\n")
cat(sprintf("  findMotifsGenome.pl %s/K119ub_Only.bed mm10 homer_out/K119ub_Only_vs_Stable -bg %s/Stable.bed -size given\n",
            bed_dir, bed_dir))
cat("\n  # Each class vs Stable background\n")
for (cls in setdiff(CLASS_ORDER, "Stable")) {
  cat(sprintf("  findMotifsGenome.pl %s/%s.bed mm10 homer_out/%s_vs_Stable -bg %s/Stable.bed -size given\n",
              bed_dir, cls, cls, bed_dir))
}
cat("\n  # T3_high vs T1_low background (gradient test)\n")
cat(sprintf("  findMotifsGenome.pl %s/K119ub_Only_T3_high.bed mm10 homer_out/T3_high_vs_T1_low -bg %s/K119ub_Only_T1_low.bed -size given\n",
            bed_dir, bed_dir))


# =============================================================================
# ACCUMULATE PART G/H STATISTICS
# =============================================================================

# Part G stats
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "Spearman_K119ub_tertile_vs_loopLogFC",
  statistic = cor_tert_loop$estimate,
  p_value = cor_tert_loop$p.value,
  comparison = "K119ub_tertile_rank_vs_loop_logFC"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "Spearman_K119ub_tertile_vs_deltaABC",
  statistic = cor_tert_abc$estimate,
  p_value = cor_tert_abc$p.value,
  comparison = "K119ub_tertile_rank_vs_mean_delta_ABC"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "Spearman_K119ub_continuous_vs_loopLogFC",
  statistic = cor_cont$estimate,
  p_value = cor_cont$p.value,
  comparison = "K119ub_log2FC_vs_loop_logFC"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "KW_loop_logFC_by_K119ub_tertile",
  statistic = kw_tert_loop$statistic,
  p_value = kw_tert_loop$p.value,
  comparison = "3_tertiles"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "KW_deltaABC_by_K119ub_tertile",
  statistic = kw_tert_abc$statistic,
  p_value = kw_tert_abc$p.value,
  comparison = "3_tertiles"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "Wilcoxon_T1_vs_T3_loopLogFC",
  statistic = wt_t1_t3_loop$statistic,
  p_value = wt_t1_t3_loop$p.value,
  comparison = "T1_low_vs_T3_high"))
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "Wilcoxon_T1_vs_T3_deltaABC",
  statistic = wt_t1_t3_abc$statistic,
  p_value = wt_t1_t3_abc$p.value,
  comparison = "T1_low_vs_T3_high"))

# Part H stats
stat_tests <- rbind(stat_tests, data.frame(
  test_name = "Binomial_K119ub_concordance_vs_50pct",
  statistic = binom_k119$statistic,
  p_value = binom_k119$p.value,
  comparison = "K119ub_Only_vs_chance"))
for (cls in names(fisher_h_results)) {
  ft <- fisher_h_results[[cls]]
  stat_tests <- rbind(stat_tests, data.frame(
    test_name = "Fisher_K119ub_concordance_vs_class",
    statistic = ft$estimate,
    p_value = ft$p.value,
    comparison = sprintf("K119ub_Only_vs_%s", cls)))
}

cat(sprintf("\n  Total statistical tests collected: %d\n", nrow(stat_tests)))


# =============================================================================
# SAVE OUTPUT FILES
# =============================================================================

cat("\n=== Saving output files ===\n\n")

# 1. Per-enhancer loop metrics
loop_out <- enh_metrics[, c("chr", "start", "end", "enhancer_class",
                             "n_loops", "mean_loop_distance", "n_promoter_loops")]
write.table(loop_out, file.path(OUTPUT_DIR, "enhancer_class_loop_metrics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  enhancer_class_loop_metrics.tsv: %d rows\n", nrow(loop_out)))

# 2. Per-enhancer ABC metrics
abc_out <- enh_abc_metrics[, c("chr", "start", "end", "enhancer_class",
                                "n_abc_pairs", "mean_delta_abc",
                                "mean_delta_unnorm", "n_gained", "n_lost")]
write.table(abc_out, file.path(OUTPUT_DIR, "enhancer_class_abc_metrics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  enhancer_class_abc_metrics.tsv: %d rows\n", nrow(abc_out)))

# 3. Class-level summary
write.table(class_summary, file.path(OUTPUT_DIR, "class_level_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  class_level_summary.tsv: %d rows\n", nrow(class_summary)))

# 4. Statistical tests
write.table(stat_tests, file.path(OUTPUT_DIR, "statistical_tests.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  statistical_tests.tsv: %d rows\n", nrow(stat_tests)))

# 5. Contact decay table
write.table(decay_agg, file.path(OUTPUT_DIR, "contact_decay_by_class.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  contact_decay_by_class.tsv: %d rows\n", nrow(decay_agg)))

# 6. Promoter-loop gene logFC
if (nrow(promoter_gene_de) > 0) {
  write.table(promoter_gene_de[, c("gene_name", "enhancer_class", "log2FC", "padj")],
              file.path(OUTPUT_DIR, "promoter_loop_gene_logfc.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  promoter_loop_gene_logfc.tsv: %d rows\n", nrow(promoter_gene_de)))
}

# 7. K119ub tertile assignments
k119_tert_out <- data.frame(
  chr = k119_enh$chr, start = k119_enh$start, end = k119_enh$end,
  tertile = k119_enh$tertile, log2FC_K119ub = k119_enh$log2FC_K119ub
)
write.table(k119_tert_out, file.path(OUTPUT_DIR, "k119ub_tertile_assignments.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  k119ub_tertile_assignments.tsv: %d rows\n", nrow(k119_tert_out)))

# 8. K119ub tertile loop summary
tert_loop_summary <- data.frame(tertile = c("T1_low", "T2_mid", "T3_high"))
for (i in 1:3) {
  trt <- tert_loop_summary$tertile[i]
  sub_loop <- k119_loop[k119_loop$tertile == trt, ]
  sub_abc_t <- k119_abc[k119_abc$tertile == trt, ]
  tert_loop_summary$n_enhancers[i] <- sum(k119_enh$tertile == trt)
  tert_loop_summary$n_loop_pairs[i] <- nrow(sub_loop)
  tert_loop_summary$median_loop_logFC[i] <- median(sub_loop$loop_logFC)
  tert_loop_summary$n_abc_enhancers[i] <- nrow(sub_abc_t)
  tert_loop_summary$median_delta_abc[i] <- median(sub_abc_t$mean_delta_abc)
}
write.table(tert_loop_summary, file.path(OUTPUT_DIR, "k119ub_tertile_loop_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  k119ub_tertile_loop_summary.tsv: %d rows\n", nrow(tert_loop_summary)))

# 9. ABC-RNA concordance by class
write.table(conc_by_class, file.path(OUTPUT_DIR, "abc_rnaseq_concordance_by_class.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  abc_rnaseq_concordance_by_class.tsv: %d rows\n", nrow(conc_by_class)))

# 10. K119ub_Only target gene DE
write.table(k119_gene_logfc, file.path(OUTPUT_DIR, "k119ub_only_target_gene_de.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  k119ub_only_target_gene_de.tsv: %d rows\n", nrow(k119_gene_logfc)))

# =============================================================================
# DONE
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 11 COMPLETE\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("Plots: 18 panels saved (PDF + SVG + JPG in subfolders)\n")
cat("Tables: 10 TSV files + 7 BED files\n")
cat("================================================================================\n")
