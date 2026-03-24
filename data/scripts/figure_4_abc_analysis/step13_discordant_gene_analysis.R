# abc/scripts/step13_discordant_gene_analysis.R
# Exploratory analysis of discordant genes — those where sign(dABC) opposes
# sign(log2FC). Understanding why ~35% of dysregulated genes are discordant
# is biologically important: are these measurement artifacts, indirect
# regulation, or conflicting multi-enhancer signals?
#
# Inputs (relative to abc/ working directory):
#   results/gene_level_summary.tsv                              — 13,588 genes
#   results/delta_abc_all_pairs.tsv                             — 180K E-G pairs
#   results/enhancer_subset_analysis/enhancer_class_abc_metrics.tsv — 21.5K classified enhancers
#
# Outputs:
#   results/figures/discordant_analysis/    — multi-panel diagnostic figures
#   results/discordant_gene_characteristics.tsv — per-gene stats
#
# Usage:
#   cd abc && Rscript scripts/step13_discordant_gene_analysis.R

# =============================================================================
# CONFIGURATION
# =============================================================================

GENE_SUMMARY_FILE <- "results/gene_level_summary.tsv"
ABC_PAIRS_FILE    <- "results/delta_abc_all_pairs.tsv"
ENH_CLASS_FILE    <- "results/enhancer_subset_analysis/enhancer_class_abc_metrics.tsv"
OUTPUT_DIR        <- "results/figures/discordant_analysis"

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

CONCORDANCE_COLORS <- c(Concordant = "#4DAF4A", Discordant = "#E41A1C")

cat("================================================================================\n")
cat("STEP 13: DISCORDANT GENE ANALYSIS\n")
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
  library(dplyr)
  library(tidyr)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating input files...\n")
stopifnot(file.exists(GENE_SUMMARY_FILE))
cat(sprintf("  [OK] Gene summary: %s\n", GENE_SUMMARY_FILE))
stopifnot(file.exists(ABC_PAIRS_FILE))
cat(sprintf("  [OK] ABC pairs: %s\n", ABC_PAIRS_FILE))
stopifnot(file.exists(ENH_CLASS_FILE))
cat(sprintf("  [OK] Enhancer classes: %s\n", ENH_CLASS_FILE))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nOutput directory: %s\n\n", OUTPUT_DIR))

# =============================================================================
# HELPERS
# =============================================================================

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

UTIL_PATH <- "../scripts/utils/multi_format_output.R"
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

fmt_p <- function(p) {
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

genes <- read.table(GENE_SUMMARY_FILE, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE)
cat(sprintf("  Genes: %d\n", nrow(genes)))

abc <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d\n", nrow(abc)))

enh_class <- read.table(ENH_CLASS_FILE, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
cat(sprintf("  Classified enhancers: %d\n", nrow(enh_class)))

# =============================================================================
# IDENTIFY CONCORDANT vs DISCORDANT GENES
# =============================================================================

cat("\n--- Identifying concordant vs discordant genes ---\n")

# Dysregulated genes have both |dABC| > 0.01 and DE (padj < 0.05, |log2FC| > 0.5)
# Column stored as string "True"/"False" from Python output
dysreg <- genes[genes$dysregulated == "True", ]
cat(sprintf("  Dysregulated genes: %d\n", nrow(dysreg)))

stopifnot(nrow(dysreg) > 0)

# Concordant: sign(max_delta_abc) == sign(log2FC)
# dABC > 0 means KO gained E-P linkage → expect upregulation (log2FC > 0)
dysreg$concordance <- ifelse(
  sign(dysreg$max_delta_abc) == sign(dysreg$log2FC),
  "Concordant", "Discordant"
)

n_conc <- sum(dysreg$concordance == "Concordant")
n_disc <- sum(dysreg$concordance == "Discordant")
cat(sprintf("  Concordant: %d (%.1f%%)\n", n_conc,
            100 * n_conc / nrow(dysreg)))
cat(sprintf("  Discordant: %d (%.1f%%)\n", n_disc,
            100 * n_disc / nrow(dysreg)))

# =============================================================================
# ANALYSIS 1: MULTI-ENHANCER CONFLICT
# =============================================================================

cat("\n--- Analysis 1: Multi-enhancer conflict ---\n")

# For each dysregulated gene, count enhancers with positive vs negative dABC
dysreg_genes <- dysreg$TargetGene

# Get max_delta_abc direction from gene summary (not in abc table)
gene_max_dir <- setNames(sign(dysreg$max_delta_abc), dysreg$TargetGene)

enh_per_gene <- abc %>%
  filter(TargetGene %in% dysreg_genes) %>%
  group_by(TargetGene) %>%
  summarise(
    n_enh_total = n(),
    n_enh_positive_dabc = sum(delta_ABC > 0),
    n_enh_negative_dabc = sum(delta_ABC < 0),
    n_enh_zero_dabc = sum(delta_ABC == 0),
    .groups = "drop"
  ) %>%
  mutate(
    max_dir = gene_max_dir[TargetGene],
    frac_agree_max = ifelse(
      max_dir == 0, NA_real_,
      ifelse(max_dir > 0,
             n_enh_positive_dabc / n_enh_total,
             n_enh_negative_dabc / n_enh_total)
    )
  )

# Merge concordance status
enh_per_gene <- left_join(enh_per_gene,
                          dysreg[, c("TargetGene", "concordance")],
                          by = "TargetGene")

cat(sprintf("  Median frac_agree_max (Concordant): %.3f\n",
            median(enh_per_gene$frac_agree_max[
              enh_per_gene$concordance == "Concordant"], na.rm = TRUE)))
cat(sprintf("  Median frac_agree_max (Discordant): %.3f\n",
            median(enh_per_gene$frac_agree_max[
              enh_per_gene$concordance == "Discordant"], na.rm = TRUE)))

wt_agree <- wilcox.test(
  frac_agree_max ~ concordance, data = enh_per_gene
)
cat(sprintf("  Wilcoxon test: %s\n", fmt_p(wt_agree$p.value)))

p_conflict <- ggplot(enh_per_gene,
                     aes(x = concordance, y = frac_agree_max,
                         fill = concordance)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "Enhancer Agreement with Dominant dABC",
    subtitle = sprintf("Wilcoxon %s", fmt_p(wt_agree$p.value)),
    x = NULL,
    y = "Fraction of enhancers\nagreeing with max dABC direction"
  ) +
  theme_pub

# =============================================================================
# ANALYSIS 2: |dABC| MAGNITUDE
# =============================================================================

cat("\n--- Analysis 2: dABC magnitude ---\n")

wt_dabc <- wilcox.test(
  abs(max_delta_abc) ~ concordance, data = dysreg
)
cat(sprintf("  Median |max_delta_abc| Concordant: %.4f\n",
            median(abs(dysreg$max_delta_abc[
              dysreg$concordance == "Concordant"]))))
cat(sprintf("  Median |max_delta_abc| Discordant: %.4f\n",
            median(abs(dysreg$max_delta_abc[
              dysreg$concordance == "Discordant"]))))
cat(sprintf("  Wilcoxon test: %s\n", fmt_p(wt_dabc$p.value)))

p_dabc_mag <- ggplot(dysreg,
                     aes(x = concordance, y = abs(max_delta_abc),
                         fill = concordance)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  scale_y_log10() +
  labs(
    title = "|max dABC| Distribution",
    subtitle = sprintf("Wilcoxon %s", fmt_p(wt_dabc$p.value)),
    x = NULL,
    y = "|max dABC| (log10 scale)"
  ) +
  theme_pub

# =============================================================================
# ANALYSIS 3: RNA-seq EFFECT SIZE
# =============================================================================

cat("\n--- Analysis 3: RNA-seq effect size ---\n")

wt_lfc <- wilcox.test(
  abs(log2FC) ~ concordance, data = dysreg
)
cat(sprintf("  Median |log2FC| Concordant: %.3f\n",
            median(abs(dysreg$log2FC[dysreg$concordance == "Concordant"]))))
cat(sprintf("  Median |log2FC| Discordant: %.3f\n",
            median(abs(dysreg$log2FC[dysreg$concordance == "Discordant"]))))
cat(sprintf("  Wilcoxon test: %s\n", fmt_p(wt_lfc$p.value)))

p_lfc_mag <- ggplot(dysreg,
                    aes(x = concordance, y = abs(log2FC),
                        fill = concordance)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "|log2FC| Distribution",
    subtitle = sprintf("Wilcoxon %s", fmt_p(wt_lfc$p.value)),
    x = NULL,
    y = "|log2FC|"
  ) +
  theme_pub

# padj comparison
wt_padj <- wilcox.test(
  padj ~ concordance, data = dysreg
)
cat(sprintf("  Median padj Concordant: %.4f\n",
            median(dysreg$padj[dysreg$concordance == "Concordant"],
                   na.rm = TRUE)))
cat(sprintf("  Median padj Discordant: %.4f\n",
            median(dysreg$padj[dysreg$concordance == "Discordant"],
                   na.rm = TRUE)))

p_padj <- ggplot(dysreg,
                 aes(x = concordance, y = -log10(padj),
                     fill = concordance)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "-log10(padj) Distribution",
    subtitle = sprintf("Wilcoxon %s", fmt_p(wt_padj$p.value)),
    x = NULL,
    y = "-log10(padj)"
  ) +
  theme_pub

# =============================================================================
# ANALYSIS 4: ENHANCER CLASS ENRICHMENT
# =============================================================================

cat("\n--- Analysis 4: Enhancer class enrichment ---\n")

# Join enhancer classes via overlap (coordinates differ between files)
abc_gr <- GRanges(seqnames = abc$chr,
                  ranges = IRanges(start = abc$start, end = abc$end))
class_gr <- GRanges(seqnames = enh_class$chr,
                    ranges = IRanges(start = enh_class$start,
                                    end = enh_class$end))
hits <- findOverlaps(abc_gr, class_gr)
first_hits <- hits[!duplicated(queryHits(hits))]

abc$enh_class_label <- NA_character_
abc$enh_class_label[queryHits(first_hits)] <-
  enh_class$enhancer_class[subjectHits(first_hits)]

# For each dysregulated gene, find its strongest enhancer's class
strongest_enh <- abc %>%
  filter(TargetGene %in% dysreg_genes) %>%
  group_by(TargetGene) %>%
  slice_max(abs(delta_ABC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(TargetGene, enh_class_label)

dysreg_with_class <- left_join(dysreg, strongest_enh, by = "TargetGene")
dysreg_with_class$enh_class_label[is.na(dysreg_with_class$enh_class_label)] <-
  "Unclassified"

# Count table for Fisher's exact test
class_conc_tab <- table(
  dysreg_with_class$enh_class_label,
  dysreg_with_class$concordance
)

cat("  Concordance by strongest enhancer class:\n")
print(class_conc_tab)

# Fisher's test (overall association)
if (nrow(class_conc_tab) >= 2 && ncol(class_conc_tab) >= 2) {
  fisher_res <- fisher.test(class_conc_tab, simulate.p.value = TRUE, B = 10000)
  cat(sprintf("  Fisher's exact test: %s\n", fmt_p(fisher_res$p.value)))
} else {
  fisher_res <- list(p.value = NA)
  cat("  Fisher's test: insufficient data\n")
}

# Bar plot of class distribution by concordance
class_conc_df <- as.data.frame(class_conc_tab)
names(class_conc_df) <- c("enhancer_class", "concordance", "count")

# Compute proportions within each concordance group
class_conc_df <- class_conc_df %>%
  group_by(concordance) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Reorder for display (classified classes first)
class_conc_df$enhancer_class <- factor(
  class_conc_df$enhancer_class,
  levels = c(CLASS_ORDER, "Unclassified")
)

all_class_colors <- c(CLASS_COLORS, Unclassified = "grey70")

p_class_enrich <- ggplot(class_conc_df,
                         aes(x = concordance, y = prop,
                             fill = enhancer_class)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = all_class_colors, name = "Strongest Enhancer Class") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Enhancer Class of Strongest E-G Pair",
    subtitle = sprintf("Fisher's %s", fmt_p(fisher_res$p.value)),
    x = NULL,
    y = "Proportion"
  ) +
  theme_pub

# =============================================================================
# ANALYSIS 5: DISTANCE DISTRIBUTION
# =============================================================================

cat("\n--- Analysis 5: Distance distribution ---\n")

# top_enh_distance is the distance of the strongest enhancer
wt_dist <- wilcox.test(
  top_enh_distance ~ concordance, data = dysreg
)
cat(sprintf("  Median distance Concordant: %s bp\n",
            format(median(dysreg$top_enh_distance[
              dysreg$concordance == "Concordant"]), big.mark = ",")))
cat(sprintf("  Median distance Discordant: %s bp\n",
            format(median(dysreg$top_enh_distance[
              dysreg$concordance == "Discordant"]), big.mark = ",")))
cat(sprintf("  Wilcoxon test: %s\n", fmt_p(wt_dist$p.value)))

p_distance <- ggplot(dysreg,
                     aes(x = concordance, y = top_enh_distance / 1e6,
                         fill = concordance)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "Distance to Strongest Enhancer",
    subtitle = sprintf("Wilcoxon %s", fmt_p(wt_dist$p.value)),
    x = NULL,
    y = "Distance (Mb)"
  ) +
  theme_pub

# =============================================================================
# ANALYSIS 6: NUMBER OF ENHANCERS
# =============================================================================

cat("\n--- Analysis 6: Number of enhancers ---\n")

wt_nenh <- wilcox.test(
  n_enhancers ~ concordance, data = dysreg
)
cat(sprintf("  Median n_enhancers Concordant: %d\n",
            median(dysreg$n_enhancers[dysreg$concordance == "Concordant"])))
cat(sprintf("  Median n_enhancers Discordant: %d\n",
            median(dysreg$n_enhancers[dysreg$concordance == "Discordant"])))
cat(sprintf("  Wilcoxon test: %s\n", fmt_p(wt_nenh$p.value)))

p_nenh <- ggplot(dysreg,
                 aes(x = concordance, y = n_enhancers,
                     fill = concordance)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  scale_y_log10() +
  labs(
    title = "Number of Enhancers per Gene",
    subtitle = sprintf("Wilcoxon %s", fmt_p(wt_nenh$p.value)),
    x = NULL,
    y = "n enhancers (log10 scale)"
  ) +
  theme_pub

# =============================================================================
# ANALYSIS 7: ENHANCER WIDTH (peak size of strongest enhancer)
# =============================================================================

cat("\n--- Analysis 7: Enhancer width ---\n")

# Peak width of the strongest enhancer per gene
dysreg$enh_width <- dysreg$top_enh_end - dysreg$top_enh_start

wt_width <- wilcox.test(enh_width ~ concordance, data = dysreg)

conc_width <- dysreg$enh_width[dysreg$concordance == "Concordant"]
disc_width <- dysreg$enh_width[dysreg$concordance == "Discordant"]

cat(sprintf("  Median width Concordant: %d bp (%.1f%% exactly 500bp)\n",
            median(conc_width), 100 * mean(conc_width == 500)))
cat(sprintf("  Median width Discordant: %d bp (%.1f%% exactly 500bp)\n",
            median(disc_width), 100 * mean(disc_width == 500)))
cat(sprintf("  Pct > 500bp — Concordant: %.1f%%, Discordant: %.1f%%\n",
            100 * mean(conc_width > 500), 100 * mean(disc_width > 500)))
cat(sprintf("  Wilcoxon test: %s\n", fmt_p(wt_width$p.value)))

p_enh_width <- ggplot(dysreg,
                      aes(x = concordance, y = enh_width,
                          fill = concordance)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.1, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "Strongest Enhancer Width",
    subtitle = sprintf("Wilcoxon %s", fmt_p(wt_width$p.value)),
    x = NULL,
    y = "Peak width (bp)"
  ) +
  theme_pub

# =============================================================================
# COMPOSITE FIGURE
# =============================================================================

cat("\n--- Assembling composite figure ---\n")

p_composite <- (p_conflict | p_dabc_mag | p_lfc_mag) /
               (p_padj | p_class_enrich | p_distance) /
               (p_nenh | p_enh_width | plot_spacer()) +
  plot_annotation(
    title = "Discordant Gene Characterization",
    subtitle = sprintf(
      "%d concordant (%.0f%%) vs %d discordant (%.0f%%) dysregulated genes",
      n_conc, 100 * n_conc / nrow(dysreg),
      n_disc, 100 * n_disc / nrow(dysreg)
    ),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(p_composite, "01_discordant_composite", w = 14, h = 14)

# Also save individual panels
save_plot(p_conflict, "02_enhancer_agreement", w = 5, h = 5)
save_plot(p_dabc_mag, "03_dabc_magnitude", w = 5, h = 5)
save_plot(p_lfc_mag, "04_log2fc_magnitude", w = 5, h = 5)
save_plot(p_class_enrich, "05_class_enrichment", w = 6, h = 5)
save_plot(p_distance, "06_distance", w = 5, h = 5)
save_plot(p_nenh, "07_n_enhancers", w = 5, h = 5)
save_plot(p_enh_width, "15_enhancer_width", w = 5, h = 5)

# =============================================================================
# SCATTER: max_delta_abc vs log2FC colored by concordance
# =============================================================================

cat("\n--- Scatter: dABC vs log2FC ---\n")

p_scatter <- ggplot(dysreg,
                    aes(x = max_delta_abc, y = log2FC,
                        color = concordance)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = CONCORDANCE_COLORS, name = "Concordance") +
  labs(
    title = "max dABC vs log2FC for Dysregulated Genes",
    subtitle = sprintf("n = %d | Concordant: %d, Discordant: %d",
                       nrow(dysreg), n_conc, n_disc),
    x = "max dABC (strongest enhancer)",
    y = "RNA-seq log2FC (KO/WT)"
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p_scatter, "08_dabc_vs_log2fc_scatter", w = 7, h = 6)

# =============================================================================
# OUTPUT TABLE
# =============================================================================

cat("\n--- Saving gene characteristics table ---\n")

# Merge enhancer conflict info with dysregulated gene table
out_df <- left_join(
  dysreg[, c("TargetGene", "max_delta_abc", "log2FC", "padj",
             "n_enhancers", "top_enh_distance", "enh_width", "concordance",
             "n_gained", "n_lost", "sum_delta_abc", "sum_abc_wt", "sum_abc_ko")],
  enh_per_gene[, c("TargetGene", "n_enh_positive_dabc",
                    "n_enh_negative_dabc", "frac_agree_max")],
  by = "TargetGene"
)

# Add strongest enhancer class
out_df <- left_join(out_df, strongest_enh, by = "TargetGene")
names(out_df)[names(out_df) == "enh_class_label"] <- "strongest_enh_class"

out_tsv <- "results/discordant_gene_characteristics.tsv"
write.table(out_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: %s (%d rows)\n", out_tsv, nrow(out_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 13 COMPLETE — DISCORDANT GENE ANALYSIS SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("Dysregulated genes: %d\n", nrow(dysreg)))
cat(sprintf("  Concordant: %d (%.1f%%)\n", n_conc, 100 * n_conc / nrow(dysreg)))
cat(sprintf("  Discordant: %d (%.1f%%)\n", n_disc, 100 * n_disc / nrow(dysreg)))

cat("\nKey comparisons (Concordant vs Discordant):\n")
cat(sprintf("  |max dABC|:      Wilcoxon %s\n", fmt_p(wt_dabc$p.value)))
cat(sprintf("  |log2FC|:        Wilcoxon %s\n", fmt_p(wt_lfc$p.value)))
cat(sprintf("  Enhancer agree:  Wilcoxon %s\n", fmt_p(wt_agree$p.value)))
cat(sprintf("  Distance:        Wilcoxon %s\n", fmt_p(wt_dist$p.value)))
cat(sprintf("  n_enhancers:     Wilcoxon %s\n", fmt_p(wt_nenh$p.value)))
cat(sprintf("  Enhancer width:  Wilcoxon %s\n", fmt_p(wt_width$p.value)))
cat(sprintf("  Class enrich:    Fisher's %s\n", fmt_p(fisher_res$p.value)))

cat(sprintf("\nOutputs:\n"))
cat(sprintf("  Figures: %s\n", OUTPUT_DIR))
cat(sprintf("  Table:   %s\n", out_tsv))
cat("================================================================================\n")
