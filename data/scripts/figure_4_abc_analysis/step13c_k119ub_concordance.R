# abc/scripts/step13c_k119ub_concordance.R
#
# Test whether K119ub signal at enhancers differs between enhancers linked
# to concordant vs discordant genes. Concordant = sign(dABC) matches sign(log2FC);
# discordant = opposite. Enhancers may link to multiple dysregulated genes, so
# we classify each enhancer by the consistency of its target concordance.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(dplyr)
  library(tidyr)
})

cat("=== Step 13c: K119ub Signal by Concordance Status ===\n")
cat(sprintf("  Started: %s\n\n", Sys.time()))

# =============================================================================
# CONFIGURATION
# =============================================================================

K119UB_FILE     <- "data/upstream/chip_peaks/k119ub_enhancer_signal.tsv"  # Original: results/k119ub_abc_enhancer_merged.tsv
DISCORDANT_FILE <- "data/tsvs/figure_4_abc_analysis/4B_discordant_gene_characteristics.tsv"  # Original: results/discordant_gene_characteristics.tsv
ABC_PAIRS_FILE  <- "data/tsvs/figure_4_abc_analysis/4A_delta_abc_all_pairs.tsv"  # Original: results/delta_abc_all_pairs.tsv
ENH_CLASS_FILE  <- "data/tsvs/figure_4_abc_analysis/4D_enhancer_class_abc_metrics.tsv"  # Original: results/enhancer_subset_analysis/enhancer_class_abc_metrics.tsv
OUTPUT_DIR      <- "data/plots/figure_4_abc_analysis"  # Original: results/figures/discordant_analysis
TSV_DIR         <- "data/tsvs/figure_4_abc_analysis"  # Original: results/figures/discordant_analysis (TSV outputs)
UTIL_PATH       <- "data/scripts/_shared/multi_format_output.R"  # Original: ../scripts/utils/multi_format_output.R

stopifnot(file.exists(K119UB_FILE))
stopifnot(file.exists(DISCORDANT_FILE))
stopifnot(file.exists(ABC_PAIRS_FILE))
stopifnot(file.exists(ENH_CLASS_FILE))
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TSV_DIR, recursive = TRUE, showWarnings = FALSE)

# Shared aesthetics
theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

CONCORDANCE_COLORS <- c(
  "Concordant-linked" = "#4DAF4A",
  "Discordant-linked" = "#E41A1C",
  "Mixed"             = "#984EA3"
)

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

k119ub <- read.table(K119UB_FILE, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
cat(sprintf("  K119ub enhancers: %d\n", nrow(k119ub)))

disc <- read.table(DISCORDANT_FILE, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)
cat(sprintf("  Dysregulated genes: %d (Concordant=%d, Discordant=%d)\n",
            nrow(disc),
            sum(disc$concordance == "Concordant"),
            sum(disc$concordance == "Discordant")))

abc <- read.table(ABC_PAIRS_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE)
cat(sprintf("  ABC pairs: %d\n", nrow(abc)))

enh_class <- read.table(ENH_CLASS_FILE, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)
cat(sprintf("  Classified enhancers: %d\n", nrow(enh_class)))

# =============================================================================
# MAP ENHANCERS → CONCORDANCE STATUS
# =============================================================================

cat("\nMapping enhancers to concordance status...\n")

# Build enhancer key in ABC pairs (same format as k119ub file)
abc$enh_key <- paste0(abc$chr, ":", abc$start, "-", abc$end)

# Subset ABC pairs to dysregulated genes only
abc_dysreg <- abc[abc$TargetGene %in% disc$TargetGene, ]
cat(sprintf("  ABC pairs linking to dysregulated genes: %d\n", nrow(abc_dysreg)))

# Add concordance status from gene table
abc_dysreg <- merge(abc_dysreg,
                    disc[, c("TargetGene", "concordance")],
                    by = "TargetGene")

# Per enhancer: classify by consistency of linked gene concordance
enh_concordance <- abc_dysreg %>%
  group_by(enh_key) %>%
  summarise(
    n_targets = n_distinct(TargetGene),
    n_concordant = sum(concordance == "Concordant"),
    n_discordant = sum(concordance == "Discordant"),
    .groups = "drop"
  ) %>%
  mutate(
    linkage_group = case_when(
      n_discordant == 0 ~ "Concordant-linked",
      n_concordant == 0 ~ "Discordant-linked",
      TRUE              ~ "Mixed"
    )
  )

cat("  Enhancer linkage groups:\n")
for (grp in c("Concordant-linked", "Discordant-linked", "Mixed")) {
  cat(sprintf("    %s: %d\n", grp,
              sum(enh_concordance$linkage_group == grp)))
}

# =============================================================================
# JOIN K119UB DATA
# =============================================================================

cat("\nJoining K119ub signal...\n")

# Inner join: enhancers with both concordance classification and K119ub data
merged <- merge(enh_concordance, k119ub, by = "enh_key")
cat(sprintf("  Enhancers with K119ub + concordance: %d\n", nrow(merged)))

# Verify no NA Fold values
stopifnot(sum(is.na(merged$Fold)) == 0)

# =============================================================================
# JOIN ENHANCER CLASSES (overlap-based, same as step12)
# =============================================================================

cat("\nJoining enhancer classes via overlap...\n")

merged_gr <- GRanges(seqnames = merged$chr,
                     ranges = IRanges(start = merged$start, end = merged$end))
class_gr <- GRanges(seqnames = enh_class$chr,
                    ranges = IRanges(start = enh_class$start,
                                    end = enh_class$end))

hits <- findOverlaps(merged_gr, class_gr)
first_hits <- hits[!duplicated(queryHits(hits))]

merged$enhancer_class <- NA_character_
merged$enhancer_class[queryHits(first_hits)] <-
  enh_class$enhancer_class[subjectHits(first_hits)]

n_classified <- sum(!is.na(merged$enhancer_class))
cat(sprintf("  Classified: %d (%.1f%%)\n", n_classified,
            100 * n_classified / nrow(merged)))

merged$enhancer_class[is.na(merged$enhancer_class)] <- "Unclassified"
merged$enhancer_class <- factor(merged$enhancer_class,
                                levels = c(CLASS_ORDER, "Unclassified"))

merged$linkage_group <- factor(merged$linkage_group,
                               levels = c("Concordant-linked",
                                          "Discordant-linked", "Mixed"))

# =============================================================================
# PLOT 11: K119ub FOLD BY CONCORDANCE — BOXPLOT
# =============================================================================

cat("\n--- Plot 11: K119ub Fold by concordance linkage ---\n")

# Wilcoxon test: concordant-linked vs discordant-linked
conc_fold <- merged$Fold[merged$linkage_group == "Concordant-linked"]
disc_fold <- merged$Fold[merged$linkage_group == "Discordant-linked"]
wt_fold <- wilcox.test(conc_fold, disc_fold)
cat(sprintf("  Wilcoxon (Concordant vs Discordant): %s\n",
            fmt_p(wt_fold$p.value)))
cat(sprintf("  Median Fold — Concordant: %.4f | Discordant: %.4f | Mixed: %.4f\n",
            median(conc_fold),
            median(disc_fold),
            median(merged$Fold[merged$linkage_group == "Mixed"])))

p11 <- ggplot(merged, aes(x = linkage_group, y = Fold,
                          fill = linkage_group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "K119ub Change at Enhancers by Gene Concordance",
    subtitle = sprintf("Wilcoxon (Conc vs Disc): %s", fmt_p(wt_fold$p.value)),
    x = "Enhancer Linkage Group",
    y = "K119ub log2(Fold Change)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(size = 10))

save_plot(p11, "4B_k119ub_by_concordance_boxplot", w = 6, h = 5.5)  # Original: 11_k119ub_by_concordance_boxplot

# =============================================================================
# PLOT 12: K119ub SIGNIFICANCE RATE — BAR PLOT
# =============================================================================

cat("\n--- Plot 12: K119ub significance rate by concordance ---\n")

sig_rate <- merged %>%
  group_by(linkage_group) %>%
  summarise(
    n_total = n(),
    n_sig = sum(FDR < 0.05),
    pct_sig = 100 * n_sig / n_total,
    .groups = "drop"
  )

cat("  Significance rates (FDR < 0.05):\n")
for (i in seq_len(nrow(sig_rate))) {
  cat(sprintf("    %s: %d / %d (%.1f%%)\n",
              sig_rate$linkage_group[i], sig_rate$n_sig[i],
              sig_rate$n_total[i], sig_rate$pct_sig[i]))
}

# Fisher's exact: sig vs not-sig × concordant vs discordant
conc_disc_only <- merged[merged$linkage_group != "Mixed", ]
cont_sig <- table(
  conc_disc_only$linkage_group,
  conc_disc_only$FDR < 0.05
)
fisher_sig <- fisher.test(cont_sig)
cat(sprintf("  Fisher's exact (sig rate ~ linkage): %s, OR=%.2f\n",
            fmt_p(fisher_sig$p.value), fisher_sig$estimate))

p12 <- ggplot(sig_rate, aes(x = linkage_group, y = pct_sig,
                            fill = linkage_group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)",
                                pct_sig, n_sig, n_total)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "K119ub Significant Change Rate by Concordance",
    subtitle = sprintf("Fisher (Conc vs Disc): %s, OR = %.2f",
                       fmt_p(fisher_sig$p.value), fisher_sig$estimate),
    x = "Enhancer Linkage Group",
    y = "% Enhancers with FDR < 0.05"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(size = 10))

save_plot(p12, "4B_k119ub_significance_rate", w = 6, h = 5.5)  # Original: 12_k119ub_significance_rate

# =============================================================================
# PLOT 13: K119ub vs DELTA ACTIVITY SCATTER
# =============================================================================

cat("\n--- Plot 13: K119ub vs delta_activity scatter ---\n")

# Per-group Spearman correlations
group_cors <- merged %>%
  group_by(linkage_group) %>%
  summarise(
    rho = cor(Fold, delta_activity, method = "spearman",
              use = "pairwise.complete.obs"),
    n = n(),
    .groups = "drop"
  )
cat("  Spearman(K119ub Fold, delta_activity):\n")
for (i in seq_len(nrow(group_cors))) {
  cat(sprintf("    %s: rho=%.3f (n=%d)\n",
              group_cors$linkage_group[i],
              group_cors$rho[i], group_cors$n[i]))
}

# Clip outliers for axis limits
clip_x <- quantile(abs(merged$delta_activity), 0.995, na.rm = TRUE)
clip_y <- quantile(abs(merged$Fold), 0.995, na.rm = TRUE)

# Subtitle with all group rhos
rho_subtitle <- paste(
  sprintf("%s: rho=%.3f", group_cors$linkage_group, group_cors$rho),
  collapse = " | "
)

p13 <- ggplot(merged, aes(x = delta_activity, y = Fold,
                          color = linkage_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_color_manual(values = CONCORDANCE_COLORS, name = "Linkage Group") +
  coord_cartesian(xlim = c(-clip_x, clip_x), ylim = c(-clip_y, clip_y)) +
  labs(
    title = "K119ub Change vs Enhancer Activity Change",
    subtitle = rho_subtitle,
    x = expression(Delta * "Activity (KO - WT)"),
    y = "K119ub log2(Fold Change)"
  ) +
  theme_pub +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

save_plot(p13, "4B_k119ub_vs_delta_activity_scatter", w = 7, h = 6)  # Original: 13_k119ub_vs_delta_activity_scatter

# =============================================================================
# PLOT 14: K119ub FOLD × CONCORDANCE, FACETED BY ENHANCER CLASS
# =============================================================================

cat("\n--- Plot 14: K119ub Fold by concordance × enhancer class ---\n")

# Classified subset only
merged_cls <- merged[merged$enhancer_class != "Unclassified", ]
merged_cls$enhancer_class <- droplevels(merged_cls$enhancer_class)
cat(sprintf("  Classified enhancers for faceted plot: %d\n", nrow(merged_cls)))

# Per-class Wilcoxon tests
class_tests <- merged_cls %>%
  filter(linkage_group != "Mixed") %>%
  group_by(enhancer_class) %>%
  summarise(
    n_conc = sum(linkage_group == "Concordant-linked"),
    n_disc = sum(linkage_group == "Discordant-linked"),
    p_value = if (n_conc >= 3 & n_disc >= 3) {
      wilcox.test(
        Fold[linkage_group == "Concordant-linked"],
        Fold[linkage_group == "Discordant-linked"]
      )$p.value
    } else {
      NA_real_
    },
    median_conc = median(Fold[linkage_group == "Concordant-linked"]),
    median_disc = median(Fold[linkage_group == "Discordant-linked"]),
    .groups = "drop"
  )

cat("  Per-class Wilcoxon (Concordant vs Discordant):\n")
for (i in seq_len(nrow(class_tests))) {
  cat(sprintf("    %s: p=%s, median_conc=%.4f, median_disc=%.4f (n=%d/%d)\n",
              class_tests$enhancer_class[i],
              if (!is.na(class_tests$p_value[i])) fmt_p(class_tests$p_value[i]) else "NA (too few)",
              class_tests$median_conc[i], class_tests$median_disc[i],
              class_tests$n_conc[i], class_tests$n_disc[i]))
}

# Facet labels with p-values
facet_labels <- class_tests %>%
  mutate(label = sprintf("%s\n%s",
                         enhancer_class,
                         sapply(p_value, function(pv) {
                           if (is.na(pv)) "n/a" else fmt_p(pv)
                         })))
facet_named <- setNames(facet_labels$label, facet_labels$enhancer_class)

p14 <- ggplot(merged_cls,
              aes(x = linkage_group, y = Fold, fill = linkage_group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ enhancer_class, ncol = 2,
             labeller = labeller(enhancer_class = facet_named)) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "K119ub Change by Concordance Status Within Enhancer Classes",
    x = "Enhancer Linkage Group",
    y = "K119ub log2(Fold Change)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
        strip.text = element_text(size = 9))

save_plot(p14, "4B_k119ub_by_concordance_x_class", w = 9, h = 8)  # Original: 14_k119ub_by_concordance_x_class

# =============================================================================
# SAVE SUMMARY TABLE
# =============================================================================

cat("\nSaving summary table...\n")

summary_tbl <- merged %>%
  group_by(linkage_group) %>%
  summarise(
    n_enhancers = n(),
    median_Fold = median(Fold),
    mean_Fold = mean(Fold),
    pct_sig_FDR05 = 100 * sum(FDR < 0.05) / n(),
    median_delta_activity = median(delta_activity),
    median_mean_delta_abc = median(mean_delta_abc),
    .groups = "drop"
  )

write.table(summary_tbl,
            file.path(TSV_DIR, "4B_k119ub_concordance_summary.tsv"),  # Original: OUTPUT_DIR/k119ub_concordance_summary.tsv
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=== Summary ===\n")
cat(sprintf("  Enhancers analyzed: %d\n", nrow(merged)))
cat(sprintf("    Concordant-linked: %d | Discordant-linked: %d | Mixed: %d\n",
            sum(merged$linkage_group == "Concordant-linked"),
            sum(merged$linkage_group == "Discordant-linked"),
            sum(merged$linkage_group == "Mixed")))
cat(sprintf("  K119ub Fold Wilcoxon (Conc vs Disc): %s\n",
            fmt_p(wt_fold$p.value)))
cat(sprintf("  K119ub sig rate Fisher OR: %.2f (%s)\n",
            fisher_sig$estimate, fmt_p(fisher_sig$p.value)))
cat(sprintf("  Plots saved: 11, 12, 13, 14 in %s\n", OUTPUT_DIR))

cat(sprintf("\nCompleted: %s\n", Sys.time()))
cat("=== Step 13c complete ===\n")
