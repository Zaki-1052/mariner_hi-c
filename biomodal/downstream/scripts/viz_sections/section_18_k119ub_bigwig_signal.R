# biomodal/downstream/scripts/viz_sections/section_18_k119ub_bigwig_signal.R
# Section 18: K119ub Continuous Bigwig Signal Analysis
#
# Replaces binary peak-count approach (Section 17) with continuous signal
# quantification from merged bigwig files. Computes mean K119ub signal per gene,
# log2(mut/ctrl) ratios, and correlates with methylation changes.
#
# Requires: data/k119ub_gene_signal.tsv (from preprocess_k119ub_bigwig.R on HPC)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_18_k119ub_bigwig_signal.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 18: K119UB BIGWIG SIGNAL ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 18: K119UB CONTINUOUS BIGWIG SIGNAL ANALYSIS\n")
cat("================================================================================\n\n")

# =============================================================================
# LOAD PREPROCESSED BIGWIG DATA
# =============================================================================

k119ub_path <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")
if (!file.exists(k119ub_path)) {
  stop(paste0(
    "Preprocessed K119ub signal file not found: ", k119ub_path, "\n",
    "Run preprocess_k119ub_bigwig.R on HPC first, then scp the output:\n",
    "  scp expanse:/.../k119ub_gene_signal.tsv biomodal/downstream/data/"
  ))
}

cat("Loading preprocessed K119ub signal data...\n")
k119ub <- read.table(k119ub_path, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Loaded %d genes with bigwig signal\n", nrow(k119ub)))
cat(sprintf("  Gene body signal classes: %s\n",
            paste(names(table(k119ub$gb_signal_class)),
                  table(k119ub$gb_signal_class), sep = "=", collapse = ", ")))

# =============================================================================
# DERIVE METHYLATION GENE GROUPS (same as Section 17)
# =============================================================================

cat("\nPreparing methylation gene groups...\n")

mc_sig <- mc_dmr %>% dplyr::filter(significant)
hmc_sig <- hmc_dmr %>% dplyr::filter(significant)

mc_up_genes   <- mc_sig$gene[mc_sig$mod_difference > 0]
mc_down_genes <- mc_sig$gene[mc_sig$mod_difference < 0]
hmc_down_genes <- hmc_sig$gene[hmc_sig$mod_difference < 0]
hmc_up_genes   <- hmc_sig$gene[hmc_sig$mod_difference > 0]

all_dmr_genes <- unique(c(mc_up_genes, mc_down_genes, hmc_down_genes, hmc_up_genes))

cat(sprintf("  mC Up:    %d genes\n", length(mc_up_genes)))
cat(sprintf("  mC Down:  %d genes\n", length(mc_down_genes)))
cat(sprintf("  hmC Down: %d genes\n", length(hmc_down_genes)))
cat(sprintf("  hmC Up:   %d genes\n", length(hmc_up_genes)))
cat(sprintf("  All DMR:  %d unique genes\n", length(all_dmr_genes)))

# =============================================================================
# JOIN K119UB SIGNAL TO METHYLATION GROUPS
# =============================================================================

cat("\nJoining K119ub signal to methylation groups...\n")

build_group_df <- function(genes, group_label, k119ub_df) {
  df <- data.frame(gene = genes, met_group = group_label, stringsAsFactors = FALSE)
  merged <- df %>%
    left_join(k119ub_df, by = c("gene" = "symbol"))
  n_matched <- sum(!is.na(merged$gb_log2fc) | !is.na(merged$gb_ctrl_signal))
  cat(sprintf("  %-10s: %d / %d matched to K119ub data (%.1f%%)\n",
              group_label, n_matched, nrow(merged), 100 * n_matched / nrow(merged)))
  merged
}

mc_up_k <- build_group_df(mc_up_genes, "mC Up", k119ub)
mc_down_k <- build_group_df(mc_down_genes, "mC Down", k119ub)
hmc_down_k <- build_group_df(hmc_down_genes, "hmC Down", k119ub)
hmc_up_k <- build_group_df(hmc_up_genes, "hmC Up", k119ub)
all_dmr_k <- build_group_df(all_dmr_genes, "All DMR Genes", k119ub)

all_groups <- bind_rows(mc_up_k, mc_down_k, hmc_down_k, hmc_up_k)
all_groups$met_group <- factor(all_groups$met_group,
                                levels = c("mC Up", "mC Down", "hmC Down", "hmC Up"))

# =============================================================================
# FIGURE 18a: K119ub Signal Distribution (Ctrl vs Mut per Group)
# =============================================================================

cat("\n--- FIGURE 18a: K119ub Signal Distribution ---\n\n")

# Pivot to long format for paired violin: ctrl vs mut
# Deduplicate first (gene symbols can map to multiple entrez IDs)
signal_long <- all_groups %>%
  dplyr::filter(!is.na(gb_ctrl_signal)) %>%
  distinct(gene, met_group, .keep_all = TRUE) %>%
  dplyr::select(gene, met_group, gb_ctrl_signal, gb_mut_signal) %>%
  pivot_longer(cols = c(gb_ctrl_signal, gb_mut_signal),
               names_to = "condition", values_to = "signal") %>%
  mutate(
    condition = ifelse(condition == "gb_ctrl_signal", "Control", "Mutant"),
    condition = factor(condition, levels = c("Control", "Mutant"))
  )

# Wilcoxon signed-rank test per group (paired ctrl vs mut)
# Use all_groups directly (already has both columns) — deduplicate to one row per gene
wilcox_paired <- all_groups %>%
  dplyr::filter(!is.na(gb_ctrl_signal) & !is.na(gb_mut_signal) &
                gb_ctrl_signal > 0 & gb_mut_signal > 0) %>%
  distinct(gene, met_group, .keep_all = TRUE) %>%
  group_by(met_group) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(gb_mut_signal, gb_ctrl_signal, paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(p_label = ifelse(is.na(p_value), "NA",
                          ifelse(p_value < 0.001, sprintf("p=%.1e", p_value),
                                 sprintf("p=%.3f", p_value))))

# Y position for annotation
max_signal <- quantile(signal_long$signal[signal_long$signal > 0], 0.99, na.rm = TRUE)

p_18a <- ggplot(signal_long %>% dplyr::filter(signal > 0),
                aes(x = met_group, y = signal, fill = condition)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.6,
              scale = "width", linewidth = 0.3) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.15,
               outlier.shape = NA, alpha = 0.8, linewidth = 0.4) +
  # Wilcoxon annotations
  geom_text(data = wilcox_paired,
            aes(x = met_group, y = max_signal * 1.1, label = p_label),
            inherit.aes = FALSE, size = 2.8, color = "grey30", fontface = "italic") +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  scale_y_log10(labels = scales::label_number()) +
  labs(
    title = "H2AK119ub Signal at DMR Gene Bodies (Ctrl vs Mut)",
    subtitle = "Continuous bigwig signal; genes with signal > 0; Wilcoxon signed-rank test",
    x = "Methylation Group",
    y = "Mean K119ub Signal (log10 scale)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_18a,
                        file.path(OUTPUT_DIR, "18a_k119ub_signal_distribution"),
                        width = 12, height = 8)

cat("  Paired Wilcoxon signed-rank (ctrl vs mut within group):\n")
for (i in 1:nrow(wilcox_paired)) {
  cat(sprintf("    %-10s: n=%d, %s\n",
              wilcox_paired$met_group[i], wilcox_paired$n[i], wilcox_paired$p_label[i]))
}

# =============================================================================
# FIGURE 18b: log2FC Distribution per Methylation Group
# =============================================================================

cat("\n--- FIGURE 18b: K119ub log2FC Distributions ---\n\n")

# Include only quantifiable genes for log2FC plots
fc_data <- all_groups %>%
  dplyr::filter(gb_signal_class == "quantifiable") %>%
  dplyr::select(gene, met_group, gb_log2fc)

# Background: all DMR genes
fc_bg <- all_dmr_k %>%
  dplyr::filter(gb_signal_class == "quantifiable") %>%
  dplyr::select(gene, gb_log2fc) %>%
  mutate(met_group = "All DMR Genes")

fc_combined <- bind_rows(
  fc_data %>% mutate(met_group = as.character(met_group)),
  fc_bg
)
fc_combined$met_group <- factor(fc_combined$met_group,
  levels = c("mC Up", "mC Down", "hmC Down", "hmC Up", "All DMR Genes"))

# Group colors: methylation colors + grey for background
group_colors <- c(
  "mC Up"         = "#D7191C",  # hypermethylated red
  "mC Down"       = "#2C7BB6",  # hypomethylated blue
  "hmC Down"      = "#377EB8",  # 5hmC blue
  "hmC Up"        = "#E41A1C",  # 5hmC red
  "All DMR Genes" = "grey60"
)

# Statistics per group
fc_stats <- fc_combined %>%
  group_by(met_group) %>%
  summarise(
    n = n(),
    median_fc = median(gb_log2fc, na.rm = TRUE),
    mean_fc = mean(gb_log2fc, na.rm = TRUE),
    # One-sample Wilcoxon: log2FC != 0
    p_vs_zero = tryCatch(
      wilcox.test(gb_log2fc, mu = 0)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# Wilcoxon rank-sum: each group vs background
bg_fc <- fc_bg$gb_log2fc
fc_stats$p_vs_bg <- NA_real_
for (i in 1:nrow(fc_stats)) {
  grp <- as.character(fc_stats$met_group[i])
  if (grp == "All DMR Genes") next
  grp_fc <- fc_combined$gb_log2fc[fc_combined$met_group == grp]
  fc_stats$p_vs_bg[i] <- tryCatch(
    wilcox.test(grp_fc, bg_fc)$p.value,
    error = function(e) NA_real_
  )
}

# Annotation labels
fc_stats$label <- sprintf("med=%+.3f\nn=%d\np(0)=%s",
  fc_stats$median_fc, fc_stats$n,
  ifelse(fc_stats$p_vs_zero < 0.001,
         sprintf("%.1e", fc_stats$p_vs_zero),
         sprintf("%.3f", fc_stats$p_vs_zero)))

# Y range for annotations
fc_ymax <- max(abs(fc_combined$gb_log2fc), na.rm = TRUE) * 1.15

p_18b <- ggplot(fc_combined, aes(x = met_group, y = gb_log2fc, fill = met_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_violin(alpha = 0.5, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, linewidth = 0.4) +
  # Median + stats annotation
  geom_text(data = fc_stats,
            aes(x = met_group, y = fc_ymax, label = label),
            inherit.aes = FALSE, size = 2.5, lineheight = 0.85,
            color = "grey20", vjust = 1) +
  scale_fill_manual(values = group_colors, guide = "none") +
  coord_cartesian(ylim = c(-fc_ymax, fc_ymax)) +
  labs(
    title = "K119ub log2 Fold Change at DMR Gene Bodies",
    subtitle = "Continuous signal; quantifiable genes only; dashed line = no change",
    x = "Methylation Group",
    y = "K119ub log2(Mutant / Control)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_18b,
                        file.path(OUTPUT_DIR, "18b_k119ub_log2fc_distributions"),
                        width = 10, height = 7)

cat("  log2FC statistics per group:\n")
cat(sprintf("  %-14s %6s %9s %9s %12s %12s\n",
            "Group", "n", "median", "mean", "p(vs 0)", "p(vs BG)"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in 1:nrow(fc_stats)) {
  cat(sprintf("  %-14s %6d %+9.4f %+9.4f %12s %12s\n",
              fc_stats$met_group[i], fc_stats$n[i],
              fc_stats$median_fc[i], fc_stats$mean_fc[i],
              ifelse(is.na(fc_stats$p_vs_zero[i]), "NA",
                     sprintf("%.2e", fc_stats$p_vs_zero[i])),
              ifelse(is.na(fc_stats$p_vs_bg[i]), "---",
                     sprintf("%.2e", fc_stats$p_vs_bg[i]))))
}

# =============================================================================
# FIGURE 18c: Methylation Change vs K119ub Signal Change (Scatter)
# =============================================================================

cat("\n--- FIGURE 18c: Methylation vs K119ub Scatter ---\n\n")

# Build scatter data: merge DMR mod_difference with K119ub log2FC
build_scatter <- function(dmr_df, modality_label, k119ub_df) {
  dmr_df %>%
    dplyr::filter(significant) %>%
    dplyr::select(gene, mod_difference) %>%
    inner_join(
      k119ub_df %>%
        dplyr::filter(gb_signal_class == "quantifiable") %>%
        dplyr::select(symbol, gb_log2fc),
      by = c("gene" = "symbol")
    ) %>%
    mutate(modality = modality_label)
}

scatter_mc  <- build_scatter(mc_dmr, "5mC", k119ub)
scatter_hmc <- build_scatter(hmc_dmr, "5hmC", k119ub)
scatter_all <- bind_rows(scatter_mc, scatter_hmc)
scatter_all$modality <- factor(scatter_all$modality, levels = c("5mC", "5hmC"))

cat(sprintf("  5mC genes with quantifiable K119ub signal: %d\n", nrow(scatter_mc)))
cat(sprintf("  5hmC genes with quantifiable K119ub signal: %d\n", nrow(scatter_hmc)))

# Spearman correlations per modality
cor_results <- scatter_all %>%
  group_by(modality) %>%
  summarise(
    n = n(),
    rho = cor(mod_difference, gb_log2fc, method = "spearman"),
    p_value = cor.test(mod_difference, gb_log2fc, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("rho = %+.3f\np = %s\nn = %d",
                          rho,
                          ifelse(p_value < 0.001, sprintf("%.1e", p_value),
                                 sprintf("%.3f", p_value)),
                          n))

# Label KEY_GENES
scatter_all$is_key <- scatter_all$gene %in% KEY_GENES

p_18c <- ggplot(scatter_all, aes(x = mod_difference, y = gb_log2fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.15, size = 0.8, color = "grey40") +
  geom_smooth(method = "lm", color = "#E41A1C", fill = "#E41A1C",
              alpha = 0.15, linewidth = 0.8) +
  # Key gene labels
  {if (any(scatter_all$is_key))
    ggrepel::geom_text_repel(
      data = scatter_all %>% dplyr::filter(is_key),
      aes(label = gene),
      size = 2.5, fontface = "italic", color = "black",
      max.overlaps = 15, nudge_y = 0.1,
      segment.size = 0.2, segment.color = "grey50"
    )
  } +
  # Correlation annotation
  geom_text(data = cor_results,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.3, size = 3.2, color = "black",
            lineheight = 0.9) +
  facet_wrap(~ modality, scales = "free_x") +
  labs(
    title = "Methylation Change vs K119ub Signal Change",
    subtitle = "Gene body mod_difference vs K119ub log2FC; Spearman correlation; quantifiable genes",
    x = "Methylation Difference (Mutant - Control)",
    y = "K119ub log2(Mutant / Control)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_18c,
                        file.path(OUTPUT_DIR, "18c_methylation_vs_k119ub_scatter"),
                        width = 14, height = 7)

cat("  Spearman correlations:\n")
for (i in 1:nrow(cor_results)) {
  cat(sprintf("    %s: rho = %+.4f, p = %s, n = %d\n",
              cor_results$modality[i], cor_results$rho[i],
              ifelse(cor_results$p_value[i] < 0.001,
                     sprintf("%.2e", cor_results$p_value[i]),
                     sprintf("%.4f", cor_results$p_value[i])),
              cor_results$n[i]))
}

# =============================================================================
# FIGURE 18d: Gene Waterfall (Ranked by K119ub Change)
# =============================================================================

cat("\n--- FIGURE 18d: K119ub Waterfall ---\n\n")

# All genes with quantifiable signal, ranked by gb_log2fc
waterfall_base <- k119ub %>%
  dplyr::filter(gb_signal_class == "quantifiable") %>%
  dplyr::select(symbol, gb_log2fc) %>%
  arrange(desc(gb_log2fc)) %>%
  mutate(rank = row_number())

# Tag DMR group membership (a gene can appear in multiple groups)
waterfall_base$dmr_group <- "None"
waterfall_base$dmr_group[waterfall_base$symbol %in% mc_up_genes]   <- "mC Up"
waterfall_base$dmr_group[waterfall_base$symbol %in% mc_down_genes] <- "mC Down"
waterfall_base$dmr_group[waterfall_base$symbol %in% hmc_down_genes] <- "hmC Down"
waterfall_base$dmr_group[waterfall_base$symbol %in% hmc_up_genes]  <- "hmC Up"

# For genes in multiple groups, prioritize: mC Up > hmC Down > mC Down > hmC Up
# (re-assign in reverse priority so highest priority overwrites)
waterfall_base$dmr_group[waterfall_base$symbol %in% hmc_up_genes]   <- "hmC Up"
waterfall_base$dmr_group[waterfall_base$symbol %in% mc_down_genes]  <- "mC Down"
waterfall_base$dmr_group[waterfall_base$symbol %in% hmc_down_genes] <- "hmC Down"
waterfall_base$dmr_group[waterfall_base$symbol %in% mc_up_genes]    <- "mC Up"

waterfall_base$dmr_group <- factor(waterfall_base$dmr_group,
  levels = c("mC Up", "mC Down", "hmC Down", "hmC Up", "None"))

waterfall_colors <- c(
  "mC Up"    = "#D7191C",
  "mC Down"  = "#2C7BB6",
  "hmC Down" = "#377EB8",
  "hmC Up"   = "#E41A1C",
  "None"     = "grey85"
)

# Compute where DMR genes cluster
n_total <- nrow(waterfall_base)
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  ranks <- waterfall_base$rank[waterfall_base$dmr_group == grp]
  if (length(ranks) > 0) {
    cat(sprintf("  %s: n=%d, median rank=%d / %d (percentile=%.1f%%)\n",
                grp, length(ranks), median(ranks), n_total,
                100 * median(ranks) / n_total))
  }
}

# Plot: draw "None" first (background), then DMR groups on top
waterfall_dmr <- waterfall_base %>% dplyr::filter(dmr_group != "None")
waterfall_none <- waterfall_base %>% dplyr::filter(dmr_group == "None")

p_18d <- ggplot() +
  # Background genes
  geom_point(data = waterfall_none,
             aes(x = rank, y = gb_log2fc),
             color = "grey85", size = 0.3, alpha = 0.5) +
  # DMR genes on top
  geom_point(data = waterfall_dmr,
             aes(x = rank, y = gb_log2fc, color = dmr_group),
             size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = waterfall_colors, name = "DMR Group",
                     breaks = c("mC Up", "mC Down", "hmC Down", "hmC Up")) +
  labs(
    title = "Gene Waterfall: All Genes Ranked by K119ub Change",
    subtitle = "Grey = non-DMR genes; colored = significant DMR genes",
    x = sprintf("Gene Rank (1 = highest K119ub gain, %d = highest loss)", n_total),
    y = "K119ub log2(Mutant / Control)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_18d,
                        file.path(OUTPUT_DIR, "18d_k119ub_waterfall"),
                        width = 12, height = 6)

# =============================================================================
# STATISTICAL TESTS (console output)
# =============================================================================

cat("\n--- STATISTICAL TESTS ---\n\n")

# Cliff's delta (non-parametric effect size)
cliffs_delta <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)
  # Count dominance pairs
  more <- sum(outer(x, y, ">"))
  less <- sum(outer(x, y, "<"))
  (more - less) / (n_x * n_y)
}

cat("1. Wilcoxon rank-sum: each group's log2FC vs All DMR Genes\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  grp_fc <- fc_combined$gb_log2fc[fc_combined$met_group == grp]
  w <- wilcox.test(grp_fc, bg_fc)
  cat(sprintf("   %-10s: W=%.0f, p=%s\n", grp, w$statistic,
              ifelse(w$p.value < 0.001, sprintf("%.2e", w$p.value),
                     sprintf("%.4f", w$p.value))))
}

cat("\n2. One-sample Wilcoxon: each group's log2FC != 0\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up", "All DMR Genes")) {
  grp_fc <- fc_combined$gb_log2fc[fc_combined$met_group == grp]
  w <- wilcox.test(grp_fc, mu = 0)
  cat(sprintf("   %-14s: V=%.0f, p=%s, median=%+.4f\n", grp, w$statistic,
              ifelse(w$p.value < 0.001, sprintf("%.2e", w$p.value),
                     sprintf("%.4f", w$p.value)),
              median(grp_fc)))
}

cat("\n3. Spearman correlation: mod_difference vs gb_log2fc\n")
for (mod in c("5mC", "5hmC")) {
  sub <- scatter_all %>% dplyr::filter(modality == mod)
  ct <- cor.test(sub$mod_difference, sub$gb_log2fc, method = "spearman")
  cat(sprintf("   %s: rho=%+.4f, S=%.0f, p=%s, n=%d\n",
              mod, ct$estimate, ct$statistic,
              ifelse(ct$p.value < 0.001, sprintf("%.2e", ct$p.value),
                     sprintf("%.4f", ct$p.value)),
              nrow(sub)))
}

cat("\n4. Cliff's delta (effect size): each group vs All DMR Genes\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  grp_fc <- fc_combined$gb_log2fc[fc_combined$met_group == grp]
  d <- cliffs_delta(grp_fc, bg_fc)
  magnitude <- if (abs(d) < 0.147) "negligible"
               else if (abs(d) < 0.33) "small"
               else if (abs(d) < 0.474) "medium"
               else "large"
  cat(sprintf("   %-10s: delta=%+.4f (%s)\n", grp, d, magnitude))
}

# =============================================================================
# SAVE SUMMARY TABLE
# =============================================================================

cat("\n--- Saving summary table ---\n\n")

summary_table <- fc_stats %>%
  dplyr::select(met_group, n, median_fc, mean_fc, p_vs_zero, p_vs_bg)

# Add Cliff's delta
summary_table$cliffs_delta <- NA_real_
for (i in 1:nrow(summary_table)) {
  grp <- as.character(summary_table$met_group[i])
  if (grp == "All DMR Genes") next
  grp_fc <- fc_combined$gb_log2fc[fc_combined$met_group == grp]
  summary_table$cliffs_delta[i] <- cliffs_delta(grp_fc, bg_fc)
}

write.table(summary_table,
            file.path(TABLES_DIR, "k119ub_bigwig_signal_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: k119ub_bigwig_signal_summary.tsv\n")

# =============================================================================
# SECTION 18 SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 18 SUMMARY: K119UB CONTINUOUS SIGNAL ANALYSIS\n")
cat("================================================================================\n\n")

cat("Figures generated:\n")
cat("  18a: K119ub signal distribution (ctrl vs mut per group)\n")
cat("  18b: K119ub log2FC distributions per methylation group\n")
cat("  18c: Methylation change vs K119ub signal change scatter\n")
cat("  18d: Gene waterfall ranked by K119ub change\n")

cat("\nKey results:\n")
for (i in 1:nrow(fc_stats)) {
  if (as.character(fc_stats$met_group[i]) != "All DMR Genes") {
    cat(sprintf("  %-10s: median log2FC = %+.4f (n=%d)\n",
                fc_stats$met_group[i], fc_stats$median_fc[i], fc_stats$n[i]))
  }
}
cat(sprintf("  Background: median log2FC = %+.4f (n=%d)\n",
            fc_stats$median_fc[fc_stats$met_group == "All DMR Genes"],
            fc_stats$n[fc_stats$met_group == "All DMR Genes"]))

cat("\nCorrelations (mod_difference vs K119ub log2FC):\n")
for (i in 1:nrow(cor_results)) {
  cat(sprintf("  %s: rho = %+.4f\n", cor_results$modality[i], cor_results$rho[i]))
}

cat("\nSection 18 complete.\n\n")
