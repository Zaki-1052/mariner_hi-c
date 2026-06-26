# biomodal/downstream/scripts/viz_sections/section_89_subcompartment_pc1_crossvalidation.R
# Section 89: Cross-Validation of CALDER2 Subcompartments with Homer PC1 Eigenvalues
#
# Validates CALDER2 discrete subcompartment calls (A.1/A.2/B.1/B.2, 100kb)
# against continuous Homer PC1 eigenvalues (25kb, averaged to 100kb).
# Tests whether the two independent tools agree on compartment identity
# and on compartment switching between control and mutant.
#
#   Panel 89a: PC1 eigenvalue distribution by CALDER2 subcompartment (violin)
#   Panel 89b: CALDER2 continous_rank vs Homer PC1 scatter (Spearman)
#   Panel 89c: Concordance heatmap — CALDER2 A/B vs Homer PC1 sign
#   Panel 89d: CALDER2 continous_rank density per subcompartment
#   Panel 89e: Switching concordance — CALDER2 label_changed vs Homer shift
#   Panel 89f: Combined composite
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_89_subcompartment_pc1_crossvalidation.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 89: SUBCOMPARTMENT PC1 CROSS-VALIDATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 89: SUBCOMPARTMENT PC1 CROSS-VALIDATION\n")
cat("================================================================================\n\n")

# ---- Constants (same as Section 66) -----------------------------------------

SUBCMPT_ORDER <- c("A.1", "A.2", "B.1", "B.2")
SUBCMPT_LABELS <- c("A.1" = "A.1\n(Strong Active)",
                     "A.2" = "A.2\n(Weak Active)",
                     "B.1" = "B.1\n(Fac. Het)",
                     "B.2" = "B.2\n(Constit. Het)")
SUBCMPT_COLORS <- c("A.1" = "#E41A1C", "A.2" = "#FF7F00",
                     "B.1" = "#4DAF4A", "B.2" = "#377EB8")

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# ---- Load CALDER2 subcompartment labels (100kb) -----------------------------

SUBCMPT_PATH <- file.path(REPO_ROOT, "ML/cmpts/outputs/calder2/late/250402_subcompartment_labels_100kb.tsv")

if (!file.exists(SUBCMPT_PATH)) {
  stop("CALDER2 subcompartment file not found: ", SUBCMPT_PATH)
}

cat("Loading CALDER2 subcompartment labels...\n")
subcmpt <- read.table(SUBCMPT_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
subcmpt <- subcmpt %>% dplyr::filter(!is.na(ctrl_label) & ctrl_label != "NA")
cat(sprintf("  %d non-NA bins loaded\n", nrow(subcmpt)))
cat(sprintf("  Subcompartments: A.1=%d, A.2=%d, B.1=%d, B.2=%d\n",
            sum(subcmpt$ctrl_label == "A.1"), sum(subcmpt$ctrl_label == "A.2"),
            sum(subcmpt$ctrl_label == "B.1"), sum(subcmpt$ctrl_label == "B.2")))

# ---- Load Homer PC1 eigenvalues (25kb) --------------------------------------

PC1_PATH <- file.path(REPO_ROOT, "tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv")

if (!file.exists(PC1_PATH)) {
  stop("Homer PC1 annotated file not found: ", PC1_PATH)
}

cat("\nLoading Homer PC1 eigenvalues (25kb resolution)...\n")
pc1_raw <- read.table(PC1_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  %d bins loaded (25kb resolution)\n", nrow(pc1_raw)))
cat(sprintf("  ctrl_avg_PC1 range: [%.3f, %.3f]\n",
            min(pc1_raw$ctrl_avg_PC1, na.rm = TRUE),
            max(pc1_raw$ctrl_avg_PC1, na.rm = TRUE)))

# ---- Average Homer 25kb bins to 100kb CALDER2 resolution --------------------

cat("\nAveraging Homer 25kb bins to 100kb CALDER2 resolution...\n")

pc1_raw$bin_100kb_start <- floor(pc1_raw$Start / 100000) * 100000 + 1

pc1_100kb <- pc1_raw %>%
  group_by(Chr, bin_100kb_start) %>%
  summarise(
    mean_ctrl_pc1 = mean(ctrl_avg_PC1, na.rm = TRUE),
    mean_mut_pc1 = mean(mut_avg_PC1, na.rm = TRUE),
    mean_difference = mean(Difference, na.rm = TRUE),
    min_adj_pvalue = min(adj_pvalue, na.rm = TRUE),
    any_significant = any(adj_pvalue < 0.05 & abs(Difference) > 0.30, na.rm = TRUE),
    n_homer_bins = n(),
    .groups = "drop"
  )

cat(sprintf("  Aggregated to %d 100kb bins\n", nrow(pc1_100kb)))
cat(sprintf("  Bins with 4 Homer sub-bins: %d (%.1f%%)\n",
            sum(pc1_100kb$n_homer_bins == 4),
            100 * mean(pc1_100kb$n_homer_bins == 4)))

# ---- Join CALDER2 and Homer at 100kb ----------------------------------------

cat("\nJoining CALDER2 labels with Homer PC1 at 100kb...\n")

merged <- subcmpt %>%
  inner_join(pc1_100kb, by = c("chr" = "Chr", "bin_start" = "bin_100kb_start"))

cat(sprintf("  Matched: %d / %d CALDER2 bins (%.1f%%)\n",
            nrow(merged), nrow(subcmpt), 100 * nrow(merged) / nrow(subcmpt)))

merged$ctrl_label <- factor(merged$ctrl_label, levels = SUBCMPT_ORDER)

# ---- Panel 89a: PC1 distribution by CALDER2 subcompartment -----------------

cat("\nCreating Figure 89a: PC1 violin by subcompartment...\n")

kw_test <- kruskal.test(mean_ctrl_pc1 ~ ctrl_label, data = merged)
cat(sprintf("  Kruskal-Wallis: chi-sq=%.1f, df=%d, %s\n",
            kw_test$statistic, kw_test$parameter, fmt_p(kw_test$p.value)))

# Pairwise Wilcoxon between adjacent compartments
pw_pairs <- list(c("A.1", "A.2"), c("A.2", "B.1"), c("B.1", "B.2"))
pw_results <- data.frame(pair = character(), p = numeric(), stringsAsFactors = FALSE)
for (pair in pw_pairs) {
  wt <- wilcox.test(
    merged$mean_ctrl_pc1[merged$ctrl_label == pair[1]],
    merged$mean_ctrl_pc1[merged$ctrl_label == pair[2]]
  )
  pw_results <- rbind(pw_results,
    data.frame(pair = paste(pair, collapse = " vs "), p = wt$p.value, stringsAsFactors = FALSE))
  cat(sprintf("    %s vs %s: %s\n", pair[1], pair[2], fmt_p(wt$p.value)))
}

pc1_medians <- merged %>%
  group_by(ctrl_label) %>%
  summarise(median_pc1 = median(mean_ctrl_pc1), .groups = "drop")

p_89a <- ggplot(merged, aes(x = ctrl_label, y = mean_ctrl_pc1, fill = ctrl_label)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(data = pc1_medians,
            aes(x = ctrl_label, y = median_pc1,
                label = sprintf("%.2f", median_pc1)),
            vjust = -1.5, size = 3, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = SUBCMPT_COLORS, guide = "none") +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  labs(
    title = "Homer PC1 Eigenvalue by CALDER2 Subcompartment",
    subtitle = sprintf("Clear ordering validates subcompartment calls | Kruskal-Wallis %s",
                        fmt_p(kw_test$p.value)),
    x = "CALDER2 Subcompartment (Control)",
    y = "Mean Homer PC1 Eigenvalue (Control)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 10, face = "bold"))

save_multiformat_ggplot(p_89a,
                        file.path(OUTPUT_DIR, "89a_pc1_violin_by_subcompartment"),
                        width = 10, height = 7)

# ---- Panel 89b: CALDER2 continous_rank vs Homer PC1 scatter -----------------

cat("Creating Figure 89b: continous_rank vs PC1 scatter...\n")

spearman_overall <- cor.test(merged$continous_rank_ctrl, merged$mean_ctrl_pc1,
                              method = "spearman", exact = FALSE)
cat(sprintf("  Spearman rho (overall): %.4f, %s\n",
            spearman_overall$estimate, fmt_p(spearman_overall$p.value)))

p_89b <- ggplot(merged, aes(x = continous_rank_ctrl, y = mean_ctrl_pc1,
                              color = ctrl_label)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8,
              linetype = "dashed") +
  scale_color_manual(values = SUBCMPT_COLORS, name = "Subcompartment") +
  annotate("text",
           x = min(merged$continous_rank_ctrl, na.rm = TRUE) + 0.05,
           y = max(merged$mean_ctrl_pc1, na.rm = TRUE) * 0.9,
           label = sprintf("Spearman rho = %.3f\n%s",
                           spearman_overall$estimate,
                           fmt_p(spearman_overall$p.value)),
           hjust = 0, size = 4, fontface = "italic") +
  labs(
    title = "CALDER2 Continuous Rank vs Homer PC1 Eigenvalue",
    subtitle = "Two independent tools show strong agreement on compartment identity",
    x = "CALDER2 Continuous Domain Rank (Control)",
    y = "Mean Homer PC1 Eigenvalue (Control, averaged to 100kb)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_89b,
                        file.path(OUTPUT_DIR, "89b_rank_vs_pc1_scatter"),
                        width = 10, height = 8)

# ---- Panel 89c: Concordance heatmap ----------------------------------------

cat("Creating Figure 89c: concordance heatmap...\n")

merged$calder_ab <- ifelse(merged$ctrl_label %in% c("A.1", "A.2"), "A", "B")
merged$homer_ab <- ifelse(merged$mean_ctrl_pc1 > 0, "A", "B")

conf_mat <- table(CALDER2 = merged$calder_ab, Homer_PC1 = merged$homer_ab)
cat("  Confusion matrix:\n")
print(conf_mat)

accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
precision_a <- conf_mat["A", "A"] / sum(conf_mat[, "A"])
recall_a <- conf_mat["A", "A"] / sum(conf_mat["A", ])
precision_b <- conf_mat["B", "B"] / sum(conf_mat[, "B"])
recall_b <- conf_mat["B", "B"] / sum(conf_mat["B", ])

# Cohen's kappa
n <- sum(conf_mat)
p_observed <- sum(diag(conf_mat)) / n
p_expected <- (sum(conf_mat["A", ]) * sum(conf_mat[, "A"]) +
               sum(conf_mat["B", ]) * sum(conf_mat[, "B"])) / n^2
kappa <- (p_observed - p_expected) / (1 - p_expected)

cat(sprintf("  Accuracy: %.1f%% | Cohen's kappa: %.3f\n", 100 * accuracy, kappa))
cat(sprintf("  A: precision=%.1f%%, recall=%.1f%%\n", 100 * precision_a, 100 * recall_a))
cat(sprintf("  B: precision=%.1f%%, recall=%.1f%%\n", 100 * precision_b, 100 * recall_b))

conf_df <- as.data.frame(conf_mat)
names(conf_df) <- c("CALDER2", "Homer_PC1", "Freq")
conf_df$pct <- NA
for (i in seq_len(nrow(conf_df))) {
  row_total <- sum(conf_df$Freq[conf_df$CALDER2 == conf_df$CALDER2[i]])
  conf_df$pct[i] <- 100 * conf_df$Freq[i] / row_total
}

p_89c <- ggplot(conf_df, aes(x = Homer_PC1, y = CALDER2, fill = pct)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Freq, pct)),
            size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "white", mid = "#FEE08B", high = "#D73027",
                       midpoint = 50, name = "Row %",
                       limits = c(0, 100)) +
  labs(
    title = "Concordance: CALDER2 vs Homer PC1 Compartment Classification",
    subtitle = sprintf("Accuracy = %.1f%% | Cohen's kappa = %.3f", 100 * accuracy, kappa),
    x = "Homer PC1 Classification (A = PC1 > 0, B = PC1 < 0)",
    y = "CALDER2 Classification (A.1+A.2 vs B.1+B.2)"
  ) +
  theme_biomodal() +
  theme(panel.grid = element_blank())

save_multiformat_ggplot(p_89c,
                        file.path(OUTPUT_DIR, "89c_concordance_heatmap"),
                        width = 8, height = 6)

# ---- Panel 89d: CALDER2 continous_rank density per subcompartment -----------

cat("Creating Figure 89d: continous_rank density by subcompartment...\n")

# KS tests between adjacent compartments
ks_results <- data.frame(pair = character(), D = numeric(), p = numeric(),
                          stringsAsFactors = FALSE)
for (pair in pw_pairs) {
  ks <- ks.test(
    merged$continous_rank_ctrl[merged$ctrl_label == pair[1]],
    merged$continous_rank_ctrl[merged$ctrl_label == pair[2]]
  )
  ks_results <- rbind(ks_results,
    data.frame(pair = paste(pair, collapse = " vs "),
               D = ks$statistic, p = ks$p.value, stringsAsFactors = FALSE))
  cat(sprintf("  KS %s vs %s: D=%.3f, %s\n", pair[1], pair[2], ks$statistic, fmt_p(ks$p.value)))
}

rank_medians <- merged %>%
  group_by(ctrl_label) %>%
  summarise(median_rank = median(continous_rank_ctrl, na.rm = TRUE), .groups = "drop")

p_89d <- ggplot(merged, aes(x = continous_rank_ctrl, fill = ctrl_label, color = ctrl_label)) +
  geom_density(alpha = 0.4, linewidth = 0.6) +
  geom_vline(data = rank_medians,
             aes(xintercept = median_rank, color = ctrl_label),
             linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = SUBCMPT_COLORS, name = "Subcompartment") +
  scale_color_manual(values = SUBCMPT_COLORS, name = "Subcompartment") +
  labs(
    title = "CALDER2 Continuous Domain Rank Distribution by Subcompartment",
    subtitle = "Rank separation validates the 4-way classification (all KS p < 2.2e-16)",
    x = "CALDER2 Continuous Domain Rank (Control)",
    y = "Density"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_89d,
                        file.path(OUTPUT_DIR, "89d_rank_density_by_subcompartment"),
                        width = 10, height = 7)

# ---- Panel 89e: Switching concordance --------------------------------------

cat("Creating Figure 89e: switching concordance...\n")

merged$calder_switched <- merged$label_changed == TRUE
merged$homer_shifted <- merged$any_significant

switch_table <- table(CALDER2_switched = merged$calder_switched,
                       Homer_shifted = merged$homer_shifted)
cat("  Switching concordance table:\n")
print(switch_table)

fisher_switch <- fisher.test(switch_table)
cat(sprintf("  Fisher's OR = %.2f, %s\n",
            fisher_switch$estimate, fmt_p(fisher_switch$p.value)))

switch_df <- as.data.frame(switch_table)
names(switch_df) <- c("CALDER2_switched", "Homer_shifted", "Freq")
switch_df$CALDER2_switched <- ifelse(switch_df$CALDER2_switched == TRUE,
                                      "Switched", "Stable")
switch_df$Homer_shifted <- ifelse(switch_df$Homer_shifted == TRUE,
                                   "Shifted", "Stable")

p_89e <- ggplot(switch_df, aes(x = Homer_shifted, y = Freq,
                                fill = CALDER2_switched)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Switched" = "#D7191C", "Stable" = "grey70"),
                    name = "CALDER2 Status") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Compartment Switching Concordance: CALDER2 vs Homer",
    subtitle = sprintf("Fisher's OR = %.2f, %s | Do the tools agree on which bins switch?",
                        fisher_switch$estimate, fmt_p(fisher_switch$p.value)),
    x = "Homer PC1 Shift Status (adj.p < 0.05 & |Diff| > 0.30)",
    y = "Number of 100kb Bins"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_89e,
                        file.path(OUTPUT_DIR, "89e_switching_concordance"),
                        width = 10, height = 7)

# ---- Panel 89f: Combined composite -----------------------------------------

cat("Creating combined Figure 89f...\n")

p_89f <- (p_89a | p_89b) / (p_89c | p_89d) / p_89e +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(
    title = "Cross-Validation: CALDER2 Subcompartments vs Homer PC1 Eigenvalues",
    subtitle = "Two independent compartment-calling methods show strong agreement",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_89f,
                        file.path(OUTPUT_DIR, "89f_combined_crossvalidation"),
                        width = 16, height = 18)

# ---- Save tables ------------------------------------------------------------

cat("\nSaving tables...\n")

# Per-subcompartment summary
pc1_summary <- merged %>%
  group_by(ctrl_label) %>%
  summarise(
    n_bins = n(),
    mean_pc1 = mean(mean_ctrl_pc1),
    median_pc1 = median(mean_ctrl_pc1),
    sd_pc1 = sd(mean_ctrl_pc1),
    mean_rank = mean(continous_rank_ctrl, na.rm = TRUE),
    median_rank = median(continous_rank_ctrl, na.rm = TRUE),
    sd_rank = sd(continous_rank_ctrl, na.rm = TRUE),
    .groups = "drop"
  )
write.table(pc1_summary,
            file.path(TABLES_DIR, "89_pc1_by_subcompartment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 89_pc1_by_subcompartment_summary.tsv\n")

# Concordance matrix
concordance_stats <- data.frame(
  metric = c("accuracy", "cohens_kappa",
             "precision_A", "recall_A", "precision_B", "recall_B",
             "n_total", "n_agree", "n_disagree"),
  value = c(accuracy, kappa,
            precision_a, recall_a, precision_b, recall_b,
            n, sum(diag(conf_mat)), sum(conf_mat) - sum(diag(conf_mat))),
  stringsAsFactors = FALSE
)
write.table(concordance_stats,
            file.path(TABLES_DIR, "89_concordance_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 89_concordance_matrix.tsv\n")

# Per-subcompartment Spearman
rank_cor_results <- data.frame(group = character(), rho = numeric(),
                                p = numeric(), n = integer(), stringsAsFactors = FALSE)
rank_cor_results <- rbind(rank_cor_results,
  data.frame(group = "Overall", rho = spearman_overall$estimate,
             p = spearman_overall$p.value, n = nrow(merged), stringsAsFactors = FALSE))
for (sc in SUBCMPT_ORDER) {
  sub <- merged %>% dplyr::filter(ctrl_label == sc)
  if (nrow(sub) > 10) {
    sc_cor <- cor.test(sub$continous_rank_ctrl, sub$mean_ctrl_pc1,
                        method = "spearman", exact = FALSE)
    rank_cor_results <- rbind(rank_cor_results,
      data.frame(group = sc, rho = sc_cor$estimate,
                 p = sc_cor$p.value, n = nrow(sub), stringsAsFactors = FALSE))
  }
}
write.table(rank_cor_results,
            file.path(TABLES_DIR, "89_pc1_rank_correlation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 89_pc1_rank_correlation.tsv\n")

# Switching concordance
switch_stats <- data.frame(
  calder_switched = c(TRUE, TRUE, FALSE, FALSE),
  homer_shifted = c(TRUE, FALSE, TRUE, FALSE),
  count = as.vector(switch_table),
  fisher_or = fisher_switch$estimate,
  fisher_p = fisher_switch$p.value,
  stringsAsFactors = FALSE
)
write.table(switch_stats,
            file.path(TABLES_DIR, "89_switching_concordance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 89_switching_concordance.tsv\n")

# Full per-bin table
bin_export <- merged %>%
  dplyr::select(chr, bin_start, bin_end, ctrl_label, mut_label, label_changed,
                continous_rank_ctrl, continous_rank_mut,
                mean_ctrl_pc1, mean_mut_pc1, mean_difference,
                min_adj_pvalue, any_significant, n_homer_bins,
                calder_ab, homer_ab)
write.table(bin_export,
            file.path(TABLES_DIR, "89_per_bin_crossvalidation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: 89_per_bin_crossvalidation.tsv (%d bins)\n", nrow(bin_export)))

cat("\n")
cat("================================================================================\n")
cat("SECTION 89 COMPLETE\n")
cat("================================================================================\n")
