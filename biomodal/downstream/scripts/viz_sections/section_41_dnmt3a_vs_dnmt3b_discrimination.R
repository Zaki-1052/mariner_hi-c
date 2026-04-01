# biomodal/downstream/scripts/viz_sections/section_41_dnmt3a_vs_dnmt3b_discrimination.R
# DNMT3A vs DNMT3B Mechanistic Discrimination (Triple Integration)
#
# Integrates H2AK119ub, H3K36me2/me3, and methylation data to test whether
# gene-body hypermethylation is driven by:
#   a) DNMT3A recruited via H2AK119ub (Gretarsson 2024)
#   b) DNMT3A recruited via H3K36me2 (NSD-dependent)
#   c) DNMT3B recruited via H3K36me3 (SETD2-dependent)
#   d) Convergent pathways
#
# Figures: 41a-41h (8 total)
#
# Usage:
#   cd downstream/
#   Rscript scripts/viz_sections/section_41_dnmt3a_vs_dnmt3b_discrimination.R

source("scripts/viz_sections/_shared_config.R")
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(pROC)
})

SECTION_DIR <- "41_dnmt3a_vs_dnmt3b"
FDR_THRESHOLD <- 0.05

dir.create(file.path(OUTPUT_DIR, SECTION_DIR), recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("================================================================================\n")
cat("SECTION 41: DNMT3A vs DNMT3B Mechanistic Discrimination\n")
cat("================================================================================\n\n")

# =============================================================================
# DATA LOADING — Comprehensive gene-level profile
# =============================================================================

COMBINED_TABLE <- file.path(TABLES_DIR, "h3k36_combined_gene_profile.tsv")
RATIO_TABLE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")

stopifnot(file.exists(RATIO_TABLE))
ratio_data <- read.table(RATIO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Load combined H3K36 profile (from section 40) or build from parts
if (file.exists(COMBINED_TABLE)) {
  full_profile <- read.table(COMBINED_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("  Combined profile loaded: %d genes, %d columns\n", nrow(full_profile), ncol(full_profile)))
} else {
  cat("  Combined table not found — building from individual tables...\n")

  ME3_TABLE <- file.path(TABLES_DIR, "h3k36me3_gene_level_summary.tsv")
  ME2_TABLE <- file.path(TABLES_DIR, "h3k36me2_gene_level_summary.tsv")
  stopifnot(file.exists(ME3_TABLE))
  stopifnot(file.exists(ME2_TABLE))

  me3_gene <- read.table(ME3_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  me2_gene <- read.table(ME2_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  full_profile <- ratio_data %>%
    left_join(me3_gene %>% dplyr::select(gene, me3_fold, me3_fdr, me3_n_peaks), by = "gene") %>%
    left_join(me2_gene %>% dplyr::select(gene, me2_fold, me2_fdr, me2_n_peaks), by = "gene")
}

# Ensure K119ub FDR data is present (combined table may have fold but not fdr)
if (!"k119ub_fdr" %in% colnames(full_profile)) {
  multi_mark_table <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
  stopifnot(file.exists(multi_mark_table))
  multi_marks <- read.table(multi_mark_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  full_profile <- full_profile %>%
    left_join(multi_marks %>% dplyr::select(gene, any_of(c("atac_fold", "k27ac_fold",
                                                            "k27me3_fold", "k119ub_fold",
                                                            "atac_fdr", "k27ac_fdr",
                                                            "k27me3_fdr", "k119ub_fdr"))),
              by = "gene")
}

# Define key binary classifications
full_profile <- full_profile %>%
  mutate(
    hyper_dmr = mc_sig & mc_diff > 0,
    hypo_dmr  = mc_sig & mc_diff <= 0,
    k119ub_gained = !is.na(k119ub_fold) & !is.na(k119ub_fdr) & k119ub_fdr < FDR_THRESHOLD & k119ub_fold > 0,
    me3_gained = !is.na(me3_fold) & !is.na(me3_fdr) & me3_fdr < FDR_THRESHOLD & me3_fold > 0,
    me3_lost   = !is.na(me3_fold) & !is.na(me3_fdr) & me3_fdr < FDR_THRESHOLD & me3_fold < 0,
    me2_gained = !is.na(me2_fold) & !is.na(me2_fdr) & me2_fdr < FDR_THRESHOLD & me2_fold > 0,
    me2_lost   = !is.na(me2_fold) & !is.na(me2_fdr) & me2_fdr < FDR_THRESHOLD & me2_fold < 0,
    coordinated = mc_sig & hmc_sig & mc_diff > 0 & hmc_diff < 0
  )

n_hyper <- sum(full_profile$hyper_dmr)
n_k119ub <- sum(full_profile$k119ub_gained, na.rm = TRUE)
cat(sprintf("  Hypermethylated genes: %d\n", n_hyper))
cat(sprintf("  K119ub gained genes: %d\n", n_k119ub))
cat(sprintf("  me3 gained genes: %d | me3 lost genes: %d\n",
            sum(full_profile$me3_gained), sum(full_profile$me3_lost)))

# =============================================================================
# FIGURE 41a: CRITICAL TEST — me3 Fold at Hypermethylated Genes WITH K119ub Gain
# =============================================================================

cat("\n--- FIGURE 41a: me3 Fold at Hyper Genes with/without K119ub ---\n\n")

violin_41a <- full_profile %>%
  dplyr::filter(!is.na(me3_fold)) %>%
  mutate(
    strat_group = case_when(
      hyper_dmr & k119ub_gained ~ "Hyper + K119ub\nGained",
      hyper_dmr & !k119ub_gained ~ "Hyper + No\nK119ub Gain",
      !hyper_dmr ~ "Non-Hyper"
    ),
    strat_group = factor(strat_group, levels = c("Hyper + K119ub\nGained",
                                                  "Hyper + No\nK119ub Gain", "Non-Hyper"))
  )

group_n <- violin_41a %>% count(strat_group)
cat("  Group sizes:\n")
for (i in seq_len(nrow(group_n))) {
  cat(sprintf("    %s: %d\n", gsub("\n", " ", group_n$strat_group[i]), group_n$n[i]))
}

# Pairwise Wilcoxon tests
groups <- levels(violin_41a$strat_group)
wilcox_pairs <- list()
for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    vals_i <- violin_41a$me3_fold[violin_41a$strat_group == groups[i]]
    vals_j <- violin_41a$me3_fold[violin_41a$strat_group == groups[j]]
    if (length(vals_i) >= 3 & length(vals_j) >= 3) {
      wt <- wilcox.test(vals_i, vals_j)
      pair_name <- paste(gsub("\n", " ", groups[i]), "vs", gsub("\n", " ", groups[j]))
      wilcox_pairs[[pair_name]] <- wt$p.value
      cat(sprintf("    %s: p = %.2e\n", pair_name, wt$p.value))
    }
  }
}

p_41a <- ggplot(violin_41a, aes(x = strat_group, y = me3_fold, fill = strat_group)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("#B2182B", "#E6AB02", "grey80")) +
  geom_text(data = group_n, aes(x = strat_group, y = -Inf, label = sprintf("n=%d", n)),
            inherit.aes = FALSE, vjust = -0.5, size = 3, color = "grey40") +
  labs(
    title = "H3K36me3 Fold at Hypermethylated Genes: K119ub Stratification",
    subtitle = paste0("If me3 LOST at K119ub-gained genes: DNMT3A (not DNMT3B) is primary effector\n",
                      "If me3 MAINTAINED/GAINED: both DNMT3A and DNMT3B may contribute"),
    x = "", y = "H3K36me3 Log2 Fold Change"
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 9, lineheight = 1.1))

save_multiformat_ggplot(p_41a, file.path(OUTPUT_DIR, SECTION_DIR, "41a_me3_at_hyper_k119ub_stratified"),
                        width = 9, height = 7)

# =============================================================================
# FIGURE 41b: me2 + me3 at Hyper Genes WITHOUT K119ub Gain
# =============================================================================

cat("\n--- FIGURE 41b: me2/me3 at Hyper Genes WITHOUT K119ub ---\n\n")

independent_df <- full_profile %>%
  dplyr::filter(!is.na(me3_fold) & !is.na(me2_fold)) %>%
  mutate(
    group = case_when(
      hyper_dmr & !k119ub_gained ~ "Hyper\n(No K119ub)",
      !hyper_dmr ~ "Non-Hyper"
    )
  ) %>%
  dplyr::filter(!is.na(group))

# Pivot for dual violin
indep_long <- independent_df %>%
  dplyr::select(gene, group, me2_fold, me3_fold) %>%
  pivot_longer(cols = c(me2_fold, me3_fold), names_to = "mark", values_to = "fold") %>%
  mutate(mark = recode(mark, "me2_fold" = "H3K36me2", "me3_fold" = "H3K36me3"))

# Wilcoxon tests
for (mark in c("H3K36me2", "H3K36me3")) {
  hyper_vals <- indep_long$fold[indep_long$group == "Hyper\n(No K119ub)" & indep_long$mark == mark]
  nonhyper_vals <- indep_long$fold[indep_long$group == "Non-Hyper" & indep_long$mark == mark]
  if (length(hyper_vals) >= 3 & length(nonhyper_vals) >= 3) {
    wt <- wilcox.test(hyper_vals, nonhyper_vals)
    cat(sprintf("  %s: Hyper(no K119ub) vs Non-Hyper p = %.2e\n", mark, wt$p.value))
  }
}

p_41b <- ggplot(indep_long, aes(x = group, y = fold, fill = mark)) +
  geom_violin(alpha = 0.5, scale = "width", position = position_dodge(width = 0.7)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7,
               position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("H3K36me2" = "#E6AB02", "H3K36me3" = "#D95F02"),
                    name = "Mark") +
  labs(
    title = "H3K36me2/me3 at Hypermethylated Genes WITHOUT K119ub Gain",
    subtitle = "Tests whether H3K36-mediated DNMT recruitment operates independently of ubiquitin",
    x = "", y = "Log2 Fold Change"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_41b, file.path(OUTPUT_DIR, SECTION_DIR, "41b_me2_me3_independent_pathway"),
                        width = 9, height = 7)

# =============================================================================
# FIGURE 41c: 6-Mark Logistic Regression Forest Plot
# =============================================================================

cat("\n--- FIGURE 41c: 6-Mark Logistic Regression ---\n\n")

model_data <- full_profile %>%
  mutate(hyper_binary = as.integer(hyper_dmr)) %>%
  dplyr::select(gene, hyper_binary, me2_fold, me3_fold,
                any_of(c("atac_fold", "k27ac_fold", "k27me3_fold", "k119ub_fold")))

# Check which predictors we have
predictor_cols <- intersect(c("atac_fold", "k27ac_fold", "k27me3_fold", "k119ub_fold",
                               "me2_fold", "me3_fold"), colnames(model_data))

model_data_cc <- model_data %>%
  dplyr::filter(complete.cases(across(all_of(c("hyper_binary", predictor_cols)))))

cat(sprintf("  Complete cases for %d-mark model: %d genes\n",
            length(predictor_cols), nrow(model_data_cc)))

if (nrow(model_data_cc) >= 50 & length(predictor_cols) >= 4) {
  # Null model
  null_model <- glm(hyper_binary ~ 1, data = model_data_cc, family = binomial)
  null_ll <- as.numeric(logLik(null_model))

  mcfadden <- function(model) { 1 - as.numeric(logLik(model)) / null_ll }

  extract_or <- function(model) {
    cc <- confint.default(model)
    coefs <- coef(model)
    data.frame(
      term = names(coefs),
      estimate = coefs,
      or = exp(coefs),
      or_lower = exp(cc[, 1]),
      or_upper = exp(cc[, 2]),
      p_value = summary(model)$coefficients[, 4],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  # 4-mark model (existing marks only)
  base_predictors <- intersect(c("atac_fold", "k27ac_fold", "k27me3_fold", "k119ub_fold"),
                                predictor_cols)
  if (length(base_predictors) == 4) {
    formula_4 <- as.formula(paste("hyper_binary ~", paste(base_predictors, collapse = " + ")))
    model_4 <- glm(formula_4, data = model_data_cc, family = binomial)
    roc_4 <- roc(model_data_cc$hyper_binary, predict(model_4, type = "response"), quiet = TRUE)
    auc_4 <- auc(roc_4)
    ci_4 <- ci.auc(roc_4, quiet = TRUE)
  }

  # Full model (with me2 + me3)
  formula_full <- as.formula(paste("hyper_binary ~", paste(predictor_cols, collapse = " + ")))
  model_full <- glm(formula_full, data = model_data_cc, family = binomial)
  roc_full <- roc(model_data_cc$hyper_binary, predict(model_full, type = "response"), quiet = TRUE)
  auc_full <- auc(roc_full)
  ci_full <- ci.auc(roc_full, quiet = TRUE)

  cat(sprintf("  4-mark AUC: %.3f [%.3f, %.3f]\n", auc_4, ci_4[1], ci_4[3]))
  cat(sprintf("  %d-mark AUC: %.3f [%.3f, %.3f]\n", length(predictor_cols), auc_full, ci_full[1], ci_full[3]))
  cat(sprintf("  AUC improvement: %+.3f\n", auc_full - auc_4))

  # Forest plot
  display_names <- c(
    "atac_fold" = "ATAC-seq", "k27ac_fold" = "H3K27ac", "k27me3_fold" = "H3K27me3",
    "k119ub_fold" = "H2AK119ub", "me2_fold" = "H3K36me2", "me3_fold" = "H3K36me3"
  )

  or_df <- extract_or(model_full) %>%
    dplyr::filter(term != "(Intercept)") %>%
    mutate(
      display_name = display_names[term],
      sig_label = ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                  ifelse(p_value < 0.05, "*", "ns"))),
      is_new = term %in% c("me2_fold", "me3_fold"),
      display_name = factor(display_name, levels = rev(display_names[predictor_cols]))
    )

  p_41c <- ggplot(or_df, aes(x = or, y = display_name)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbarh(aes(xmin = or_lower, xmax = or_upper), height = 0.2, linewidth = 0.8) +
    geom_point(aes(color = is_new), size = 3.5) +
    geom_text(aes(label = sprintf("OR=%.2f %s", or, sig_label)),
              hjust = -0.2, size = 3.5) +
    scale_color_manual(values = c("FALSE" = "#E41A1C", "TRUE" = "#D95F02"),
                       labels = c("Existing marks", "H3K36me2/me3 (new)"),
                       name = "") +
    scale_x_log10() +
    labs(
      title = sprintf("%d-Mark Logistic Regression: Predicting mC Hypermethylation", length(predictor_cols)),
      subtitle = sprintf("N=%s genes | 4-mark AUC=%.3f | %d-mark AUC=%.3f (\u0394=%+.3f)",
                          format(nrow(model_data_cc), big.mark = ","),
                          auc_4, length(predictor_cols), auc_full, auc_full - auc_4),
      x = "Odds Ratio (log scale)", y = ""
    ) +
    theme_biomodal() +
    theme(axis.text.y = element_text(size = 12, face = "bold"))

  save_multiformat_ggplot(p_41c, file.path(OUTPUT_DIR, SECTION_DIR, "41c_logistic_regression_forest"),
                          width = 11, height = 7)

  # Export model coefficients
  write.table(or_df, file.path(TABLES_DIR, "extended_logistic_regression_6mark.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  cat("  Insufficient data for logistic regression\n")
}

# =============================================================================
# FIGURE 41d: Stratified Bar — Hypermethylation Rate by Pathway Group
# =============================================================================

cat("\n--- FIGURE 41d: Hypermethylation Rate by Pathway Group ---\n\n")

pathway_df <- full_profile %>%
  mutate(
    pathway_group = case_when(
      k119ub_gained & (me3_gained | me2_gained) ~ "Convergent\n(K119ub + H3K36)",
      k119ub_gained & !me3_gained & !me2_gained ~ "K119ub Only",
      !k119ub_gained & (me3_gained | me2_gained) ~ "H3K36 Only",
      TRUE ~ "Neither"
    ),
    pathway_group = factor(pathway_group,
                           levels = c("K119ub Only", "H3K36 Only",
                                      "Convergent\n(K119ub + H3K36)", "Neither"))
  )

pathway_summary <- pathway_df %>%
  group_by(pathway_group) %>%
  summarise(
    n_total = n(),
    n_hyper = sum(hyper_dmr),
    pct_hyper = 100 * n_hyper / n_total,
    .groups = "drop"
  )

cat("  Pathway group summary:\n")
for (i in seq_len(nrow(pathway_summary))) {
  cat(sprintf("    %s: %d/%d hyper (%.1f%%)\n",
              gsub("\n", " ", pathway_summary$pathway_group[i]),
              pathway_summary$n_hyper[i], pathway_summary$n_total[i], pathway_summary$pct_hyper[i]))
}

# Chi-square test
chi_table <- pathway_df %>%
  group_by(pathway_group) %>%
  summarise(hyper = sum(hyper_dmr), not_hyper = sum(!hyper_dmr), .groups = "drop")
chi_mat <- as.matrix(chi_table[, c("hyper", "not_hyper")])
rownames(chi_mat) <- chi_table$pathway_group
chi_result <- chisq.test(chi_mat)

p_41d <- ggplot(pathway_summary, aes(x = pathway_group, y = pct_hyper, fill = pathway_group)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_hyper, n_hyper, n_total)),
            vjust = -0.3, size = 3.5, lineheight = 0.85) +
  scale_fill_manual(values = c("K119ub Only" = "#756BB1", "H3K36 Only" = "#D95F02",
                                "Convergent\n(K119ub + H3K36)" = "#984EA3", "Neither" = "grey70")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Hypermethylation Rate by DNMT Recruitment Pathway",
    subtitle = sprintf("Chi-square p = %.2e | Convergent = K119ub gained AND me2/me3 gained",
                       chi_result$p.value),
    x = "Pathway Group", y = "% Hypermethylated"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_41d, file.path(OUTPUT_DIR, SECTION_DIR, "41d_pathway_hypermethylation_rate"),
                        width = 10, height = 7)

# =============================================================================
# FIGURE 41e: 3x2 Decision Matrix Heatmap
# =============================================================================

cat("\n--- FIGURE 41e: me3 Status x K119ub Status Decision Matrix ---\n\n")

decision_df <- full_profile %>%
  dplyr::filter(hyper_dmr & !is.na(me3_fold) & !is.na(k119ub_fold)) %>%
  mutate(
    me3_status = case_when(
      me3_gained ~ "me3 Gained",
      me3_lost   ~ "me3 Lost",
      TRUE       ~ "me3 Unchanged"
    ),
    k119ub_status = ifelse(k119ub_gained, "K119ub Gained", "K119ub Not Gained"),
    me3_status = factor(me3_status, levels = c("me3 Gained", "me3 Unchanged", "me3 Lost")),
    k119ub_status = factor(k119ub_status, levels = c("K119ub Gained", "K119ub Not Gained"))
  )

decision_table <- table(decision_df$me3_status, decision_df$k119ub_status)
total_n <- sum(decision_table)
row_totals <- rowSums(decision_table)
col_totals <- colSums(decision_table)
expected <- outer(row_totals, col_totals) / total_n
enrichment <- decision_table / expected

# Build heatmap data
dm_df <- as.data.frame(as.table(decision_table))
colnames(dm_df) <- c("me3_status", "k119ub_status", "Observed")
exp_df <- as.data.frame(as.table(round(expected, 1)))
colnames(exp_df) <- c("me3_status", "k119ub_status", "Expected")
enr_df <- as.data.frame(as.table(round(enrichment, 2)))
colnames(enr_df) <- c("me3_status", "k119ub_status", "Enrichment")

decision_hm <- dm_df %>%
  left_join(exp_df, by = c("me3_status", "k119ub_status")) %>%
  left_join(enr_df, by = c("me3_status", "k119ub_status")) %>%
  mutate(
    label = sprintf("n=%d\nExp=%.0f\nO/E=%.2f", Observed, Expected, Enrichment),
    mechanism = case_when(
      me3_status == "me3 Lost" & k119ub_status == "K119ub Gained" ~ "DNMT3A\n(K119ub)",
      me3_status == "me3 Gained" & k119ub_status == "K119ub Not Gained" ~ "DNMT3B\n(me3)",
      me3_status == "me3 Gained" & k119ub_status == "K119ub Gained" ~ "Both",
      TRUE ~ ""
    )
  )

p_41e <- ggplot(decision_hm, aes(x = k119ub_status, y = me3_status, fill = Enrichment)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = label), size = 3.5, lineheight = 1.1) +
  geom_text(aes(label = mechanism), y = as.numeric(decision_hm$me3_status) - 0.35,
            size = 2.5, fontface = "bold", color = "grey20") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 1, name = "O/E") +
  labs(
    title = "Decision Matrix: me3 Status x K119ub Status (Hypermethylated Genes Only)",
    subtitle = sprintf("n = %d hypermethylated genes with both marks | Mechanism annotations shown",
                       nrow(decision_df)),
    x = "H2AK119ub Status", y = "H3K36me3 Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  )

save_multiformat_ggplot(p_41e, file.path(OUTPUT_DIR, SECTION_DIR, "41e_decision_matrix_heatmap"),
                        width = 9, height = 7)

# =============================================================================
# FIGURE 41f: Triple Scatter — K119ub vs me3 Fold, Colored by mC
# =============================================================================

cat("\n--- FIGURE 41f: K119ub vs me3 Scatter ---\n\n")

triple_scatter <- full_profile %>%
  dplyr::filter(!is.na(k119ub_fold) & !is.na(me3_fold)) %>%
  mutate(
    mc_label = case_when(
      hyper_dmr ~ "mC Hyper",
      hypo_dmr  ~ "mC Hypo",
      TRUE      ~ "Non-DMR"
    )
  )

cor_triple <- cor.test(triple_scatter$k119ub_fold, triple_scatter$me3_fold, method = "spearman")

p_41f <- ggplot(triple_scatter, aes(x = k119ub_fold, y = me3_fold, color = mc_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.3, size = 1) +
  scale_color_manual(values = c("mC Hyper" = "#D7191C", "mC Hypo" = "#2C7BB6", "Non-DMR" = "grey70"),
                     name = "mC Status") +
  # Quadrant labels
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.3, size = 3,
           label = "K119ub\u2191 + me3\u2193\n= Pure DNMT3A", color = "grey30") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3,
           label = "K119ub\u2191 + me3\u2191\n= Both pathways", color = "grey30") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3,
           label = "K119ub\u2193 + me3\u2191\n= DNMT3B only", color = "grey30") +
  labs(
    title = "H2AK119ub vs H3K36me3 Fold Change",
    subtitle = sprintf("Spearman rho = %+.3f, p = %.2e | Colored by mC DMR status",
                       cor_triple$estimate, cor_triple$p.value),
    x = "H2AK119ub Log2FC", y = "H3K36me3 Log2FC"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_41f, file.path(OUTPUT_DIR, SECTION_DIR, "41f_k119ub_vs_me3_scatter"),
                        width = 10, height = 9)

# =============================================================================
# FIGURE 41g: Pathway Attribution Pie Chart
# =============================================================================

cat("\n--- FIGURE 41g: Pathway Attribution ---\n\n")

attr_df <- full_profile %>%
  dplyr::filter(hyper_dmr) %>%
  mutate(
    pathway = case_when(
      k119ub_gained & me3_gained  ~ "Convergent\n(K119ub + me3)",
      k119ub_gained & me2_gained  ~ "Convergent\n(K119ub + me2)",
      k119ub_gained               ~ "DNMT3A via\nH2AK119ub",
      me3_gained                  ~ "DNMT3B via\nH3K36me3",
      me2_gained                  ~ "DNMT3A via\nH3K36me2",
      TRUE                        ~ "Unknown\nMechanism"
    )
  )

pathway_counts <- attr_df %>%
  count(pathway) %>%
  arrange(desc(n)) %>%
  mutate(pct = 100 * n / sum(n),
         label = sprintf("%s\n%d (%.1f%%)", pathway, n, pct))

cat("  Pathway attribution for hypermethylated genes:\n")
for (i in seq_len(nrow(pathway_counts))) {
  cat(sprintf("    %s: %d (%.1f%%)\n", gsub("\n", " ", pathway_counts$pathway[i]),
              pathway_counts$n[i], pathway_counts$pct[i]))
}

pathway_colors <- c(
  "DNMT3A via\nH2AK119ub" = "#756BB1",
  "DNMT3B via\nH3K36me3"  = "#D95F02",
  "DNMT3A via\nH3K36me2"  = "#E6AB02",
  "Convergent\n(K119ub + me3)" = "#984EA3",
  "Convergent\n(K119ub + me2)" = "#BC80BD",
  "Unknown\nMechanism"    = "#999999"
)

p_41g <- ggplot(pathway_counts, aes(x = "", y = n, fill = pathway)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_manual(values = pathway_colors, name = "Pathway") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 2.8, lineheight = 0.85) +
  labs(
    title = "Pathway Attribution: Hypermethylated Gene Bodies",
    subtitle = sprintf("n = %d hypermethylated genes classified by DNMT recruitment evidence", nrow(attr_df))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.position = "right"
  )

save_multiformat_ggplot(p_41g, file.path(OUTPUT_DIR, SECTION_DIR, "41g_pathway_attribution_pie"),
                        width = 10, height = 8)

# =============================================================================
# FIGURE 41h: Top-50 Hypermethylated Genes Summary Heatmap
# =============================================================================

cat("\n--- FIGURE 41h: Top-50 Summary Heatmap ---\n\n")

top50 <- full_profile %>%
  dplyr::filter(hyper_dmr) %>%
  arrange(desc(mc_diff)) %>%
  slice_head(n = 50)

# Build matrix for heatmap
hm_cols <- c("mc_diff", "hmc_diff")
hm_labels <- c("5mC diff", "5hmC diff")

# Add available fold-change columns
optional_cols <- c("k119ub_fold", "me3_fold", "me2_fold", "k27me3_fold", "k27ac_fold", "atac_fold")
optional_labels <- c("H2AK119ub", "H3K36me3", "H3K36me2", "H3K27me3", "H3K27ac", "ATAC-seq")

for (i in seq_along(optional_cols)) {
  if (optional_cols[i] %in% colnames(top50)) {
    hm_cols <- c(hm_cols, optional_cols[i])
    hm_labels <- c(hm_labels, optional_labels[i])
  }
}

hm_mat <- as.matrix(top50[, hm_cols])
rownames(hm_mat) <- top50$gene
colnames(hm_mat) <- hm_labels

# Row annotation: pathway
row_anno_df <- data.frame(
  Pathway = attr_df$pathway[match(top50$gene, attr_df$gene)],
  row.names = top50$gene
)
row_anno_df$Pathway[is.na(row_anno_df$Pathway)] <- "Unknown\nMechanism"

# Scale columns independently for visualization
hm_mat_scaled <- apply(hm_mat, 2, function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
})
rownames(hm_mat_scaled) <- top50$gene

hm_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

anno_colors <- list(Pathway = pathway_colors[unique(row_anno_df$Pathway)])

hm_call <- bquote(pheatmap(
  .(hm_mat_scaled),
  color = .(hm_colors),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_method = "complete",
  fontsize_row = 7,
  fontsize_col = 10,
  fontsize = 10,
  main = "Top 50 Hypermethylated Genes: Multi-Mark Profile (Z-scored)",
  annotation_row = .(row_anno_df),
  annotation_colors = .(anno_colors),
  border_color = "grey90",
  na_col = "grey95"
))

save_multiformat_pheatmap(
  hm_call,
  file.path(OUTPUT_DIR, SECTION_DIR, "41h_top50_summary_heatmap"),
  width = 10, height = 14
)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting Tables ---\n\n")

# Pathway attribution per gene
attr_export <- attr_df %>%
  dplyr::select(gene, mc_diff, hmc_diff, k119ub_gained, me3_gained, me3_lost,
                me2_gained, me2_lost, coordinated, pathway,
                any_of(c("me3_fold", "me2_fold", "k119ub_fold")))

write.table(attr_export, file.path(TABLES_DIR, "dnmt3a_vs_dnmt3b_pathway_attribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_vs_dnmt3b_pathway_attribution.tsv (%d rows)\n", nrow(attr_export)))

cat("\n================================================================================\n")
cat("SECTION 41 COMPLETE\n")
cat("================================================================================\n")
