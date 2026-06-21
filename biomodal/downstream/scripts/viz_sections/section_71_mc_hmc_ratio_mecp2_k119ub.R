# biomodal/downstream/scripts/viz_sections/section_71_mc_hmc_ratio_mecp2_k119ub.R
# Section 71: 5mC/5hmC Change Ratio vs MeCP2 Binding at H2AK119ub Loci
# Standalone script - sources shared config for all dependencies and data
#
# Computes the ratio of 5mC gain to 5hmC loss per gene (mc_diff / -hmc_diff),
# then asks: does MeCP2 increasingly bind where the ratio is more 5mC-skewed,
# or is it recruited equally to all genes with ratio > 1?
#
#   Panel 71a: Scatter (x = K119ub log2FC, y = MeCP2 fold, color = log2 ratio)
#   Panel 71b: MeCP2 fold by ratio quintile (dose-response: gradient vs threshold)
#   Panel 71c: Direct ratio vs MeCP2 scatter with OLS + loess
#   Panel 71d: Nested model comparison (K119ub alone vs ratio alone vs both vs interaction)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_71_mc_hmc_ratio_mecp2_k119ub.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 71: 5mC/5hmC RATIO vs MeCP2 BINDING AT H2AK119ub LOCI
# =============================================================================

cat("================================================================================\n")
cat("SECTION 71: 5mC/5hmC RATIO vs MeCP2 BINDING AT H2AK119ub LOCI\n")
cat("================================================================================\n\n")

# ---- Load data ---------------------------------------------------------------

DIFFBIND_GENE_PATH <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
K119UB_SIGNAL_PATH <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")

stopifnot(
  "diffbind_gene_level_all_marks.tsv not found (run section_33 first)" =
    file.exists(DIFFBIND_GENE_PATH),
  "k119ub_gene_signal.tsv not found" =
    file.exists(K119UB_SIGNAL_PATH),
  "MeCP2 annotated file not found" =
    file.exists(MECP2_FILES$annotated)
)

cat("Loading gene-level data...\n")
diffbind_genes <- read.table(DIFFBIND_GENE_PATH, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")
cat(sprintf("  DiffBind gene table: %d genes\n", nrow(diffbind_genes)))

cat("Loading K119ub gene body signal...\n")
k119ub_signal <- read.table(K119UB_SIGNAL_PATH, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "")
cat(sprintf("  K119ub gene signal: %d genes\n", nrow(k119ub_signal)))

cat("Loading and aggregating MeCP2 from full annotated file...\n")
mecp2_raw <- read.table(MECP2_FILES$annotated, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, fill = TRUE, quote = "")
mecp2_raw$Fold <- as.numeric(mecp2_raw$Fold)

mecp2_gene <- mecp2_raw %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(
    mecp2_mean_fold = mean(Fold, na.rm = TRUE),
    .groups = "drop"
  )
cat(sprintf("  MeCP2: %d unique genes with peaks\n", nrow(mecp2_gene)))

# ---- Build master table ------------------------------------------------------

cat("\nBuilding master gene-level table...\n")

master <- diffbind_genes %>%
  dplyr::select(gene, chr, start, end, mc_diff, hmc_diff) %>%
  left_join(mecp2_gene, by = c("gene" = "SYMBOL")) %>%
  left_join(
    k119ub_signal %>% dplyr::select(symbol, gb_log2fc),
    by = c("gene" = "symbol")
  ) %>%
  dplyr::filter(
    is.finite(mc_diff) & is.finite(hmc_diff) &
    is.finite(mecp2_mean_fold) & is.finite(gb_log2fc)
  )

cat(sprintf("  Master table: %d genes with mc_diff, hmc_diff, MeCP2, and K119ub data\n",
            nrow(master)))

# ---- Compute 5mC/5hmC change ratio -------------------------------------------

cat("\nComputing 5mC/5hmC change ratio...\n")

MIN_HMC_DIFF <- 0.001
LOG2_RATIO_CAP <- 3

master <- master %>%
  dplyr::filter(abs(hmc_diff) > MIN_HMC_DIFF) %>%
  dplyr::mutate(
    meth_ratio = mc_diff / (-hmc_diff),
    log2_ratio = log2(meth_ratio),
    log2_ratio_capped = pmax(pmin(log2_ratio, LOG2_RATIO_CAP), -LOG2_RATIO_CAP)
  ) %>%
  dplyr::filter(is.finite(log2_ratio))

cat(sprintf("  After filtering (|hmc_diff| > %.3f, finite ratio): %d genes\n",
            MIN_HMC_DIFF, nrow(master)))

n_positive <- sum(master$meth_ratio > 0)
n_above_1 <- sum(master$meth_ratio > 1)
n_below_1 <- sum(master$meth_ratio > 0 & master$meth_ratio < 1)
n_negative <- sum(master$meth_ratio < 0)

cat(sprintf("  Positive ratio (canonical mC up/hmC down): %d (%.1f%%)\n",
            n_positive, 100 * n_positive / nrow(master)))
cat(sprintf("    Ratio > 1 (mC gain > hmC loss): %d\n", n_above_1))
cat(sprintf("    Ratio < 1 (hmC loss > mC gain): %d\n", n_below_1))
cat(sprintf("  Negative ratio (non-canonical): %d (%.1f%%)\n",
            n_negative, 100 * n_negative / nrow(master)))
cat(sprintf("  Log2(ratio) range: [%.2f, %.2f], median = %.3f\n",
            min(master$log2_ratio), max(master$log2_ratio),
            median(master$log2_ratio)))

# ---- Statistics: ratio vs MeCP2 correlation ----------------------------------

cat("\nCorrelation tests...\n")

rho_ratio_mecp2 <- cor.test(master$log2_ratio, master$mecp2_mean_fold,
                            method = "spearman")
cat(sprintf("  Spearman (log2_ratio vs MeCP2 fold): rho = %.3f, p = %.2e\n",
            rho_ratio_mecp2$estimate, rho_ratio_mecp2$p.value))

rho_k119ub_mecp2 <- cor.test(master$gb_log2fc, master$mecp2_mean_fold,
                             method = "spearman")
cat(sprintf("  Spearman (K119ub vs MeCP2 fold):     rho = %.3f, p = %.2e\n",
            rho_k119ub_mecp2$estimate, rho_k119ub_mecp2$p.value))

rho_k119ub_ratio <- cor.test(master$gb_log2fc, master$log2_ratio,
                             method = "spearman")
cat(sprintf("  Spearman (K119ub vs log2_ratio):     rho = %.3f, p = %.2e\n",
            rho_k119ub_ratio$estimate, rho_k119ub_ratio$p.value))

# ---- Panel 71a: Main scatter (K119ub x, MeCP2 y, color = ratio) -------------

cat("\nCreating Figure 71a: K119ub vs MeCP2 colored by 5mC/5hmC ratio...\n")

SECTION_DIR <- file.path(OUTPUT_DIR, "71_mc_hmc_ratio_mecp2_k119ub")

x_lim <- quantile(master$gb_log2fc, c(0.005, 0.995), na.rm = TRUE)
y_lim <- quantile(master$mecp2_mean_fold, c(0.005, 0.995), na.rm = TRUE)

rho_label <- ifelse(rho_ratio_mecp2$p.value < 2.2e-16,
                    sprintf("Ratio~MeCP2: rho = %.3f, p < 2.2e-16",
                            rho_ratio_mecp2$estimate),
                    sprintf("Ratio~MeCP2: rho = %.3f, p = %.1e",
                            rho_ratio_mecp2$estimate, rho_ratio_mecp2$p.value))

label_df <- master %>%
  dplyr::filter(gene %in% KEY_GENES) %>%
  dplyr::filter(is.finite(gb_log2fc) & is.finite(mecp2_mean_fold))

p_71a <- ggplot(master, aes(x = gb_log2fc, y = mecp2_mean_fold,
                            color = log2_ratio_capped)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(alpha = 0.35, size = 0.7) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-LOG2_RATIO_CAP, LOG2_RATIO_CAP),
    name = "log2(5mC/5hmC\nchange ratio)",
    guide = guide_colorbar(barwidth = 1.2, barheight = 8)
  ) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  annotate("text", x = x_lim[1] + diff(x_lim) * 0.02,
           y = y_lim[2] - diff(y_lim) * 0.02,
           label = rho_label, hjust = 0, vjust = 1, size = 3.5) +
  annotate("text", x = x_lim[1] + diff(x_lim) * 0.02,
           y = y_lim[2] - diff(y_lim) * 0.08,
           label = sprintf("n = %s genes", format(nrow(master), big.mark = ",")),
           hjust = 0, vjust = 1, size = 3.5) +
  labs(
    title = "MeCP2 Binding vs H2AK119ub, Colored by 5mC/5hmC Ratio",
    subtitle = "Red = mC gain dominates | Blue = hmC loss dominates | White = stoichiometric",
    x = "H2AK119ub Gene Body log2FC (Mut/Ctrl)",
    y = "MeCP2 Mean Fold (log2)"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

if (nrow(label_df) > 0) {
  p_71a <- p_71a + geom_text_repel(
    data = label_df,
    aes(x = gb_log2fc, y = mecp2_mean_fold, label = gene),
    inherit.aes = FALSE, size = 2.5, max.overlaps = 15,
    segment.color = "grey50", segment.size = 0.3,
    fontface = "italic", color = "black"
  )
}

save_multiformat_ggplot(p_71a,
                        file.path(SECTION_DIR, "71a_k119ub_mecp2_ratio_scatter"),
                        width = 10, height = 8)

# ---- Panel 71b: MeCP2 dose-response by ratio quintile -----------------------

cat("Creating Figure 71b: MeCP2 fold by ratio quintile...\n")

master$ratio_quintile <- cut(master$log2_ratio,
                             breaks = quantile(master$log2_ratio, probs = seq(0, 1, 0.2)),
                             include.lowest = TRUE, labels = FALSE)

quintile_summary <- master %>%
  group_by(ratio_quintile) %>%
  summarise(
    n = n(),
    median_log2_ratio = median(log2_ratio),
    median_mecp2 = median(mecp2_mean_fold),
    mean_mecp2 = mean(mecp2_mean_fold),
    q25_mecp2 = quantile(mecp2_mean_fold, 0.25),
    q75_mecp2 = quantile(mecp2_mean_fold, 0.75),
    .groups = "drop"
  )

quintile_labels <- master %>%
  group_by(ratio_quintile) %>%
  summarise(
    range_label = sprintf("[%.1f, %.1f]",
                          min(log2_ratio), max(log2_ratio)),
    n = n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(x_label = sprintf("Q%d\n%s\nn=%d", ratio_quintile, range_label, n))

master <- master %>%
  left_join(quintile_labels %>% dplyr::select(ratio_quintile, x_label),
            by = "ratio_quintile")

kw_test <- kruskal.test(mecp2_mean_fold ~ ratio_quintile, data = master)
cat(sprintf("  Kruskal-Wallis across quintiles: H = %.1f, df = %d, p = %.2e\n",
            kw_test$statistic, kw_test$parameter, kw_test$p.value))

pairwise_results <- data.frame(
  comparison = character(), p = numeric(), stringsAsFactors = FALSE
)
for (q in 1:4) {
  v1 <- master$mecp2_mean_fold[master$ratio_quintile == q]
  v2 <- master$mecp2_mean_fold[master$ratio_quintile == q + 1]
  wt <- wilcox.test(v1, v2)
  pairwise_results <- rbind(pairwise_results, data.frame(
    comparison = sprintf("Q%d vs Q%d", q, q + 1),
    p = wt$p.value, stringsAsFactors = FALSE
  ))
  cat(sprintf("  Q%d vs Q%d: Wilcoxon p = %.2e\n", q, q + 1, wt$p.value))
}

quintile_fill_vals <- quintile_summary$median_log2_ratio

p_71b <- ggplot(master, aes(x = factor(ratio_quintile), y = mecp2_mean_fold)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(aes(fill = factor(ratio_quintile)), alpha = 0.6,
              scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(
    values = setNames(
      colorRampPalette(c("#2166AC", "grey90", "#B2182B"))(5),
      as.character(1:5)
    ),
    guide = "none"
  ) +
  scale_x_discrete(labels = quintile_labels$x_label) +
  coord_cartesian(
    ylim = quantile(master$mecp2_mean_fold, c(0.005, 0.995), na.rm = TRUE)
  ) +
  annotate("text", x = 3, y = Inf,
           label = sprintf("Kruskal-Wallis p = %.1e", kw_test$p.value),
           vjust = 1.5, size = 3.5, fontface = "italic") +
  labs(
    title = "MeCP2 Fold by 5mC/5hmC Ratio Quintile",
    subtitle = "Gradient: MeCP2 increases across quintiles | Threshold: flat above Q3",
    x = "log2(5mC/5hmC Ratio) Quintile", y = "MeCP2 Mean Fold (log2)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 8))

save_multiformat_ggplot(p_71b,
                        file.path(SECTION_DIR, "71b_mecp2_by_ratio_quintile"),
                        width = 9, height = 7)

# ---- Panel 71c: Direct ratio vs MeCP2 scatter with OLS ----------------------

cat("Creating Figure 71c: direct ratio vs MeCP2 scatter...\n")

ols_model <- lm(mecp2_mean_fold ~ log2_ratio, data = master)
ols_summary <- summary(ols_model)
ols_slope <- coef(ols_model)["log2_ratio"]
ols_r2 <- ols_summary$r.squared
ols_ci <- confint(ols_model, "log2_ratio", level = 0.95)

cat(sprintf("  OLS: MeCP2 ~ log2_ratio\n"))
cat(sprintf("    slope = %.4f [%.4f, %.4f], R2 = %.4f\n",
            ols_slope, ols_ci[1], ols_ci[2], ols_r2))

rho_direct_label <- ifelse(
  rho_ratio_mecp2$p.value < 2.2e-16,
  sprintf("rho = %.3f, p < 2.2e-16", rho_ratio_mecp2$estimate),
  sprintf("rho = %.3f, p = %.1e", rho_ratio_mecp2$estimate, rho_ratio_mecp2$p.value)
)

p_71c <- ggplot(master, aes(x = log2_ratio, y = mecp2_mean_fold)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(aes(color = log2_ratio_capped), alpha = 0.25, size = 0.6) +
  geom_smooth(method = "loess", color = "black", linewidth = 0.8, se = TRUE,
              alpha = 0.2) +
  geom_smooth(method = "lm", color = "#E41A1C", linewidth = 0.6, se = FALSE,
              linetype = "dashed") +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-LOG2_RATIO_CAP, LOG2_RATIO_CAP),
    guide = "none"
  ) +
  coord_cartesian(
    xlim = c(-LOG2_RATIO_CAP, LOG2_RATIO_CAP),
    ylim = quantile(master$mecp2_mean_fold, c(0.005, 0.995), na.rm = TRUE)
  ) +
  annotate("text", x = -LOG2_RATIO_CAP + 0.1, y = Inf,
           label = sprintf("Spearman %s\nOLS slope = %.4f, R2 = %.4f",
                           rho_direct_label, ols_slope, ols_r2),
           hjust = 0, vjust = 1.3, size = 3.5) +
  labs(
    title = "5mC/5hmC Ratio vs MeCP2 Fold (Direct)",
    subtitle = "Black = loess | Red dashed = OLS | Is there a continuous gradient?",
    x = "log2(5mC/5hmC Change Ratio)",
    y = "MeCP2 Mean Fold (log2)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_71c,
                        file.path(SECTION_DIR, "71c_ratio_vs_mecp2_direct"),
                        width = 9, height = 7)

# ---- Panel 71d: Nested model comparison (secondary validation) ---------------

cat("\nCreating Figure 71d: nested model comparison...\n")

model_df <- master %>%
  dplyr::select(mecp2_mean_fold, gb_log2fc, log2_ratio) %>%
  dplyr::filter(complete.cases(.))

cat(sprintf("  Complete cases for models: %d genes\n", nrow(model_df)))

m1 <- lm(mecp2_mean_fold ~ gb_log2fc, data = model_df)
m2 <- lm(mecp2_mean_fold ~ log2_ratio, data = model_df)
m3 <- lm(mecp2_mean_fold ~ gb_log2fc + log2_ratio, data = model_df)
m4 <- lm(mecp2_mean_fold ~ gb_log2fc * log2_ratio, data = model_df)

collect_model_stats <- function(model, name) {
  s <- summary(model)
  f_p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  data.frame(
    model = name,
    n = nrow(model$model),
    r_squared = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    aic = AIC(model),
    f_stat = s$fstatistic[1],
    f_p = f_p,
    stringsAsFactors = FALSE
  )
}

model_comparison <- rbind(
  collect_model_stats(m1, "K119ub only"),
  collect_model_stats(m2, "Ratio only"),
  collect_model_stats(m3, "K119ub + Ratio"),
  collect_model_stats(m4, "K119ub * Ratio")
)

cat("\n  Model comparison:\n")
for (i in seq_len(nrow(model_comparison))) {
  row <- model_comparison[i, ]
  cat(sprintf("    %-20s R2=%.4f  adj_R2=%.4f  AIC=%.0f  p=%.2e\n",
              row$model, row$r_squared, row$adj_r_squared,
              row$aic, row$f_p))
}

# Variance partitioning
r2_k119ub <- summary(m1)$r.squared
r2_ratio <- summary(m2)$r.squared
r2_full <- summary(m3)$r.squared

vp <- data.frame(
  component = c("K119ub unique", "Shared", "Ratio unique", "Unexplained"),
  fraction = c(
    r2_full - r2_ratio,
    r2_k119ub + r2_ratio - r2_full,
    r2_full - r2_k119ub,
    1 - r2_full
  ),
  stringsAsFactors = FALSE
)
vp$pct <- sprintf("%.1f%%", 100 * vp$fraction)

cat("\n  Variance partition:\n")
for (i in seq_len(nrow(vp))) {
  cat(sprintf("    %-18s %.4f (%s)\n", vp$component[i], vp$fraction[i], vp$pct[i]))
}

# Interaction term from Model 4
m4_coefs <- coef(summary(m4))
interaction_row <- m4_coefs["gb_log2fc:log2_ratio", , drop = FALSE]
interaction_est <- interaction_row[1, "Estimate"]
interaction_p <- interaction_row[1, "Pr(>|t|)"]

cat(sprintf("\n  Interaction test (K119ub x ratio):\n"))
cat(sprintf("    estimate = %+.4f, p = %.2e\n", interaction_est, interaction_p))
cat(sprintf("    %s\n", ifelse(interaction_p < 0.05,
            "=> K119ub and ratio interact (synergy)",
            "=> No significant interaction (additive effects)")))

# Standardized betas for additive model (Model 3)
model_df_z <- model_df
model_df_z$gb_log2fc <- as.numeric(scale(model_df_z$gb_log2fc))
model_df_z$log2_ratio <- as.numeric(scale(model_df_z$log2_ratio))
m3_z <- lm(mecp2_mean_fold ~ gb_log2fc + log2_ratio, data = model_df_z)
m3_z_coefs <- coef(summary(m3_z))
m3_z_ci <- confint(m3_z)

std_betas <- data.frame(
  predictor = c("K119ub log2FC", "log2(mC/hmC ratio)"),
  beta = m3_z_coefs[c("gb_log2fc", "log2_ratio"), "Estimate"],
  se = m3_z_coefs[c("gb_log2fc", "log2_ratio"), "Std. Error"],
  ci_lo = m3_z_ci[c("gb_log2fc", "log2_ratio"), 1],
  ci_hi = m3_z_ci[c("gb_log2fc", "log2_ratio"), 2],
  p = m3_z_coefs[c("gb_log2fc", "log2_ratio"), "Pr(>|t|)"],
  stringsAsFactors = FALSE
)

cat("\n  Standardized betas (additive model):\n")
for (i in seq_len(nrow(std_betas))) {
  cat(sprintf("    %-25s beta=%+.4f [%+.4f, %+.4f]  p=%.2e\n",
              std_betas$predictor[i], std_betas$beta[i],
              std_betas$ci_lo[i], std_betas$ci_hi[i], std_betas$p[i]))
}

# Visualization: R2 bar chart
model_comparison$model <- factor(model_comparison$model,
                                 levels = model_comparison$model)

p_71d <- ggplot(model_comparison, aes(x = model, y = r_squared)) +
  geom_bar(stat = "identity", width = 0.6,
           fill = c("#377EB8", "#E41A1C", "#984EA3", "#FF7F00"),
           alpha = 0.8) +
  geom_text(aes(label = sprintf("R2=%.4f\nAIC=%.0f", r_squared, aic)),
            vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Nested Model Comparison: Predicting MeCP2 Fold",
    subtitle = sprintf("Interaction (K119ub x Ratio): est=%+.4f, p=%.2e | n=%s",
                        interaction_est, interaction_p,
                        format(nrow(model_df), big.mark = ",")),
    x = "", y = expression(R^2)
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 10))

save_multiformat_ggplot(p_71d,
                        file.path(SECTION_DIR, "71d_nested_model_comparison"),
                        width = 9, height = 7)

# ---- Combined figure ---------------------------------------------------------

cat("\nCreating combined Figure 71...\n")

p_71_combined <- (p_71a | p_71b) / (p_71c | p_71d) +
  plot_annotation(
    title = "5mC/5hmC Ratio vs MeCP2 at H2AK119ub Loci",
    subtitle = "Does MeCP2 binding scale with the degree of 5mC-skew (gradient), or respond equally above ratio=1 (threshold)?",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_71_combined,
                        file.path(SECTION_DIR, "71_mc_hmc_ratio_mecp2_k119ub"),
                        width = 18, height = 14)

# ---- Save tables -------------------------------------------------------------

cat("\nSaving tables...\n")

# Per-gene table
gene_table <- master %>%
  dplyr::select(gene, chr, start, end, mc_diff, hmc_diff,
                meth_ratio, log2_ratio, gb_log2fc, mecp2_mean_fold,
                ratio_quintile) %>%
  arrange(desc(log2_ratio))

write.table(gene_table,
            file.path(TABLES_DIR, "71_ratio_mecp2_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 71_ratio_mecp2_summary.tsv\n")

# Quintile dose-response
quintile_table <- quintile_summary %>%
  left_join(pairwise_results %>%
              dplyr::mutate(ratio_quintile = as.integer(gsub("Q(\\d) vs Q\\d", "\\1",
                                                             comparison))),
            by = "ratio_quintile") %>%
  dplyr::mutate(wilcox_p_vs_next = p) %>%
  dplyr::select(-comparison, -p)

write.table(quintile_table,
            file.path(TABLES_DIR, "71_quintile_dose_response.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 71_quintile_dose_response.tsv\n")

# Model comparison
write.table(model_comparison,
            file.path(TABLES_DIR, "71_model_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 71_model_comparison.tsv\n")

# Variance partition
write.table(vp,
            file.path(TABLES_DIR, "71_variance_partition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 71_variance_partition.tsv\n")

# Statistics summary
stats_table <- data.frame(
  test = c("Spearman (log2_ratio vs MeCP2)",
           "Spearman (K119ub vs MeCP2)",
           "Spearman (K119ub vs ratio)",
           "Kruskal-Wallis (quintiles)",
           "OLS (MeCP2 ~ ratio)",
           "Interaction (K119ub x ratio)"),
  statistic = c(sprintf("rho = %.3f", rho_ratio_mecp2$estimate),
                sprintf("rho = %.3f", rho_k119ub_mecp2$estimate),
                sprintf("rho = %.3f", rho_k119ub_ratio$estimate),
                sprintf("H = %.1f", kw_test$statistic),
                sprintf("slope = %.4f, R2 = %.4f", ols_slope, ols_r2),
                sprintf("est = %+.4f", interaction_est)),
  p_value = c(rho_ratio_mecp2$p.value,
              rho_k119ub_mecp2$p.value,
              rho_k119ub_ratio$p.value,
              kw_test$p.value,
              ols_summary$coefficients["log2_ratio", "Pr(>|t|)"],
              interaction_p),
  n = nrow(master),
  stringsAsFactors = FALSE
)

write.table(stats_table,
            file.path(TABLES_DIR, "71_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 71_statistics.tsv\n")

# Standardized betas
write.table(std_betas,
            file.path(TABLES_DIR, "71_standardized_betas.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 71_standardized_betas.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 71 COMPLETE\n")
cat("================================================================================\n")
