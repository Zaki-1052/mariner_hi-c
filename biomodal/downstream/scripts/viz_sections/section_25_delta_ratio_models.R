# biomodal/downstream/scripts/viz_sections/section_25_delta_ratio_models.R
# Section 25: Delta-Ratio Linear Model Refits (TODO 14d)
# Standalone script - sources shared config for all dependencies and data
#
# Refits the logistic regression models from Sections 23 and 24 using
# delta_ratio (continuous 5hmC/(5mC+5hmC) change) as response instead of
# binary DMR status. This tests whether the same features that predict
# *which* genes become DMRs also predict *how much* TET activity changes.
#
# Sign convention: delta_ratio is more negative = more TET impairment.
# Features that positively predicted binary DMR status should have
# negative linear coefficients (same biology, flipped sign).
#
# Data sources (all pre-existing TSVs, joined by gene):
#   1. demethylation_ratio_all_genes.tsv (Section 22) — delta_ratio, ratio_ctrl, ratio_mut
#   2. baseline_hmc_predictor_all_genes.tsv (Section 23) — wt_hmc, gb_ctrl_signal
#   3. dnmt3a_feature_matrix.tsv (Section 24) — k119ub, atac_count, cpg_density, etc.
#
# Figures:
#   25a: Dose-response scatter — wt_hmc (x) vs delta_ratio (y)
#   25b: Binary vs continuous model comparison (Section 23) — AUC vs R²
#   25c: Feature importance comparison (Section 24) — logistic |beta| vs linear |beta|
#   25d: Binary vs continuous model comparison (Section 24) — AUC vs R²
#   25e: Residual diagnostics — 4-panel
#   25f: Predicted vs observed delta_ratio scatter (Full model)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_25_delta_ratio_models.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 25 CONFIGURATION
# =============================================================================

# Input tables (exported by Sections 22, 23, 24)
RATIO_TABLE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
HMC_TABLE <- file.path(TABLES_DIR, "baseline_hmc_predictor_all_genes.tsv")
DNMT3A_TABLE <- file.path(TABLES_DIR, "dnmt3a_feature_matrix.tsv")

# Original logistic model comparison tables
HMC_MODEL_TABLE <- file.path(TABLES_DIR, "baseline_hmc_model_comparison.tsv")
DNMT3A_MODEL_TABLE <- file.path(TABLES_DIR, "dnmt3a_model_comparison.tsv")
DNMT3A_IMPORTANCE_TABLE <- file.path(TABLES_DIR, "dnmt3a_feature_importance.tsv")

# Helper: format p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 25: DELTA-RATIO LINEAR MODEL REFITS (TODO 14d)\n")
cat("========================================================================\n\n")

cat("Validating inputs...\n")
stopifnot("demethylation_ratio_all_genes.tsv not found" = file.exists(RATIO_TABLE))
stopifnot("baseline_hmc_predictor_all_genes.tsv not found" = file.exists(HMC_TABLE))
stopifnot("dnmt3a_feature_matrix.tsv not found" = file.exists(DNMT3A_TABLE))
stopifnot("baseline_hmc_model_comparison.tsv not found" = file.exists(HMC_MODEL_TABLE))
stopifnot("dnmt3a_model_comparison.tsv not found" = file.exists(DNMT3A_MODEL_TABLE))
stopifnot("dnmt3a_feature_importance.tsv not found" = file.exists(DNMT3A_IMPORTANCE_TABLE))
cat("  All inputs validated.\n\n")

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("--- Step 1: Loading data ---\n")

ratio_df <- read.table(RATIO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Demethylation ratio: %d genes\n", nrow(ratio_df)))

hmc_df <- read.table(HMC_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Baseline hmC predictor: %d genes\n", nrow(hmc_df)))

dnmt3a_df <- read.table(DNMT3A_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  DNMT3A feature matrix: %d genes\n", nrow(dnmt3a_df)))

# Load original logistic model summaries for side-by-side comparison
hmc_models_orig <- read.table(HMC_MODEL_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dnmt3a_models_orig <- read.table(DNMT3A_MODEL_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dnmt3a_importance_orig <- read.table(DNMT3A_IMPORTANCE_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# =============================================================================
# STEP 2: BUILD SECTION 23 REFIT DATA (delta_ratio ~ wt_hmc, gb_ctrl_signal)
# =============================================================================

cat("\n--- Step 2: Building Section 23 refit data ---\n")

# Join ratio data with Section 23 predictors
hmc_merged <- dplyr::inner_join(
  ratio_df %>% dplyr::select(gene, delta_ratio, chromatin_state),
  hmc_df %>% dplyr::select(gene, wt_hmc, gb_ctrl_signal),
  by = "gene"
) %>%
  dplyr::filter(complete.cases(.))

cat(sprintf("  Merged genes (ratio + 5hmC + K119ub): %d\n", nrow(hmc_merged)))
cat(sprintf("  Delta-ratio range: [%.4f, %.4f], median: %.4f\n",
            min(hmc_merged$delta_ratio), max(hmc_merged$delta_ratio),
            median(hmc_merged$delta_ratio)))

# =============================================================================
# STEP 3: FIT SECTION 23 LINEAR MODELS
# =============================================================================

cat("\n--- Step 3: Fitting Section 23 linear models ---\n")

# Model A refit: delta_ratio ~ wt_hmc
cat("  Model A: delta_ratio ~ wt_hmc\n")
lm_a <- lm(delta_ratio ~ wt_hmc, data = hmc_merged)
lm_a_summary <- summary(lm_a)

# Model B refit: delta_ratio ~ gb_ctrl_signal
cat("  Model B: delta_ratio ~ gb_ctrl_signal\n")
lm_b <- lm(delta_ratio ~ gb_ctrl_signal, data = hmc_merged)
lm_b_summary <- summary(lm_b)

# Model C refit: delta_ratio ~ wt_hmc + gb_ctrl_signal
cat("  Model C: delta_ratio ~ wt_hmc + gb_ctrl_signal\n")
lm_c <- lm(delta_ratio ~ wt_hmc + gb_ctrl_signal, data = hmc_merged)
lm_c_summary <- summary(lm_c)

# Standardized versions (z-scored predictors)
hmc_merged_z <- hmc_merged %>%
  dplyr::mutate(
    wt_hmc_z = as.numeric(scale(wt_hmc)),
    gb_ctrl_signal_z = as.numeric(scale(gb_ctrl_signal))
  )

lm_a_z <- lm(delta_ratio ~ wt_hmc_z, data = hmc_merged_z)
lm_b_z <- lm(delta_ratio ~ gb_ctrl_signal_z, data = hmc_merged_z)
lm_c_z <- lm(delta_ratio ~ wt_hmc_z + gb_ctrl_signal_z, data = hmc_merged_z)

# Print summaries
cat(sprintf("\n  Model A (Baseline 5hmC):\n"))
cat(sprintf("    R² = %.4f, adj R² = %.4f\n", lm_a_summary$r.squared, lm_a_summary$adj.r.squared))
cat(sprintf("    F(%d,%d) = %.1f, %s\n", lm_a_summary$fstatistic[2], lm_a_summary$fstatistic[3],
            lm_a_summary$fstatistic[1], fmt_p(pf(lm_a_summary$fstatistic[1],
            lm_a_summary$fstatistic[2], lm_a_summary$fstatistic[3], lower.tail = FALSE))))
cat(sprintf("    wt_hmc: coef = %.6f, %s\n",
            coef(lm_a)["wt_hmc"], fmt_p(lm_a_summary$coefficients["wt_hmc", 4])))

cat(sprintf("\n  Model B (K119ub signal):\n"))
cat(sprintf("    R² = %.4f, adj R² = %.4f\n", lm_b_summary$r.squared, lm_b_summary$adj.r.squared))
cat(sprintf("    F(%d,%d) = %.1f, %s\n", lm_b_summary$fstatistic[2], lm_b_summary$fstatistic[3],
            lm_b_summary$fstatistic[1], fmt_p(pf(lm_b_summary$fstatistic[1],
            lm_b_summary$fstatistic[2], lm_b_summary$fstatistic[3], lower.tail = FALSE))))
cat(sprintf("    gb_ctrl_signal: coef = %.6f, %s\n",
            coef(lm_b)["gb_ctrl_signal"], fmt_p(lm_b_summary$coefficients["gb_ctrl_signal", 4])))

cat(sprintf("\n  Model C (Combined):\n"))
cat(sprintf("    R² = %.4f, adj R² = %.4f\n", lm_c_summary$r.squared, lm_c_summary$adj.r.squared))
cat(sprintf("    F(%d,%d) = %.1f, %s\n", lm_c_summary$fstatistic[2], lm_c_summary$fstatistic[3],
            lm_c_summary$fstatistic[1], fmt_p(pf(lm_c_summary$fstatistic[1],
            lm_c_summary$fstatistic[2], lm_c_summary$fstatistic[3], lower.tail = FALSE))))

cat(sprintf("\n  Standardized coefficients (Model C):\n"))
cat(sprintf("    wt_hmc_z:         beta = %.4f\n", coef(lm_c_z)["wt_hmc_z"]))
cat(sprintf("    gb_ctrl_signal_z: beta = %.4f\n", coef(lm_c_z)["gb_ctrl_signal_z"]))

# =============================================================================
# STEP 4: BUILD SECTION 24 REFIT DATA (delta_ratio ~ 7 DNMT3A features)
# =============================================================================

cat("\n--- Step 4: Building Section 24 refit data ---\n")

predictor_cols <- c("k119ub", "atac_count", "cpg_density", "baseline_mc",
                    "baseline_hmc", "log_gene_length", "log_expression")

dnmt3a_merged <- dplyr::inner_join(
  ratio_df %>% dplyr::select(gene, delta_ratio, chromatin_state),
  dnmt3a_df %>% dplyr::select(gene, all_of(predictor_cols)),
  by = "gene"
) %>%
  dplyr::filter(complete.cases(.))

cat(sprintf("  Merged genes (ratio + DNMT3A features): %d\n", nrow(dnmt3a_merged)))

# =============================================================================
# STEP 5: FIT SECTION 24 LINEAR MODELS
# =============================================================================

cat("\n--- Step 5: Fitting Section 24 linear models ---\n")

# Full model
cat("  Full: delta_ratio ~ all 7 features\n")
lm_full <- lm(delta_ratio ~ k119ub + atac_count + cpg_density + baseline_mc +
                baseline_hmc + log_gene_length + log_expression,
              data = dnmt3a_merged)
lm_full_summary <- summary(lm_full)

# DNMT3A recruitment
cat("  DNMT3A recruitment: delta_ratio ~ k119ub + atac_count + cpg_density\n")
lm_dnmt3a <- lm(delta_ratio ~ k119ub + atac_count + cpg_density,
                data = dnmt3a_merged)
lm_dnmt3a_summary <- summary(lm_dnmt3a)

# TET impediment
cat("  TET impediment: delta_ratio ~ baseline_hmc + atac_count\n")
lm_tet <- lm(delta_ratio ~ baseline_hmc + atac_count, data = dnmt3a_merged)
lm_tet_summary <- summary(lm_tet)

# K119ub only
cat("  K119ub only: delta_ratio ~ k119ub\n")
lm_k119ub <- lm(delta_ratio ~ k119ub, data = dnmt3a_merged)
lm_k119ub_summary <- summary(lm_k119ub)

# Standardized versions
dnmt3a_merged_z <- dnmt3a_merged
for (col in predictor_cols) {
  dnmt3a_merged_z[[paste0(col, "_z")]] <- as.numeric(scale(dnmt3a_merged_z[[col]]))
}

lm_full_z <- lm(delta_ratio ~ k119ub_z + atac_count_z + cpg_density_z +
                  baseline_mc_z + baseline_hmc_z + log_gene_length_z + log_expression_z,
                data = dnmt3a_merged_z)
lm_dnmt3a_z <- lm(delta_ratio ~ k119ub_z + atac_count_z + cpg_density_z,
                   data = dnmt3a_merged_z)
lm_tet_z <- lm(delta_ratio ~ baseline_hmc_z + atac_count_z, data = dnmt3a_merged_z)
lm_k119ub_z <- lm(delta_ratio ~ k119ub_z, data = dnmt3a_merged_z)

# Collect all Section 24 linear models
s24_lm_list <- list(
  Full = lm_full, `DNMT3A recruitment` = lm_dnmt3a,
  `TET impediment` = lm_tet, `K119ub only` = lm_k119ub
)
s24_lm_z_list <- list(
  Full = lm_full_z, `DNMT3A recruitment` = lm_dnmt3a_z,
  `TET impediment` = lm_tet_z, `K119ub only` = lm_k119ub_z
)
s24_lm_summaries <- list(
  Full = lm_full_summary, `DNMT3A recruitment` = lm_dnmt3a_summary,
  `TET impediment` = lm_tet_summary, `K119ub only` = lm_k119ub_summary
)

# Print Section 24 summaries
cat("\n  Section 24 Linear Model Comparison:\n")
for (nm in names(s24_lm_list)) {
  s <- s24_lm_summaries[[nm]]
  f_p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  cat(sprintf("    %-22s R²=%.4f  adj_R²=%.4f  F=%.1f  %s\n",
              paste0(nm, ":"), s$r.squared, s$adj.r.squared,
              s$fstatistic[1], fmt_p(f_p)))
}

# Standardized coefficients (full model)
cat(sprintf("\n  Standardized coefficients (Full linear model):\n"))
full_z_coefs <- coef(lm_full_z)
full_z_se <- summary(lm_full_z)$coefficients[, 2]
full_z_p <- summary(lm_full_z)$coefficients[, 4]
for (term in names(full_z_coefs)[-1]) {
  cat(sprintf("    %-25s beta=%7.4f  %s\n",
              gsub("_z$", "", term), full_z_coefs[term], fmt_p(full_z_p[term])))
}

# =============================================================================
# FIGURE 25a: DOSE-RESPONSE SCATTER (wt_hmc vs delta_ratio)
# =============================================================================

cat("\n--- Figure 25a: Dose-response scatter ---\n")

# Spearman correlation
rho_dose <- cor.test(hmc_merged$wt_hmc, hmc_merged$delta_ratio, method = "spearman")
cat(sprintf("  Spearman rho (wt_hmc vs delta_ratio): %.3f (%s)\n",
            rho_dose$estimate, fmt_p(rho_dose$p.value)))

# Merge significance from ratio_df
hmc_plot_data <- hmc_merged %>%
  dplyr::left_join(
    ratio_df %>% dplyr::select(gene, hmc_sig),
    by = "gene"
  ) %>%
  dplyr::mutate(
    is_key_gene = gene %in% KEY_GENES,
    significant = ifelse(is.na(hmc_sig), FALSE, hmc_sig)
  )

# Top genes to label
top_label <- hmc_plot_data %>%
  dplyr::arrange(delta_ratio) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::bind_rows(
    hmc_plot_data %>%
      dplyr::arrange(dplyr::desc(delta_ratio)) %>%
      dplyr::slice_head(n = 5)
  ) %>%
  dplyr::bind_rows(
    hmc_plot_data %>% dplyr::filter(is_key_gene)
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

p_25a <- ggplot(hmc_plot_data, aes(x = wt_hmc, y = delta_ratio)) +
  geom_point(aes(color = significant), alpha = 0.12, size = 0.6) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE,
              fill = "grey80", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"),
                     labels = c("TRUE" = "hmC DMR", "FALSE" = "Not significant"),
                     name = "hmC DMR (q<0.05)") +
  geom_text_repel(data = top_label, aes(label = gene), size = 2.8,
                  max.overlaps = 25, segment.alpha = 0.5,
                  color = "black", fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.8,
           label = sprintf("Spearman rho = %.3f\n%s\nLinear R\u00B2 = %.4f\nN = %s genes",
                           rho_dose$estimate, fmt_p(rho_dose$p.value),
                           lm_a_summary$r.squared,
                           format(nrow(hmc_merged), big.mark = ","))) +
  labs(
    title = "Dose-Response: Baseline WT 5hmC vs Demethylation Efficiency Change",
    subtitle = "Higher baseline 5hmC \u2192 more negative delta_ratio (greater TET impairment)",
    x = "Baseline WT 5hmC (regional fraction)",
    y = "Delta-ratio (KO \u2212 WT demethylation efficiency)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_25a, file.path(OUTPUT_DIR, "25a_dose_response_delta_ratio"), 10, 8)

# =============================================================================
# FIGURE 25b: BINARY VS CONTINUOUS MODEL COMPARISON (SECTION 23)
# =============================================================================

cat("--- Figure 25b: Binary vs continuous model comparison (Section 23) ---\n")

# Build comparison data frame
s23_comparison <- data.frame(
  model = rep(c("A: Baseline 5hmC", "B: K119ub signal", "C: Combined"), 2),
  framework = rep(c("Logistic (AUC)", "Linear (R\u00B2)"), each = 3),
  metric = c(
    hmc_models_orig$auc[1], hmc_models_orig$auc[2], hmc_models_orig$auc[3],
    lm_a_summary$r.squared, lm_b_summary$r.squared, lm_c_summary$r.squared
  ),
  stringsAsFactors = FALSE
)
s23_comparison$model <- factor(s23_comparison$model,
                                levels = c("A: Baseline 5hmC", "B: K119ub signal", "C: Combined"))

s23_colors <- c("Logistic (AUC)" = "#377EB8", "Linear (R\u00B2)" = "#E41A1C")

p_25b <- ggplot(s23_comparison, aes(x = model, y = metric, fill = framework)) +
  geom_col(position = position_dodge(width = 0.7), alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", metric)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = s23_colors, name = "Framework") +
  labs(
    title = "Section 23 Refit: Binary vs Continuous Response",
    subtitle = "Logistic (AUC for DMR status) vs Linear (R\u00B2 for delta_ratio)",
    x = NULL,
    y = "Model Performance (AUC or R\u00B2)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom")

save_multiformat_ggplot(p_25b, file.path(OUTPUT_DIR, "25b_binary_vs_continuous_s23"), 10, 7)

# =============================================================================
# FIGURE 25c: FEATURE IMPORTANCE COMPARISON (SECTION 24)
# =============================================================================

cat("--- Figure 25c: Feature importance comparison (Section 24) ---\n")

# Display names for features
display_names <- c(
  k119ub = "H2AK119ub", atac_count = "ATAC peaks",
  cpg_density = "CpG density", baseline_mc = "Baseline 5mC",
  baseline_hmc = "Baseline 5hmC", log_gene_length = "Gene length (log10)",
  log_expression = "Expression (log10)"
)

# Get logistic betas from original importance table
logistic_betas <- dnmt3a_importance_orig %>%
  dplyr::select(feature, lr_standardized_beta, lr_rank)

# Get linear betas from current full model
linear_z_coefs <- coef(lm_full_z)[-1]
linear_betas <- data.frame(
  feature = gsub("_z$", "", names(linear_z_coefs)),
  linear_standardized_beta = as.numeric(linear_z_coefs),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    linear_abs_beta = abs(linear_standardized_beta),
    linear_rank = rank(-linear_abs_beta)
  )

# Merge
importance_cmp <- dplyr::inner_join(logistic_betas, linear_betas, by = "feature") %>%
  dplyr::mutate(display_name = display_names[feature])

# Rank correlation between logistic and linear
rank_rho <- cor(importance_cmp$lr_rank, importance_cmp$linear_rank, method = "spearman")
cat(sprintf("  Logistic vs Linear rank correlation (Spearman): %.3f\n", rank_rho))

# Scatter plot: logistic |beta| vs linear |beta|
p_25c <- ggplot(importance_cmp,
                aes(x = abs(lr_standardized_beta), y = linear_abs_beta)) +
  geom_point(size = 4, color = "#333333") +
  geom_text_repel(aes(label = display_name), size = 3.5, max.overlaps = 10) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.8,
           label = sprintf("Spearman rho = %.3f\n(rank correlation of\nfeature importance)",
                           rank_rho)) +
  labs(
    title = "Feature Importance: Logistic vs Linear Models (Section 24)",
    subtitle = "Standardized |beta| from full models; same biology, different response",
    x = "|Standardized beta| (Logistic: hyper-DMR)",
    y = "|Standardized beta| (Linear: delta_ratio)"
  ) +
  theme_biomodal() +
  coord_equal()

save_multiformat_ggplot(p_25c, file.path(OUTPUT_DIR, "25c_feature_importance_comparison"), 9, 9)

# =============================================================================
# FIGURE 25d: BINARY VS CONTINUOUS MODEL COMPARISON (SECTION 24)
# =============================================================================

cat("--- Figure 25d: Binary vs continuous model comparison (Section 24) ---\n")

# Match model names between original and current
s24_model_names <- c("Full", "DNMT3A recruitment", "TET impediment", "K119ub only")

s24_comparison <- data.frame(
  model = rep(s24_model_names, 2),
  framework = rep(c("Logistic (AUC)", "Linear (R\u00B2)"), each = 4),
  metric = c(
    dnmt3a_models_orig$auc[match(s24_model_names, dnmt3a_models_orig$model)],
    sapply(s24_lm_summaries[s24_model_names], function(s) s$r.squared)
  ),
  stringsAsFactors = FALSE
)
s24_comparison$model <- factor(s24_comparison$model, levels = s24_model_names)

s24_colors <- c("Logistic (AUC)" = "#377EB8", "Linear (R\u00B2)" = "#E41A1C")

p_25d <- ggplot(s24_comparison, aes(x = model, y = metric, fill = framework)) +
  geom_col(position = position_dodge(width = 0.7), alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", metric)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = s24_colors, name = "Framework") +
  labs(
    title = "Section 24 Refit: Binary vs Continuous Response",
    subtitle = "Logistic (AUC for hyper-DMR) vs Linear (R\u00B2 for delta_ratio)",
    x = NULL,
    y = "Model Performance (AUC or R\u00B2)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom")

save_multiformat_ggplot(p_25d, file.path(OUTPUT_DIR, "25d_binary_vs_continuous_s24"), 11, 7)

# =============================================================================
# FIGURE 25e: RESIDUAL DIAGNOSTICS (FULL LINEAR MODEL)
# =============================================================================

cat("--- Figure 25e: Residual diagnostics ---\n")

resid_data <- data.frame(
  fitted = fitted(lm_full),
  residuals = residuals(lm_full),
  std_residuals = rstandard(lm_full),
  leverage = hatvalues(lm_full),
  cooks = cooks.distance(lm_full)
)

# Panel 1: Residuals vs Fitted
p_25e_1 <- ggplot(resid_data, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.1, size = 0.5, color = "grey40") +
  geom_smooth(method = "loess", color = "#E41A1C", linewidth = 0.8, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_biomodal(base_size = 10)

# Panel 2: QQ plot
p_25e_2 <- ggplot(resid_data, aes(sample = std_residuals)) +
  geom_qq(alpha = 0.1, size = 0.5, color = "grey40") +
  geom_qq_line(color = "#E41A1C", linewidth = 0.8) +
  labs(title = "Normal Q-Q", x = "Theoretical quantiles", y = "Standardized residuals") +
  theme_biomodal(base_size = 10)

# Panel 3: Scale-Location
p_25e_3 <- ggplot(resid_data, aes(x = fitted, y = sqrt(abs(std_residuals)))) +
  geom_point(alpha = 0.1, size = 0.5, color = "grey40") +
  geom_smooth(method = "loess", color = "#E41A1C", linewidth = 0.8, se = FALSE) +
  labs(title = "Scale-Location", x = "Fitted values",
       y = expression(sqrt("|Standardized residuals|"))) +
  theme_biomodal(base_size = 10)

# Panel 4: Residuals vs Leverage
p_25e_4 <- ggplot(resid_data, aes(x = leverage, y = std_residuals)) +
  geom_point(aes(size = cooks), alpha = 0.15, color = "grey40") +
  geom_smooth(method = "loess", color = "#E41A1C", linewidth = 0.8, se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_size_continuous(range = c(0.3, 3), name = "Cook's D") +
  labs(title = "Residuals vs Leverage", x = "Leverage", y = "Standardized residuals") +
  theme_biomodal(base_size = 10) +
  theme(legend.position = "bottom")

p_25e <- (p_25e_1 | p_25e_2) / (p_25e_3 | p_25e_4) +
  plot_annotation(
    title = "Residual Diagnostics: Full Linear Model (delta_ratio ~ 7 features)",
    subtitle = sprintf("N=%s genes, R\u00B2=%.4f",
                       format(nrow(dnmt3a_merged), big.mark = ","),
                       lm_full_summary$r.squared),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

save_multiformat_ggplot(p_25e, file.path(OUTPUT_DIR, "25e_residual_diagnostics"), 12, 10)

# =============================================================================
# FIGURE 25f: PREDICTED VS OBSERVED (FULL MODEL, COLORED BY CHROMATIN STATE)
# =============================================================================

cat("--- Figure 25f: Predicted vs observed ---\n")

dnmt3a_merged$predicted_delta_ratio <- predict(lm_full, newdata = dnmt3a_merged)

# Pearson correlation
pred_cor <- cor.test(dnmt3a_merged$delta_ratio, dnmt3a_merged$predicted_delta_ratio)
cat(sprintf("  Pearson r (predicted vs observed): %.4f (%s)\n",
            pred_cor$estimate, fmt_p(pred_cor$p.value)))

p_25f <- ggplot(dnmt3a_merged, aes(x = predicted_delta_ratio, y = delta_ratio)) +
  geom_point(aes(color = chromatin_state), alpha = 0.2, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_smooth(method = "lm", color = "#E41A1C", linewidth = 0.8, se = TRUE, fill = "#E41A1C22") +
  scale_color_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.8,
           label = sprintf("Pearson r = %.3f\nR\u00B2 = %.4f\n%s\nN = %s",
                           pred_cor$estimate, lm_full_summary$r.squared,
                           fmt_p(pred_cor$p.value),
                           format(nrow(dnmt3a_merged), big.mark = ","))) +
  labs(
    title = "Predicted vs Observed Delta-Ratio (Full Model)",
    subtitle = "Linear regression with 7 features; colored by chromatin state",
    x = "Predicted delta_ratio",
    y = "Observed delta_ratio"
  ) +
  theme_biomodal() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

save_multiformat_ggplot(p_25f, file.path(OUTPUT_DIR, "25f_predicted_vs_observed"), 10, 9)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

# Table 1: Linear model comparison (Sections 23 + 24 refits)
extract_lm_stats <- function(model, model_name) {
  s <- summary(model)
  f_p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  data.frame(
    model = model_name,
    r_squared = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    f_statistic = s$fstatistic[1],
    f_df1 = s$fstatistic[2],
    f_df2 = s$fstatistic[3],
    f_pvalue = f_p,
    n_genes = nobs(model),
    residual_se = s$sigma,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

lm_comparison <- rbind(
  extract_lm_stats(lm_a, "S23_A_Baseline_5hmC"),
  extract_lm_stats(lm_b, "S23_B_K119ub_signal"),
  extract_lm_stats(lm_c, "S23_C_Combined"),
  extract_lm_stats(lm_full, "S24_Full"),
  extract_lm_stats(lm_dnmt3a, "S24_DNMT3A_recruitment"),
  extract_lm_stats(lm_tet, "S24_TET_impediment"),
  extract_lm_stats(lm_k119ub, "S24_K119ub_only")
)

write.table(lm_comparison,
            file.path(TABLES_DIR, "delta_ratio_linear_model_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved delta_ratio_linear_model_comparison.tsv (%d rows)\n", nrow(lm_comparison)))

# Table 2: Full coefficient table with standardized betas
build_lm_coef_table <- function(model, model_z, model_name) {
  s <- summary(model)$coefficients
  sz <- summary(model_z)$coefficients

  raw_terms <- rownames(s)
  result <- data.frame(
    model = model_name,
    term = raw_terms,
    estimate = s[, 1],
    std_error = s[, 2],
    t_value = s[, 3],
    p_value = s[, 4],
    standardized_beta = NA_real_,
    standardized_se = NA_real_,
    stringsAsFactors = FALSE, row.names = NULL
  )

  for (i in seq_len(nrow(result))) {
    z_name <- paste0(result$term[i], "_z")
    if (z_name %in% rownames(sz)) {
      result$standardized_beta[i] <- sz[z_name, 1]
      result$standardized_se[i] <- sz[z_name, 2]
    } else if (result$term[i] == "(Intercept)" && "(Intercept)" %in% rownames(sz)) {
      result$standardized_beta[i] <- sz["(Intercept)", 1]
      result$standardized_se[i] <- sz["(Intercept)", 2]
    }
  }

  return(result)
}

coef_table <- rbind(
  build_lm_coef_table(lm_a, lm_a_z, "S23_A_Baseline_5hmC"),
  build_lm_coef_table(lm_b, lm_b_z, "S23_B_K119ub_signal"),
  build_lm_coef_table(lm_c, lm_c_z, "S23_C_Combined"),
  build_lm_coef_table(lm_full, lm_full_z, "S24_Full"),
  build_lm_coef_table(lm_dnmt3a, lm_dnmt3a_z, "S24_DNMT3A_recruitment"),
  build_lm_coef_table(lm_tet, lm_tet_z, "S24_TET_impediment"),
  build_lm_coef_table(lm_k119ub, lm_k119ub_z, "S24_K119ub_only")
)

write.table(coef_table,
            file.path(TABLES_DIR, "delta_ratio_linear_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved delta_ratio_linear_coefficients.tsv (%d rows)\n", nrow(coef_table)))

# Table 3: Binary vs continuous side-by-side comparison
binary_vs_continuous <- rbind(
  # Section 23 models
  data.frame(
    section = "S23",
    model = c("A: Baseline 5hmC", "B: K119ub signal", "C: Combined"),
    logistic_auc = hmc_models_orig$auc,
    logistic_mcfadden_r2 = hmc_models_orig$mcfadden_r2,
    linear_r2 = c(lm_a_summary$r.squared, lm_b_summary$r.squared, lm_c_summary$r.squared),
    linear_adj_r2 = c(lm_a_summary$adj.r.squared, lm_b_summary$adj.r.squared, lm_c_summary$adj.r.squared),
    stringsAsFactors = FALSE
  ),
  # Section 24 models
  data.frame(
    section = "S24",
    model = s24_model_names,
    logistic_auc = dnmt3a_models_orig$auc[match(s24_model_names, dnmt3a_models_orig$model)],
    logistic_mcfadden_r2 = dnmt3a_models_orig$mcfadden_r2[match(s24_model_names, dnmt3a_models_orig$model)],
    linear_r2 = sapply(s24_lm_summaries[s24_model_names], function(s) s$r.squared),
    linear_adj_r2 = sapply(s24_lm_summaries[s24_model_names], function(s) s$adj.r.squared),
    stringsAsFactors = FALSE, row.names = NULL
  )
)

write.table(binary_vs_continuous,
            file.path(TABLES_DIR, "delta_ratio_binary_vs_continuous.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved delta_ratio_binary_vs_continuous.tsv (%d rows)\n", nrow(binary_vs_continuous)))

# Table 4: Feature importance comparison (logistic vs linear betas)
importance_export <- importance_cmp %>%
  dplyr::select(feature, display_name,
                lr_standardized_beta, lr_rank,
                linear_standardized_beta, linear_abs_beta, linear_rank) %>%
  dplyr::arrange(linear_rank)

write.table(importance_export,
            file.path(TABLES_DIR, "delta_ratio_feature_importance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved delta_ratio_feature_importance.tsv (%d rows)\n", nrow(importance_export)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 25 SUMMARY: DELTA-RATIO LINEAR MODEL REFITS\n")
cat("========================================================================\n\n")

cat("Section 23 Refits (delta_ratio ~ baseline 5hmC / K119ub):\n")
cat(sprintf("  Model A (Baseline 5hmC): R\u00B2=%.4f (logistic AUC=%.3f)\n",
            lm_a_summary$r.squared, hmc_models_orig$auc[1]))
cat(sprintf("  Model B (K119ub):        R\u00B2=%.4f (logistic AUC=%.3f)\n",
            lm_b_summary$r.squared, hmc_models_orig$auc[2]))
cat(sprintf("  Model C (Combined):      R\u00B2=%.4f (logistic AUC=%.3f)\n",
            lm_c_summary$r.squared, hmc_models_orig$auc[3]))

cat(sprintf("\nSection 24 Refits (delta_ratio ~ DNMT3A features):\n"))
for (nm in s24_model_names) {
  s <- s24_lm_summaries[[nm]]
  orig_auc <- dnmt3a_models_orig$auc[match(nm, dnmt3a_models_orig$model)]
  cat(sprintf("  %-22s R\u00B2=%.4f (logistic AUC=%.3f)\n",
              paste0(nm, ":"), s$r.squared, orig_auc))
}

cat(sprintf("\nFeature importance rank correlation (logistic vs linear): rho=%.3f\n",
            rank_rho))

cat(sprintf("\nSign convention check (Full model, standardized betas):\n"))
for (term in names(full_z_coefs)[-1]) {
  feat_name <- gsub("_z$", "", term)
  orig_beta <- dnmt3a_importance_orig$lr_standardized_beta[
    dnmt3a_importance_orig$feature == feat_name]
  if (length(orig_beta) > 0) {
    cat(sprintf("  %-20s logistic=%+.3f  linear=%+.3f  %s\n",
                display_names[feat_name], orig_beta, full_z_coefs[term],
                ifelse(sign(orig_beta) == -sign(full_z_coefs[term]),
                       "CONSISTENT (flipped)", "SAME SIGN")))
  }
}

cat(sprintf("\nNote: R\u00B2 values are much lower than AUC values by design.\n"))
cat("  AUC measures classification ability (0-1 scale, 0.5 = chance)\n")
cat("  R\u00B2 measures variance explained (0-1 scale, 0 = no variance explained)\n")
cat("  These are NOT directly comparable. Both can be informative.\n")

cat("\n--- Output files ---\n")
cat(sprintf("  Figures: %s/25{a-f}_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/delta_ratio_*.tsv (4 files)\n", TABLES_DIR))

cat("\n=== Section 25 complete ===\n")
