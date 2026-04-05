# biomodal/downstream/scripts/viz_sections/section_23_baseline_hmc_predictor.R
# Section 23: Baseline 5hmC as Predictor of DMR Susceptibility
# Standalone script - sources shared config for all dependencies and data
#
# Tests the "substrate availability" hypothesis: do genes with higher
# wildtype 5hmC levels have higher probability of becoming hmC DMRs?
# Compares predictive power of baseline 5hmC vs K119ub signal, and
# combined model. A stronger 5hmC effect confirms TET-mediated
# demethylation block interpretation.
#
# Data sources:
#   1. hmC feature extraction (per-sample regional fractions) — WT 5hmC levels
#   2. hmC DMR BED (from _shared_config.R) — binary DMR status + effect sizes
#   3. K119ub gene signal TSV — WT gene body K119ub signal
#
# Figures:
#   23a: ROC curves (3 logistic regression models overlaid)
#   23b: Predicted probability curve (Model A: baseline 5hmC)
#   23c: Model comparison (AIC + AUC bar/dot chart)
#   23d: Dose-response scatter (baseline 5hmC vs hmC mod_difference)
#   23e: Dose-response by chromatin state (faceted)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_23_baseline_hmc_predictor.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(pROC)
})

# =============================================================================
# SECTION 23 CONFIGURATION
# =============================================================================

# Feature extraction file path (per-sample hmC regional fractions, from _shared_config.R)
HMC_EXTRACT_FILE <- EXTRACT_PATHS$hmc_regional_frac

# K119ub gene signal file
K119UB_SIGNAL_FILE <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")

# Model colors
MODEL_COLORS <- c(
  "Model A: Baseline 5hmC"  = "#377EB8",
  "Model B: K119ub signal"  = "#756BB1",
  "Model C: Combined"       = "#E41A1C"
)

# Section output subdirectory
SECTION_DIR <- file.path(OUTPUT_DIR, "23_baseline_hmc_predictor")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

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
cat("SECTION 23: BASELINE 5hmC AS PREDICTOR OF DMR SUSCEPTIBILITY\n")
cat("========================================================================\n\n")

cat("Validating inputs...\n")
stopifnot("hmc_dmr not loaded from shared config" = exists("hmc_dmr") && !is.null(hmc_dmr))
stopifnot("hmC extract file not found" = file.exists(HMC_EXTRACT_FILE))
stopifnot("K119ub gene signal file not found" = file.exists(K119UB_SIGNAL_FILE))
cat("  All inputs validated.\n\n")

# =============================================================================
# STEP 1: LOAD & COMPUTE BASELINE WT 5hmC
# =============================================================================

cat("--- Step 1: Loading baseline WT 5hmC from feature extraction ---\n")

hmc_extract <- read.table(gzfile(HMC_EXTRACT_FILE), header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, na.strings = c("", "NA"))
cat(sprintf("  hmC extract: %d genes x %d columns\n", nrow(hmc_extract), ncol(hmc_extract)))

# Rename sample columns by position
# Column order: Chromosome, Start, End, ctrl-M, mut-M, ctrl-F, mut-F, Annotation, Name
names(hmc_extract)[4:7] <- c("hmc_ctrl_M", "hmc_mut_M", "hmc_ctrl_F", "hmc_mut_F")

hmc_baseline <- hmc_extract %>%
  dplyr::select(gene = Name, hmc_ctrl_M, hmc_ctrl_F) %>%
  dplyr::mutate(wt_hmc = (hmc_ctrl_M + hmc_ctrl_F) / 2) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::filter(complete.cases(wt_hmc))

cat(sprintf("  Baseline WT 5hmC computed for %d genes\n", nrow(hmc_baseline)))
cat(sprintf("  Median WT 5hmC: %.4f, Mean: %.4f, Range: [%.4f, %.4f]\n",
            median(hmc_baseline$wt_hmc), mean(hmc_baseline$wt_hmc),
            min(hmc_baseline$wt_hmc), max(hmc_baseline$wt_hmc)))

# =============================================================================
# STEP 2: MERGE WITH DMR STATUS
# =============================================================================

cat("\n--- Step 2: Merging with hmC DMR status ---\n")

hmc_dmr_for_join <- hmc_dmr %>%
  dplyr::select(gene, chr, start, end,
                dmr_qvalue, mod_difference, significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

merged <- dplyr::inner_join(hmc_baseline, hmc_dmr_for_join, by = "gene")
cat(sprintf("  Merged: %d genes with both WT 5hmC and DMR status\n", nrow(merged)))
cat(sprintf("  DMR+ (q<0.05): %d genes (%.1f%%)\n",
            sum(merged$significant), 100 * mean(merged$significant)))
cat(sprintf("  DMR-: %d genes\n", sum(!merged$significant)))

# =============================================================================
# STEP 3: LOAD K119ub SIGNAL
# =============================================================================

cat("\n--- Step 3: Loading K119ub gene signal ---\n")

k119ub_raw <- read.table(K119UB_SIGNAL_FILE, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
cat(sprintf("  K119ub signal: %d genes\n", nrow(k119ub_raw)))

k119ub <- k119ub_raw %>%
  dplyr::select(gene = symbol, gb_ctrl_signal) %>%
  dplyr::filter(is.finite(gb_ctrl_signal)) %>%
  dplyr::distinct(gene, .keep_all = TRUE)
cat(sprintf("  After filtering: %d genes with finite WT K119ub signal\n", nrow(k119ub)))

# Combined dataset for Model C
combined <- dplyr::inner_join(merged, k119ub, by = "gene")
cat(sprintf("  Combined (5hmC + K119ub + DMR): %d genes\n", nrow(combined)))
cat(sprintf("  DMR+ in combined set: %d genes (%.1f%%)\n",
            sum(combined$significant), 100 * mean(combined$significant)))

# =============================================================================
# STEP 4: FIT LOGISTIC REGRESSION MODELS
# =============================================================================

cat("\n--- Step 4: Fitting logistic regression models ---\n")

# Use combined dataset for all models so they're comparable on same N
model_data <- combined

# Model A: significant ~ wt_hmc
cat("  Model A: significant ~ wt_hmc\n")
model_a <- glm(significant ~ wt_hmc, data = model_data, family = binomial)
coef_a <- summary(model_a)$coefficients

# Model B: significant ~ gb_ctrl_signal
cat("  Model B: significant ~ gb_ctrl_signal\n")
model_b <- glm(significant ~ gb_ctrl_signal, data = model_data, family = binomial)
coef_b <- summary(model_b)$coefficients

# Model C: significant ~ wt_hmc + gb_ctrl_signal
cat("  Model C: significant ~ wt_hmc + gb_ctrl_signal\n")
model_c <- glm(significant ~ wt_hmc + gb_ctrl_signal, data = model_data, family = binomial)
coef_c <- summary(model_c)$coefficients

# Standardized versions (z-scored predictors) for comparable coefficients
model_data_z <- model_data %>%
  dplyr::mutate(
    wt_hmc_z = as.numeric(scale(wt_hmc)),
    gb_ctrl_signal_z = as.numeric(scale(gb_ctrl_signal))
  )

model_a_z <- glm(significant ~ wt_hmc_z, data = model_data_z, family = binomial)
model_b_z <- glm(significant ~ gb_ctrl_signal_z, data = model_data_z, family = binomial)
model_c_z <- glm(significant ~ wt_hmc_z + gb_ctrl_signal_z, data = model_data_z, family = binomial)

# McFadden pseudo-R²
null_ll <- logLik(glm(significant ~ 1, data = model_data, family = binomial))
mcfadden <- function(model) {
  1 - as.numeric(logLik(model)) / as.numeric(null_ll)
}

r2_a <- mcfadden(model_a)
r2_b <- mcfadden(model_b)
r2_c <- mcfadden(model_c)

# ROC/AUC
roc_a <- roc(model_data$significant, predict(model_a, type = "response"), quiet = TRUE)
roc_b <- roc(model_data$significant, predict(model_b, type = "response"), quiet = TRUE)
roc_c <- roc(model_data$significant, predict(model_c, type = "response"), quiet = TRUE)

ci_a <- ci.auc(roc_a, quiet = TRUE)
ci_b <- ci.auc(roc_b, quiet = TRUE)
ci_c <- ci.auc(roc_c, quiet = TRUE)

# Extract odds ratios with 95% CI
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

or_a <- extract_or(model_a)
or_b <- extract_or(model_b)
or_c <- extract_or(model_c)

# Print summaries
cat(sprintf("\n  Model A (Baseline 5hmC):\n"))
cat(sprintf("    OR = %.2f [%.2f, %.2f], %s\n",
            or_a$or[2], or_a$or_lower[2], or_a$or_upper[2], fmt_p(or_a$p_value[2])))
cat(sprintf("    AIC = %.1f, McFadden R² = %.4f\n", AIC(model_a), r2_a))
cat(sprintf("    AUC = %.3f [%.3f, %.3f]\n", auc(roc_a), ci_a[1], ci_a[3]))

cat(sprintf("\n  Model B (K119ub signal):\n"))
cat(sprintf("    OR = %.2f [%.2f, %.2f], %s\n",
            or_b$or[2], or_b$or_lower[2], or_b$or_upper[2], fmt_p(or_b$p_value[2])))
cat(sprintf("    AIC = %.1f, McFadden R² = %.4f\n", AIC(model_b), r2_b))
cat(sprintf("    AUC = %.3f [%.3f, %.3f]\n", auc(roc_b), ci_b[1], ci_b[3]))

cat(sprintf("\n  Model C (Combined):\n"))
cat(sprintf("    wt_hmc:        OR = %.2f [%.2f, %.2f], %s\n",
            or_c$or[2], or_c$or_lower[2], or_c$or_upper[2], fmt_p(or_c$p_value[2])))
cat(sprintf("    gb_ctrl_signal: OR = %.2f [%.2f, %.2f], %s\n",
            or_c$or[3], or_c$or_lower[3], or_c$or_upper[3], fmt_p(or_c$p_value[3])))
cat(sprintf("    AIC = %.1f, McFadden R² = %.4f\n", AIC(model_c), r2_c))
cat(sprintf("    AUC = %.3f [%.3f, %.3f]\n", auc(roc_c), ci_c[1], ci_c[3]))

# Standardized coefficient comparison
cat(sprintf("\n  Standardized coefficients (Model C, z-scored predictors):\n"))
cat(sprintf("    wt_hmc_z:        beta=%.3f, OR=%.2f\n",
            coef(model_c_z)["wt_hmc_z"], exp(coef(model_c_z)["wt_hmc_z"])))
cat(sprintf("    gb_ctrl_signal_z: beta=%.3f, OR=%.2f\n",
            coef(model_c_z)["gb_ctrl_signal_z"], exp(coef(model_c_z)["gb_ctrl_signal_z"])))

# =============================================================================
# STEP 5: CHROMATIN STATE ANNOTATION (for 23e)
# =============================================================================

cat("\n--- Step 5: Annotating chromatin states ---\n")

chip_peaks <- list(
  ctcf     = load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF"),
  h3k27ac  = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac"),
  h3k27me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  h3k4me1  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1"),
  h3k4me3  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3"),
  bivalent = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
)
stopifnot("Missing ChIP peak files" = all(!sapply(chip_peaks, is.null)))

# Convert merged data to GRanges
merged_gr <- GRanges(
  seqnames = merged$chr,
  ranges = IRanges(start = merged$start, end = merged$end)
)

# Compute TSS distances
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gene_gr <- genes(txdb)
tss_gr <- resize(gene_gr, width = 1, fix = "start")
nearest_hits <- distanceToNearest(merged_gr, tss_gr)
distance_to_tss <- rep(NA_real_, length(merged_gr))
distance_to_tss[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

# Compute overlaps and classify
chip_overlaps <- compute_chip_overlaps(merged_gr, chip_peaks)
merged$chromatin_state <- classify_chromatin_state(chip_overlaps, distance_to_tss, TSS_THRESHOLD)

cat("  Chromatin state breakdown:\n")
state_counts <- table(merged$chromatin_state)
for (s in CHROMATIN_STATE_ORDER) {
  if (s %in% names(state_counts)) {
    cat(sprintf("    %s: %d genes\n", s, state_counts[s]))
  }
}

# =============================================================================
# FIGURE 23a: ROC CURVES (3 MODELS OVERLAID)
# =============================================================================

cat("\n--- Figure 23a: ROC curves ---\n")

# Build ROC data frames for ggplot
roc_to_df <- function(roc_obj, model_name) {
  data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities,
    model = model_name,
    stringsAsFactors = FALSE
  )
}

roc_df <- rbind(
  roc_to_df(roc_a, "Model A: Baseline 5hmC"),
  roc_to_df(roc_b, "Model B: K119ub signal"),
  roc_to_df(roc_c, "Model C: Combined")
)
roc_df$model <- factor(roc_df$model, levels = names(MODEL_COLORS))

# Legend labels with AUC
legend_labels <- c(
  sprintf("Model A: Baseline 5hmC (AUC=%.3f)", auc(roc_a)),
  sprintf("Model B: K119ub signal (AUC=%.3f)", auc(roc_b)),
  sprintf("Model C: Combined (AUC=%.3f)", auc(roc_c))
)

# Annotation text for AIC + R²
annot_text <- sprintf(
  "Model A: AIC=%.0f, R\u00B2=%.4f\nModel B: AIC=%.0f, R\u00B2=%.4f\nModel C: AIC=%.0f, R\u00B2=%.4f\nN=%s genes (%d DMR+)",
  AIC(model_a), r2_a,
  AIC(model_b), r2_b,
  AIC(model_c), r2_c,
  format(nrow(model_data), big.mark = ","), sum(model_data$significant)
)

p_23a <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_line(linewidth = 1.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = MODEL_COLORS, labels = legend_labels, name = NULL) +
  annotate("text", x = 0.95, y = 0.15, hjust = 1, vjust = 0, size = 3.2,
           label = annot_text) +
  labs(
    title = "ROC Curves: Predicting hmC DMR Status",
    subtitle = "Logistic regression models comparing baseline 5hmC vs K119ub as predictors",
    x = "1 \u2212 Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_biomodal() +
  theme(legend.position = c(0.65, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 9)) +
  coord_equal()

save_multiformat_ggplot(p_23a, file.path(OUTPUT_DIR, "23a_roc_curves"), 9, 9)

# =============================================================================
# FIGURE 23b: PREDICTED PROBABILITY CURVE (MODEL A)
# =============================================================================

cat("--- Figure 23b: Predicted probability curve ---\n")

# Add predicted probabilities
model_data$pred_prob_a <- predict(model_a, type = "response")

# OR and CI annotation for Model A predictor
or_a_pred <- or_a[2, ]
or_label <- sprintf("OR = %.2f [%.2f, %.2f]\n%s\nAUC = %.3f",
                    or_a_pred$or, or_a_pred$or_lower, or_a_pred$or_upper,
                    fmt_p(or_a_pred$p_value), auc(roc_a))

p_23b <- ggplot(model_data, aes(x = wt_hmc, y = as.numeric(significant))) +
  geom_point(aes(color = significant), alpha = 0.08, size = 0.8) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              color = MODEL_COLORS["Model A: Baseline 5hmC"],
              linewidth = 1.2, se = TRUE, fill = "#377EB844") +
  geom_rug(sides = "b", alpha = 0.02) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"),
                     labels = c("TRUE" = "DMR", "FALSE" = "Not DMR"),
                     name = "hmC DMR Status") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = scales::percent) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = or_label) +
  labs(
    title = "Predicted Probability of hmC DMR by Baseline WT 5hmC",
    subtitle = "Logistic regression: higher baseline 5hmC \u2192 higher DMR probability",
    x = "Baseline WT 5hmC (regional fraction, average of ctrl-M and ctrl-F)",
    y = "P(hmC DMR)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_23b, file.path(OUTPUT_DIR, "23b_predicted_probability"), 10, 7)

# =============================================================================
# FIGURE 23c: MODEL COMPARISON (AIC + AUC)
# =============================================================================

cat("--- Figure 23c: Model comparison ---\n")

comparison_df <- data.frame(
  model = factor(c("A: Baseline 5hmC", "B: K119ub signal", "C: Combined"),
                 levels = c("A: Baseline 5hmC", "B: K119ub signal", "C: Combined")),
  aic = c(AIC(model_a), AIC(model_b), AIC(model_c)),
  auc = c(auc(roc_a), auc(roc_b), auc(roc_c)),
  auc_lower = c(ci_a[1], ci_b[1], ci_c[1]),
  auc_upper = c(ci_a[3], ci_b[3], ci_c[3]),
  r2 = c(r2_a, r2_b, r2_c),
  color = c("#377EB8", "#756BB1", "#E41A1C"),
  stringsAsFactors = FALSE
)

# Left panel: AIC bars (lower = better)
p_23c_left <- ggplot(comparison_df, aes(x = model, y = aic, fill = model)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f", aic)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = setNames(comparison_df$color, comparison_df$model)) +
  labs(title = "AIC (lower = better)", x = NULL, y = "AIC") +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

# Right panel: AUC point + CI whiskers
p_23c_right <- ggplot(comparison_df, aes(x = model, y = auc, color = model)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = auc_lower, ymax = auc_upper), width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -1.5, size = 3.5) +
  scale_color_manual(values = setNames(comparison_df$color, comparison_df$model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  labs(title = "AUC (higher = better)", x = NULL, y = "AUC (95% CI)") +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

p_23c <- (p_23c_left | p_23c_right) +
  plot_annotation(
    title = "Logistic Regression Model Comparison",
    subtitle = sprintf("Predicting hmC DMR status (N=%s genes, %d DMR+)",
                       format(nrow(model_data), big.mark = ","),
                       sum(model_data$significant)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

save_multiformat_ggplot(p_23c, file.path(OUTPUT_DIR, "23c_model_comparison"), 12, 7)

# =============================================================================
# FIGURE 23d: DOSE-RESPONSE SCATTER (BASELINE 5hmC vs mod_difference)
# =============================================================================

cat("--- Figure 23d: Dose-response scatter ---\n")

# Spearman correlation (all genes)
rho_dose <- cor.test(merged$wt_hmc, merged$mod_difference, method = "spearman")
cat(sprintf("  Spearman rho (wt_hmc vs mod_difference): %.3f (%s)\n",
            rho_dose$estimate, fmt_p(rho_dose$p.value)))

# Identify top genes to label (most extreme mod_difference among DMRs)
merged <- merged %>%
  dplyr::mutate(is_key_gene = gene %in% KEY_GENES)

top_label <- merged %>%
  dplyr::filter(significant) %>%
  dplyr::arrange(mod_difference) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::bind_rows(
    merged %>%
      dplyr::filter(significant) %>%
      dplyr::arrange(dplyr::desc(mod_difference)) %>%
      dplyr::slice_head(n = 5)
  ) %>%
  dplyr::bind_rows(
    merged %>% dplyr::filter(is_key_gene)
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

p_23d <- ggplot(merged, aes(x = wt_hmc, y = mod_difference)) +
  geom_point(aes(color = significant), alpha = 0.15, size = 0.6) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE,
              fill = "grey80", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"),
                     labels = c("TRUE" = "Significant DMR", "FALSE" = "Not significant"),
                     name = "hmC DMR (q<0.05)") +
  geom_text_repel(data = top_label, aes(label = gene), size = 2.8,
                  max.overlaps = 25, segment.alpha = 0.5,
                  color = "black", fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.8,
           label = sprintf("Spearman rho = %.3f\n%s\nN = %s genes",
                           rho_dose$estimate, fmt_p(rho_dose$p.value),
                           format(nrow(merged), big.mark = ","))) +
  labs(
    title = "Dose-Response: Baseline WT 5hmC vs hmC Change in BAP1-KO",
    subtitle = "Higher baseline 5hmC \u2192 more hmC loss (negative mod_difference)",
    x = "Baseline WT 5hmC (regional fraction)",
    y = "hmC mod_difference (KO \u2212 WT)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_23d, file.path(OUTPUT_DIR, "23d_dose_response_scatter"), 10, 8)

# =============================================================================
# FIGURE 23e: DOSE-RESPONSE BY CHROMATIN STATE (FACETED)
# =============================================================================

cat("--- Figure 23e: Dose-response by chromatin state ---\n")

# Per-state Spearman correlations
state_rhos <- merged %>%
  dplyr::group_by(chromatin_state) %>%
  dplyr::summarise(
    n = dplyr::n(),
    rho = tryCatch(
      cor.test(wt_hmc, mod_difference, method = "spearman")$estimate,
      error = function(e) NA_real_
    ),
    p_value = tryCatch(
      cor.test(wt_hmc, mod_difference, method = "spearman")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = sprintf("rho=%.3f\nn=%s", rho, format(n, big.mark = ","))
  )

cat("  Per-chromatin-state Spearman rho:\n")
for (i in 1:nrow(state_rhos)) {
  cat(sprintf("    %s: rho=%.3f (n=%d, %s)\n",
              state_rhos$chromatin_state[i], state_rhos$rho[i],
              state_rhos$n[i], fmt_p(state_rhos$p_value[i])))
}

p_23e <- ggplot(merged, aes(x = wt_hmc, y = mod_difference)) +
  geom_point(aes(color = significant), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.7, se = TRUE,
              fill = "grey80", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"),
                     labels = c("TRUE" = "DMR", "FALSE" = "Not sig."),
                     name = "hmC DMR") +
  geom_text(data = state_rhos,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.3, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~chromatin_state, scales = "free_y", ncol = 4) +
  labs(
    title = "Dose-Response by Chromatin State",
    subtitle = "Baseline WT 5hmC vs hmC change, faceted by chromatin environment",
    x = "Baseline WT 5hmC (regional fraction)",
    y = "hmC mod_difference (KO \u2212 WT)"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9))

save_multiformat_ggplot(p_23e, file.path(OUTPUT_DIR, "23e_dose_response_by_chromatin"), 14, 10)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

# Table 1: Full merged data with predictions
export_all <- merged %>%
  dplyr::left_join(
    k119ub %>% dplyr::rename(gb_ctrl_signal = gb_ctrl_signal),
    by = "gene"
  ) %>%
  dplyr::mutate(
    pred_prob_a = predict(
      glm(significant ~ wt_hmc, data = merged, family = binomial),
      newdata = ., type = "response"
    )
  ) %>%
  dplyr::select(gene, chr, start, end,
                wt_hmc, hmc_ctrl_M, hmc_ctrl_F,
                gb_ctrl_signal,
                dmr_qvalue, mod_difference, significant,
                chromatin_state, pred_prob_a) %>%
  dplyr::arrange(dplyr::desc(wt_hmc))

write.table(export_all,
            file.path(TABLES_DIR, "baseline_hmc_predictor_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved baseline_hmc_predictor_all_genes.tsv (%d rows)\n", nrow(export_all)))

# Table 2: Model comparison summary
model_comparison <- data.frame(
  model = c("A: Baseline 5hmC", "B: K119ub signal", "C: Combined"),
  n_genes = rep(nrow(model_data), 3),
  n_dmr = rep(sum(model_data$significant), 3),
  aic = c(AIC(model_a), AIC(model_b), AIC(model_c)),
  auc = c(as.numeric(auc(roc_a)), as.numeric(auc(roc_b)), as.numeric(auc(roc_c))),
  auc_ci_lower = c(ci_a[1], ci_b[1], ci_c[1]),
  auc_ci_upper = c(ci_a[3], ci_b[3], ci_c[3]),
  mcfadden_r2 = c(r2_a, r2_b, r2_c),
  stringsAsFactors = FALSE
)

write.table(model_comparison,
            file.path(TABLES_DIR, "baseline_hmc_model_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved baseline_hmc_model_comparison.tsv (%d rows)\n", nrow(model_comparison)))

# Table 3: Detailed coefficient table (raw + standardized)
build_coef_table <- function(model, model_z, model_name) {
  raw_or <- extract_or(model)
  z_coefs <- coef(model_z)
  z_se <- summary(model_z)$coefficients[, 2]

  raw_or$model <- model_name
  raw_or$standardized_beta <- z_coefs[raw_or$term]
  raw_or$standardized_se <- z_se[raw_or$term]

  # Match z-scored term names back to raw terms
  term_map <- c("(Intercept)" = "(Intercept)",
                "wt_hmc" = "wt_hmc_z",
                "gb_ctrl_signal" = "gb_ctrl_signal_z")
  for (i in seq_len(nrow(raw_or))) {
    z_name <- term_map[raw_or$term[i]]
    if (!is.na(z_name) && z_name %in% names(z_coefs)) {
      raw_or$standardized_beta[i] <- z_coefs[z_name]
      raw_or$standardized_se[i] <- z_se[z_name]
    }
  }

  return(raw_or)
}

coef_table <- rbind(
  build_coef_table(model_a, model_a_z, "A: Baseline 5hmC"),
  build_coef_table(model_b, model_b_z, "B: K119ub signal"),
  build_coef_table(model_c, model_c_z, "C: Combined")
)

write.table(coef_table,
            file.path(TABLES_DIR, "baseline_hmc_model_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved baseline_hmc_model_coefficients.tsv (%d rows)\n", nrow(coef_table)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 23 SUMMARY: BASELINE 5hmC AS PREDICTOR OF DMR SUSCEPTIBILITY\n")
cat("========================================================================\n\n")

cat(sprintf("Total genes with baseline WT 5hmC: %d\n", nrow(hmc_baseline)))
cat(sprintf("Merged with hmC DMR status: %d genes\n", nrow(merged)))
cat(sprintf("Combined with K119ub signal: %d genes\n", nrow(combined)))
cat(sprintf("DMR+ (q<0.05): %d genes (%.1f%%)\n",
            sum(merged$significant), 100 * mean(merged$significant)))

cat(sprintf("\nLogistic Regression Model Comparison:\n"))
cat(sprintf("  %-25s AIC=%-9.0f AUC=%.3f [%.3f, %.3f] R\u00B2=%.4f\n",
            "Model A (Baseline 5hmC):", AIC(model_a), auc(roc_a), ci_a[1], ci_a[3], r2_a))
cat(sprintf("  %-25s AIC=%-9.0f AUC=%.3f [%.3f, %.3f] R\u00B2=%.4f\n",
            "Model B (K119ub):", AIC(model_b), auc(roc_b), ci_b[1], ci_b[3], r2_b))
cat(sprintf("  %-25s AIC=%-9.0f AUC=%.3f [%.3f, %.3f] R\u00B2=%.4f\n",
            "Model C (Combined):", AIC(model_c), auc(roc_c), ci_c[1], ci_c[3], r2_c))

cat(sprintf("\nStandardized coefficients (Model C):\n"))
cat(sprintf("  wt_hmc_z:         beta=%.3f (OR=%.2f per SD)\n",
            coef(model_c_z)["wt_hmc_z"], exp(coef(model_c_z)["wt_hmc_z"])))
cat(sprintf("  gb_ctrl_signal_z: beta=%.3f (OR=%.2f per SD)\n",
            coef(model_c_z)["gb_ctrl_signal_z"], exp(coef(model_c_z)["gb_ctrl_signal_z"])))

cat(sprintf("\nDose-response (all genes):\n"))
cat(sprintf("  Spearman rho (wt_hmc vs mod_difference): %.3f (%s)\n",
            rho_dose$estimate, fmt_p(rho_dose$p.value)))

cat(sprintf("\nInterpretation:\n"))
if (auc(roc_a) > auc(roc_b)) {
  cat("  Baseline 5hmC is a STRONGER predictor of DMR status than K119ub signal.\n")
  cat("  This supports the substrate availability hypothesis: genes with more 5hmC\n")
  cat("  have more substrate for TET impairment to affect in BAP1-KO.\n")
} else {
  cat("  K119ub signal is a STRONGER predictor of DMR status than baseline 5hmC.\n")
  cat("  This suggests chromatin-level PRC1 occupancy is more predictive than\n")
  cat("  substrate abundance alone.\n")
}

if (as.numeric(rho_dose$estimate) < 0) {
  cat(sprintf("  Negative dose-response (rho=%.3f) confirms: higher baseline 5hmC\n",
              rho_dose$estimate))
  cat("  leads to more 5hmC loss in BAP1-KO, consistent with TET impairment.\n")
}

cat("\n--- Output files ---\n")
cat(sprintf("  Figures: %s/23{a-e}_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/baseline_hmc_*.tsv (3 files)\n", TABLES_DIR))

cat("\n=== Section 23 complete ===\n")
