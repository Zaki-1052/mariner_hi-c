# biomodal/downstream/scripts/viz_sections/section_24_dnmt3a_prediction.R
# Section 24: DNMT3A Binding Prediction — Dual-Mechanism Hypothesis Testing
# Standalone script - sources shared config for all dependencies and data
#
# Tests whether gene body hypermethylation (mC-up DMRs) in BAP1-KO can be
# computationally predicted from features related to the DNMT3A-UDR recruitment
# model (Chen et al. 2024). The UDR domain makes a bidentate interaction with
# H2AK119ub and the nucleosome acidic patch, providing a structural basis for
# PRC1-dependent de novo methylation.
#
# Dual-mechanism hypothesis:
#   1. DNMT3A recruitment — UDR-K119ub interaction recruits DNMT3A to gene bodies
#   2. TET impediment    — K119ub blocks TET access, reducing 5hmC turnover
#
# Compares predictive power of "DNMT3A recruitment features" vs "TET impediment
# features" for mC hypermethylation using logistic regression and random forest.
#
# Data sources:
#   1. mC DMR BED (from _shared_config.R) — outcome + baseline mC + gene length + CpG count
#   2. K119ub gene signal TSV — WT gene body K119ub signal
#   3. ATAC consensus ctrl peaks — gene-level accessibility (via ChIPseeker)
#   4. hmC feature extraction (per-sample regional fractions) — baseline 5hmC
#   5. RNA-seq Excel — baseMean expression level
#   6. ChIP-seq peaks (from _shared_config.R) — chromatin state classification
#
# Figures:
#   24a: Feature correlation heatmap (Spearman, 7 predictors)
#   24b: ROC curves (5 logistic regression models overlaid)
#   24c: Feature importance (2-panel: LR standardized betas + RF importance)
#   24d: Model comparison (AIC + AUC bar/dot chart, 5 models)
#   24e: Predicted probability by K119ub signal (sigmoid curve)
#   24f: Predicted probability by chromatin state (faceted violin/box)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_24_dnmt3a_prediction.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(pROC)
  library(ChIPseeker)
  library(readxl)
  library(randomForest)
})

# =============================================================================
# SECTION 24 CONFIGURATION
# =============================================================================

# Feature extraction paths (same as section 23)
EXTRACT_DIR <- file.path(BASE_DIR, "modality/outputs/run-3/outputs_CG/Results",
                         "gencode.vM25.mouse.genes.annotation/Extract_20260221_185106")
HMC_EXTRACT_FILE <- file.path(EXTRACT_DIR, "Extract_hmc_regional-frac_20260221_185106.tsv.gz")

# K119ub gene signal file (same as section 23)
K119UB_SIGNAL_FILE <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")

# RNA-seq file (same as section 20)
RNA_SEQ_FILE <- "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"

# Model colors (5 models)
MODEL_COLORS <- c(
  "Full"              = "#E41A1C",
  "DNMT3A recruitment" = "#756BB1",
  "TET impediment"    = "#377EB8",
  "K119ub only"       = "#4DAF4A",
  "Stepwise"          = "#FF7F00"
)

# Section output directory
SECTION_DIR <- file.path(OUTPUT_DIR, "24_dnmt3a_prediction")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: format p-value (same as section 23)
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
cat("SECTION 24: DNMT3A BINDING PREDICTION — DUAL-MECHANISM HYPOTHESIS\n")
cat("========================================================================\n\n")

cat("Validating inputs...\n")
stopifnot("mc_dmr not loaded from shared config" = exists("mc_dmr") && !is.null(mc_dmr))
stopifnot("hmc_dmr not loaded from shared config" = exists("hmc_dmr") && !is.null(hmc_dmr))
stopifnot("hmC extract file not found" = file.exists(HMC_EXTRACT_FILE))
stopifnot("K119ub gene signal file not found" = file.exists(K119UB_SIGNAL_FILE))
stopifnot("RNA-seq Excel file not found" = file.exists(RNA_SEQ_FILE))
stopifnot("ATAC consensus ctrl file not found" = file.exists(ATAC_FILES$consensus_ctrl))
cat("  All inputs validated.\n\n")

# =============================================================================
# STEP 1: BUILD FEATURE MATRIX
# =============================================================================

cat("--- Step 1: Building feature matrix ---\n\n")

# --- 1a: mC DMR features (outcome + baseline mC + gene geometry) ---
cat("  1a: Extracting mC DMR features...\n")

mc_features <- mc_dmr %>%
  dplyr::select(gene, chr, start, end,
                baseline_mc = mean_mod_group1,
                num_contexts, mod_difference, dmr_qvalue, significant) %>%
  dplyr::mutate(
    gene_length = end - start,
    cpg_density = num_contexts / gene_length,
    hyper_dmr = significant & mod_difference > 0
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

cat(sprintf("    %d genes from mC DMR data\n", nrow(mc_features)))
cat(sprintf("    Hyper-DMR (mC-up, q<0.05): %d genes (%.1f%%)\n",
            sum(mc_features$hyper_dmr), 100 * mean(mc_features$hyper_dmr)))

# --- 1b: K119ub signal ---
cat("  1b: Loading K119ub gene signal...\n")

k119ub_raw <- read.table(K119UB_SIGNAL_FILE, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
k119ub <- k119ub_raw %>%
  dplyr::select(gene = symbol, gb_ctrl_signal) %>%
  dplyr::filter(is.finite(gb_ctrl_signal)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::rename(k119ub = gb_ctrl_signal)

cat(sprintf("    %d genes with finite WT K119ub signal\n", nrow(k119ub)))

# --- 1c: ATAC accessibility (consensus ctrl peaks per gene) ---
cat("  1c: Annotating ATAC consensus ctrl peaks to genes...\n")

atac_ctrl_gr <- load_chip_peaks(ATAC_FILES$consensus_ctrl, "ATAC consensus ctrl")
stopifnot("ATAC consensus ctrl peaks failed to load" = !is.null(atac_ctrl_gr))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
atac_anno <- suppressMessages(annotatePeak(
  atac_ctrl_gr, TxDb = txdb, annoDb = "org.Mm.eg.db", level = "gene"
))
atac_anno_df <- as.data.frame(atac_anno)

atac_per_gene <- atac_anno_df %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(atac_count = dplyr::n(), .groups = "drop") %>%
  dplyr::rename(gene = SYMBOL)

cat(sprintf("    %d genes with at least one ATAC peak\n", nrow(atac_per_gene)))

# --- 1d: Baseline 5hmC (average of ctrl-M and ctrl-F) ---
cat("  1d: Loading baseline WT 5hmC from feature extraction...\n")

hmc_extract <- read.table(gzfile(HMC_EXTRACT_FILE), header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, na.strings = c("", "NA"))
names(hmc_extract)[4:7] <- c("hmc_ctrl_M", "hmc_mut_M", "hmc_ctrl_F", "hmc_mut_F")

hmc_baseline <- hmc_extract %>%
  dplyr::select(gene = Name, hmc_ctrl_M, hmc_ctrl_F) %>%
  dplyr::mutate(baseline_hmc = (hmc_ctrl_M + hmc_ctrl_F) / 2) %>%
  dplyr::select(gene, baseline_hmc) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::filter(complete.cases(baseline_hmc))

cat(sprintf("    %d genes with baseline WT 5hmC\n", nrow(hmc_baseline)))

# --- 1e: RNA-seq expression (baseMean) ---
cat("  1e: Loading RNA-seq expression levels...\n")

rna_all <- read_excel(RNA_SEQ_FILE) %>%
  dplyr::select(gene = ensembl_gene_id, baseMean) %>%
  dplyr::filter(!is.na(gene) & gene != "" & !is.na(baseMean)) %>%
  dplyr::mutate(log_expression = log10(baseMean + 1)) %>%
  dplyr::select(gene, log_expression) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

cat(sprintf("    %d genes with RNA-seq baseMean\n", nrow(rna_all)))

# --- 1f: Inner join cascade ---
cat("\n  1f: Merging all features...\n")

feature_matrix <- mc_features %>%
  dplyr::inner_join(k119ub, by = "gene") %>%
  dplyr::inner_join(atac_per_gene, by = "gene") %>%
  dplyr::inner_join(hmc_baseline, by = "gene") %>%
  dplyr::inner_join(rna_all, by = "gene")

# Log-transform gene_length
feature_matrix <- feature_matrix %>%
  dplyr::mutate(log_gene_length = log10(gene_length))

cat(sprintf("    Final feature matrix: %d genes with all 7 features\n", nrow(feature_matrix)))
cat(sprintf("    Hyper-DMR (outcome): %d genes (%.1f%%)\n",
            sum(feature_matrix$hyper_dmr),
            100 * mean(feature_matrix$hyper_dmr)))

# Filter to complete cases for modeling
model_data <- feature_matrix %>%
  dplyr::select(gene, chr, start, end,
                hyper_dmr, k119ub, atac_count, cpg_density,
                baseline_mc, baseline_hmc, log_gene_length, log_expression) %>%
  dplyr::filter(complete.cases(.))

cat(sprintf("    Complete cases for modeling: %d genes\n", nrow(model_data)))
cat(sprintf("    Hyper-DMR rate: %d / %d (%.1f%%)\n",
            sum(model_data$hyper_dmr), nrow(model_data),
            100 * mean(model_data$hyper_dmr)))

# =============================================================================
# STEP 2: FIT LOGISTIC REGRESSION MODELS
# =============================================================================

cat("\n--- Step 2: Fitting logistic regression models ---\n\n")

# McFadden pseudo-R²
null_model <- glm(hyper_dmr ~ 1, data = model_data, family = binomial)
null_ll <- as.numeric(logLik(null_model))
mcfadden <- function(model) {
  1 - as.numeric(logLik(model)) / null_ll
}

# Extract OR with 95% CI
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

# --- Model 1: Full ---
cat("  Model 1: Full (all 7 features)\n")
model_full <- glm(
  hyper_dmr ~ k119ub + atac_count + cpg_density + baseline_mc +
    baseline_hmc + log_gene_length + log_expression,
  data = model_data, family = binomial
)

# --- Model 2: DNMT3A recruitment ---
cat("  Model 2: DNMT3A recruitment (k119ub + atac_count + cpg_density)\n")
model_dnmt3a <- glm(
  hyper_dmr ~ k119ub + atac_count + cpg_density,
  data = model_data, family = binomial
)

# --- Model 3: TET impediment ---
cat("  Model 3: TET impediment (baseline_hmc + atac_count)\n")
model_tet <- glm(
  hyper_dmr ~ baseline_hmc + atac_count,
  data = model_data, family = binomial
)

# --- Model 4: K119ub only ---
cat("  Model 4: K119ub only\n")
model_k119ub <- glm(
  hyper_dmr ~ k119ub,
  data = model_data, family = binomial
)

# --- Model 5: Stepwise ---
cat("  Model 5: Stepwise (AIC-based bidirectional selection)\n")
model_step <- step(model_full, direction = "both", trace = 0)
cat(sprintf("    Selected terms: %s\n",
            paste(names(coef(model_step))[-1], collapse = " + ")))

# --- Standardized versions (z-scored predictors) ---
predictor_cols <- c("k119ub", "atac_count", "cpg_density", "baseline_mc",
                    "baseline_hmc", "log_gene_length", "log_expression")

model_data_z <- model_data
for (col in predictor_cols) {
  model_data_z[[paste0(col, "_z")]] <- as.numeric(scale(model_data_z[[col]]))
}

model_full_z <- glm(
  hyper_dmr ~ k119ub_z + atac_count_z + cpg_density_z + baseline_mc_z +
    baseline_hmc_z + log_gene_length_z + log_expression_z,
  data = model_data_z, family = binomial
)

model_dnmt3a_z <- glm(
  hyper_dmr ~ k119ub_z + atac_count_z + cpg_density_z,
  data = model_data_z, family = binomial
)

model_tet_z <- glm(
  hyper_dmr ~ baseline_hmc_z + atac_count_z,
  data = model_data_z, family = binomial
)

model_k119ub_z <- glm(
  hyper_dmr ~ k119ub_z,
  data = model_data_z, family = binomial
)

# Rebuild stepwise formula with z-scored terms
step_terms <- names(coef(model_step))[-1]
step_formula_z <- as.formula(paste("hyper_dmr ~",
                                   paste0(step_terms, "_z", collapse = " + ")))
model_step_z <- glm(step_formula_z, data = model_data_z, family = binomial)

# --- ROC / AUC ---
roc_full   <- roc(model_data$hyper_dmr, predict(model_full, type = "response"), quiet = TRUE)
roc_dnmt3a <- roc(model_data$hyper_dmr, predict(model_dnmt3a, type = "response"), quiet = TRUE)
roc_tet    <- roc(model_data$hyper_dmr, predict(model_tet, type = "response"), quiet = TRUE)
roc_k119ub <- roc(model_data$hyper_dmr, predict(model_k119ub, type = "response"), quiet = TRUE)
roc_step   <- roc(model_data$hyper_dmr, predict(model_step, type = "response"), quiet = TRUE)

ci_full   <- ci.auc(roc_full, quiet = TRUE)
ci_dnmt3a <- ci.auc(roc_dnmt3a, quiet = TRUE)
ci_tet    <- ci.auc(roc_tet, quiet = TRUE)
ci_k119ub <- ci.auc(roc_k119ub, quiet = TRUE)
ci_step   <- ci.auc(roc_step, quiet = TRUE)

# AIC and R²
models_list <- list(
  Full = model_full, `DNMT3A recruitment` = model_dnmt3a,
  `TET impediment` = model_tet, `K119ub only` = model_k119ub,
  Stepwise = model_step
)
models_z_list <- list(
  Full = model_full_z, `DNMT3A recruitment` = model_dnmt3a_z,
  `TET impediment` = model_tet_z, `K119ub only` = model_k119ub_z,
  Stepwise = model_step_z
)
roc_list <- list(
  Full = roc_full, `DNMT3A recruitment` = roc_dnmt3a,
  `TET impediment` = roc_tet, `K119ub only` = roc_k119ub,
  Stepwise = roc_step
)
ci_list <- list(
  Full = ci_full, `DNMT3A recruitment` = ci_dnmt3a,
  `TET impediment` = ci_tet, `K119ub only` = ci_k119ub,
  Stepwise = ci_step
)

# Print model summaries
cat("\n  Model comparison:\n")
for (nm in names(models_list)) {
  mod <- models_list[[nm]]
  roc_obj <- roc_list[[nm]]
  ci_obj <- ci_list[[nm]]
  cat(sprintf("    %-22s AIC=%8.1f  AUC=%.3f [%.3f, %.3f]  R\u00B2=%.4f\n",
              paste0(nm, ":"), AIC(mod), auc(roc_obj), ci_obj[1], ci_obj[3],
              mcfadden(mod)))
}

# VIF check on full model
cat("\n  VIF check (full model):\n")
vif_vals <- tryCatch({
  # Compute VIF manually to avoid car package dependency
  predictor_names <- names(coef(model_full))[-1]
  vif_result <- numeric(length(predictor_names))
  names(vif_result) <- predictor_names
  for (i in seq_along(predictor_names)) {
    pred_i <- predictor_names[i]
    other_preds <- predictor_names[-i]
    formula_vif <- as.formula(paste(pred_i, "~", paste(other_preds, collapse = " + ")))
    r2_i <- summary(lm(formula_vif, data = model_data))$r.squared
    vif_result[i] <- 1 / (1 - r2_i)
  }
  vif_result
}, error = function(e) {
  cat(sprintf("    VIF computation failed: %s\n", e$message))
  NULL
})

if (!is.null(vif_vals)) {
  for (i in seq_along(vif_vals)) {
    flag <- ifelse(vif_vals[i] > 5, " *** HIGH", "")
    cat(sprintf("    %-20s VIF = %.2f%s\n", names(vif_vals)[i], vif_vals[i], flag))
  }
}

# DeLong test: DNMT3A recruitment vs TET impediment
cat("\n  DeLong test (DNMT3A recruitment vs TET impediment):\n")
delong_result <- roc.test(roc_dnmt3a, roc_tet, method = "delong")
cat(sprintf("    AUC difference: %.4f (DNMT3A=%.3f vs TET=%.3f)\n",
            as.numeric(auc(roc_dnmt3a)) - as.numeric(auc(roc_tet)),
            auc(roc_dnmt3a), auc(roc_tet)))
cat(sprintf("    %s\n", fmt_p(delong_result$p.value)))

# Likelihood ratio test: DNMT3A vs K119ub-only (nested)
lr_stat <- -2 * (as.numeric(logLik(model_k119ub)) - as.numeric(logLik(model_dnmt3a)))
lr_df <- length(coef(model_dnmt3a)) - length(coef(model_k119ub))
lr_pval <- pchisq(lr_stat, df = lr_df, lower.tail = FALSE)
cat(sprintf("\n  Likelihood ratio test (DNMT3A recruitment vs K119ub only):\n"))
cat(sprintf("    LR statistic = %.2f, df = %d, %s\n", lr_stat, lr_df, fmt_p(lr_pval)))

# Print standardized coefficients for full model
cat("\n  Standardized coefficients (Full model, z-scored):\n")
full_z_coefs <- coef(model_full_z)
full_z_se <- summary(model_full_z)$coefficients[, 2]
full_z_p <- summary(model_full_z)$coefficients[, 4]
for (term in names(full_z_coefs)[-1]) {
  cat(sprintf("    %-25s beta=%7.3f  OR=%5.2f  %s\n",
              gsub("_z$", "", term),
              full_z_coefs[term], exp(full_z_coefs[term]),
              fmt_p(full_z_p[term])))
}

# =============================================================================
# STEP 3: RANDOM FOREST FEATURE IMPORTANCE
# =============================================================================

cat("\n--- Step 3: Random forest feature importance ---\n")

set.seed(42)
rf_data <- model_data %>%
  dplyr::select(hyper_dmr, all_of(predictor_cols)) %>%
  dplyr::mutate(hyper_dmr = factor(hyper_dmr, levels = c(FALSE, TRUE)))

rf_model <- randomForest(
  hyper_dmr ~ ., data = rf_data, ntree = 500,
  importance = TRUE, classwt = c("FALSE" = 1, "TRUE" = sum(!model_data$hyper_dmr) / sum(model_data$hyper_dmr))
)

rf_importance <- as.data.frame(importance(rf_model))
rf_importance$feature <- rownames(rf_importance)
rf_importance <- rf_importance %>%
  dplyr::arrange(dplyr::desc(MeanDecreaseAccuracy))

cat("  Random forest variable importance (sorted by MeanDecreaseAccuracy):\n")
for (i in seq_len(nrow(rf_importance))) {
  cat(sprintf("    %-20s Accuracy=%.2f  Gini=%.2f\n",
              rf_importance$feature[i],
              rf_importance$MeanDecreaseAccuracy[i],
              rf_importance$MeanDecreaseGini[i]))
}

cat(sprintf("  RF OOB error rate: %.1f%%\n", 100 * rf_model$err.rate[500, "OOB"]))

# =============================================================================
# STEP 4: CHROMATIN STATE ANNOTATION
# =============================================================================

cat("\n--- Step 4: Annotating chromatin states ---\n")

chip_peaks <- list(
  ctcf     = load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF"),
  h3k27ac  = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac"),
  h3k27me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  h3k4me1  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1"),
  h3k4me3  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3"),
  bivalent = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
)
stopifnot("Missing ChIP peak files" = all(!sapply(chip_peaks, is.null)))

model_gr <- GRanges(
  seqnames = model_data$chr,
  ranges = IRanges(start = model_data$start, end = model_data$end)
)

gene_gr <- genes(txdb)
tss_gr <- resize(gene_gr, width = 1, fix = "start")
nearest_hits <- distanceToNearest(model_gr, tss_gr)
distance_to_tss <- rep(NA_real_, length(model_gr))
distance_to_tss[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

chip_overlaps <- compute_chip_overlaps(model_gr, chip_peaks)
model_data$chromatin_state <- classify_chromatin_state(chip_overlaps, distance_to_tss, TSS_THRESHOLD)

cat("  Chromatin state breakdown:\n")
state_counts <- table(model_data$chromatin_state)
for (s in CHROMATIN_STATE_ORDER) {
  if (s %in% names(state_counts)) {
    n_hyper <- sum(model_data$hyper_dmr[model_data$chromatin_state == s])
    cat(sprintf("    %-20s %5d genes (%d hyper-DMR, %.1f%%)\n",
                s, state_counts[s], n_hyper,
                100 * n_hyper / state_counts[s]))
  }
}

# Add full-model predicted probabilities
model_data$pred_prob_full <- predict(model_full, newdata = model_data, type = "response")

# =============================================================================
# FIGURE 24a: FEATURE CORRELATION HEATMAP
# =============================================================================

cat("\n--- Figure 24a: Feature correlation heatmap ---\n")

cor_matrix <- cor(
  model_data[, predictor_cols],
  method = "spearman", use = "complete.obs"
)

# Clean feature names for display
display_names <- c(
  k119ub = "H2AK119ub", atac_count = "ATAC peaks",
  cpg_density = "CpG density", baseline_mc = "Baseline 5mC",
  baseline_hmc = "Baseline 5hmC", log_gene_length = "Gene length (log10)",
  log_expression = "Expression (log10)"
)
rownames(cor_matrix) <- display_names[rownames(cor_matrix)]
colnames(cor_matrix) <- display_names[colnames(cor_matrix)]

# Save correlation heatmap via save_multiformat_pheatmap (produces pdf+svg+jpg)
save_multiformat_pheatmap(
  quote(pheatmap(
    cor_matrix,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = sprintf("Feature Correlation Matrix (Spearman, N=%s genes)",
                   format(nrow(model_data), big.mark = ",")),
    fontsize = 11,
    fontsize_number = 9,
    border_color = "grey90",
    clustering_method = "ward.D2"
  )),
  file.path(OUTPUT_DIR, "24a_feature_correlation"),
  width = 9, height = 8
)

# =============================================================================
# FIGURE 24b: ROC CURVES (5 MODELS OVERLAID)
# =============================================================================

cat("\n--- Figure 24b: ROC curves ---\n")

roc_to_df <- function(roc_obj, model_name) {
  data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities,
    model = model_name,
    stringsAsFactors = FALSE
  )
}

roc_df <- rbind(
  roc_to_df(roc_full, "Full"),
  roc_to_df(roc_dnmt3a, "DNMT3A recruitment"),
  roc_to_df(roc_tet, "TET impediment"),
  roc_to_df(roc_k119ub, "K119ub only"),
  roc_to_df(roc_step, "Stepwise")
)
roc_df$model <- factor(roc_df$model, levels = names(MODEL_COLORS))

legend_labels <- sprintf(
  "%s (AUC=%.3f)",
  names(MODEL_COLORS),
  c(auc(roc_full), auc(roc_dnmt3a), auc(roc_tet), auc(roc_k119ub), auc(roc_step))
)

annot_text <- sprintf(
  paste0(
    "DeLong: DNMT3A vs TET %s\n",
    "N=%s genes (%d hyper-DMR)"
  ),
  fmt_p(delong_result$p.value),
  format(nrow(model_data), big.mark = ","), sum(model_data$hyper_dmr)
)

p_24b <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_line(linewidth = 1.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = MODEL_COLORS, labels = legend_labels, name = NULL) +
  annotate("text", x = 0.95, y = 0.15, hjust = 1, vjust = 0, size = 3.2,
           label = annot_text) +
  labs(
    title = "ROC Curves: Predicting mC Hypermethylation (hyper-DMR)",
    subtitle = "DNMT3A recruitment vs TET impediment vs full model",
    x = "1 \u2212 Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_biomodal() +
  theme(legend.position = c(0.65, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 9)) +
  coord_equal()

save_multiformat_ggplot(p_24b, file.path(OUTPUT_DIR, "24b_roc_curves"), 10, 10)

# =============================================================================
# FIGURE 24c: FEATURE IMPORTANCE (2-PANEL)
# =============================================================================

cat("--- Figure 24c: Feature importance ---\n")

# Left panel: LR standardized betas (full model)
lr_betas <- data.frame(
  feature = gsub("_z$", "", names(full_z_coefs)[-1]),
  beta = as.numeric(full_z_coefs[-1]),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    display_name = display_names[feature],
    sign = ifelse(beta > 0, "Positive", "Negative"),
    abs_beta = abs(beta)
  ) %>%
  dplyr::arrange(abs_beta)

lr_betas$display_name <- factor(lr_betas$display_name, levels = lr_betas$display_name)

p_24c_left <- ggplot(lr_betas, aes(x = abs_beta, y = display_name, fill = sign)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", beta), x = abs_beta + max(abs_beta) * 0.03),
            hjust = 0, size = 3.2) +
  scale_fill_manual(values = c("Positive" = "#D7191C", "Negative" = "#2C7BB6"),
                    name = "Direction") +
  labs(title = "Logistic Regression\n(standardized betas)",
       x = "|Standardized beta|", y = NULL) +
  theme_biomodal() +
  theme(legend.position = "bottom")

# Right panel: RF variable importance (MeanDecreaseAccuracy)
rf_plot_data <- rf_importance %>%
  dplyr::mutate(
    display_name = display_names[feature],
    display_name = factor(display_name, levels = display_name[order(MeanDecreaseAccuracy)])
  )

p_24c_right <- ggplot(rf_plot_data, aes(x = MeanDecreaseAccuracy, y = display_name)) +
  geom_col(fill = "#4DAF4A", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", MeanDecreaseAccuracy),
                x = MeanDecreaseAccuracy + max(MeanDecreaseAccuracy) * 0.03),
            hjust = 0, size = 3.2) +
  labs(title = "Random Forest\n(MeanDecreaseAccuracy)",
       x = "Mean Decrease in Accuracy", y = NULL) +
  theme_biomodal()

p_24c <- (p_24c_left | p_24c_right) +
  plot_annotation(
    title = "Feature Importance: Predicting mC Hypermethylation",
    subtitle = sprintf("N=%s genes, %d hyper-DMR (%.1f%%)",
                       format(nrow(model_data), big.mark = ","),
                       sum(model_data$hyper_dmr),
                       100 * mean(model_data$hyper_dmr)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

save_multiformat_ggplot(p_24c, file.path(OUTPUT_DIR, "24c_feature_importance"), 14, 8)

# =============================================================================
# FIGURE 24d: MODEL COMPARISON (AIC + AUC)
# =============================================================================

cat("--- Figure 24d: Model comparison ---\n")

comparison_df <- data.frame(
  model = factor(names(MODEL_COLORS), levels = names(MODEL_COLORS)),
  aic = sapply(models_list, AIC),
  auc = sapply(roc_list, function(r) as.numeric(auc(r))),
  auc_lower = sapply(ci_list, function(c) c[1]),
  auc_upper = sapply(ci_list, function(c) c[3]),
  r2 = sapply(models_list, mcfadden),
  color = as.character(MODEL_COLORS),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Left panel: AIC bars
p_24d_left <- ggplot(comparison_df, aes(x = model, y = aic, fill = model)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f", aic)), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = MODEL_COLORS) +
  labs(title = "AIC (lower = better)", x = NULL, y = "AIC") +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

# Right panel: AUC points + CI
delong_label <- sprintf("DeLong p=%s\n(DNMT3A vs TET)", fmt_p(delong_result$p.value))

p_24d_right <- ggplot(comparison_df, aes(x = model, y = auc, color = model)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = auc_lower, ymax = auc_upper), width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -1.5, size = 3.2) +
  scale_color_manual(values = MODEL_COLORS) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 4.5, y = min(comparison_df$auc_lower) - 0.01,
           label = delong_label, size = 3, hjust = 0.5) +
  labs(title = "AUC (higher = better)", x = NULL, y = "AUC (95% CI)") +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

p_24d <- (p_24d_left | p_24d_right) +
  plot_annotation(
    title = "Model Comparison: Predicting mC Hypermethylation",
    subtitle = sprintf("5 logistic regression models (N=%s genes, %d hyper-DMR)",
                       format(nrow(model_data), big.mark = ","),
                       sum(model_data$hyper_dmr)),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

save_multiformat_ggplot(p_24d, file.path(OUTPUT_DIR, "24d_model_comparison"), 14, 7)

# =============================================================================
# FIGURE 24e: PREDICTED PROBABILITY BY K119ub SIGNAL
# =============================================================================

cat("--- Figure 24e: Predicted probability by K119ub signal ---\n")

or_k119ub_info <- extract_or(model_k119ub)
or_label <- sprintf("OR = %.2f [%.2f, %.2f]\n%s\nAUC = %.3f",
                    or_k119ub_info$or[2], or_k119ub_info$or_lower[2],
                    or_k119ub_info$or_upper[2],
                    fmt_p(or_k119ub_info$p_value[2]), auc(roc_k119ub))

p_24e <- ggplot(model_data, aes(x = k119ub, y = as.numeric(hyper_dmr))) +
  geom_point(aes(color = hyper_dmr), alpha = 0.08, size = 0.8) +
  geom_smooth(method = "glm", method.args = list(family = binomial),
              color = MODEL_COLORS["K119ub only"],
              linewidth = 1.2, se = TRUE, fill = "#4DAF4A44") +
  geom_rug(sides = "b", alpha = 0.02) +
  scale_color_manual(values = c("TRUE" = "#D7191C", "FALSE" = "grey70"),
                     labels = c("TRUE" = "Hyper-DMR", "FALSE" = "Not hyper-DMR"),
                     name = "mC DMR Status") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     labels = scales::percent) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = or_label) +
  labs(
    title = "Predicted Probability of mC Hypermethylation by K119ub Signal",
    subtitle = "Logistic regression: higher K119ub \u2192 higher hypermethylation probability (DNMT3A-UDR model)",
    x = "WT H2AK119ub Gene Body Signal",
    y = "P(mC Hyper-DMR)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_24e, file.path(OUTPUT_DIR, "24e_k119ub_predicted_probability"), 10, 7)

# =============================================================================
# FIGURE 24f: PREDICTED PROBABILITY BY CHROMATIN STATE
# =============================================================================

cat("--- Figure 24f: Predicted probability by chromatin state ---\n")

# Per-state AUC
state_auc <- model_data %>%
  dplyr::group_by(chromatin_state) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_hyper = sum(hyper_dmr),
    auc_val = tryCatch({
      if (sum(hyper_dmr) >= 5 & sum(!hyper_dmr) >= 5) {
        as.numeric(auc(roc(hyper_dmr, pred_prob_full, quiet = TRUE)))
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = ifelse(is.na(auc_val),
                   sprintf("n=%s\n(%d hyper)", format(n, big.mark = ","), n_hyper),
                   sprintf("AUC=%.2f\nn=%s (%d hyper)",
                           auc_val, format(n, big.mark = ","), n_hyper))
  )

cat("  Per-chromatin-state AUC (full model):\n")
for (i in seq_len(nrow(state_auc))) {
  auc_str <- ifelse(is.na(state_auc$auc_val[i]), "N/A (too few)",
                    sprintf("%.3f", state_auc$auc_val[i]))
  cat(sprintf("    %-20s AUC=%s  (n=%d, %d hyper)\n",
              state_auc$chromatin_state[i], auc_str,
              state_auc$n[i], state_auc$n_hyper[i]))
}

p_24f <- ggplot(model_data, aes(x = hyper_dmr, y = pred_prob_full, fill = hyper_dmr)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  facet_wrap(~chromatin_state, ncol = 4, scales = "free_y") +
  geom_text(data = state_auc,
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.3, size = 2.6, inherit.aes = FALSE) +
  scale_fill_manual(values = c("TRUE" = "#D7191C", "FALSE" = "#2C7BB6"),
                    labels = c("TRUE" = "Hyper-DMR", "FALSE" = "Not hyper"),
                    name = "mC Status") +
  scale_x_discrete(labels = c("TRUE" = "Hyper", "FALSE" = "Not")) +
  labs(
    title = "Full Model Predictions by Chromatin State",
    subtitle = "Predicted P(hyper-DMR) from full logistic model, split by actual status",
    x = "Actual mC Hyper-DMR Status",
    y = "Predicted P(mC Hyper-DMR)"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9))

save_multiformat_ggplot(p_24f, file.path(OUTPUT_DIR, "24f_chromatin_state_predictions"), 14, 10)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

# Table 1: Full per-gene feature matrix + predicted probabilities
export_features <- feature_matrix %>%
  dplyr::left_join(
    model_data %>% dplyr::select(gene, chromatin_state, pred_prob_full),
    by = "gene"
  ) %>%
  dplyr::select(
    gene, chr, start, end,
    hyper_dmr, baseline_mc, baseline_hmc, k119ub, atac_count,
    cpg_density, log_gene_length, log_expression,
    chromatin_state, pred_prob_full,
    mod_difference, dmr_qvalue, significant
  ) %>%
  dplyr::arrange(dplyr::desc(pred_prob_full))

write.table(export_features,
            file.path(TABLES_DIR, "dnmt3a_feature_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_feature_matrix.tsv (%d rows)\n", nrow(export_features)))

# Table 2: Model comparison summary
model_comparison <- data.frame(
  model = names(models_list),
  n_genes = rep(nrow(model_data), length(models_list)),
  n_hyper_dmr = rep(sum(model_data$hyper_dmr), length(models_list)),
  n_predictors = sapply(models_list, function(m) length(coef(m)) - 1),
  aic = sapply(models_list, AIC),
  auc = sapply(roc_list, function(r) as.numeric(auc(r))),
  auc_ci_lower = sapply(ci_list, function(c) c[1]),
  auc_ci_upper = sapply(ci_list, function(c) c[3]),
  mcfadden_r2 = sapply(models_list, mcfadden),
  stringsAsFactors = FALSE,
  row.names = NULL
)

write.table(model_comparison,
            file.path(TABLES_DIR, "dnmt3a_model_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_model_comparison.tsv (%d rows)\n", nrow(model_comparison)))

# Table 3: Coefficient table (raw + standardized, all models)
build_coef_table <- function(model, model_z, model_name, predictor_map) {
  raw_or <- extract_or(model)
  z_coefs <- coef(model_z)
  z_se <- summary(model_z)$coefficients[, 2]

  raw_or$model <- model_name
  raw_or$standardized_beta <- NA_real_
  raw_or$standardized_se <- NA_real_

  for (i in seq_len(nrow(raw_or))) {
    z_name <- paste0(raw_or$term[i], "_z")
    if (z_name %in% names(z_coefs)) {
      raw_or$standardized_beta[i] <- z_coefs[z_name]
      raw_or$standardized_se[i] <- z_se[z_name]
    } else if (raw_or$term[i] == "(Intercept)" && "(Intercept)" %in% names(z_coefs)) {
      raw_or$standardized_beta[i] <- z_coefs["(Intercept)"]
      raw_or$standardized_se[i] <- z_se["(Intercept)"]
    }
  }

  return(raw_or)
}

coef_table <- rbind(
  build_coef_table(model_full, model_full_z, "Full", predictor_cols),
  build_coef_table(model_dnmt3a, model_dnmt3a_z, "DNMT3A recruitment", predictor_cols),
  build_coef_table(model_tet, model_tet_z, "TET impediment", predictor_cols),
  build_coef_table(model_k119ub, model_k119ub_z, "K119ub only", predictor_cols),
  build_coef_table(model_step, model_step_z, "Stepwise", predictor_cols)
)

write.table(coef_table,
            file.path(TABLES_DIR, "dnmt3a_model_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_model_coefficients.tsv (%d rows)\n", nrow(coef_table)))

# Table 4: Feature importance comparison (LR + RF side-by-side)
importance_table <- data.frame(
  feature = predictor_cols,
  display_name = display_names[predictor_cols],
  stringsAsFactors = FALSE
)

# LR standardized betas from full model
lr_beta_vals <- full_z_coefs[-1]
importance_table$lr_standardized_beta <- lr_beta_vals[paste0(importance_table$feature, "_z")]
importance_table$lr_abs_beta <- abs(importance_table$lr_standardized_beta)
importance_table$lr_rank <- rank(-importance_table$lr_abs_beta)

# RF importance
rf_imp_ordered <- rf_importance[match(importance_table$feature, rf_importance$feature), ]
importance_table$rf_mean_decrease_accuracy <- rf_imp_ordered$MeanDecreaseAccuracy
importance_table$rf_mean_decrease_gini <- rf_imp_ordered$MeanDecreaseGini
importance_table$rf_rank <- rank(-importance_table$rf_mean_decrease_accuracy)

# Rank correlation between LR and RF
rank_cor <- cor(importance_table$lr_rank, importance_table$rf_rank, method = "spearman")
cat(sprintf("  LR vs RF rank correlation (Spearman): %.3f\n", rank_cor))

write.table(importance_table,
            file.path(TABLES_DIR, "dnmt3a_feature_importance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_feature_importance.tsv (%d rows)\n", nrow(importance_table)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 24 SUMMARY: DNMT3A BINDING PREDICTION\n")
cat("========================================================================\n\n")

cat(sprintf("Feature matrix: %d genes with all 7 features\n", nrow(model_data)))
cat(sprintf("Outcome: %d hyper-DMR (%.1f%%), %d non-hyper\n",
            sum(model_data$hyper_dmr), 100 * mean(model_data$hyper_dmr),
            sum(!model_data$hyper_dmr)))

cat(sprintf("\nLogistic Regression Model Comparison:\n"))
for (nm in names(models_list)) {
  mod <- models_list[[nm]]
  roc_obj <- roc_list[[nm]]
  ci_obj <- ci_list[[nm]]
  cat(sprintf("  %-22s AIC=%8.1f  AUC=%.3f [%.3f, %.3f]  R\u00B2=%.4f\n",
              paste0(nm, ":"), AIC(mod), auc(roc_obj), ci_obj[1], ci_obj[3],
              mcfadden(mod)))
}

cat(sprintf("\nDeLong test (DNMT3A recruitment vs TET impediment): %s\n",
            fmt_p(delong_result$p.value)))

cat(sprintf("\nTop 3 LR standardized predictors (full model):\n"))
lr_sorted <- sort(abs(full_z_coefs[-1]), decreasing = TRUE)
for (i in 1:min(3, length(lr_sorted))) {
  term <- gsub("_z$", "", names(lr_sorted)[i])
  cat(sprintf("  %d. %-20s |beta|=%.3f (beta=%.3f)\n",
              i, display_names[term],
              lr_sorted[i], full_z_coefs[names(lr_sorted)[i]]))
}

cat(sprintf("\nTop 3 RF predictors (MeanDecreaseAccuracy):\n"))
for (i in 1:min(3, nrow(rf_importance))) {
  cat(sprintf("  %d. %-20s Accuracy=%.2f\n",
              i, display_names[rf_importance$feature[i]],
              rf_importance$MeanDecreaseAccuracy[i]))
}

cat(sprintf("\nLR vs RF rank correlation: %.3f\n", rank_cor))

# Biological interpretation
cat(sprintf("\nInterpretation:\n"))
dnmt3a_auc <- as.numeric(auc(roc_dnmt3a))
tet_auc <- as.numeric(auc(roc_tet))

if (dnmt3a_auc > tet_auc & delong_result$p.value < 0.05) {
  cat("  DNMT3A recruitment model SIGNIFICANTLY outperforms TET impediment.\n")
  cat("  Supports DIRECT DNMT3A-UDR recruitment as primary mechanism.\n")
} else if (tet_auc > dnmt3a_auc & delong_result$p.value < 0.05) {
  cat("  TET impediment model SIGNIFICANTLY outperforms DNMT3A recruitment.\n")
  cat("  Supports TET blockade as primary mechanism for hypermethylation.\n")
} else if (delong_result$p.value >= 0.05) {
  cat("  No significant difference between DNMT3A recruitment and TET impediment models.\n")
  if (abs(dnmt3a_auc - tet_auc) < 0.01) {
    cat("  Both mechanisms contribute comparably — dual mechanism supported.\n")
  } else {
    better <- ifelse(dnmt3a_auc > tet_auc, "DNMT3A recruitment", "TET impediment")
    cat(sprintf("  %s trends better but not significantly (may need more power).\n", better))
  }
}

# K119ub interpretation
k119ub_beta <- full_z_coefs["k119ub_z"]
if (!is.na(k119ub_beta) && k119ub_beta > 0) {
  cat(sprintf("  K119ub is a POSITIVE predictor (beta=%.3f): genes with more K119ub\n", k119ub_beta))
  cat("  are more likely to gain methylation — consistent with DNMT3A-UDR recruitment.\n")
} else if (!is.na(k119ub_beta)) {
  cat(sprintf("  K119ub is a NEGATIVE predictor (beta=%.3f): unexpected direction.\n", k119ub_beta))
  cat("  May reflect indirect effects or confounding with gene expression.\n")
}

cat("\n--- Output files ---\n")
cat(sprintf("  Figures: %s/24{a-f}_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/dnmt3a_*.tsv (4 files)\n", TABLES_DIR))

cat("\n=== Section 24 complete ===\n")
