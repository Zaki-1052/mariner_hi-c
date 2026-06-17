# biomodal/downstream/scripts/viz_sections/section_24_dnmt3a_prediction.R
# Section 24: DNMT3A Binding Prediction — Dual-Mechanism Hypothesis Testing
# Standalone script - sources shared config for all dependencies and data
#
# Tests whether gene body hypermethylation (mC-up DMRs) in BAP1-KO can be
# computationally predicted from features related to two competing mechanisms:
# (1) DNMT3A-UDR recruitment (Chen et al. 2024) and (2) TET impediment.
#
# Result: TET impediment significantly outperforms DNMT3A recruitment
# (DeLong p < 2.2e-16). K119ub is a NEGATIVE predictor of hyper-DMR (OR < 1),
# opposite to the DNMT3A-UDR prediction. Baseline 5hmC (TET substrate) is the
# strongest predictor. An exclusive model comparison (removing shared atac_count)
# confirms this: TET-exclusive (baseline_hmc only) > DNMT3A-exclusive (k119ub +
# cpg_density).
#
# Dual-mechanism hypothesis tested:
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
#   24g: 10-fold cross-validated AUC comparison (boxplot + in-sample overlay)
#   24h: Stratified model comparison — all genes vs non-promoter subset
#   24i: K119ub x baseline 5hmC interaction (tertile-stratified sigmoid + forest plot)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_24_dnmt3a_prediction.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(pROC)
  library(ChIPseeker)
  library(readxl)
  library(randomForest)
  library(caret)
})

# =============================================================================
# SECTION 24 CONFIGURATION
# =============================================================================

# Feature extraction paths (from _shared_config.R)
HMC_EXTRACT_FILE <- EXTRACT_PATHS$hmc_regional_frac

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
# STEP 2b: EXCLUSIVE MECHANISTIC MODEL COMPARISON (no shared atac_count)
# =============================================================================

cat("\n--- Step 2b: Exclusive models (no shared features) ---\n\n")

# DNMT3A-exclusive: only K119ub-pathway features (no atac_count)
cat("  DNMT3A exclusive: hyper_dmr ~ k119ub + cpg_density\n")
model_dnmt3a_excl <- glm(
  hyper_dmr ~ k119ub + cpg_density,
  data = model_data, family = binomial
)

# TET-exclusive: only 5hmC substrate (no atac_count)
cat("  TET exclusive: hyper_dmr ~ baseline_hmc\n")
model_tet_excl <- glm(
  hyper_dmr ~ baseline_hmc,
  data = model_data, family = binomial
)

# Z-scored exclusive models
model_dnmt3a_excl_z <- glm(
  hyper_dmr ~ k119ub_z + cpg_density_z,
  data = model_data_z, family = binomial
)
model_tet_excl_z <- glm(
  hyper_dmr ~ baseline_hmc_z,
  data = model_data_z, family = binomial
)

# ROC/AUC
roc_dnmt3a_excl <- roc(model_data$hyper_dmr,
                        predict(model_dnmt3a_excl, type = "response"), quiet = TRUE)
roc_tet_excl    <- roc(model_data$hyper_dmr,
                        predict(model_tet_excl, type = "response"), quiet = TRUE)
ci_dnmt3a_excl  <- ci.auc(roc_dnmt3a_excl, quiet = TRUE)
ci_tet_excl     <- ci.auc(roc_tet_excl, quiet = TRUE)

# DeLong test: exclusive models
delong_excl <- roc.test(roc_dnmt3a_excl, roc_tet_excl, method = "delong")

cat(sprintf("\n  Exclusive model comparison:\n"))
cat(sprintf("    %-25s AIC=%8.1f  AUC=%.3f [%.3f, %.3f]  R²=%.4f\n",
            "DNMT3A excl (k119+cpg):", AIC(model_dnmt3a_excl),
            auc(roc_dnmt3a_excl), ci_dnmt3a_excl[1], ci_dnmt3a_excl[3],
            mcfadden(model_dnmt3a_excl)))
cat(sprintf("    %-25s AIC=%8.1f  AUC=%.3f [%.3f, %.3f]  R²=%.4f\n",
            "TET excl (baseline_hmc):", AIC(model_tet_excl),
            auc(roc_tet_excl), ci_tet_excl[1], ci_tet_excl[3],
            mcfadden(model_tet_excl)))
cat(sprintf("    DeLong (exclusive): AUC diff = %.4f, %s\n",
            as.numeric(auc(roc_dnmt3a_excl)) - as.numeric(auc(roc_tet_excl)),
            fmt_p(delong_excl$p.value)))
cat(sprintf("    (Shared models: DNMT3A=%.3f vs TET=%.3f, %s)\n",
            auc(roc_dnmt3a), auc(roc_tet), fmt_p(delong_result$p.value)))

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
  annotate("text", x = 0.05, y = 0.45, hjust = 0, vjust = 1, size = 3.2,
           label = annot_text) +
  labs(
    title = "ROC Curves: Predicting mC Hypermethylation (hyper-DMR)",
    subtitle = "DNMT3A recruitment vs TET impediment vs full model",
    x = "1 \u2212 Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    caption = paste0("Full: K119ub + ATAC + CpG dens. + 5mC + 5hmC + length + expr | ",
                     "DNMT3A: K119ub + ATAC + CpG dens. | TET: 5hmC + ATAC | ",
                     "K119ub: K119ub only | Step: AIC-selected")
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
delong_label <- sprintf("DeLong %s\n(DNMT3A vs TET)", fmt_p(delong_result$p.value))

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
    caption = paste0("Full: K119ub + ATAC + CpG dens. + 5mC + 5hmC + length + expr | ",
                     "DNMT3A: K119ub + ATAC + CpG dens. | TET: 5hmC + ATAC | K119ub: K119ub only | Step: AIC-selected"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 8.5, color = "grey40")
    )
  )

save_multiformat_ggplot(p_24d, file.path(OUTPUT_DIR, "24d_model_comparison"), 14, 8)

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
    subtitle = ifelse(or_k119ub_info$or[2] > 1,
      "Logistic regression: higher K119ub \u2192 higher hypermethylation probability (DNMT3A-UDR model)",
      "Logistic regression: higher K119ub \u2192 lower hypermethylation probability (inverse of DNMT3A-UDR prediction)"),
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
# STEP 6: 10-FOLD CROSS-VALIDATION
# =============================================================================

cat("\n--- Step 6: 10-fold cross-validation ---\n")

set.seed(42)

# Stratified folds preserving hyper_dmr class balance
cv_folds <- createFolds(
  factor(model_data$hyper_dmr),
  k = 10, list = TRUE, returnTrain = FALSE
)

# Model formulas (reuse from Step 2)
cv_formulas <- list(
  Full = hyper_dmr ~ k119ub + atac_count + cpg_density +
    baseline_mc + baseline_hmc + log_gene_length + log_expression,
  `DNMT3A recruitment` = hyper_dmr ~ k119ub + atac_count + cpg_density,
  `TET impediment` = hyper_dmr ~ baseline_hmc + atac_count,
  `K119ub only` = hyper_dmr ~ k119ub,
  Stepwise = formula(model_step)
)

# Cross-validate each model
cv_results_list <- list()
for (model_name in names(cv_formulas)) {
  fold_aucs <- numeric(10)
  for (fold_i in seq_along(cv_folds)) {
    test_idx <- cv_folds[[fold_i]]
    train_data <- model_data[-test_idx, ]
    test_data <- model_data[test_idx, ]

    fit <- glm(cv_formulas[[model_name]],
               data = train_data, family = binomial)
    preds <- predict(fit, newdata = test_data, type = "response")
    fold_aucs[fold_i] <- tryCatch(
      as.numeric(auc(roc(test_data$hyper_dmr, preds, quiet = TRUE))),
      error = function(e) NA_real_
    )
  }
  cv_results_list[[model_name]] <- fold_aucs
}

# Build results data frame
cv_long <- do.call(rbind, lapply(names(cv_results_list), function(nm) {
  data.frame(
    model = nm,
    fold = seq_along(cv_results_list[[nm]]),
    auc = cv_results_list[[nm]],
    stringsAsFactors = FALSE
  )
}))
cv_long$model <- factor(cv_long$model, levels = names(MODEL_COLORS))

# Summary statistics
cv_summary <- cv_long %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    cv_mean_auc = mean(auc, na.rm = TRUE),
    cv_sd_auc = sd(auc, na.rm = TRUE),
    cv_min_auc = min(auc, na.rm = TRUE),
    cv_max_auc = max(auc, na.rm = TRUE),
    .groups = "drop"
  )

# Add in-sample AUC for comparison
cv_summary$in_sample_auc <- sapply(
  as.character(cv_summary$model),
  function(nm) as.numeric(auc(roc_list[[nm]]))
)
cv_summary$optimism <- cv_summary$in_sample_auc - cv_summary$cv_mean_auc

cat("  10-fold CV results:\n")
for (i in seq_len(nrow(cv_summary))) {
  cat(sprintf("    %-22s CV AUC=%.3f +/- %.3f  (in-sample=%.3f, optimism=%.3f)\n",
              paste0(cv_summary$model[i], ":"),
              cv_summary$cv_mean_auc[i], cv_summary$cv_sd_auc[i],
              cv_summary$in_sample_auc[i], cv_summary$optimism[i]))
}

# --- Figure 24g: CV AUC boxplot ---
cat("\n--- Figure 24g: Cross-validated AUC comparison ---\n")

# In-sample AUC overlay points
insample_df <- data.frame(
  model = factor(names(MODEL_COLORS), levels = names(MODEL_COLORS)),
  auc = sapply(roc_list, function(r) as.numeric(auc(r)))
)

# Annotation labels
cv_label_df <- cv_summary %>%
  dplyr::mutate(
    label = sprintf("%.3f +/- %.3f", cv_mean_auc, cv_sd_auc)
  )

p_24g <- ggplot(cv_long, aes(x = model, y = auc, fill = model)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = model), width = 0.15, size = 1.5, alpha = 0.7) +
  geom_point(data = insample_df, aes(x = model, y = auc),
             shape = 18, size = 4, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_text(data = cv_label_df,
            aes(x = model, y = cv_min_auc - 0.015, label = label),
            size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = MODEL_COLORS) +
  scale_color_manual(values = MODEL_COLORS) +
  labs(
    title = "10-Fold Cross-Validated AUC",
    subtitle = paste0(
      "Boxes = CV fold AUCs; black diamonds = in-sample AUC\n",
      sprintf("N=%s genes, 10 stratified folds",
              format(nrow(model_data), big.mark = ","))
    ),
    x = NULL, y = "AUC",
    caption = paste0("Full: K119ub + ATAC + CpG dens. + 5mC + 5hmC + length + expr | ",
                     "DNMT3A: K119ub + ATAC + CpG dens. | TET: 5hmC + ATAC | K119ub: K119ub only | Step: AIC-selected")
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

save_multiformat_ggplot(
  p_24g, file.path(OUTPUT_DIR, "24g_cv_auc_comparison"), 11, 8
)

# Export CV table
cv_export <- cv_summary %>%
  dplyr::select(model, in_sample_auc, cv_mean_auc, cv_sd_auc,
                cv_min_auc, cv_max_auc, optimism)
write.table(cv_export,
            file.path(TABLES_DIR, "dnmt3a_cv_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_cv_results.tsv (%d rows)\n", nrow(cv_export)))

# =============================================================================
# STEP 7: STRATIFIED ANALYSIS — NON-PROMOTER SUBSET
# =============================================================================

cat("\n--- Step 7: Stratified analysis (non-promoter subset) ---\n")

nonprom <- model_data %>%
  dplyr::filter(chromatin_state != "Active_Promoter")

cat(sprintf("  Non-promoter subset: %d genes (%d hyper-DMR, %.1f%%)\n",
            nrow(nonprom), sum(nonprom$hyper_dmr),
            100 * mean(nonprom$hyper_dmr)))

# Refit all 5 models on non-promoter subset
np_null_ll <- as.numeric(logLik(
  glm(hyper_dmr ~ 1, data = nonprom, family = binomial)
))
np_mcfadden <- function(model) {
  1 - as.numeric(logLik(model)) / np_null_ll
}

np_models <- list()
np_rocs <- list()
np_cis <- list()
for (nm in names(cv_formulas)) {
  np_models[[nm]] <- glm(cv_formulas[[nm]],
                          data = nonprom, family = binomial)
  np_rocs[[nm]] <- roc(nonprom$hyper_dmr,
                        predict(np_models[[nm]], type = "response"),
                        quiet = TRUE)
  np_cis[[nm]] <- ci.auc(np_rocs[[nm]], quiet = TRUE)
}

# DeLong test on non-promoter subset
np_delong <- roc.test(
  np_rocs[["DNMT3A recruitment"]],
  np_rocs[["TET impediment"]],
  method = "delong"
)

cat("  Non-promoter model comparison:\n")
for (nm in names(np_models)) {
  cat(sprintf("    %-22s AUC=%.3f [%.3f, %.3f]  AIC=%.1f  R2=%.4f\n",
              paste0(nm, ":"),
              auc(np_rocs[[nm]]), np_cis[[nm]][1], np_cis[[nm]][3],
              AIC(np_models[[nm]]), np_mcfadden(np_models[[nm]])))
}
cat(sprintf("  DeLong (non-promoter, DNMT3A vs TET): %s\n",
            fmt_p(np_delong$p.value)))

# Non-promoter K119ub coefficient direction
np_full_z_data <- nonprom %>%
  dplyr::mutate(dplyr::across(
    all_of(predictor_cols),
    ~ as.numeric(scale(.x)),
    .names = "{.col}_z"
  ))
np_full_z <- glm(
  hyper_dmr ~ k119ub_z + atac_count_z + cpg_density_z +
    baseline_mc_z + baseline_hmc_z + log_gene_length_z + log_expression_z,
  data = np_full_z_data, family = binomial
)
np_k119ub_beta <- coef(np_full_z)["k119ub_z"]
cat(sprintf("  Non-promoter K119ub standardized beta: %.3f (vs All: %.3f)\n",
            np_k119ub_beta, full_z_coefs["k119ub_z"]))

# --- Figure 24h: Stratified model comparison ---
cat("\n--- Figure 24h: Stratified model comparison ---\n")

# Build combined comparison data
strat_df <- rbind(
  data.frame(
    model = factor(names(MODEL_COLORS), levels = names(MODEL_COLORS)),
    stratum = "All genes",
    auc = sapply(roc_list, function(r) as.numeric(auc(r))),
    auc_lower = sapply(ci_list, function(c) c[1]),
    auc_upper = sapply(ci_list, function(c) c[3]),
    stringsAsFactors = FALSE, row.names = NULL
  ),
  data.frame(
    model = factor(names(MODEL_COLORS), levels = names(MODEL_COLORS)),
    stratum = "Non-promoter",
    auc = sapply(np_rocs, function(r) as.numeric(auc(r))),
    auc_lower = sapply(np_cis, function(c) c[1]),
    auc_upper = sapply(np_cis, function(c) c[3]),
    stringsAsFactors = FALSE, row.names = NULL
  )
)
strat_df$stratum <- factor(strat_df$stratum,
                           levels = c("All genes", "Non-promoter"))

# Left panel: AUC comparison
p_24h_left <- ggplot(strat_df,
                     aes(x = model, y = auc,
                         color = stratum, group = stratum)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = auc_lower, ymax = auc_upper),
                width = 0.2, linewidth = 0.7,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("All genes" = "#333333", "Non-promoter" = "#E66100"),
    name = "Subset"
  ) +
  labs(title = "AUC by Subset", x = NULL, y = "AUC (95% CI)") +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        legend.position = "bottom")

# Right panel: top feature importance comparison
np_z_coefs <- coef(np_full_z)[-1]
all_z_coefs <- full_z_coefs[-1]

importance_cmp <- data.frame(
  feature = gsub("_z$", "", names(all_z_coefs)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    display_name = display_names[feature],
    all_beta = as.numeric(all_z_coefs),
    np_beta = as.numeric(np_z_coefs[paste0(feature, "_z")])
  ) %>%
  tidyr::pivot_longer(
    cols = c(all_beta, np_beta),
    names_to = "stratum", values_to = "beta"
  ) %>%
  dplyr::mutate(
    stratum = ifelse(stratum == "all_beta", "All genes", "Non-promoter"),
    stratum = factor(stratum, levels = c("All genes", "Non-promoter"))
  )

p_24h_right <- ggplot(importance_cmp,
                      aes(x = beta, y = reorder(display_name, abs(beta)),
                          fill = stratum)) +
  geom_col(position = position_dodge(width = 0.7),
           alpha = 0.8, width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(
    values = c("All genes" = "#333333", "Non-promoter" = "#E66100"),
    name = "Subset"
  ) +
  labs(title = "Standardized Betas (Full Model)",
       x = "Standardized beta", y = NULL) +
  theme_biomodal() +
  theme(legend.position = "bottom")

np_delong_label <- sprintf(
  "Non-promoter DeLong: %s\nAll genes: N=%s, Non-promoter: N=%s",
  fmt_p(np_delong$p.value),
  format(nrow(model_data), big.mark = ","),
  format(nrow(nonprom), big.mark = ",")
)

p_24h <- (p_24h_left | p_24h_right) +
  plot_annotation(
    title = "Model Performance: All Genes vs Non-Promoter Subset",
    subtitle = np_delong_label,
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

save_multiformat_ggplot(
  p_24h, file.path(OUTPUT_DIR, "24h_stratified_comparison"), 15, 9
)

# Export stratified comparison table
strat_export <- data.frame(
  model = rep(names(models_list), 2),
  stratum = rep(c("All genes", "Non-promoter"), each = 5),
  n_genes = c(rep(nrow(model_data), 5), rep(nrow(nonprom), 5)),
  n_hyper = c(rep(sum(model_data$hyper_dmr), 5),
              rep(sum(nonprom$hyper_dmr), 5)),
  auc = c(sapply(roc_list, function(r) as.numeric(auc(r))),
          sapply(np_rocs, function(r) as.numeric(auc(r)))),
  auc_ci_lower = c(sapply(ci_list, function(c) c[1]),
                   sapply(np_cis, function(c) c[1])),
  auc_ci_upper = c(sapply(ci_list, function(c) c[3]),
                   sapply(np_cis, function(c) c[3])),
  aic = c(sapply(models_list, AIC),
          sapply(np_models, AIC)),
  mcfadden_r2 = c(sapply(models_list, mcfadden),
                  sapply(np_models, np_mcfadden)),
  stringsAsFactors = FALSE, row.names = NULL
)

write.table(strat_export,
            file.path(TABLES_DIR, "dnmt3a_stratified_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_stratified_comparison.tsv (%d rows)\n",
            nrow(strat_export)))

# =============================================================================
# STEP 8: K119ub x BASELINE 5hmC INTERACTION
# =============================================================================

cat("\n--- Step 8: K119ub x baseline 5hmC interaction ---\n")

# Fit interaction model
model_interact <- glm(
  hyper_dmr ~ k119ub * baseline_hmc + atac_count + cpg_density +
    baseline_mc + log_gene_length + log_expression,
  data = model_data, family = binomial
)

# Standardized interaction model
model_interact_z <- glm(
  hyper_dmr ~ k119ub_z * baseline_hmc_z + atac_count_z + cpg_density_z +
    baseline_mc_z + log_gene_length_z + log_expression_z,
  data = model_data_z, family = binomial
)

# Interaction term statistics
interact_coefs <- summary(model_interact)$coefficients
interact_z_coefs <- summary(model_interact_z)$coefficients

interact_term <- "k119ub:baseline_hmc"
interact_z_term <- "k119ub_z:baseline_hmc_z"

interact_or <- exp(coef(model_interact)[interact_term])
interact_ci <- exp(confint.default(model_interact)[interact_term, ])
interact_p <- interact_coefs[interact_term, 4]

interact_z_beta <- coef(model_interact_z)[interact_z_term]
interact_z_p <- interact_z_coefs[interact_z_term, 4]

cat(sprintf("  Interaction term (k119ub:baseline_hmc):\n"))
cat(sprintf("    Raw OR = %.3f [%.3f, %.3f], %s\n",
            interact_or, interact_ci[1], interact_ci[2],
            fmt_p(interact_p)))
cat(sprintf("    Standardized beta = %.3f, %s\n",
            interact_z_beta, fmt_p(interact_z_p)))

# ROC/AUC for interaction model
roc_interact <- roc(model_data$hyper_dmr,
                    predict(model_interact, type = "response"),
                    quiet = TRUE)
ci_interact <- ci.auc(roc_interact, quiet = TRUE)
cat(sprintf("    AUC = %.3f [%.3f, %.3f]  AIC = %.1f\n",
            auc(roc_interact), ci_interact[1], ci_interact[3],
            AIC(model_interact)))

# LR test: interaction model vs full model (nested)
lr_interact <- -2 * (as.numeric(logLik(model_full)) -
                      as.numeric(logLik(model_interact)))
lr_interact_p <- pchisq(lr_interact, df = 1, lower.tail = FALSE)
cat(sprintf("    LR test vs Full model: chi2=%.2f, df=1, %s\n",
            lr_interact, fmt_p(lr_interact_p)))

# Tertile-stratified K119ub effect
model_data$hmc_tertile <- cut(
  model_data$baseline_hmc,
  breaks = quantile(model_data$baseline_hmc, probs = c(0, 1/3, 2/3, 1)),
  labels = c("Low 5hmC", "Mid 5hmC", "High 5hmC"),
  include.lowest = TRUE
)

tertile_results <- list()
for (tert in levels(model_data$hmc_tertile)) {
  tert_data <- model_data %>% dplyr::filter(hmc_tertile == tert)
  tert_fit <- glm(hyper_dmr ~ k119ub, data = tert_data, family = binomial)
  tert_or_info <- extract_or(tert_fit)
  tert_roc <- roc(tert_data$hyper_dmr,
                  predict(tert_fit, type = "response"), quiet = TRUE)
  tertile_results[[tert]] <- data.frame(
    tertile = tert,
    n = nrow(tert_data),
    n_hyper = sum(tert_data$hyper_dmr),
    k119ub_or = tert_or_info$or[2],
    k119ub_or_lower = tert_or_info$or_lower[2],
    k119ub_or_upper = tert_or_info$or_upper[2],
    k119ub_p = tert_or_info$p_value[2],
    auc = as.numeric(auc(tert_roc)),
    stringsAsFactors = FALSE
  )
}
tertile_df <- do.call(rbind, tertile_results)
rownames(tertile_df) <- NULL

cat("\n  Per-tertile K119ub effect:\n")
for (i in seq_len(nrow(tertile_df))) {
  cat(sprintf("    %-12s OR=%.3f [%.3f, %.3f] %s (n=%d, AUC=%.3f)\n",
              paste0(tertile_df$tertile[i], ":"),
              tertile_df$k119ub_or[i],
              tertile_df$k119ub_or_lower[i],
              tertile_df$k119ub_or_upper[i],
              fmt_p(tertile_df$k119ub_p[i]),
              tertile_df$n[i],
              tertile_df$auc[i]))
}

# --- Figure 24i: Interaction visualization ---
cat("\n--- Figure 24i: K119ub x baseline 5hmC interaction ---\n")

# Left panel: predicted probability curves by 5hmC tertile
# Generate prediction grid
k119ub_range <- seq(
  min(model_data$k119ub), max(model_data$k119ub), length.out = 200
)

# Median values for other predictors
median_vals <- model_data %>%
  dplyr::summarise(dplyr::across(
    c(atac_count, cpg_density, baseline_mc,
      log_gene_length, log_expression),
    median
  ))

# Tertile median 5hmC values
hmc_medians <- model_data %>%
  dplyr::group_by(hmc_tertile) %>%
  dplyr::summarise(baseline_hmc = median(baseline_hmc), .groups = "drop")

pred_grid <- expand.grid(
  k119ub = k119ub_range,
  hmc_tertile = levels(model_data$hmc_tertile),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(hmc_medians, by = "hmc_tertile") %>%
  dplyr::mutate(
    atac_count = median_vals$atac_count,
    cpg_density = median_vals$cpg_density,
    baseline_mc = median_vals$baseline_mc,
    log_gene_length = median_vals$log_gene_length,
    log_expression = median_vals$log_expression
  )

pred_grid$pred_prob <- predict(model_interact,
                               newdata = pred_grid, type = "response")

tertile_colors <- c(
  "Low 5hmC" = "#2C7BB6", "Mid 5hmC" = "#ABD9E9", "High 5hmC" = "#D7191C"
)

interact_annot <- sprintf(
  "Interaction: %s\nbeta(z) = %.3f",
  fmt_p(interact_p), interact_z_beta
)

p_24i_left <- ggplot(pred_grid,
                     aes(x = k119ub, y = pred_prob,
                         color = hmc_tertile)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = tertile_colors,
                     name = "Baseline 5hmC Tertile") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           size = 3.5, label = interact_annot) +
  labs(
    title = "Predicted P(hyper-DMR) by K119ub",
    subtitle = "Interaction model, at tertile-median 5hmC",
    x = "WT H2AK119ub Gene Body Signal",
    y = "P(mC Hyper-DMR)"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

# Right panel: forest plot of per-tertile K119ub OR
tertile_df$tertile <- factor(tertile_df$tertile,
                             levels = rev(levels(model_data$hmc_tertile)))
tertile_df$or_label <- sapply(seq_len(nrow(tertile_df)), function(i) {
  sprintf("OR=%.2f\n%s", tertile_df$k119ub_or[i], fmt_p(tertile_df$k119ub_p[i]))
})

p_24i_right <- ggplot(tertile_df,
                      aes(x = k119ub_or, y = tertile)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = k119ub_or_lower, xmax = k119ub_or_upper),
                width = 0.2, linewidth = 0.8, orientation = "y") +
  geom_point(size = 4, aes(color = tertile)) +
  geom_text(aes(label = or_label, x = k119ub_or_upper),
            hjust = -0.15, vjust = 0.5, size = 2.8) +
  scale_color_manual(values = tertile_colors) +
  scale_x_log10() +
  labs(
    title = "K119ub OR by 5hmC Tertile",
    subtitle = "Univariate K119ub effect within each tertile",
    x = "Odds Ratio (log scale)", y = NULL
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

p_24i <- (p_24i_left | p_24i_right) +
  plot_annotation(
    title = "K119ub Effect Modulated by Baseline 5hmC Level",
    subtitle = sprintf(
      "Interaction model AUC=%.3f (Full=%.3f); LR test %s",
      auc(roc_interact), auc(roc_full), fmt_p(lr_interact_p)
    ),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

save_multiformat_ggplot(
  p_24i, file.path(OUTPUT_DIR, "24i_interaction_k119ub_hmc"), 15, 8
)

# =============================================================================
# FIGURE 24j: EXCLUSIVE MODEL COMPARISON (no shared atac_count)
# =============================================================================

cat("\n--- Figure 24j: Exclusive mechanistic model comparison ---\n")

# Build ROC data for exclusive + shared pairs
excl_roc_df <- rbind(
  roc_to_df(roc_dnmt3a, "DNMT3A (shared)"),
  roc_to_df(roc_tet, "TET (shared)"),
  roc_to_df(roc_dnmt3a_excl, "DNMT3A (exclusive)"),
  roc_to_df(roc_tet_excl, "TET (exclusive)")
)

EXCL_COLORS <- c(
  "DNMT3A (shared)"    = "#756BB1",
  "TET (shared)"       = "#377EB8",
  "DNMT3A (exclusive)" = "#BCBDDC",
  "TET (exclusive)"    = "#9ECAE1"
)
excl_roc_df$model <- factor(excl_roc_df$model, levels = names(EXCL_COLORS))

excl_legend <- sprintf(
  "%s (AUC=%.3f)",
  names(EXCL_COLORS),
  c(auc(roc_dnmt3a), auc(roc_tet), auc(roc_dnmt3a_excl), auc(roc_tet_excl))
)

# Left panel: ROC overlay
p_24j_left <- ggplot(excl_roc_df, aes(x = 1 - specificity, y = sensitivity, color = model)) +
  geom_line(aes(linetype = ifelse(grepl("exclusive", model), "dashed", "solid")),
            linewidth = 1.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = EXCL_COLORS, labels = excl_legend, name = NULL) +
  scale_linetype_identity() +
  annotate("text", x = 0.05, y = 0.45, hjust = 0, vjust = 1, size = 3.2,
           label = sprintf("Shared DeLong: %s\nExclusive DeLong: %s\nN=%s genes",
                           fmt_p(delong_result$p.value),
                           fmt_p(delong_excl$p.value),
                           format(nrow(model_data), big.mark = ","))) +
  labs(
    title = "ROC: Shared vs Exclusive Features",
    x = "1 − Specificity", y = "Sensitivity"
  ) +
  theme_biomodal() +
  theme(legend.position = c(0.65, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 8)) +
  coord_equal()

# Right panel: AUC point + CI comparison
excl_comparison <- data.frame(
  model = factor(names(EXCL_COLORS), levels = names(EXCL_COLORS)),
  auc = c(as.numeric(auc(roc_dnmt3a)), as.numeric(auc(roc_tet)),
          as.numeric(auc(roc_dnmt3a_excl)), as.numeric(auc(roc_tet_excl))),
  auc_lower = c(ci_dnmt3a[1], ci_tet[1], ci_dnmt3a_excl[1], ci_tet_excl[1]),
  auc_upper = c(ci_dnmt3a[3], ci_tet[3], ci_dnmt3a_excl[3], ci_tet_excl[3]),
  type = c("Shared", "Shared", "Exclusive", "Exclusive"),
  stringsAsFactors = FALSE
)

p_24j_right <- ggplot(excl_comparison, aes(x = model, y = auc, color = model)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = auc_lower, ymax = auc_upper), width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -1.5, size = 3.2) +
  scale_color_manual(values = EXCL_COLORS) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  labs(title = "AUC Comparison", x = NULL, y = "AUC (95% CI)") +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

p_24j <- (p_24j_left | p_24j_right) +
  plot_annotation(
    title = "Exclusive Mechanistic Comparison (No Shared Features)",
    subtitle = sprintf(
      "Removing shared atac_count: DNMT3A (k119ub+cpg) vs TET (baseline_hmc); DeLong %s",
      fmt_p(delong_excl$p.value)),
    caption = paste0("Shared: DNMT3A = K119ub + ATAC + CpG dens. | TET = 5hmC + ATAC\n",
                     "Exclusive: DNMT3A = K119ub + CpG dens. (no ATAC) | TET = Baseline 5hmC only (no ATAC)"),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0.5, size = 9, color = "grey40")
    )
  )

save_multiformat_ggplot(p_24j, file.path(OUTPUT_DIR, "24j_exclusive_model_comparison"), 16, 9)

# Export exclusive model comparison table
excl_export <- data.frame(
  model = c("DNMT3A (shared)", "TET (shared)", "DNMT3A (exclusive)", "TET (exclusive)"),
  features = c("k119ub + atac_count + cpg_density", "baseline_hmc + atac_count",
               "k119ub + cpg_density", "baseline_hmc"),
  auc = excl_comparison$auc,
  auc_ci_lower = excl_comparison$auc_lower,
  auc_ci_upper = excl_comparison$auc_upper,
  aic = c(AIC(model_dnmt3a), AIC(model_tet), AIC(model_dnmt3a_excl), AIC(model_tet_excl)),
  mcfadden_r2 = c(mcfadden(model_dnmt3a), mcfadden(model_tet),
                  mcfadden(model_dnmt3a_excl), mcfadden(model_tet_excl)),
  delong_p = c(delong_result$p.value, delong_result$p.value,
               delong_excl$p.value, delong_excl$p.value),
  stringsAsFactors = FALSE
)

write.table(excl_export,
            file.path(TABLES_DIR, "dnmt3a_exclusive_model_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_exclusive_model_comparison.tsv (%d rows)\n", nrow(excl_export)))

# Export interaction results table
interact_export <- rbind(
  data.frame(
    term = "Interaction model summary",
    value = NA_character_,
    aic = AIC(model_interact),
    auc = as.numeric(auc(roc_interact)),
    lr_test_p = lr_interact_p,
    stringsAsFactors = FALSE
  ),
  data.frame(
    term = "k119ub:baseline_hmc (raw)",
    value = sprintf("OR=%.3f [%.3f, %.3f], %s",
                    interact_or, interact_ci[1], interact_ci[2],
                    fmt_p(interact_p)),
    aic = NA_real_, auc = NA_real_, lr_test_p = NA_real_,
    stringsAsFactors = FALSE
  ),
  data.frame(
    term = "k119ub_z:baseline_hmc_z (standardized)",
    value = sprintf("beta=%.3f, %s", interact_z_beta, fmt_p(interact_z_p)),
    aic = NA_real_, auc = NA_real_, lr_test_p = NA_real_,
    stringsAsFactors = FALSE
  )
)

# Append tertile rows
for (i in seq_len(nrow(tertile_df))) {
  interact_export <- rbind(interact_export, data.frame(
    term = sprintf("K119ub OR in %s tertile", tertile_df$tertile[i]),
    value = sprintf("OR=%.3f [%.3f, %.3f], %s, n=%d",
                    tertile_df$k119ub_or[i],
                    tertile_df$k119ub_or_lower[i],
                    tertile_df$k119ub_or_upper[i],
                    fmt_p(tertile_df$k119ub_p[i]),
                    tertile_df$n[i]),
    aic = NA_real_,
    auc = tertile_df$auc[i],
    lr_test_p = NA_real_,
    stringsAsFactors = FALSE
  ))
}

write.table(interact_export,
            file.path(TABLES_DIR, "dnmt3a_interaction_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved dnmt3a_interaction_results.tsv (%d rows)\n",
            nrow(interact_export)))

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

cat("\n--- Enhancement results ---\n")

cat(sprintf("\n10-Fold Cross-Validation:\n"))
for (i in seq_len(nrow(cv_summary))) {
  cat(sprintf("  %-22s CV=%.3f +/- %.3f (optimism=%.3f)\n",
              paste0(cv_summary$model[i], ":"),
              cv_summary$cv_mean_auc[i], cv_summary$cv_sd_auc[i],
              cv_summary$optimism[i]))
}

cat(sprintf("\nNon-Promoter Stratification:\n"))
cat(sprintf("  N=%d genes (excl. %d Active_Promoter)\n",
            nrow(nonprom), nrow(model_data) - nrow(nonprom)))
cat(sprintf("  K119ub beta: All=%.3f, Non-promoter=%.3f\n",
            full_z_coefs["k119ub_z"], np_k119ub_beta))
cat(sprintf("  DeLong (non-prom, DNMT3A vs TET): %s\n",
            fmt_p(np_delong$p.value)))

cat(sprintf("\nK119ub x 5hmC Interaction:\n"))
cat(sprintf("  Interaction OR=%.3f, %s\n", interact_or, fmt_p(interact_p)))
cat(sprintf("  Interaction beta(z)=%.3f, %s\n",
            interact_z_beta, fmt_p(interact_z_p)))
for (i in seq_len(nrow(tertile_df))) {
  cat(sprintf("  %s: K119ub OR=%.3f\n",
              tertile_df$tertile[i], tertile_df$k119ub_or[i]))
}

if (interact_p < 0.05 && interact_z_beta > 0) {
  cat("  SIGNIFICANT POSITIVE interaction: K119ub promotes hyper-DMR\n")
  cat("  specifically at high-5hmC genes. Supports dual mechanism.\n")
} else if (interact_p < 0.05 && interact_z_beta < 0) {
  cat("  SIGNIFICANT NEGATIVE interaction: K119ub effect is even more\n")
  cat("  negative at high-5hmC genes. Argues against dual mechanism.\n")
} else {
  cat("  Non-significant interaction: K119ub effect does not depend\n")
  cat("  on baseline 5hmC level. No evidence for dual mechanism.\n")
}

cat("\n--- Output files ---\n")
cat(sprintf("  Figures: %s/24{a-i}_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/dnmt3a_*.tsv (7 files)\n", TABLES_DIR))

cat("\n=== Section 24 complete ===\n")
