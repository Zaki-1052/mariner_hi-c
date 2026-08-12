# biomodal/downstream/scripts/viz_sections/section_62_multifeature_chromatin_regression.R
# Section 62: Multi-Feature Chromatin Regression at MeCP2 Peaks
#
# Core test of Jai's hypothesis: can chromatin features explain MeCP2 binding
# better than CG methylation alone? Extracts continuous BigWig signal for 7
# histone marks + CG mC/hmC at all 202,650 MeCP2 DiffBind peaks, then fits
# nested regression models comparing CG-only vs Chromatin-only vs Full.
#
# Figures:
#   62a: R² bar comparing binding-level + differential models
#   62b: Variance partition stacked bar (CG-unique | Shared | Chromatin-unique)
#   62c: LASSO coefficient path plot
#   62d: Standardized coefficient forest plot (full model)
#   62e: Chromatin-state stratified R² comparison
#   62f: Section 56 residual vs K27me3 scatter
#   62g: Composite
#
# Input:  peaks/mecp2/MeCP2_annotated.txt, 56_mecp2_peak_annotated.tsv,
#         HISTONE_BIGWIGS + METHYLATION_BIGWIGS from _shared_config.R
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_62_multifeature_chromatin_regression.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(patchwork)
})

# =============================================================================
# SECTION 62 CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 62: MULTI-FEATURE CHROMATIN REGRESSION AT MeCP2 PEAKS\n")
cat("  Can chromatin features explain MeCP2 binding better than methylation?\n")
cat("================================================================================\n\n")

SEC62_DIR <- file.path(OUTPUT_DIR, "62_multifeature_chromatin_regression")
dir.create(SEC62_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNAL_CACHE <- file.path(TABLES_DIR, "62_mecp2_peak_chromatin_signal.tsv")
SEC56_PEAKS  <- file.path(TABLES_DIR, "56_mecp2_peak_annotated.tsv")

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC62_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# BIGWIG SIGNAL EXTRACTION (checkpointed)
# =============================================================================

extract_signal_at_peaks <- function(bw_path, peak_gr, label) {
  cat(sprintf("  %s (%d peaks)...", label, length(peak_gr)))
  bw <- rtracklayer::import.bw(bw_path, which = peak_gr)
  cov <- coverage(bw, weight = "score")
  ranges_by_chr <- split(ranges(peak_gr), seqnames(peak_gr))
  shared_chrs <- intersect(names(cov), names(ranges_by_chr))
  v <- Views(cov[shared_chrs], ranges_by_chr[shared_chrs])
  means <- viewMeans(v)
  result <- rep(NA_real_, length(peak_gr))
  chr_names <- as.character(seqnames(peak_gr))
  for (chr in names(means)) {
    idx <- which(chr_names == chr)
    result[idx] <- as.numeric(means[[chr]])
  }
  nonzero <- sum(result > 0, na.rm = TRUE)
  cat(sprintf(" median=%.4f, non-zero=%d/%d\n",
              median(result, na.rm = TRUE), nonzero, length(result)))
  result
}

if (file.exists(SIGNAL_CACHE)) {
  cat("--- Loading cached signal from", SIGNAL_CACHE, "---\n")
  peaks <- read.table(SIGNAL_CACHE, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "")
  cat(sprintf("  Loaded %d peaks × %d columns\n", nrow(peaks), ncol(peaks)))

} else {
  cat("--- Extracting BigWig signal at MeCP2 peaks (this will take ~30 min) ---\n\n")

  peaks <- load_diffbind_flex(MECP2_FILES$annotated, "MeCP2")
  peak_gr <- GRanges(seqnames = peaks$Chr,
                     ranges = IRanges(start = peaks$Start, end = peaks$End))
  cat(sprintf("  %d MeCP2 DiffBind peaks loaded\n\n", length(peak_gr)))

  # MeCP2 signal
  cat("  MeCP2 BigWigs:\n")
  peaks$mecp2_ctrl <- extract_signal_at_peaks(MECP2_FILES$ctrl_bw, peak_gr, "MeCP2_ctrl")
  peaks$mecp2_mut  <- extract_signal_at_peaks(MECP2_FILES$mut_bw, peak_gr, "MeCP2_mut")

  # CG methylation
  cat("  CG methylation BigWigs:\n")
  peaks$cg_mc_ctrl  <- extract_signal_at_peaks(METHYLATION_BIGWIGS$cg_mc_ctrl, peak_gr, "CG_mC_ctrl")
  peaks$cg_mc_mut   <- extract_signal_at_peaks(METHYLATION_BIGWIGS$cg_mc_mut, peak_gr, "CG_mC_mut")
  peaks$cg_hmc_ctrl <- extract_signal_at_peaks(METHYLATION_BIGWIGS$cg_hmc_ctrl, peak_gr, "CG_hmC_ctrl")
  peaks$cg_hmc_mut  <- extract_signal_at_peaks(METHYLATION_BIGWIGS$cg_hmc_mut, peak_gr, "CG_hmC_mut")

  # Histone marks
  cat("  Histone mark BigWigs:\n")
  mark_names <- c("k27me3", "k27ac", "k119ub", "atac", "k4me3", "k36me3", "k27me1")
  for (mk in mark_names) {
    ctrl_key <- paste0(mk, "_ctrl")
    mut_key  <- paste0(mk, "_mut")
    peaks[[ctrl_key]] <- extract_signal_at_peaks(HISTONE_BIGWIGS[[ctrl_key]], peak_gr, paste0(mk, "_ctrl"))
    peaks[[mut_key]]  <- extract_signal_at_peaks(HISTONE_BIGWIGS[[mut_key]], peak_gr, paste0(mk, "_mut"))
  }

  # Compute log2FC for each mark
  cat("\n  Computing log2FC per mark...\n")
  all_marks <- c("mecp2", "cg_mc", "cg_hmc", mark_names)
  for (mk in all_marks) {
    ctrl_col <- paste0(mk, "_ctrl")
    mut_col  <- paste0(mk, "_mut")
    fc_col   <- paste0(mk, "_log2fc")

    ctrl_vals <- peaks[[ctrl_col]]
    mut_vals  <- peaks[[mut_col]]
    all_nonzero <- c(ctrl_vals[ctrl_vals > 0], mut_vals[mut_vals > 0])
    pseudo <- if (length(all_nonzero) > 0) {
      max(quantile(all_nonzero, 0.05, na.rm = TRUE), 1e-6)
    } else { 1e-6 }

    peaks[[fc_col]] <- log2((mut_vals + pseudo) / (ctrl_vals + pseudo))
    cat(sprintf("    %s: pseudo=%.6f, median_log2fc=%.4f\n",
                mk, pseudo, median(peaks[[fc_col]], na.rm = TRUE)))
  }

  cat(sprintf("\n  Saving signal cache: %s\n", SIGNAL_CACHE))
  write.table(peaks, SIGNAL_CACHE, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Cached %d peaks × %d columns\n", nrow(peaks), ncol(peaks)))
}

# =============================================================================
# PREPARE REGRESSION DATA
# =============================================================================

cat("\n--- Preparing regression data ---\n")

# Log2-transform ctrl signals for binding-level models
ctrl_marks <- c("mecp2_ctrl", "cg_mc_ctrl", "cg_hmc_ctrl",
                "k27me3_ctrl", "k27ac_ctrl", "k119ub_ctrl", "atac_ctrl",
                "k4me3_ctrl", "k36me3_ctrl", "k27me1_ctrl")

for (col in ctrl_marks) {
  log_col <- paste0(col, "_log")
  vals <- peaks[[col]]
  all_nonzero <- vals[!is.na(vals) & vals > 0]
  pseudo <- if (length(all_nonzero) > 0) {
    max(quantile(all_nonzero, 0.05, na.rm = TRUE), 1e-6)
  } else { 1e-6 }
  peaks[[log_col]] <- log2(vals + pseudo)
}

# Classify MeCP2 status
peaks$mecp2_class <- dplyr::case_when(
  peaks$FDR < 0.05 & peaks$Fold > 0 ~ "MeCP2 Up",
  peaks$FDR < 0.05 & peaks$Fold < 0 ~ "MeCP2 Down",
  TRUE ~ "Not Significant"
)

cat(sprintf("  MeCP2 Up: %d | Down: %d | NS: %d\n",
            sum(peaks$mecp2_class == "MeCP2 Up"),
            sum(peaks$mecp2_class == "MeCP2 Down"),
            sum(peaks$mecp2_class == "Not Significant")))

# =============================================================================
# 62a: NESTED REGRESSION MODELS — BINDING LEVEL (WHERE MeCP2 binds)
# =============================================================================

cat("\n--- 62a: Binding-level regression models ---\n")

binding_predictors_cg <- c("cg_mc_ctrl_log", "cg_hmc_ctrl_log")
binding_predictors_chr <- c("k27me3_ctrl_log", "k27ac_ctrl_log", "k119ub_ctrl_log",
                            "atac_ctrl_log", "k4me3_ctrl_log", "k36me3_ctrl_log",
                            "k27me1_ctrl_log")
binding_predictors_full <- c(binding_predictors_cg, binding_predictors_chr)

binding_valid <- complete.cases(peaks[, c("mecp2_ctrl_log", binding_predictors_full)])
cat(sprintf("  Complete cases for binding models: %d / %d peaks\n",
            sum(binding_valid), nrow(peaks)))

binding_df <- peaks[binding_valid, ]

bind_m1 <- lm(as.formula(paste("mecp2_ctrl_log ~",
              paste(binding_predictors_cg, collapse = " + "))), data = binding_df)
bind_m2 <- lm(as.formula(paste("mecp2_ctrl_log ~",
              paste(binding_predictors_chr, collapse = " + "))), data = binding_df)
bind_m3 <- lm(as.formula(paste("mecp2_ctrl_log ~",
              paste(binding_predictors_full, collapse = " + "))), data = binding_df)

# --- Differential models (WHERE MeCP2 changes) ---

cat("\n--- 62c: Differential regression models ---\n")

diff_predictors_cg <- c("cg_mc_log2fc")
diff_predictors_chr <- c("k27me3_log2fc", "k27ac_log2fc", "k119ub_log2fc",
                         "atac_log2fc", "k4me3_log2fc", "k36me3_log2fc",
                         "k27me1_log2fc")
diff_predictors_full <- c("cg_mc_log2fc", "cg_hmc_log2fc", diff_predictors_chr)

diff_valid <- complete.cases(peaks[, c("Fold", diff_predictors_full)])
cat(sprintf("  Complete cases for differential models: %d / %d peaks\n",
            sum(diff_valid), nrow(peaks)))

diff_df <- peaks[diff_valid, ]

diff_m1 <- lm(as.formula(paste("Fold ~",
              paste(diff_predictors_cg, collapse = " + "))), data = diff_df)
diff_m2 <- lm(as.formula(paste("Fold ~",
              paste(diff_predictors_chr, collapse = " + "))), data = diff_df)
diff_m3 <- lm(as.formula(paste("Fold ~",
              paste(diff_predictors_full, collapse = " + "))), data = diff_df)

# Collect model summaries
collect_model_summary <- function(model, name, type) {
  s <- summary(model)
  f_p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  data.frame(
    type = type, model = name, n = nrow(model$model),
    r_squared = s$r.squared, adj_r_squared = s$adj.r.squared,
    f_stat = s$fstatistic[1], f_p = f_p,
    aic = AIC(model), bic = BIC(model),
    stringsAsFactors = FALSE
  )
}

model_summary <- rbind(
  collect_model_summary(bind_m1, "CG only", "Binding level"),
  collect_model_summary(bind_m2, "Chromatin only", "Binding level"),
  collect_model_summary(bind_m3, "Full (CG + Chromatin)", "Binding level"),
  collect_model_summary(diff_m1, "CG only", "Differential"),
  collect_model_summary(diff_m2, "Chromatin only", "Differential"),
  collect_model_summary(diff_m3, "Full (CG + Chromatin)", "Differential")
)

cat("\n  Model comparison:\n")
for (i in seq_len(nrow(model_summary))) {
  r <- model_summary[i, ]
  cat(sprintf("    [%s] %-22s R²=%.4f  adj_R²=%.4f  AIC=%.0f\n",
              r$type, r$model, r$r_squared, r$adj_r_squared, r$aic))
}

# Plot 62a: R² bar
model_summary$label <- paste0(model_summary$type, "\n", model_summary$model)
model_summary$label <- factor(model_summary$label,
  levels = c("Binding level\nCG only", "Binding level\nChromatin only",
             "Binding level\nFull (CG + Chromatin)",
             "Differential\nCG only", "Differential\nChromatin only",
             "Differential\nFull (CG + Chromatin)"))

fill_colors <- rep(c("#D6604D", "#4393C3", "#762A83"), 2)

p_62a <- ggplot(model_summary, aes(x = label, y = r_squared, fill = label)) +
  geom_col(alpha = 0.85, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.4f", r_squared)), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "MeCP2 regression: CG methylation vs chromatin marks",
    subtitle = "Left: predicting WHERE MeCP2 binds | Right: predicting WHERE MeCP2 changes",
    x = NULL, y = expression(R^2)
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 9, lineheight = 0.9))

save_plot(p_62a, "62a_model_r2_comparison", w = 13, h = 7)

# =============================================================================
# 62b: VARIANCE PARTITIONING
# =============================================================================

cat("\n--- 62b: Variance partitioning ---\n")

partition_variance <- function(r2_cg, r2_chr, r2_full, label) {
  cg_unique  <- r2_full - r2_chr
  chr_unique <- r2_full - r2_cg
  shared     <- r2_cg + r2_chr - r2_full
  unexplained <- 1 - r2_full
  data.frame(
    type = label,
    component = c("CG unique", "Shared", "Chromatin unique", "Unexplained"),
    fraction = c(cg_unique, shared, chr_unique, unexplained),
    stringsAsFactors = FALSE
  )
}

vp <- rbind(
  partition_variance(summary(bind_m1)$r.squared, summary(bind_m2)$r.squared,
                     summary(bind_m3)$r.squared, "Binding level"),
  partition_variance(summary(diff_m1)$r.squared, summary(diff_m2)$r.squared,
                     summary(diff_m3)$r.squared, "Differential")
)

cat("  Variance partition:\n")
for (i in seq_len(nrow(vp))) {
  cat(sprintf("    [%s] %-18s %.4f (%.1f%%)\n",
              vp$type[i], vp$component[i], vp$fraction[i], 100 * vp$fraction[i]))
}

vp$component <- factor(vp$component,
  levels = c("Unexplained", "Chromatin unique", "Shared", "CG unique"))
vp_colors <- c("CG unique" = "#D6604D", "Shared" = "#FDB863",
               "Chromatin unique" = "#4393C3", "Unexplained" = "grey85")

p_62b <- ggplot(vp, aes(x = type, y = fraction, fill = component)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(fraction > 0.02,
                sprintf("%.1f%%", 100 * fraction), "")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = vp_colors, name = "Component") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Variance partitioning: CG methylation vs chromatin marks",
    subtitle = "How much MeCP2 variance is uniquely explained by each?",
    x = NULL, y = "Fraction of total variance"
  ) +
  theme_biomodal() +
  coord_flip()

save_plot(p_62b, "62b_variance_partition", w = 11, h = 5)

# =============================================================================
# 62c: LASSO VARIABLE SELECTION
# =============================================================================

cat("\n--- 62c: LASSO variable selection ---\n")

# Binding-level LASSO
x_bind <- as.matrix(binding_df[, binding_predictors_full])
y_bind <- binding_df$mecp2_ctrl_log

set.seed(42)
cv_bind <- cv.glmnet(x_bind, y_bind, alpha = 1, nfolds = 10)
lasso_bind <- glmnet(x_bind, y_bind, alpha = 1)

coef_1se <- as.matrix(coef(cv_bind, s = "lambda.1se"))
coef_min <- as.matrix(coef(cv_bind, s = "lambda.min"))

selected_1se <- rownames(coef_1se)[coef_1se[, 1] != 0 & rownames(coef_1se) != "(Intercept)"]
selected_min <- rownames(coef_min)[coef_min[, 1] != 0 & rownames(coef_min) != "(Intercept)"]

cat(sprintf("  Binding-level LASSO:\n"))
cat(sprintf("    lambda.1se: %d features selected: %s\n",
            length(selected_1se), paste(selected_1se, collapse = ", ")))
cat(sprintf("    lambda.min: %d features selected: %s\n",
            length(selected_min), paste(selected_min, collapse = ", ")))

# Differential LASSO
x_diff <- as.matrix(diff_df[, diff_predictors_full])
y_diff <- diff_df$Fold

set.seed(42)
cv_diff <- cv.glmnet(x_diff, y_diff, alpha = 1, nfolds = 10)

coef_diff_1se <- as.matrix(coef(cv_diff, s = "lambda.1se"))
selected_diff <- rownames(coef_diff_1se)[coef_diff_1se[, 1] != 0 &
                                          rownames(coef_diff_1se) != "(Intercept)"]
cat(sprintf("  Differential LASSO (lambda.1se): %s\n",
            paste(selected_diff, collapse = ", ")))

# LASSO coefficient path plot (binding-level)
p_62c <- {
  lambda_seq <- lasso_bind$lambda
  coef_mat <- as.matrix(coef(lasso_bind))[-1, ]
  path_df <- do.call(rbind, lapply(seq_len(ncol(coef_mat)), function(j) {
    data.frame(
      predictor = rownames(coef_mat),
      log_lambda = log(lambda_seq[j]),
      coefficient = coef_mat[, j],
      stringsAsFactors = FALSE
    )
  }))

  ggplot(path_df, aes(x = log_lambda, y = coefficient, color = predictor)) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = log(cv_bind$lambda.1se), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = log(cv_bind$lambda.min), linetype = "dotted", color = "grey40") +
    labs(
      title = "LASSO coefficient path (binding-level model)",
      subtitle = "Dashed = λ.1se | Dotted = λ.min",
      x = "log(λ)", y = "Coefficient", color = "Predictor"
    ) +
    theme_biomodal() +
    theme(legend.position = "right")
}

save_plot(p_62c, "62c_lasso_path", w = 12, h = 7)

# Save LASSO results
lasso_df <- data.frame(
  predictor = rownames(coef_1se)[-1],
  binding_coef_1se = coef_1se[-1, 1],
  binding_coef_min = coef_min[-1, 1],
  diff_coef_1se = coef_diff_1se[-1, 1],
  stringsAsFactors = FALSE
)

# =============================================================================
# 62d: STANDARDIZED COEFFICIENT FOREST PLOT
# =============================================================================

cat("\n--- 62d: Standardized coefficients ---\n")

# Scale predictors for binding full model
binding_df_z <- binding_df
for (col in binding_predictors_full) {
  binding_df_z[[col]] <- as.numeric(scale(binding_df_z[[col]]))
}

bind_m3_z <- lm(as.formula(paste("mecp2_ctrl_log ~",
                paste(binding_predictors_full, collapse = " + "))), data = binding_df_z)

std_coefs <- data.frame(summary(bind_m3_z)$coefficients, check.names = FALSE)
std_coefs$predictor <- rownames(std_coefs)
std_coefs <- std_coefs[std_coefs$predictor != "(Intercept)", ]
colnames(std_coefs) <- c("estimate", "se", "t_value", "p_value", "predictor")

ci <- confint(bind_m3_z)
std_coefs$ci_lo <- ci[std_coefs$predictor, 1]
std_coefs$ci_hi <- ci[std_coefs$predictor, 2]

std_coefs$mark_type <- ifelse(grepl("^cg_", std_coefs$predictor), "Methylation", "Chromatin")
std_coefs <- std_coefs[order(abs(std_coefs$estimate)), ]
std_coefs$predictor <- factor(std_coefs$predictor, levels = std_coefs$predictor)

cat("  Standardized betas (binding-level full model):\n")
for (i in seq_len(nrow(std_coefs))) {
  r <- std_coefs[i, ]
  cat(sprintf("    %-22s β=%+.4f [%+.4f, %+.4f]  %s\n",
              r$predictor, r$estimate, r$ci_lo, r$ci_hi, fmt_p(r$p_value)))
}

p_62d <- ggplot(std_coefs, aes(x = predictor, y = estimate, color = mark_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.7, linewidth = 0.8) +
  scale_color_manual(values = c("Methylation" = "#D6604D", "Chromatin" = "#4393C3"),
                     name = "Feature type") +
  coord_flip() +
  labs(
    title = "Standardized coefficients: full binding-level model",
    subtitle = "Which marks have the largest effect on MeCP2 binding?",
    x = NULL, y = "Standardized β (95% CI)"
  ) +
  theme_biomodal()

save_plot(p_62d, "62d_standardized_coefficients", w = 11, h = 7)

# Also do differential model standardized coefficients
diff_df_z <- diff_df
for (col in diff_predictors_full) {
  diff_df_z[[col]] <- as.numeric(scale(diff_df_z[[col]]))
}

diff_m3_z <- lm(as.formula(paste("Fold ~",
                paste(diff_predictors_full, collapse = " + "))), data = diff_df_z)

std_diff <- data.frame(summary(diff_m3_z)$coefficients, check.names = FALSE)
std_diff$predictor <- rownames(std_diff)
std_diff <- std_diff[std_diff$predictor != "(Intercept)", ]
colnames(std_diff) <- c("estimate", "se", "t_value", "p_value", "predictor")

cat("\n  Standardized betas (differential full model):\n")
for (i in seq_len(nrow(std_diff))) {
  r <- std_diff[nrow(std_diff) - i + 1, ]
  cat(sprintf("    %-22s β=%+.4f  %s\n",
              r$predictor, r$estimate, fmt_p(r$p_value)))
}

# =============================================================================
# 62e: CHROMATIN-STATE STRATIFIED R²
# =============================================================================

cat("\n--- 62e: Chromatin-state stratified R² ---\n")

# Classify peaks by chromatin state
peak_gr_all <- GRanges(seqnames = peaks$Chr,
                       ranges = IRanges(start = peaks$Start, end = peaks$End))

chip_overlaps <- data.frame(
  H3K27ac_overlap  = countOverlaps(peak_gr_all, load_chip_peaks(CHIP_PEAK_FILES$h3k27ac)) > 0,
  H3K27me3_overlap = countOverlaps(peak_gr_all, load_chip_peaks(CHIP_PEAK_FILES$h3k27me3)) > 0,
  H3K4me1_overlap  = countOverlaps(peak_gr_all, load_chip_peaks(CHIP_PEAK_FILES$h3k4me1)) > 0,
  H3K4me3_overlap  = countOverlaps(peak_gr_all, load_chip_peaks(CHIP_PEAK_FILES$h3k4me3)) > 0,
  Bivalent_overlap = countOverlaps(peak_gr_all, load_chip_peaks(CHIP_PEAK_FILES$bivalent)) > 0
)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
tss_gr <- GenomicFeatures::promoters(txdb, upstream = 0, downstream = 1)
dist_to_tss <- GenomicRanges::distanceToNearest(peak_gr_all, tss_gr)
peaks$dist_to_tss <- rep(NA_real_, nrow(peaks))
peaks$dist_to_tss[queryHits(dist_to_tss)] <- mcols(dist_to_tss)$distance

peaks$chromatin_state <- classify_chromatin_state(chip_overlaps, peaks$dist_to_tss)

state_counts <- table(peaks$chromatin_state)
cat("  Peak chromatin states:\n")
for (s in names(sort(state_counts, decreasing = TRUE))) {
  cat(sprintf("    %-20s %d\n", s, state_counts[s]))
}

# Fit CG-only and Full models within each state
state_r2 <- data.frame()
for (state in names(state_counts[state_counts >= 50])) {
  state_mask <- peaks$chromatin_state == state & diff_valid
  n <- sum(state_mask)
  if (n < 30) next

  s_df <- peaks[state_mask, ]
  m_cg <- tryCatch(lm(Fold ~ cg_mc_log2fc, data = s_df), error = function(e) NULL)
  m_full <- tryCatch(lm(as.formula(paste("Fold ~",
                     paste(diff_predictors_full, collapse = " + "))), data = s_df),
                     error = function(e) NULL)

  if (!is.null(m_cg) && !is.null(m_full)) {
    r2_cg <- summary(m_cg)$r.squared
    r2_full <- summary(m_full)$r.squared
    state_r2 <- rbind(state_r2, data.frame(
      chromatin_state = state, n = n,
      r2_cg_only = r2_cg, r2_full = r2_full,
      r2_gain = r2_full - r2_cg,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("    %-20s n=%5d  R²_CG=%.4f  R²_Full=%.4f  gain=+%.4f\n",
                state, n, r2_cg, r2_full, r2_full - r2_cg))
  }
}

# Plot 62e
if (nrow(state_r2) > 0) {
  state_long <- tidyr::pivot_longer(state_r2, cols = c(r2_cg_only, r2_full),
                                    names_to = "model", values_to = "r2")
  state_long$model <- ifelse(state_long$model == "r2_cg_only", "CG only", "Full")
  state_long$model <- factor(state_long$model, levels = c("CG only", "Full"))

  p_62e <- ggplot(state_long, aes(x = reorder(chromatin_state, -r2),
                                  y = r2, fill = model)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("CG only" = "#D6604D", "Full" = "#4393C3"),
                      name = "Model") +
    labs(
      title = "R² by chromatin state: CG-only vs Full model",
      subtitle = "Where does adding chromatin marks improve MeCP2 prediction most?",
      x = NULL, y = expression(R^2)
    ) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  save_plot(p_62e, "62e_chromatin_state_r2", w = 12, h = 7)
}

# =============================================================================
# 62f: SECTION 56 RESIDUAL vs CHROMATIN MARKS
# =============================================================================

cat("\n--- 62f: CG-corrected residual vs chromatin marks ---\n")

if (file.exists(SEC56_PEAKS)) {
  sec56 <- read.table(SEC56_PEAKS, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
  cat(sprintf("  Loaded section 56 peaks: %d rows\n", nrow(sec56)))

  peaks_with_resid <- merge(
    peaks[, c("Chr", "Start", "End", diff_predictors_chr)],
    sec56[, c("Chr", "Start", "End", "residual", "mecp2_class")],
    by = c("Chr", "Start", "End")
  )
  cat(sprintf("  Merged with residuals: %d peaks\n", nrow(peaks_with_resid)))

  resid_valid <- complete.cases(peaks_with_resid[, c("residual", diff_predictors_chr)])
  resid_df <- peaks_with_resid[resid_valid, ]

  lm_resid_chr <- lm(as.formula(paste("residual ~",
                      paste(diff_predictors_chr, collapse = " + "))), data = resid_df)
  r2_resid <- summary(lm_resid_chr)$r.squared
  cat(sprintf("  Residual ~ chromatin marks: R²=%.4f (n=%d)\n", r2_resid, nrow(resid_df)))

  cat("  Per-mark Spearman correlations with residual:\n")
  resid_cor <- data.frame()
  for (mk in diff_predictors_chr) {
    ct <- cor.test(resid_df$residual, resid_df[[mk]], method = "spearman")
    cat(sprintf("    %-18s rho=%+.4f  %s\n", mk, ct$estimate, fmt_p(ct$p.value)))
    resid_cor <- rbind(resid_cor, data.frame(
      mark = mk, rho = as.numeric(ct$estimate),
      p_value = ct$p.value, stringsAsFactors = FALSE
    ))
  }

  # Scatter: residual vs k27me3
  k27me3_rho <- resid_cor$rho[resid_cor$mark == "k27me3_log2fc"]
  k27me3_p   <- resid_cor$p_value[resid_cor$mark == "k27me3_log2fc"]

  p_62f <- ggplot(resid_df, aes(x = k27me3_log2fc, y = residual)) +
    geom_point(alpha = 0.15, size = 0.5, color = "grey50") +
    geom_smooth(method = "lm", color = "#4393C3", linewidth = 1, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(
      title = "CG-unexplained MeCP2 enrichment vs H3K27me3 change",
      subtitle = sprintf("Spearman ρ = %.4f, %s | Chromatin model R² = %.4f",
                         k27me3_rho, fmt_p(k27me3_p), r2_resid),
      x = "H3K27me3 log2FC", y = "MeCP2 residual (CG-unexplained)"
    ) +
    theme_biomodal()

  save_plot(p_62f, "62f_residual_vs_k27me3", w = 9, h = 8)
} else {
  cat("  WARNING: 56_mecp2_peak_annotated.tsv not found, skipping residual analysis\n")
}

# =============================================================================
# SAVE TABLES
# =============================================================================

cat("\n--- Saving tables ---\n")

write.table(model_summary, file.path(TABLES_DIR, "62_model_comparison_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 62_model_comparison_summary.tsv\n")

write.table(std_coefs, file.path(TABLES_DIR, "62_standardized_coefficients.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 62_standardized_coefficients.tsv\n")

write.table(lasso_df, file.path(TABLES_DIR, "62_lasso_selected_features.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 62_lasso_selected_features.tsv\n")

write.table(vp, file.path(TABLES_DIR, "62_variance_partition.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 62_variance_partition.tsv\n")

if (nrow(state_r2) > 0) {
  write.table(state_r2, file.path(TABLES_DIR, "62_chromatin_state_r2.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved: 62_chromatin_state_r2.tsv\n")
}

if (exists("resid_cor") && nrow(resid_cor) > 0) {
  write.table(resid_cor, file.path(TABLES_DIR, "62_residual_chromatin_regression.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved: 62_residual_chromatin_regression.tsv\n")
}

# =============================================================================
# COMPOSITE
# =============================================================================

cat("\n--- 62g: Composite ---\n")

p_62g <- (p_62a | p_62b) / (p_62d | p_62f) +
  plot_annotation(
    title = "Section 62: Multi-Feature Chromatin Regression at MeCP2 Peaks",
    subtitle = "Can chromatin features explain MeCP2 binding better than CG methylation?",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

save_plot(p_62g, "62g_composite", w = 20, h = 14)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 62 SUMMARY\n")
cat("================================================================================\n\n")

bind_r2_cg  <- model_summary$r_squared[model_summary$type == "Binding level" & model_summary$model == "CG only"]
bind_r2_chr <- model_summary$r_squared[model_summary$type == "Binding level" & model_summary$model == "Chromatin only"]
diff_r2_cg  <- model_summary$r_squared[model_summary$type == "Differential" & model_summary$model == "CG only"]
diff_r2_chr <- model_summary$r_squared[model_summary$type == "Differential" & model_summary$model == "Chromatin only"]

cat("  KEY RESULTS:\n")
cat(sprintf("    Binding level: Chromatin R²=%.4f %s CG R²=%.4f → %s\n",
            bind_r2_chr, ifelse(bind_r2_chr > bind_r2_cg, ">", "≤"), bind_r2_cg,
            ifelse(bind_r2_chr > bind_r2_cg,
                   "CHROMATIN explains MeCP2 binding better than methylation",
                   "CG methylation still dominant")))
cat(sprintf("    Differential:  Chromatin R²=%.4f %s CG R²=%.4f → %s\n",
            diff_r2_chr, ifelse(diff_r2_chr > diff_r2_cg, ">", "≤"), diff_r2_cg,
            ifelse(diff_r2_chr > diff_r2_cg,
                   "MeCP2 redistribution tracks CHROMATIN REMODELING",
                   "MeCP2 redistribution tracks METHYLATION CHANGE")))
cat(sprintf("    LASSO selected (binding): %s\n", paste(selected_1se, collapse = ", ")))
if (exists("r2_resid")) {
  cat(sprintf("    CG-unexplained residual → chromatin: R²=%.4f\n", r2_resid))
}

cat("\n  Plots saved to:", SEC62_DIR, "\n")
cat("  Tables saved to:", TABLES_DIR, "\n")

cat("\nSection 62 complete.\n")
